package org.DenoiseStorm.tvgf.core;
import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import org.DenoiseStorm.tvgf.math.FFTGradient;

public class TVGFCore {
    private static final double ROF_THETA=25.0;
    public static double[][] process(double[][] G, int s0WinSize, int staticWinSize, double lambda , boolean IfDenoise) {
        int h = G.length;
        int w = G[0].length;

        // --- Step 1: Calculate the initial S0 (quarter-window) ---
        double[][] S0LT = calculateInitialS0LeftTop(G, s0WinSize, s0WinSize);

        // --- Step 2: Offset Correction (Corresponding to the S0LT1 calculation in Matlab) ---
        double[][] S0 = shiftS0(S0LT, s0WinSize);
        // --- step 2.5: ROF denoise to get the guided image
        if (IfDenoise){
            G = rofDenoise(G,G[0].length,G.length,ROF_THETA);}
        // --- Step 3: Eliminate scattering by TVGF ---
        return TVGFRemoveScattering(G, S0, staticWinSize, lambda);
    }

    // --- step 2.5
    private static double[][] rofDenoise(double[][] source, int w, int h, double theta) {
        // ---1. double[][] -> float[]
        float[] pixels = new float[w * h];
        for (int y = 0; y < h; y++) {
            double[] srcRow = source[y];
            int base = y * w;
            for (int x = 0; x < w; x++) {
                double v = srcRow[x];
                // [0,1]
                if (v < 0.0) {
                    v = 0.0;
                } else if (v > 1.0) {
                    v = 1.0;
                }
                pixels[base + x] = (float) v;
            }
        }
        // --- 2. 32-bit ImagePlus ---
        FloatProcessor fp = new FloatProcessor(w, h, pixels, null);
        ImagePlus tempImp = new ImagePlus("TVGF_ROF_tmp", fp);
        // --- 3. ImageJ ROF Denoise ---
        ImagePlus beforeCurrent = WindowManager.getCurrentImage();
        IJ.run(tempImp, "ROF Denoise", "theta=" + theta);
        // --- 4. ROF ---
        ImagePlus outImp = tempImp;
        ImagePlus afterCurrent = WindowManager.getCurrentImage();
        //
        if (afterCurrent != null && afterCurrent != tempImp && afterCurrent != beforeCurrent) {
            outImp = afterCurrent;
        } else {
            //
            int[] ids = WindowManager.getIDList();
            if (ids != null) {
                for (int id : ids) {
                    ImagePlus img = WindowManager.getImage(id);
                    if (img != null && img != tempImp && img != beforeCurrent) {
                        String title = img.getTitle() == null ? "" : img.getTitle().toLowerCase();
                        if (title.contains("rof") || title.contains("denoise")) {
                            outImp = img;
                            break;
                        }
                    }
                }
            }
        }
        // ---5.  double[][] ---
        ImageProcessor outIp = outImp.getProcessor();
        if (!(outIp instanceof FloatProcessor)) {
            outIp = outIp.convertToFloat();
        }
        float[] outPixels = (float[]) outIp.getPixels();
        double[][] denoised = new double[h][w];
        for (int y = 0; y < h; y++) {
            int base = y * w;
            double[] outRow = denoised[y];
            for (int x = 0; x < w; x++) {
                double v = outPixels[base + x];
                //
                if (v < 0.0) {
                    v = 0.0;
                } else if (v > 1.0) {
                    v = 1.0;
                }
                outRow[x] = v;
            }
        }
        // ---6---
        if (outImp != null && outImp != tempImp) {
            outImp.changes = false;
            outImp.close();
        }
        tempImp.changes = false;
        tempImp.close();
        return denoised;
    }
    // --- Correspondingly, for my Matlab's "CalculateInitialS0_LeftTop", after testing, there is almost no difference between LT, LB, RT, and RB. ---
    private static double[][] calculateInitialS0LeftTop(double[][] G, int WinRow, int WinCol) {
        int h = G.length;
        int w = G[0].length;
        double[][] S0 = new double[h][w];
        int halfWinRow = WinRow - 1;
        int halfWinCol = WinCol - 1;

        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                int r_start, r_end;
                if (y + halfWinRow < h) {
                    r_start = y;
                    r_end = y + halfWinRow;
                } else {
                    r_start = y - halfWinRow;
                    r_end = y;
                    if (r_start < 0) {
                        r_start = 0;
                        r_end = h - 1;
                    }
                }

                int c_start, c_end;
                if (x + halfWinCol < w) {
                    c_start = x;
                    c_end = x + halfWinCol;
                } else {
                    c_start = x - halfWinCol;
                    c_end = x;
                    if (c_start < 0) {
                        c_start = 0;
                        c_end = w - 1;
                    }
                }

                double minVal = Double.MAX_VALUE;
                for (int ry = r_start; ry <= r_end; ry++) {
                    for (int cx = c_start; cx <= c_end; cx++) {
                        if (G[ry][cx] < minVal) {
                            minVal = G[ry][cx];
                        }
                    }
                }
                S0[y][x] = minVal;
            }
        }
        return S0;
    }

    // --- Corresponding to the S0LT1 offset correction in Matlab ---
    private static double[][] shiftS0(double[][] S0LT, int s0WinSize) {
        int h = S0LT.length;
        int w = S0LT[0].length;
        double[][] S0 = new double[h][w];

        int s_r = s0WinSize - 1;
        int s_c = s0WinSize;
        int ssc = (int) Math.round((double) s_c / 2);
        int ssr = ssc - 1;

        for (int y = ssc; y < h; y++) {
            for (int x = ssc; x < w; x++) {
                int srcY = y - ssc;
                int srcX = x - ssc;
                if (srcY < h - ssr && srcX < w - ssr) {
                    S0[y][x] = S0LT[srcY][srcX];
                }
            }
        }
        return S0;
    }

    // --- TVGFRemoveScattering=TVGF_Dehaze ---
    private static double[][] TVGFRemoveScattering(double[][] G, double[][] S0, int WinSize, double lambda) {
        int h = G.length;
        int w = G[0].length;
        int halfWin = WinSize / 2;

        // 1. Pre-compute the 1-norm of the gradient of G (corresponding to the FFT implementation in Matlab)
        double[][][] grad = FFTGradient.computeGradient(G);
        double[][] G_x = grad[0];
        double[][] G_y = grad[1];
        double[][] absGradG = new double[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                absGradG[y][x] = Math.abs(G_x[y][x]) + Math.abs(G_y[y][x]);
            }
        }

        // 2. Initialize output graph
        double[][] S = new double[h][w];
        double[][] a_map = new double[h][w];
        double[][] b_map = new double[h][w];

        // 3. Sliding window traversal
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                int r_start = Math.max(0, y - halfWin);
                int r_end = Math.min(h - 1, y + halfWin);
                int c_start = Math.max(0, x - halfWin);
                int c_end = Math.min(w - 1, x + halfWin);

                // Extract window data
                int winH = r_end - r_start + 1;
                int winW = c_end - c_start + 1;
                int currentN = winH * winW;
                double[] winG = new double[currentN];
                double[] winS0 = new double[currentN];
                double[] winGrad = new double[currentN];
                int idx = 0;
                for (int ry = r_start; ry <= r_end; ry++) {
                    for (int cx = c_start; cx <= c_end; cx++) {
                        winG[idx] = G[ry][cx];
                        winS0[idx] = S0[ry][cx];
                        winGrad[idx] = absGradG[ry][cx];
                        idx++;
                    }
                }

                // Calculate 5 statistical quantities
                double meanG = mean(winG);
                double meanS0 = mean(winS0);
                double meanGrad = mean(winGrad);

                // centering
                double[] centG = subtract(winG, meanG);
                double[] centS0 = subtract(winS0, meanS0);
                double varG = sumOfSquares(centG) / currentN;
                double covGS = sumOfProducts(centG, centS0) / currentN;

                // Soft threshold T
                double T = lambda * meanGrad;

                // a_xy
                double a_xy;
                if (varG < 1e-8) {
                    a_xy = 0;
                } else {
                    if (covGS > T) {
                        a_xy = (covGS - T) / varG;
                    } else if (covGS < -T) {
                        a_xy = (covGS + T) / varG;
                    } else {
                        a_xy = 0;
                    }
                }

                // b_xy
                double b_xy = meanS0 - a_xy * meanG;

                // S_xy
                double S_xy = a_xy * G[y][x] + b_xy;
                S_xy = Math.max(0, Math.min(1, S_xy));

                S[y][x] = S_xy;
                a_map[y][x] = a_xy;
                b_map[y][x] = b_xy;
            }
        }

        // 4. Calculate the intermediate graph V and the final graph U
        double[][] V = new double[h][w];
        double maxV = 0;
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                double denominator = Math.max(1 - S[y][x], 1e-8);
                V[y][x] = (G[y][x] - S[y][x]) / denominator;
                if (V[y][x] > maxV) maxV = V[y][x];
            }
        }

        double alpha = (maxV < 1e-8) ? 1 : 1 / maxV;
        double[][] U = new double[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                U[y][x] = alpha * V[y][x];
                U[y][x] = Math.max(0, Math.min(1, U[y][x]));
            }
        }

        return U;
    }

    // --- Auxiliary mathematical functions ---
    private static double mean(double[] arr) {
        double sum = 0;
        for (double v : arr) sum += v;
        return sum / arr.length;
    }

    private static double[] subtract(double[] arr, double val) {
        double[] res = new double[arr.length];
        for (int i = 0; i < arr.length; i++) res[i] = arr[i] - val;
        return res;
    }

    private static double sumOfSquares(double[] arr) {
        double sum = 0;
        for (double v : arr) sum += v * v;
        return sum;
    }

    private static double sumOfProducts(double[] a, double[] b) {
        double sum = 0;
        for (int i = 0; i < a.length; i++) sum += a[i] * b[i];
        return sum;
    }
}