package org.oxh.deconvolution.core;

import java.util.Arrays;

public class DeconvolutionCore {

    private static final double EPS = 1e-6;

    /**
     * Iterative deblur wrapper. Normalizes input and calls deblurCore.
     */
    public static double[][] iterativeDeblur(double[][] data, double[][] psf,
                                             int iteration, int rule) {
        int rows = data.length;
        int cols = data[0].length;

        double maxVal = 0;
        for (int j = 0; j < rows; j++) {
            for (int i = 0; i < cols; i++) {
                if (data[j][i] > maxVal) maxVal = data[j][i];
            }
        }

        double[][] normalized = new double[rows][cols];
        for (int j = 0; j < rows; j++) {
            for (int i = 0; i < cols; i++) {
                normalized[j][i] = data[j][i] / maxVal;
            }
        }

        return deblurCore(normalized, psf, iteration, rule);
    }

    /**
     * Core deconvolution – exactly replicates MATLAB's deblur_core.
     */
    private static double[][] deblurCore(double[][] data, double[][] psf,
                                         int iteration, int rule) {
        // Normalize PSF (MATLAB's kernel = kernel ./ sum(kernel(:)))
        double sumPsf = 0;
        for (int j = 0; j < psf.length; j++)
            for (int i = 0; i < psf[0].length; i++)
                sumPsf += psf[j][i];
        double[][] kernel = new double[psf.length][psf[0].length];
        for (int j = 0; j < psf.length; j++)
            for (int i = 0; i < psf[0].length; i++)
                kernel[j][i] = psf[j][i] / sumPsf;

        int dx = data.length;
        int dy = data[0].length;
        int B = Math.min(dx, dy) / 6;

        // Pad data
        double[][] paddedData = padReplicate(data, B);
        int pRows = paddedData.length;
        int pCols = paddedData[0].length;

        // --- Place PSF into upper-left of a zeros matrix (size of data) ---
        double[][] kernelLarge = new double[pRows][pCols];
        int copyRows = Math.min(kernel.length, pRows);
        int copyCols = Math.min(kernel[0].length, pCols);
        for (int j = 0; j < copyRows; j++) {
            System.arraycopy(kernel[j], 0, kernelLarge[j], 0, copyCols);
        }

        // Find peak and circshift to origin (1,1) in MATLAB -> (0,0) in Java
        int peakR = 0, peakC = 0;
        double peakVal = -1;
        for (int j = 0; j < pRows; j++) {
            for (int i = 0; i < pCols; i++) {
                if (kernelLarge[j][i] > peakVal) {
                    peakVal = kernelLarge[j][i];
                    peakR = j;
                    peakC = i;
                }
            }
        }
        kernelLarge = circshift2D(kernelLarge, -peakR, -peakC);

        // OTF and its conjugate
        double[] otf = FFTUtil.fft2(kernelLarge);
        double[] otfConj = FFTUtil.complexConj(otf);

        // Initialization (MATLAB: yk = data; xk = zeros; vk = zeros)
        double[][] yk = copy2D(paddedData);
        double[][] xk = new double[pRows][pCols];
        double[][] vk = new double[pRows][pCols];

        if (rule == 2) {
            // ========== Accelerated Landweber (Nesterov) ==========
            double t = 1.0;
            double gamma1 = 1.0;
            double[][] xk_update = null;

            // Precompute FFT of data (used multiple times)
            double[] fftData = FFTUtil.fft2(paddedData);

            for (int i = 1; i <= iteration; i++) {
                if (i == 1) {
                    // MATLAB: xk_update = data;
                    xk_update = copy2D(paddedData);
                    // MATLAB: xk = data + t * ifftn( conj(otf) .* (fftn(data) - otf .* fftn(data)) )
                    double[] otfData = FFTUtil.complexMul(otf, fftData);
                    double[] diff = complexSubtract(fftData, otfData);
                    double[] corr = FFTUtil.complexMul(otfConj, diff);
                    double[][] grad = FFTUtil.ifft2Real(corr, pRows, pCols);
                    xk = addWeighted(paddedData, grad, 1.0, t);
                    // Note: No clamp here (MATLAB does not apply max(yk,1e-6) in first iteration)
                } else {
                    // MATLAB: gamma2 = 1/2 * sqrt(4*gamma1^2 + gamma1^4) - gamma1^2
                    double gamma2 = 0.5 * Math.sqrt(4 * gamma1 * gamma1 + Math.pow(gamma1, 4)) - gamma1 * gamma1;
                    double beta = -gamma2 * (1.0 - 1.0 / gamma1);

                    // yk_update = xk + beta*(xk - xk_update)
                    double[][] yk_update = addWeighted(xk, xk_update, 1.0 + beta, -beta);

                    // yk = yk_update + t * ifftn( conj(otf) .* (fftn(data) - otf .* fftn(yk_update)) )
                    double[] ykFFT = FFTUtil.fft2(yk_update);
                    double[] otfYk = FFTUtil.complexMul(otf, ykFFT);
                    double[] diff = complexSubtract(fftData, otfYk);
                    double[] corr = FFTUtil.complexMul(otfConj, diff);
                    double[][] grad = FFTUtil.ifft2Real(corr, pRows, pCols);
                    yk = addWeighted(yk_update, grad, 1.0, t);

                    // MATLAB: yk = max(yk,1e-6);
                    clampMin(yk, EPS);

                    // MATLAB: gamma1 = gamma2; xk_update = xk; xk = yk;
                    gamma1 = gamma2;
                    xk_update = xk;
                    xk = copy2D(yk);
                }
            }
            // After loop, yk holds the last estimate (if iteration>1)
            // For iteration==1, yk is still paddedData (MATLAB behavior)
        } else {
            // ========== Accelerated Richardson-Lucy (Biggs-Andrews) ==========
            double[][] xk_update;
            double[][] vk_update;

            for (int iter = 1; iter <= iteration; iter++) {
                // MATLAB: xk_update = xk;
                xk_update = xk;

                // --- RL update (exactly inline as MATLAB) ---
                // rliter = fftn( data ./ max(ifftn(otf .* fftn(yk)), 1e-6) )
                double[] ykFFT = FFTUtil.fft2(yk);
                double[] blurredFreq = FFTUtil.complexMul(otf, ykFFT);
                double[][] blurred = FFTUtil.ifft2Real(blurredFreq, pRows, pCols);
                clampMin(blurred, EPS);

                double[][] ratio = divide2D(paddedData, blurred);
                double[] ratioFFT = FFTUtil.fft2(ratio);
                double[] corrFreq = FFTUtil.complexMul(otfConj, ratioFFT);
                double[][] corrSpatial = FFTUtil.ifft2Real(corrFreq, pRows, pCols);

                // denom = real(ifftn(fftn(ones(size(data))) .* otf)) (recalculated each iteration)
                double[][] ones = constant2D(1.0, pRows, pCols);
                double[] onesFFT = FFTUtil.fft2(ones);
                double[] denomFreq = FFTUtil.complexMul(onesFFT, otf);
                double[][] denomSpatial = FFTUtil.ifft2Real(denomFreq, pRows, pCols);
                clampMin(denomSpatial, EPS);

                // xk = yk .* corrSpatial ./ denomSpatial
                xk = multiplyDivide(yk, corrSpatial, denomSpatial);
                clampMin(xk, EPS);

                // --- Biggs-Andrews acceleration ---
                vk_update = vk;
                vk = subtractPositive(xk, yk);

                double alpha;
                if (iter == 1) {
                    alpha = 0.0;
                } else {
                    double num = dotProduct(vk_update, vk);
                    double den = dotProduct(vk_update, vk_update) + EPS;
                    alpha = num / den;
                    alpha = Math.max(Math.min(alpha, 1.0), EPS);
                }

                // yk = xk + alpha * (xk - xk_update)
                yk = addWeighted(xk, xk_update, 1.0 + alpha, -alpha);
                clampMin(yk, EPS);

                // MATLAB: yk(isnan(yk)) = 1e-6;
                for (int j = 0; j < pRows; j++) {
                    for (int i = 0; i < pCols; i++) {
                        if (Double.isNaN(yk[j][i])) yk[j][i] = EPS;
                    }
                }
            }
        }

        // Final: yk(yk < 0) = 0; and crop border
        clampMin(yk, 0.0);
        return cropBorder(yk, B);
    }

    // ================== Helper functions ==================

    private static double[] complexSubtract(double[] a, double[] b) {
        int n = a.length;
        double[] result = new double[n];
        for (int i = 0; i < n; i++) {
            result[i] = a[i] - b[i];
        }
        return result;
    }

    private static double[][] padReplicate(double[][] src, int pad) {
        int h = src.length, w = src[0].length;
        int nh = h + 2 * pad, nw = w + 2 * pad;
        double[][] dst = new double[nh][nw];
        for (int j = 0; j < h; j++)
            System.arraycopy(src[j], 0, dst[j + pad], pad, w);
        for (int j = 0; j < pad; j++) {
            System.arraycopy(dst[pad], 0, dst[j], 0, nw);
            System.arraycopy(dst[pad + h - 1], 0, dst[pad + h + j], 0, nw);
        }
        for (int j = 0; j < nh; j++) {
            double left = dst[j][pad], right = dst[j][pad + w - 1];
            for (int i = 0; i < pad; i++) {
                dst[j][i] = left;
                dst[j][pad + w + i] = right;
            }
        }
        return dst;
    }

    private static double[][] circshift2D(double[][] src, int shiftR, int shiftC) {
        int rows = src.length, cols = src[0].length;
        double[][] dst = new double[rows][cols];
        for (int j = 0; j < rows; j++) {
            int jr = ((j + shiftR) % rows + rows) % rows;
            for (int i = 0; i < cols; i++) {
                int ic = ((i + shiftC) % cols + cols) % cols;
                dst[jr][ic] = src[j][i];
            }
        }
        return dst;
    }

    private static double[][] cropBorder(double[][] src, int border) {
        int h = src.length - 2 * border;
        int w = src[0].length - 2 * border;
        double[][] dst = new double[h][w];
        for (int j = 0; j < h; j++)
            System.arraycopy(src[j + border], border, dst[j], 0, w);
        return dst;
    }

    private static double[][] copy2D(double[][] src) {
        double[][] dst = new double[src.length][];
        for (int j = 0; j < src.length; j++) dst[j] = src[j].clone();
        return dst;
    }

    private static double[][] constant2D(double val, int rows, int cols) {
        double[][] m = new double[rows][cols];
        for (int j = 0; j < rows; j++) Arrays.fill(m[j], val);
        return m;
    }

    private static void clampMin(double[][] m, double minVal) {
        for (int j = 0; j < m.length; j++)
            for (int i = 0; i < m[0].length; i++)
                if (m[j][i] < minVal) m[j][i] = minVal;
    }

    private static double[][] multiplyDivide(double[][] a, double[][] b, double[][] c) {
        int rows = a.length, cols = a[0].length;
        double[][] r = new double[rows][cols];
        for (int j = 0; j < rows; j++)
            for (int i = 0; i < cols; i++)
                r[j][i] = a[j][i] * b[j][i] / c[j][i];
        return r;
    }

    private static double[][] divide2D(double[][] a, double[][] b) {
        int rows = a.length, cols = a[0].length;
        double[][] r = new double[rows][cols];
        for (int j = 0; j < rows; j++)
            for (int i = 0; i < cols; i++)
                r[j][i] = a[j][i] / b[j][i];
        return r;
    }

    private static double[][] addWeighted(double[][] a, double[][] b, double wa, double wb) {
        int rows = a.length, cols = a[0].length;
        double[][] r = new double[rows][cols];
        for (int j = 0; j < rows; j++)
            for (int i = 0; i < cols; i++)
                r[j][i] = wa * a[j][i] + wb * b[j][i];
        return r;
    }

    private static double[][] subtractPositive(double[][] a, double[][] b) {
        int rows = a.length, cols = a[0].length;
        double[][] r = new double[rows][cols];
        for (int j = 0; j < rows; j++)
            for (int i = 0; i < cols; i++)
                r[j][i] = Math.max(a[j][i] - b[j][i], EPS);
        return r;
    }

    private static double dotProduct(double[][] a, double[][] b) {
        double sum = 0;
        for (int j = 0; j < a.length; j++)
            for (int i = 0; i < a[0].length; i++)
                sum += a[j][i] * b[j][i];
        return sum;
    }
}