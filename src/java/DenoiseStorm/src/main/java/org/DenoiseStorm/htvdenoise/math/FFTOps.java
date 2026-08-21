package org.DenoiseStorm.htvdenoise.math;

import org.jtransforms.fft.FloatFFT_2D;

/**
 * FFT and frequency domain operation class, corresponding to MATLAB's fftn/ifftn and operation_* frequency domain kernel functions
 * Implemented based on the JTransforms library, fully aligned with the FFT logic of MATLAB
 */
public class FFTOps {

    // ------------------------------ Core FFT function ------------------------------
    /**
     * 2D Fast Fourier Transform, corresponding to MATLAB's fftn
     * @param data Input real matrix of float type [h][w]
     * @return Complex FFT result of float type [2][h][w], where [0] represents the real part and [1] represents the imaginary part
     */
    public static float[][][] fftn2D(float[][] data) {
        int h = data.length;
        int w = data[0].length;

        // Convert to the one-dimensional complex number format required by JTransforms：[re0, im0, re1, im1, ...]
        float[] complexData = new float[2 * h * w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                int idx = 2 * (y * w + x);
                complexData[idx] = data[y][x];  // Real
                complexData[idx + 1] = 0.0f;    // The imaginary part is initially zero.
            }
        }

        // Perform 2D complex FFT
        FloatFFT_2D fft = new FloatFFT_2D(h, w);
        fft.complexForward(complexData);

        // Convert back to a two-dimensional complex number array [real part / imaginary part][h][w]
        float[][][] result = new float[2][h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                int idx = 2 * (y * w + x);
                result[0][y][x] = complexData[idx];
                result[1][y][x] = complexData[idx + 1];
            }
        }
        return result;
    }

    /**
     * 2D inverse Fourier transform, corresponding to MATLAB's ifftn, returns the real part
     * @param complexData Input complex array of type float[2][h][w], [0] is the real part, [1] is the imaginary part
     * @return The real matrix after inverse transformation float[h][w]
     */
    public static float[][] ifftn2D(float[][][] complexData) {
        int h = complexData[0].length;
        int w = complexData[0][0].length;

        // Convert to the one-dimensional complex number format required by JTransforms
        float[] complex1D = new float[2 * h * w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                int idx = 2 * (y * w + x);
                complex1D[idx] = complexData[0][y][x];
                complex1D[idx + 1] = complexData[1][y][x];
            }
        }

        // Perform 2D complex inverse FFT. The second parameter, true, indicates automatic scaling (divided by h*w, energy)
        FloatFFT_2D fft = new FloatFFT_2D(h, w);
        fft.complexInverse(complex1D, true);

        // Extract the real part and ignore the extremely small imaginary part caused by floating-point numerical errors.
        float[][] result = new float[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                int idx = 2 * (y * w + x);
                result[y][x] = complex1D[idx];
            }
        }
        return result;
    }

    // ------------------------------ Frequency-domain kernel generating function ------------------------------
    /**
     * Frequency domain kernel in the xx direction, corresponding to MATLAB's operation_xx
     * The kernel is [1, -2, 1], and after FFT, take the square of the amplitude (fftn(kernel) .* conj(fftn(kernel)))
     * @param w Image width
     * @param h Image height
     * @return Frequency domain kernel matrix of type float[h][w]
     */
    public static float[][] operationXX(int w, int h) {
        // Construct the kernel [1,-2,1]
        float[][] kernel = new float[1][3];
        kernel[0][0] = 1.0f;
        kernel[0][1] = -2.0f;
        kernel[0][2] = 1.0f;

        // Padding
        float[][] paddedKernel = new float[h][w];
        paddedKernel[0][0] = kernel[0][0];
        paddedKernel[0][1] = kernel[0][1];
        paddedKernel[0][2] = kernel[0][2];

        // Calculate the square of the amplitude after FFT
        float[][][] fftKernel = fftn2D(paddedKernel);
        float[][] result = new float[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                float re = fftKernel[0][y][x];
                float im = fftKernel[1][y][x];
                result[y][x] = re * re + im * im; // square of the amplitude = FFT * conj(FFT)
            }
        }
        return result;
    }

    /**
     * Frequency domain kernel in the yy direction, corresponding to MATLAB's operation_yy
     * The kernel is [1; -2; 1], and the amplitude square is taken after FFT
     * @param w Image width
     * @param h Image height
     * @return Frequency domain kernel matrix of type float[h][w]
     */
    public static float[][] operationYY(int w, int h) {
        // Construct the kernel [1;-2;1]
        float[][] kernel = new float[3][1];
        kernel[0][0] = 1.0f;
        kernel[1][0] = -2.0f;
        kernel[2][0] = 1.0f;

        // Pad to the image size
        float[][] paddedKernel = new float[h][w];
        paddedKernel[0][0] = kernel[0][0];
        paddedKernel[1][0] = kernel[1][0];
        paddedKernel[2][0] = kernel[2][0];

        // Calculate the square of the amplitude after FFT
        float[][][] fftKernel = fftn2D(paddedKernel);
        float[][] result = new float[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                float re = fftKernel[0][y][x];
                float im = fftKernel[1][y][x];
                result[y][x] = re * re + im * im;
            }
        }
        return result;
    }

    /**
     * Frequency domain kernel in the xy direction, corresponding to MATLAB's operation_xy
     * The kernel is [[1,-1],[-1,1]], and the amplitude square is taken after FFT
     * @param w Image width
     * @param h Image height
     * @return Frequency domain kernel matrix of type float[h][w]
     */
    public static float[][] operationXY(int w, int h) {
        // construct the kernel [[1,-1],[-1,1]]
        float[][] kernel = new float[2][2];
        kernel[0][0] = 1.0f;
        kernel[0][1] = -1.0f;
        kernel[1][0] = -1.0f;
        kernel[1][1] = 1.0f;

        // padding
        float[][] paddedKernel = new float[h][w];
        paddedKernel[0][0] = kernel[0][0];
        paddedKernel[0][1] = kernel[0][1];
        paddedKernel[1][0] = kernel[1][0];
        paddedKernel[1][1] = kernel[1][1];

        // Calculate the square of the amplitude after FFT
        float[][][] fftKernel = fftn2D(paddedKernel);
        float[][] result = new float[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                float re = fftKernel[0][y][x];
                float im = fftKernel[1][y][x];
                result[y][x] = re * re + im * im;
            }
        }
        return result;
    }

    // ------------------------------ Frequency domain solution of the core function ------------------------------
    /**
     * Frequency domain solution, corresponding to MATLAB's g = real(ifftn(g_update ./ denom))
     * Handling the special logic for the first iteration's denominator
     * @param gUpdate the right-hand term g_update
     * @param denom pre-calculated frequency domain denominator normalize
    //     * @param isFirstIter whether it is the first iteration
    //     * @param fidelityOverMu the fixed denominator for the first iteration fidelity/mu
     * @return the solved image g
     */
    public static float[][] solveInFourierDomain(float[][] gUpdate, float[][] denom) {
        int h = gUpdate.length;
        int w = gUpdate[0].length;

        // Perform FFT on g_update
        float[][][] fftG = fftn2D(gUpdate);
        float[][] re = fftG[0];
        float[][] im = fftG[1];

        // Frequency domain division:FFT(g_update) ./ denom
        // The denominator is the pre-calculated matrix "denom"
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                float d = denom[y][x];
                if (d < 1e-10f) d = 1e-10f; // Avoid division by zero
                re[y][x] /= d;
                im[y][x] /= d;
            }
        }

        // Inverse FFT returns the real part
        return ifftn2D(new float[][][]{re, im});
    }
    // ------------------------------ Fundamental Matrix Tool Method ------------------------------
    /**
     * deep copy of matrix
     */
    public static float[][] copy(float[][] src) {
        int h = src.length;
        int w = src[0].length;
        float[][] dst = new float[h][w];
        for (int y = 0; y < h; y++) {
            System.arraycopy(src[y], 0, dst[y], 0, w);
        }
        return dst;
    }

    /**
     * Add the two matrices(element-wise)
     */
    public static float[][] add(float[][] a, float[][] b) {
        int h = a.length;
        int w = a[0].length;
        float[][] result = new float[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                result[y][x] = a[y][x] + b[y][x];
            }
        }
        return result;
    }

    /**
     * Subtract the two matrices(element-wise)
     */
    public static float[][] subtract(float[][] a, float[][] b) {
        int h = a.length;
        int w = a[0].length;
        float[][] result = new float[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                result[y][x] = a[y][x] - b[y][x];
            }
        }
        return result;
    }

    /**
     * Matrix element-wise multiplication by a scalar
     */
    public static float[][] multiply(float[][] src, float scalar) {
        int h = src.length;
        int w = src[0].length;
        float[][] result = new float[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                result[y][x] = src[y][x] * scalar;
            }
        }
        return result;
    }

    /**
     * Calculate the Frobenius norm of the matrix, corresponding to MATLAB's norm(g(:))
     */
    public static float norm(float[][] src) {
        float sumSq = 0.0f;
        for (float[] row : src) {
            for (float val : row) {
                sumSq += val * val;
            }
        }
        return (float) Math.sqrt(sumSq);
    }
}