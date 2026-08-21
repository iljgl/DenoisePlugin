package org.oxh.deconvolution.core;

import org.jtransforms.fft.DoubleFFT_2D;

public class FFTUtil {

    /**
     * 2D FFT of a real image, output interleaved complex array.
     */
    public static double[] fft2(double[][] img) {
        int rows = img.length;
        int cols = img[0].length;
        double[] data = new double[rows * cols * 2];
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                data[2 * (r * cols + c)] = img[r][c];
            }
        }
        DoubleFFT_2D fft = new DoubleFFT_2D(rows, cols);
        fft.complexForward(data);
        return data;
    }

    /**
     * MATLAB-like fftn(X, sz): truncate/pad X to (rows x cols) then FFT.
     * Truncation takes top-left portion if X is larger than sz.
     */
    public static double[] fftn(double[][] x, int rows, int cols) {
        double[][] padded = new double[rows][cols];
        int copyRows = Math.min(x.length, rows);
        int copyCols = Math.min(x[0].length, cols);
        for (int r = 0; r < copyRows; r++) {
            System.arraycopy(x[r], 0, padded[r], 0, copyCols);
        }
        return fft2(padded);
    }

    /**
     * Inverse FFT, returning real part only (scaled).
     */
    public static double[][] ifft2Real(double[] data, int rows, int cols) {
        double[] copy = data.clone();
        DoubleFFT_2D fft = new DoubleFFT_2D(rows, cols);
        fft.complexInverse(copy, true);
        double[][] result = new double[rows][cols];
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                result[r][c] = copy[2 * (r * cols + c)];
            }
        }
        return result;
    }

    /**
     * Element-wise complex multiplication: a .* b
     */
    public static double[] complexMul(double[] a, double[] b) {
        int n = a.length / 2;
        double[] result = new double[a.length];
        for (int i = 0; i < n; i++) {
            double ar = a[2 * i], ai = a[2 * i + 1];
            double br = b[2 * i], bi = b[2 * i + 1];
            result[2 * i]     = ar * br - ai * bi;
            result[2 * i + 1] = ar * bi + ai * br;
        }
        return result;
    }

    /**
     * Complex conjugate.
     */
    public static double[] complexConj(double[] a) {
        double[] result = a.clone();
        for (int i = 0; i < a.length / 2; i++) {
            result[2 * i + 1] = -result[2 * i + 1];
        }
        return result;
    }
}