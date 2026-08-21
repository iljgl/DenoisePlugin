package org.DenoiseStorm.tvgf.math;

import org.jtransforms.fft.DoubleFFT_2D;

public class FFTGradient {

    // --- 对应 Matlab 的 FFT 梯度计算 ---
    public static double[][][] computeGradient(double[][] G) {
        int h = G.length;
        int w = G[0].length;

        // 1. 创建核 (对应 Matlab 的 kernelx 和 kernely),采用中心差分
        double[][] kernelx = new double[h][w];
        double[][] kernely = new double[h][w];
        kernelx[0][w-1] = -0.5;
        if (w > 1) kernelx[0][1] = 0.5;
        if (h > 1) kernely[1][0] = 0.5;
        kernely[h-1][0] = -0.5;

        // 2. 准备 FFT 数据 (JTransforms 需要一维数组，实部虚部交替)
        DoubleFFT_2D fft = new DoubleFFT_2D(h, w);

        double[] fftG = new double[2 * h * w];
        double[] fftKernelX = new double[2 * h * w];
        double[] fftKernelY = new double[2 * h * w];

        // 填充实部
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                fftG[2 * (y * w + x)] = G[y][x];
                fftKernelX[2 * (y * w + x)] = kernelx[y][x];
                fftKernelY[2 * (y * w + x)] = kernely[y][x];
            }
        }

        // 3. 执行 FFT
        fft.complexForward(fftG);
        fft.complexForward(fftKernelX);
        fft.complexForward(fftKernelY);

        // 4. 频域相乘
        double[] fftGx = complexMultiply(fftG, fftKernelX, h, w);
        double[] fftGy = complexMultiply(fftG, fftKernelY, h, w);

        // 5. 逆 FFT
        fft.complexInverse(fftGx, true);
        fft.complexInverse(fftGy, true);

        // 6. 提取实部
        double[][] Gx = new double[h][w];
        double[][] Gy = new double[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                Gx[y][x] = fftGx[2 * (y * w + x)];
                Gy[y][x] = fftGy[2 * (y * w + x)];
            }
        }

        return new double[][][]{Gx, Gy};
    }

    // --- 复数乘法 ---
    private static double[] complexMultiply(double[] a, double[] b, int h, int w) {
        double[] res = new double[2 * h * w];
        for (int i = 0; i < h * w; i++) {
            double reA = a[2 * i];
            double imA = a[2 * i + 1];
            double reB = b[2 * i];
            double imB = b[2 * i + 1];
            res[2 * i] = reA * reB - imA * imB;
            res[2 * i + 1] = reA * imB + imA * reB;
        }
        return res;
    }
}