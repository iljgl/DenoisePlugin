package org.DenoiseStorm.deconvolution.core;

/**
 * Generate scalar diffraction PSF based on Debye integral approximation.
 * Equivalent to Matlab's Generate_PSF + kernel functions.
 */
public class PSFGenerator {

    private static final int INTEGRATION_POINTS = 200; // Simpson integration resolution

    /**
     * Generate normalized PSF.
     * @param pixel pixel size in micrometers
     * @param lambda emission wavelength in micrometers
     * @param NA numerical aperture
     * @param z defocus distance in micrometers (0 = focal plane)
     * @param rCutoff radial cutoff in micrometers
     * @return normalized 2D PSF (sum = 1)
     */
    public static double[][] generate(double pixel, double lambda, double NA,
                                      double z, double rCutoff) {
        // Compute half-size of PSF grid in pixels
//        int nn = (int) Math.ceil(rCutoff / pixel) + 1;
        int nn = 64;
        int size = 2 * nn + 1;

        double sin2 = NA * NA / 2.0;
        double u = 8 * Math.PI * z * sin2 / lambda;

        double[][] IP = new double[size][size];

        for (int j = 0; j < size; j++) {
            double y = (j - nn) * pixel;
            for (int i = 0; i < size; i++) {
                double x = (i - nn) * pixel;
                double r = Math.sqrt(x * x + y * y);

                if (r > rCutoff) {
                    IP[j][i] = 0;
                } else {
                    IP[j][i] = integrateRadial(r, NA, lambda, u);
                }
            }
        }

        // PSF = |amplitude|^2
        double[][] psf = new double[size][size];
        double sum = 0;
        for (int j = 0; j < size; j++) {
            for (int i = 0; i < size; i++) {
                psf[j][i] = IP[j][i] * IP[j][i];
                sum += psf[j][i];
            }
        }

        // Normalize so that sum(psf) = 1
        for (int j = 0; j < size; j++) {
            for (int i = 0; i < size; i++) {
                psf[j][i] /= sum;
            }
        }

        return psf;
    }

    /**
     * Integrate h(r, p) over p in [0, 1] using Simpson's rule.
     * h(r,p) = 2 * exp(i*u*p^2/2) * besselj(0, 2*pi*r*NA/lambda * p)
     * Returns the real magnitude of the integral result.
     */
    private static double integrateRadial(double r, double NA, double lambda, double u) {
        double k = 2 * Math.PI * r * NA / lambda;
        double halfU = u / 2.0;

        int n = INTEGRATION_POINTS;
        if (n % 2 != 0) n++; // ensure even number of intervals

        double h = 1.0 / n;

        // Real and imaginary parts of the integral
        double sumReal = 0;
        double sumImag = 0;

        // Endpoints (p=0 and p=1)
        // At p=0: besselj(0,0)=1, exp(0)=1 → integrand = 2
        sumReal += 2.0;
        // At p=1
        double p1 = 1.0;
        double phase1 = halfU * p1 * p1;
        double bess1 = besselJ0(k * p1);
        sumReal += 2.0 * Math.cos(phase1) * bess1;
        sumImag += 2.0 * Math.sin(phase1) * bess1;

        // Interior points
        for (int i = 1; i < n; i++) {
            double p = i * h;
            double phase = halfU * p * p;
            double bess = besselJ0(k * p);
            double weight = (i % 2 == 0) ? 2.0 : 4.0;
            sumReal += weight * 2.0 * Math.cos(phase) * bess;
            sumImag += weight * 2.0 * Math.sin(phase) * bess;
        }

        double integralReal = sumReal * h / 3.0;
        double integralImag = sumImag * h / 3.0;

        // Return magnitude
        return Math.sqrt(integralReal * integralReal + integralImag * integralImag);
    }

    /**
     * Zero-order Bessel function of the first kind J0(x).
     * Abramowitz & Stegun approximation, accuracy ~ 1e-8.
     */
    private static double besselJ0(double x) {
        x = Math.abs(x);
        if (x < 8.0) {
            double y = x * x;
            double num = 57568490574.0 + y * (-13362590354.0 + y * (651619640.7
                    + y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
            double den = 57568490411.0 + y * (1029532985.0 + y * (9494680.718
                    + y * (59272.64853 + y * (267.8532712 + y * 1.0))));
            return num / den;
        } else {
            double z = 8.0 / x;
            double y = z * z;
            double xx = x - 0.785398164; // pi/4
            double p = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
                    + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
            double q = z * (-0.1562499995e-1 + y * (0.1430488765e-3
                    + y * (-0.6911147651e-5 + y * (0.7621095161e-6 - y * 0.934935152e-7))));
            return Math.sqrt(0.636619772 / x) * (Math.cos(xx) * p - Math.sin(xx) * q);
        }
    }
}