package org.oxh.deconvolution;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import org.oxh.deconvolution.core.DeconvolutionCore;
import org.oxh.deconvolution.core.PSFGenerator;

/**
 * ImageJ plugin for iterative deconvolution.
 * Supports both Richardson-Lucy (Poisson) and Accelerated Landweber (Gaussian).
 * Generates PSF from optical parameters and outputs both deconvolved image and PSF.
 */
public class Deconvolution_Plugin implements PlugIn {

    private static final double Z_DEFOCUS = 0.0;    // focal plane
    private static final double R_CUTOFF = 2.0;     // radial cutoff in micrometers

    private static final String ALGO_RL = "Richardson-Lucy (Poisson)";
    private static final String ALGO_LW = "Accelerated Landweber (Gaussian)";

    @Override
    public void run(String s) {
        ImagePlus imp = IJ.getImage();
        if (imp == null) {
            IJ.noImage();
            return;
        }

        // --- 1. Parameter dialog ---
        GenericDialog gd = new GenericDialog("Deconvolution Parameters");
        gd.addMessage("Iterative Deconvolution with theoretical PSF");
        gd.addNumericField("EX wavelength lambda (μm):", 0.561, 3);
        gd.addNumericField("Numerical Aperture NA:", 1.49, 2);
        gd.addNumericField("Pixel size (μm):", 0.0433, 4);
        gd.addNumericField("Iterations:", 10, 0);

        gd.addMessage("Algorithm selection:");
        gd.addRadioButtonGroup("Algorithm:",
                new String[]{ALGO_RL, ALGO_LW},
                2, 1,
                ALGO_RL);
        gd.addCheckbox("Dispaly PSF?", false);
        gd.showDialog();
        if (gd.wasCanceled()) return;

        double lambda = gd.getNextNumber();
        double NA = gd.getNextNumber();
        double pixel = gd.getNextNumber();
        int iteration = (int) gd.getNextNumber();

        String selectedAlgo = gd.getNextRadioButton();
        boolean ShowPSF = gd.getNextBoolean();
        int rule = ALGO_RL.equals(selectedAlgo) ? 1 : 2;

        String suffix = (rule == 1) ? "_DeconvPoisson" : "_DeconvGaussian";

        // --- 2. Generate PSF and display it ---
        IJ.showStatus("Generating PSF...");
        double[][] psf = PSFGenerator.generate(pixel, lambda, NA, Z_DEFOCUS, R_CUTOFF);
        if (ShowPSF){
            displayPSF(psf, lambda, NA, pixel);}

        // --- 3. Check for time series ---
        int numFrames = imp.getStackSize();
        boolean processAllFrames = false;

        if (numFrames > 1) {
            GenericDialog frameDialog = new GenericDialog("Process Frames");
            frameDialog.addMessage("This is a time-series image with " + numFrames + " frames.");
            frameDialog.addCheckbox("Process all frames?", true);
            frameDialog.showDialog();
            if (frameDialog.wasCanceled()) return;
            processAllFrames = frameDialog.getNextBoolean();
        }

        // --- 4. Process image(s) ---
        ImagePlus resultImp;

        if (processAllFrames) {
            resultImp = processAllFrames(imp, psf, iteration, rule, suffix);
        } else {
            resultImp = processSingleFrame(imp, psf, iteration, rule, suffix);
        }

        resultImp.show();
        IJ.showStatus("Deconvolution Complete");
    }

    // --- Single frame processing ---
    private ImagePlus processSingleFrame(ImagePlus imp, double[][] psf,
                                         int iteration, int rule, String suffix) {
        ImageProcessor ip = imp.getProcessor();
        int w = ip.getWidth();
        int h = ip.getHeight();

        double[][] img = loadImage(ip, w, h);
        double recordMax = getMax(img);

        IJ.showStatus("Running Deconvolution...");
        IJ.showProgress(0);
        double[][] deconv = DeconvolutionCore.iterativeDeblur(img, psf, iteration, rule);
        IJ.showProgress(1);

        double[][] normalized = normalizeToRange(deconv, recordMax);
        return createResultImage(normalized, imp.getTitle(), w, h, suffix);
    }

    // --- Stack processing ---
    private ImagePlus processAllFrames(ImagePlus imp, double[][] psf,
                                       int iteration, int rule, String suffix) {
        int w = imp.getWidth();
        int h = imp.getHeight();
        int numFrames = imp.getStackSize();
        ImageStack stack = imp.getStack();
        ImageStack resultStack = new ImageStack(w, h);

        for (int i = 1; i <= numFrames; i++) {
            IJ.showStatus("Deconvolving frame " + i + " of " + numFrames);
            IJ.showProgress(i, numFrames);

            ImageProcessor ip = stack.getProcessor(i);
            double[][] img = loadImage(ip, w, h);

            double[][] deconv = DeconvolutionCore.iterativeDeblur(img, psf, iteration, rule);
            double recordMax = getMax(img);
            double[][] normalized = normalizeToRange(deconv, recordMax);

            short[] uint16Pixels = convertToUint16(normalized, w, h);
            ShortProcessor sp = new ShortProcessor(w, h);
            sp.setPixels(uint16Pixels);
            resultStack.addSlice("Frame " + i, sp);
        }

        IJ.showStatus("Deconvolution Complete");
        IJ.showProgress(1.0);

        String originalTitle = imp.getTitle();
        String baseName = originalTitle.lastIndexOf(".") > 0
                ? originalTitle.substring(0, originalTitle.lastIndexOf("."))
                : originalTitle;
        String newTitle = baseName + suffix + ".tif";

        return new ImagePlus(newTitle, resultStack);
    }

    // --- Display PSF as separate 16-bit image ---
    private void displayPSF(double[][] psf, double lambda, double NA, double pixel) {
        int h = psf.length;
        int w = psf[0].length;

        double maxPsf = 0;
        for (int j = 0; j < h; j++) {
            for (int i = 0; i < w; i++) {
                if (psf[j][i] > maxPsf) maxPsf = psf[j][i];
            }
        }

        short[] pixels = new short[w * h];
        for (int j = 0; j < h; j++) {
            for (int i = 0; i < w; i++) {
//                int val = (int) Math.round(psf[j][i] / maxPsf * 65535);
                int val = (int) Math.round(psf[j][i] * 65535);
                pixels[j * w + i] = (short) Math.max(0, Math.min(65535, val));
            }
        }

        ShortProcessor sp = new ShortProcessor(w, h);
        sp.setPixels(pixels);

        String title = String.format("PSF_NA_%.2f_lambda_%.3fum_pixel_%.4fum.tif",
                NA, lambda, pixel);
        new ImagePlus(title, sp).show();
    }

    // --- Helper functions ---

    private double[][] loadImage(ImageProcessor ip, int w, int h) {
        double[][] img = new double[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                img[y][x] = ip.getPixel(x, y);
            }
        }
        return img;
    }

    private double getMax(double[][] img) {
        double max = 0;
        for (double[] row : img) {
            for (double v : row) {
                if (v > max) max = v;
            }
        }
        return max;
    }

    private double[][] normalizeToRange(double[][] img, double maxVal) {
        double min = Double.MAX_VALUE;
        double max = -Double.MAX_VALUE;
        int h = img.length;
        int w = img[0].length;

        for (int j = 0; j < h; j++) {
            for (int i = 0; i < w; i++) {
                if (img[j][i] < min) min = img[j][i];
                if (img[j][i] > max) max = img[j][i];
            }
        }

        double range = max - min;
        if (range < 1e-10) range = 1e-10;

        double[][] out = new double[h][w];
        for (int j = 0; j < h; j++) {
            for (int i = 0; i < w; i++) {
                out[j][i] = (img[j][i] - min) / range * maxVal;
            }
        }
        return out;
    }

    private ImagePlus createResultImage(double[][] result, String originalTitle,
                                        int w, int h, String suffix) {
        short[] uint16Pixels = convertToUint16(result, w, h);

        String baseName = originalTitle.lastIndexOf(".") > 0
                ? originalTitle.substring(0, originalTitle.lastIndexOf("."))
                : originalTitle;
        String newTitle = baseName + suffix + ".tif";

        ShortProcessor resultSp = new ShortProcessor(w, h);
        resultSp.setPixels(uint16Pixels);
        return new ImagePlus(newTitle, resultSp);
    }

    private short[] convertToUint16(double[][] result, int w, int h) {
        short[] uint16Pixels = new short[w * h];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                double val = result[y][x];
                int uint16Val = Math.max(0, Math.min(65535, (int) Math.round(val)));
                uint16Pixels[y * w + x] = (short) uint16Val;
            }
        }
        return uint16Pixels;
    }
}