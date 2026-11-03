package org.DenoiseStorm.htvdenoise;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ShortProcessor;
import ij.process.ImageProcessor;
import org.DenoiseStorm.htvdenoise.core.SparseHessianCore;

public class HTVDenoise_Plugin implements PlugIn {

    @Override
    public void run(String s) {
        ImagePlus imp = IJ.getImage();
        if (imp == null) {
            IJ.noImage();
            return;
        }

        // 1.get parameters
        GenericDialog gd = new GenericDialog("HTV Denoise Parameters");
//        gd.addMessage("Recommended range:Fidelity=150, Hess in (1,50), Paral1 in (1,50)");
        gd.addNumericField("Fidelity (lambda),keep 150:", 150, 0);
        gd.addNumericField("Hess (hessian),(1-50):", 10, 0);
        gd.addNumericField("paral1(Sparsity),(1-50):", 15, 0);
        gd.addNumericField("Max Iterations,(50-100):", 50, 0);
        gd.addCheckbox("Use Boost (faster/with fewer details):", false);
        gd.addNumericField("Rel Tolerance,(0.01-0.001):", 0.01, 4);
        gd.addCheckbox("Use ZeroPercent(The percentage of zeros in the background):", false);
        gd.addNumericField("ZeroPercent (0-1):", 0.2, 3);

        gd.showDialog();
        if (gd.wasCanceled()) return;

        float fidelity = (float) gd.getNextNumber();
        float hess = (float) gd.getNextNumber();
        float paral1 = (float) gd.getNextNumber();
        int iter = (int) gd.getNextNumber();
        boolean boost = gd.getNextBoolean();
        float rel = (float) gd.getNextNumber();
        boolean useZeroPercent = gd.getNextBoolean();
        float zeroPercent = (float) gd.getNextNumber();

        //Set mu to 1
        float mu = 1.0f;

        // 2.Check if it is a time series image
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

        // 3.process images
        ImagePlus resultImp;

        if (processAllFrames) {
            resultImp = processAllFrames(imp, fidelity, hess, paral1, mu, iter,
                    boost ? 1 : 0, rel, useZeroPercent, zeroPercent);
        } else {
            resultImp = processSingleFrame(imp, fidelity, hess, paral1, mu, iter,
                    boost ? 1 : 0, rel, useZeroPercent, zeroPercent);
        }

        // 4.Display Results
        resultImp.show();
    }

    //Process a single frame
    private ImagePlus processSingleFrame(ImagePlus imp, float fidelity, float hess, float paral1,
                                         float mu, int maxIter, int boost, float rel,
                                         boolean useZeroPercent, float zeroPercent) {
        ImageProcessor ip = imp.getProcessor();
        FloatProcessor fp = ip.convertToFloatProcessor();
        int w = fp.getWidth();
        int h = fp.getHeight();
        float[] pixels = (float[]) fp.getPixels();

        // Convert to 2D array
        float[][] img2D = new float[h][w];
        for (int y = 0; y < h; y++) {
            System.arraycopy(pixels, y * w, img2D[y], 0, w);
        }

        // Run the core algorithm (with a while loop)
        float[][] result = processWithWhileLoop(img2D, fidelity, hess, paral1, mu, maxIter,
                boost, rel, useZeroPercent, zeroPercent);

        // Switch back to 1D and display
        return createResultImage(result, imp.getTitle(), w, h);
    }

    //Process all frames
    private ImagePlus processAllFrames(ImagePlus imp, float fidelity, float hess, float paral1,
                                       float mu, int maxIter, int boost, float rel,
                                       boolean useZeroPercent, float zeroPercent) {
        int w = imp.getWidth();
        int h = imp.getHeight();
        int numFrames = imp.getStackSize();
        ImageStack stack = imp.getStack();
        ImageStack resultStack = new ImageStack(w, h);

        // processing frame-by-frame
        for (int i = 1; i <= numFrames; i++) {
            IJ.showStatus("Processing frame " + i + " of " + numFrames);
            IJ.showProgress(i, numFrames);

            ImageProcessor ip = stack.getProcessor(i);
            FloatProcessor fp = ip.convertToFloatProcessor();
            float[] pixels = (float[]) fp.getPixels();

            // 2D
            float[][] img2D = new float[h][w];
            for (int y = 0; y < h; y++) {
                System.arraycopy(pixels, y * w, img2D[y], 0, w);
            }

            // while
            float[][] result = processWithWhileLoop(img2D, fidelity, hess, paral1, mu, maxIter,
                    boost, rel, useZeroPercent, zeroPercent);

            // Switch back to 1D and add it to the result stack
            short[] uint16Pixels = convertToUint16(result, w, h);
            ShortProcessor sp = new ShortProcessor(w, h);
            sp.setPixels(uint16Pixels);
            resultStack.addSlice("Frame " + i, sp);
        }

        IJ.showStatus("HTV Denoise Complete");
        IJ.showProgress(1.0);

        // Create the resulting image
        String originalTitle = imp.getTitle();
        String baseName = originalTitle.lastIndexOf(".") > 0
                ? originalTitle.substring(0, originalTitle.lastIndexOf("."))
                : originalTitle;
        String newTitle = baseName + "_HTV.tif";

        return new ImagePlus(newTitle, resultStack);
    }

    // The processing logic with a while loop
    private float[][] processWithWhileLoop(float[][] img2D, float fidelity, float hess, float paral1,
                                           float mu, int maxIter, int boost, float rel,
                                           boolean useZeroPercent, float zeroPercent) {
        float[][] currentImg = img2D;
        boolean done = false;

        while (!done) {
            SparseHessianCore.ProcessResult result = SparseHessianCore.process(
                    currentImg, fidelity, hess, paral1, mu, maxIter, boost, rel,
                    useZeroPercent, zeroPercent
            );

            currentImg = result.image;

            // completed?
            if (useZeroPercent) {
                done = result.zeroPercentReached;
            } else {
                done = true; // The relative-change convergence criterion is completed in just one iteration.
            }
        }

        return currentImg;
    }

    // Create a single-frame result image
    private ImagePlus createResultImage(float[][] result, String originalTitle, int w, int h) {
        short[] uint16Pixels = convertToUint16(result, w, h);

        String baseName = originalTitle.lastIndexOf(".") > 0
                ? originalTitle.substring(0, originalTitle.lastIndexOf("."))
                : originalTitle;
        String newTitle = baseName + "_HTV.tif";

        ShortProcessor resultSp = new ShortProcessor(w, h);
        resultSp.setPixels(uint16Pixels);
        return new ImagePlus(newTitle, resultSp);
    }

    // convert to uint16
    private short[] convertToUint16(float[][] result, int w, int h) {
        float[] resultPixels = new float[w * h];
        for (int y = 0; y < h; y++) {
            System.arraycopy(result[y], 0, resultPixels, y * w, w);
        }

        short[] uint16Pixels = new short[w * h];
        for (int i = 0; i < resultPixels.length; i++) {
            float floatVal = resultPixels[i];
            int uint16Val = Math.max(0, Math.min(65535, Math.round(floatVal)));
            uint16Pixels[i] = (short) uint16Val;
        }
        return uint16Pixels;
    }
}