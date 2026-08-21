package org.DenoiseStorm.tvgf;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import org.DenoiseStorm.tvgf.core.TVGFCore;

public class TVGF_Plugin implements PlugIn {

    @Override
    public void run(String s) {
        ImagePlus imp = IJ.getImage();
        if (imp == null) {
            IJ.noImage();
            return;
        }

        // 1.get parameters
        GenericDialog gd = new GenericDialog("TVGF Parameters");
        gd.addMessage("U=G-S/(1-S),S=aG+b,ε(a,b)=|S-S0|^2+Lambda*|a∇S|");
        gd.addNumericField("S0 WindowSize(1%~10% of image size),:", 10, 0);
        gd.addNumericField("Statistics WindowSize(Usually equals to S0):", 10, 0);
        gd.addNumericField("Lambda(Higher,S near to S0),0.1~0.001:", 0.01, 4);
        gd.showDialog();
        if (gd.wasCanceled()) return;

        int s0WinSize = (int) gd.getNextNumber();
        int staticWinSize = (int) gd.getNextNumber();
        double lambda = gd.getNextNumber();

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

        // 3.processing
        ImagePlus resultImp;

        if (processAllFrames) {
            resultImp = processAllFrames(imp, s0WinSize, staticWinSize, lambda);
        } else {
            resultImp = processSingleFrame(imp, s0WinSize, staticWinSize, lambda);
        }

        // 4.display result
        resultImp.show();
    }

    // Single frame
    private ImagePlus processSingleFrame(ImagePlus imp, int s0WinSize, int staticWinSize, double lambda) {
        ImageProcessor ip = imp.getProcessor();
        int w = ip.getWidth();
        int h = ip.getHeight();

        // Convert to double[][] and normalize
        double[][] G = loadAndNormalize(ip, w, h);
        double recordMax = getRecordMax(ip);

        // Core
        IJ.showStatus("Running TVGF...");
        double[][] U = TVGFCore.process(G, s0WinSize, staticWinSize, lambda);

        // create the results
        return createResultImage(U, imp.getTitle(), w, h, recordMax, "_TVGF.tif");
    }

    // All frames
    private ImagePlus processAllFrames(ImagePlus imp, int s0WinSize, int staticWinSize, double lambda) {
        int w = imp.getWidth();
        int h = imp.getHeight();
        int numFrames = imp.getStackSize();
        ImageStack stack = imp.getStack();
        ImageStack resultStack = new ImageStack(w, h);

        for (int i = 1; i <= numFrames; i++) {
            IJ.showStatus("Processing frame " + i + " of " + numFrames);
            IJ.showProgress(i, numFrames);

            ImageProcessor ip = stack.getProcessor(i);
            double[][] G = loadAndNormalize(ip, w, h);
            double recordMax = getRecordMax(ip);

            double[][] U = TVGFCore.process(G, s0WinSize, staticWinSize, lambda);

            short[] uint16Pixels = convertToUint16(U, w, h, recordMax);
            ShortProcessor sp = new ShortProcessor(w, h);
            sp.setPixels(uint16Pixels);
            resultStack.addSlice("Frame " + i, sp);
        }

        IJ.showStatus("TVGF Complete");
        IJ.showProgress(1.0);

        String originalTitle = imp.getTitle();
        String baseName = originalTitle.lastIndexOf(".") > 0
                ? originalTitle.substring(0, originalTitle.lastIndexOf("."))
                : originalTitle;
        String newTitle = baseName + "_TVGF.tif";

        return new ImagePlus(newTitle, resultStack);
    }

    // --- auxiliary functions ---

    private double getRecordMax(ImageProcessor ip) {
        int w = ip.getWidth();
        int h = ip.getHeight();
        double recordMax = 0;
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                int val = ip.getPixel(x, y);
                if (val > recordMax) recordMax = val;
            }
        }
        return recordMax;
    }

    private double[][] loadAndNormalize(ImageProcessor ip, int w, int h) {
        double recordMax = getRecordMax(ip);
        double[][] G = new double[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                G[y][x] = ip.getPixel(x, y) / recordMax;
            }
        }
        return G;
    }

    private ImagePlus createResultImage(double[][] result, String originalTitle,
                                        int w, int h, double recordMax, String suffix) {
        short[] uint16Pixels = convertToUint16(result, w, h, recordMax);

        String baseName = originalTitle.lastIndexOf(".") > 0
                ? originalTitle.substring(0, originalTitle.lastIndexOf("."))
                : originalTitle;
        String newTitle = baseName + suffix;

        ShortProcessor resultSp = new ShortProcessor(w, h);
        resultSp.setPixels(uint16Pixels);
        return new ImagePlus(newTitle, resultSp);
    }

    private short[] convertToUint16(double[][] result, int w, int h, double recordMax) {
        short[] uint16Pixels = new short[w * h];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                double floatVal = result[y][x] * recordMax;
                int uint16Val = Math.max(0, Math.min(65535, (int) Math.round(floatVal)));
                uint16Pixels[y * w + x] = (short) uint16Val;
            }
        }
        return uint16Pixels;
    }
}