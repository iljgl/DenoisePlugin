package org.DenoiseStorm.GeodesicFilter;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import org.DenoiseStorm.GeodesicFilter.Core.GeodesicFilterCore;

public class GeodesicFilter_Plugin implements PlugIn {

    @Override
    public void run(String s) {
        ImagePlus imp = IJ.getImage();
        if (imp == null) {
            IJ.noImage();
            return;
        }

        // 1.get parameters
        GenericDialog gd = new GenericDialog("Geodesic Filter Parameters");
        gd.addMessage("Geodesic filtering, based on Dijkstra's shortest path algorithm");
        gd.addNumericField("Window Size (Odd number, default 7):", 7, 0);
        gd.addNumericField("Alpha (Gray Scale Weight, default 10):", 10, 2);
        gd.addNumericField("Sigma (Gaussian standard deviation, default 1):", 1, 2);
        gd.showDialog();
        if (gd.wasCanceled()) return;

        int W = (int) gd.getNextNumber();
        double alpha = gd.getNextNumber();
        double sigma = gd.getNextNumber();

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
            resultImp = processAllFrames(imp, W, alpha, sigma);
        } else {
            resultImp = processSingleFrame(imp, W, alpha, sigma);
        }

        // 4.display result
        resultImp.show();
    }

    // Single frame
    private ImagePlus processSingleFrame(ImagePlus imp, int W, double alpha, double sigma) {
        ImageProcessor ip = imp.getProcessor();
        int w = ip.getWidth();
        int h = ip.getHeight();

        double recordMax = getRecordMax(ip);
        double[][] G = loadAndNormalize(ip, w, h, recordMax);

        IJ.showStatus("Running Geodesic Filter...");
        IJ.showProgress(0);
        double[][] filtered = GeodesicFilterCore.process(G, W, alpha, sigma);
        IJ.showProgress(1);

        return createResultImage(filtered, imp.getTitle(), w, h, recordMax, "_Geodesic.tif");
    }

    // All frames
    private ImagePlus processAllFrames(ImagePlus imp, int W, double alpha, double sigma) {
        int w = imp.getWidth();
        int h = imp.getHeight();
        int numFrames = imp.getStackSize();
        ImageStack stack = imp.getStack();
        ImageStack resultStack = new ImageStack(w, h);

        for (int i = 1; i <= numFrames; i++) {
            IJ.showStatus("Processing frame " + i + " of " + numFrames);
            IJ.showProgress(i, numFrames);

            ImageProcessor ip = stack.getProcessor(i);
            double recordMax = getRecordMax(ip);
            double[][] G = loadAndNormalize(ip, w, h, recordMax);

            double[][] filtered = GeodesicFilterCore.process(G, W, alpha, sigma);

            short[] uint16Pixels = convertToUint16(filtered, w, h, recordMax);
            ShortProcessor sp = new ShortProcessor(w, h);
            sp.setPixels(uint16Pixels);
            resultStack.addSlice("Frame " + i, sp);
        }

        IJ.showStatus("Geodesic Filter Complete");
        IJ.showProgress(1.0);

        String originalTitle = imp.getTitle();
        String baseName = originalTitle.lastIndexOf(".") > 0
                ? originalTitle.substring(0, originalTitle.lastIndexOf("."))
                : originalTitle;
        String newTitle = baseName + "_Geodesic.tif";

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

    private double[][] loadAndNormalize(ImageProcessor ip, int w, int h, double recordMax) {
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