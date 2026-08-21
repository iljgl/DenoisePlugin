package org.oxh.htvdenoise.math;

/**
 * Differential operation tool class, corresponding to the forward_diff and back_diff functions in MATLAB version.
 * Supports 2D image differentiation in the y-direction (dim = 1) and x-direction (dim = 2)
 */


public class DifferentialOps {

    /**
     * Forward difference, exactly corresponding to the forward_diff function in MATLAB
     * @param data Input 2D image of type float[h][w], where h represents the number of rows (y-axis) and w represents the number of columns (x-axis)
     * @param step Difference step size. The original code is fixed at 1
     * @param dim Difference dimension: 1 = y-direction (rows), 2 = x-direction (columns), corresponding to the dim parameter in MATLAB
     * @return Difference result, of the same size as the input
     */
    public static float[][] forwardDiff(float[][] data, float step, int dim) {
        int h = data.length;
        int w = data[0].length;
        float[][] out = new float[h][w];

        if (dim == 1) {
            // Forward difference in the y direction (row)：out[y][x] = data[y+1][x] - data[y][x]，the last row padded with zeros
            for (int y = 0; y < h - 1; y++) {
                for (int x = 0; x < w; x++) {
                    out[y][x] = (data[y + 1][x] - data[y][x]) / step;
                }
            }
            // The last row remain 0
        } else if (dim == 2) {
            // Forward difference in the x direction (col)：：out[y][x] = data[y][x+1] - data[y][x]，the last col padded with zeros
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w - 1; x++) {
                    out[y][x] = (data[y][x + 1] - data[y][x]) / step;
                }
            }
            // The last col remain 0
        }
        return out;
    }

    /**
     * Backward difference, exactly corresponding to the back_diff function in MATLAB
     * @param data Input 2D image of type float[h][w]
     * @param step Difference step size, the original code is fixed at 1
     * @param dim Difference dimension: 1 = y direction (rows), 2 = x direction (columns)
     * @return Difference result, the same size as the input
     */
    public static float[][] backDiff(float[][] data, float step, int dim) {
        int h = data.length;
        int w = data[0].length;
        float[][] out = new float[h][w];

        if (dim == 1) {
            // y(row)：out[y][x] = data[y][x] - data[y-1][x]
            for (int y = 1; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    out[y][x] = (data[y][x] - data[y - 1][x]) / step;
                }
            }
            // the first row keep 0
        } else if (dim == 2) {
            // x(col)：out[y][x] = data[y][x] - data[y][x-1]
            for (int y = 0; y < h; y++) {
                for (int x = 1; x < w; x++) {
                    out[y][x] = (data[y][x] - data[y][x - 1]) / step;
                }
            }
            // the first col keep 0
        }
        return out;
    }
}