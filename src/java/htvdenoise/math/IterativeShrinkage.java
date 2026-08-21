package org.oxh.htvdenoise.math;

/**
 * Refine the contraction threshold tool class, corresponding to the 4 iter_* core iterative functions in MATLAB
 * Implement the ADMM iterative update for the Hessian regularization term and the L1 sparse regularization term
 */

//Note: In Java, the default initial value for basic data types is generally 0.False,\u0000=0 or \uffff=65535
public class IterativeShrinkage {

    /**
     * Iteration result encapsulation class, returning the updated L term and the dual variable b
     */
    public static class IterResult {
        public float[][] L;  // The linear terms generated after iteration,used to update g_update
        public float[][] b;  // The updated dual variables,used for the next iteration

        public IterResult(float[][] L, float[][] b) {
            this.L = L;
            this.b = b;
        }
    }

    // ------------------------------ Core iterative function ------------------------------
    /**
     * Second-order derivative iteration in the xx direction (Hessian term), corresponding to MATLAB's iter_xx
     * @param g The current iteration image
     * @param bxx The dual variable from the previous iteration
     * @param para The weight of the regularization term, originally fixed at 1
     * @param mu The penalty parameter, originally set to 1
     * @return The iteration result Lxx and the updated bxx
     */
    public static IterResult iterXX(float[][] g, float[][] bxx, float para, float mu) {
        int h = g.length;
        int w = g[0].length;

        // Calculate the second derivative gxx = back_diff(forward_diff(g,1,1),1,1)
        float[][] fx = DifferentialOps.forwardDiff(g, 1.0f, 1);
        float[][] gxx = DifferentialOps.backDiff(fx, 1.0f, 1);

        // soft threshold shrinkage:signdxx = sign(gxx+bxx) * max(|gxx+bxx| - 1/mu, 0)
        float[][] dxx = new float[h][w];
        float threshold = 1.0f / mu;
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                float z = gxx[y][x] + bxx[y][x];
                float absZ = Math.abs(z);
                if (absZ > threshold) {
                    dxx[y][x] = Math.signum(z) * (absZ - threshold);
                } else {
                    dxx[y][x] = 0.0f;
                }
            }
        }

        // Update the dual variables bxx = bxx + (gxx - dxx)
        float[][] newBxx = new float[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                newBxx[y][x] = bxx[y][x] + (gxx[y][x] - dxx[y][x]);
            }
        }

        // Calculate Lxx = para * back_diff(forward_diff(dxx - newBxx, 1, 1), 1, 1)
        float[][] diff = new float[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                diff[y][x] = dxx[y][x] - newBxx[y][x];
            }
        }
        float[][] fx2 = DifferentialOps.forwardDiff(diff, 1.0f, 1);
        float[][] Lxx = DifferentialOps.backDiff(fx2, 1.0f, 1);
        // Multiply by the weight(para)
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                Lxx[y][x] *= para;
            }
        }

        return new IterResult(Lxx, newBxx);
    }

    /**
     * Second-order derivative iteration in the yy direction (Hessian term), corresponding to MATLAB's iter_yy
     * @param g The current iteration image
     * @param byy The dual variable from the previous iteration
     * @param para The regularization term weight, originally fixed at 1
     * @param mu The penalty parameter, originally set to 1
     * @return The iteration result Lyy and the updated byy
     */
    public static IterResult iterYY(float[][] g, float[][] byy, float para, float mu) {
        int h = g.length;
        int w = g[0].length;

        // Calculate the second derivative gyy = back_diff(forward_diff(g,1,2),1,2)
        float[][] fx = DifferentialOps.forwardDiff(g, 1.0f, 2);
        float[][] gyy = DifferentialOps.backDiff(fx, 1.0f, 2);

        // soft threshold shrinkage
        float[][] dyy = new float[h][w];
        float threshold = 1.0f / mu;
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                float z = gyy[y][x] + byy[y][x];
                float absZ = Math.abs(z);
                if (absZ > threshold) {
                    dyy[y][x] = Math.signum(z) * (absZ - threshold);
                } else {
                    dyy[y][x] = 0.0f;
                }
            }
        }

        // Update the dual variables
        float[][] newByy = new float[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                newByy[y][x] = byy[y][x] + (gyy[y][x] - dyy[y][x]);
            }
        }

        // CalculateLyy
        float[][] diff = new float[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                diff[y][x] = dyy[y][x] - newByy[y][x];
            }
        }
        float[][] fx2 = DifferentialOps.forwardDiff(diff, 1.0f, 2);
        float[][] Lyy = DifferentialOps.backDiff(fx2, 1.0f, 2);
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                Lyy[y][x] *= para;
            }
        }

        return new IterResult(Lyy, newByy);
    }

    /**
     * Two-order derivative iteration in the xy direction (Hessian term), corresponding to MATLAB's iter_xy
     * @param g The current iteration image
     * @param bxy The dual variable from the previous iteration
     * @param para The weight of the regularization term, originally fixed at 2
     * @param mu The penalty parameter, originally set to 1
     * @return The iteration result Lxy and the updated bxy
     */
    public static IterResult iterXY(float[][] g, float[][] bxy, float para, float mu) {
        int h = g.length;
        int w = g[0].length;

        // Calculate the mixed second-order derivative gxy = forward_diff(forward_diff(g,1,1),1,2)
        float[][] fx1 = DifferentialOps.forwardDiff(g, 1.0f, 1);
        float[][] gxy = DifferentialOps.forwardDiff(fx1, 1.0f, 2);

        // soft threshold shrinkage
        float[][] dxy = new float[h][w];
        float threshold = 1.0f / mu;
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                float z = gxy[y][x] + bxy[y][x];
                float absZ = Math.abs(z);
                if (absZ > threshold) {
                    dxy[y][x] = Math.signum(z) * (absZ - threshold);
                } else {
                    dxy[y][x] = 0.0f;
                }
            }
        }

        // Update the dual variables
        float[][] newBxy = new float[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                newBxy[y][x] = bxy[y][x] + (gxy[y][x] - dxy[y][x]);
            }
        }

        // Calculate the Lxy = para * back_diff(back_diff(dxy - newBxy, 1, 2), 1, 1)
        float[][] diff = new float[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                diff[y][x] = dxy[y][x] - newBxy[y][x];
            }
        }
        float[][] bx1 = DifferentialOps.backDiff(diff, 1.0f, 2);
        float[][] Lxy = DifferentialOps.backDiff(bx1, 1.0f, 1);
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                Lxy[y][x] *= para;
            }
        }

        return new IterResult(Lxy, newBxy);
    }

    /**
     Sparse term iteration (L1 regularization), corresponding to MATLAB's iter_sparse
     * @param g The current iteration image
     * @param bsparse The dual variable from the previous iteration
     * @param para The weight of the sparse term paral1
     * @param boost Acceleration mode: 0 = off, 1 = on
     * @param mu Penalty parameter, default value in the original code is 1
     * @return Iteration result Lsparse and updated bsparse
     */
    public static IterResult iterSparse(float[][] g, float[][] bsparse, float para, int boost, float mu) {
        int h = g.length;
        int w = g[0].length;

        // dealing the boost mode
        float[][] gsparse;
        if (boost == 1) {
            gsparse = multiply(g, para);
        } else {
            gsparse = copy(g);
        }

        // soft threshold shrinkage
        float[][] d = new float[h][w];
        float threshold = 1.0f / mu;
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                float z = gsparse[y][x] + bsparse[y][x];
                float absZ = Math.abs(z);
                if (absZ > threshold) {
                    d[y][x] = Math.signum(z) * (absZ - threshold);
                } else {
                    d[y][x] = 0.0f;
                }
            }
        }

        // Update the dual variables
        float[][] newBsparse = new float[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                newBsparse[y][x] = bsparse[y][x] + (gsparse[y][x] - d[y][x]);
            }
        }

        // Calculate Lsparse = para * (d - newBsparse)
        float[][] Lsparse = new float[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                Lsparse[y][x] = para * (d[y][x] - newBsparse[y][x]);
            }
        }

        return new IterResult(Lsparse, newBsparse);
    }

    // ------------------------------ Overloaded methods (default parameters) ------------------------------
    public static IterResult iterXX(float[][] g, float[][] bxx) {
        return iterXX(g, bxx, 1.0f, 1.0f);
    }

    public static IterResult iterYY(float[][] g, float[][] byy) {
        return iterYY(g, byy, 1.0f, 1.0f);
    }

    public static IterResult iterXY(float[][] g, float[][] bxy) {
        return iterXY(g, bxy, 2.0f, 1.0f);
    }

    public static IterResult iterSparse(float[][] g, float[][] bsparse, float para) {
        return iterSparse(g, bsparse, para, 0, 1.0f);
    }

    // ------------------------------ Fundamental Matrix Tool Method ------------------------------
    /**
     * deep copy of matrix
     */
    private static float[][] copy(float[][] src) {
        int h = src.length;
        int w = src[0].length;
        float[][] dst = new float[h][w];
        for (int y = 0; y < h; y++) {
            System.arraycopy(src[y], 0, dst[y], 0, w);
        }
        return dst;
    }

    /**
     * Matrix element-wise multiplication by a scalar
     */
    private static float[][] multiply(float[][] src, float scalar) {
        int h = src.length;
        int w = src[0].length;
        float[][] dst = new float[h][w];
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                dst[y][x] = src[y][x] * scalar;
            }
        }
        return dst;
    }
}