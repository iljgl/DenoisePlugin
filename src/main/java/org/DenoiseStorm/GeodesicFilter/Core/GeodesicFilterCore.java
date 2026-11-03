package org.DenoiseStorm.GeodesicFilter.Core;

import java.util.PriorityQueue;
import java.util.Arrays;
import java.util.Comparator;

public class GeodesicFilterCore {

    // 8 Neighborhood Offset (Row, Column)
    private static final int[][] OFFSETS = {
            {-1, -1}, {-1, 0}, {-1, 1},
            {0, -1},           {0, 1},
            {1, -1},  {1, 0},  {1, 1}
    };
    // The spatial Euclidean distance corresponding to 8 directions
    private static final double[] SPATIAL_DISTS = {
            Math.sqrt(2), 1, Math.sqrt(2),
            1,                1,
            Math.sqrt(2), 1, Math.sqrt(2)
    };

    /**
     * Geodesic Filter Main Entry Point
     * @param img     Input grayscale image, value range [0, 1], Java row-first double[h][w]
     * @param W       Window size (odd number)
     * @param alpha   Gray difference weight
     * @param sigma   Gaussian kernel standard deviation
     * @return        Filtered image, value range [0, 1]
     */
    public static double[][] process(double[][] img, int W, double alpha, double sigma) {
        if (W % 2 == 0) {
            throw new IllegalArgumentException("The window size W must be an odd number.");
        }
        if (sigma <= 0 || alpha <= 0) {
            throw new IllegalArgumentException("Alpha and sigma must be greater than 0.");
        }

        int h = img.length;
        int w = img[0].length;
        int half = W / 2;

        // Symmetric Padding
        double[][] padded = padSymmetric(img, half);
        double[][] filtered = new double[h][w];

        // The row and column coordinates (row-first index -> (row, col)) of each node within the pre-computed window
        // Java row-major：idx = row * W + col
        int[] nodeRows = new int[W * W];
        int[] nodeCols = new int[W * W];
        for (int r = 0; r < W; r++) {
            for (int c = 0; c < W; c++) {
                int idx = r * W + c;
                nodeRows[idx] = r;
                nodeCols[idx] = c;
            }
        }

        //Central node index
        int centerIdx = half * W + half;

        //Traverse each pixel of the original image
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                //Extract window (extracting W x W from the padded image)
                double[][] window = extractWindow(padded, i + half, j + half, half, W);

                // --- 1.Dijkstra calculates the geodesic distances from the center to all nodes.  ---
                double[] dist = dijkstra(window, W, centerIdx, nodeRows, nodeCols, alpha);

                // --- 2.Gaussian weighting + weighted average---
                double sigmaSq2 = 2 * sigma * sigma;
                double sumWeights = 0;
                double sumWeighted = 0;

                for (int idx = 0; idx < W * W; idx++) {
                    double weight = Math.exp(-dist[idx] * dist[idx] / sigmaSq2);
                    sumWeights += weight;
                    int r = nodeRows[idx];
                    int c = nodeCols[idx];
                    sumWeighted += weight * window[r][c];
                }

                filtered[i][j] = sumWeighted / sumWeights;
            }
        }

        return filtered;
    }

    /**
     * Dijkstra's shortest path: Calculate the geodesic distances from the center to all nodes within the calculation window
     */
    private static double[] dijkstra(double[][] window, int W, int centerIdx,
                                     int[] nodeRows, int[] nodeCols, double alpha) {
        int Nnodes = W * W;
        double[] dist = new double[Nnodes];
        boolean[] visited = new boolean[Nnodes];


        //noinspection MismatchedReadAndWriteOfArray
        @SuppressWarnings("unused")//Retain the path predecessors For future expansion convenience.
        int[] prev = new int[Nnodes]; //Retain the path predecessors and be consistent with the Matlab version.

        Arrays.fill(dist, Double.POSITIVE_INFINITY);
        Arrays.fill(prev, -1);
        dist[centerIdx] = 0;

        // Priority Queue: Sorted by distance in ascending order, elements are (distance, nodeIndex)
        // PriorityQueue:java(Minimum Binary Heap,O(nlogn)) better than matlab(linear list,O(n^2))
        PriorityQueue<double[]> pq = new PriorityQueue<>(Comparator.comparingDouble(a -> a[0]));
        pq.add(new double[]{0, centerIdx});

        while (!pq.isEmpty()) {
            double[] curr = pq.poll();//Inert deletion,Allow duplicate entries,O(nlogn)
            int u = (int) curr[1];
            double d = curr[0];

            if (visited[u]) continue;
            visited[u] = true;

            int ur = nodeRows[u];
            int uc = nodeCols[u];
            double valU = window[ur][uc];

            //Traverse the 8-neighborhood
            for (int k = 0; k < 8; k++) {
                int vr = ur + OFFSETS[k][0];
                int vc = uc + OFFSETS[k][1];

                //Window boundary check
                if (vr < 0 || vr >= W || vc < 0 || vc >= W) continue;

                int v = vr * W + vc;//int v = vr * W + (vc+1)-1;
                if (visited[v]) continue;

                //Edge weight: sqrt((spatial distance)² + (alpha * gray difference)²)
                double intensityDiff = valU - window[vr][vc];
                double edgeWeight = Math.sqrt(
                        SPATIAL_DISTS[k] * SPATIAL_DISTS[k] + alpha * alpha * intensityDiff * intensityDiff
                );

                double newDist = d + edgeWeight;
                if (newDist < dist[v]) {
                    dist[v] = newDist;
                    prev[v] = u;
                    pq.add(new double[]{newDist, v});//Inert deletion,Allow duplicate entries,O(nlogn)
                }
            }
        }

        return dist;
    }

    /**
     * Mirror-filling boundary (equivalent implementation in Matlab using padarray 'symmetric'))
     * Symmetrical mirroring: Boundary pixels do not repeat
     */
    private static double[][] padSymmetric(double[][] img, int pad) {
        int h = img.length;
        int w = img[0].length;
        int newH = h + 2 * pad;
        int newW = w + 2 * pad;
        double[][] padded = new double[newH][newW];

        //Fill the central area
        for (int i = 0; i < h; i++) {
            System.arraycopy(img[i], 0, padded[i + pad], pad, w);
        }

        //Fill the top and bottom (mirror flip)
        for (int i = 0; i < pad; i++) {
            //top: pad-1-i \to i
            System.arraycopy(padded[pad + pad - 1 - i], 0, padded[i], 0, newW);
            //bottom
            System.arraycopy(padded[h + pad - 1 - i], 0, padded[h + pad + i], 0, newW);
        }

        //Fill the left and right
        for (int i = 0; i < newH; i++) {
            for (int j = 0; j < pad; j++) {
                padded[i][j] = padded[i][pad + pad - 1 - j];
                padded[i][newW - 1 - j] = padded[i][w + j];
            }
        }

        return padded;
    }

    /**
     *Extract a W×W window centered at (cr, cc) and with a radius of half from the filled image.
     */
    private static double[][] extractWindow(double[][] padded, int cr, int cc, int half, int W) {
        double[][] window = new double[W][W];
        for (int r = 0; r < W; r++) {
            System.arraycopy(padded[cr - half + r], cc - half, window[r], 0, W);
        }
        return window;
    }
}
