///////////////////////////////////////////////////////////////////////////
//                                                                       //
// Program file name: Nuclear.java                                       //
//                                                                       //
// Tao Pang 2006                                                         //
//                                                                       //
// Last modified: January 18, 2006                                       //
//                                                                       //
// (1) This Java program is part of the book, "An Introduction to        //
//     Computational Physics, 2nd Edition," written by Tao Pang and      //
//     published by Cambridge University Press on January 19, 2006.      //
//                                                                       //
// (2) No warranties, express or implied, are made for this program.     //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

// A program to study the time-dependent temperature
// field around a nuclear waste rod in a 2D model.

import java.lang.*;

public class Nuclear {
    final static int nx = 100, n = 5, nt = 1000, mt = 1000;

    public static void main(String argv[]) {
        double d[] = new double[nx];
        double e[] = new double[nx];
        double c[] = new double[nx];
        double b[] = new double[nx];
        double p[] = new double[nx];
        double s[][] = new double[nt + 1][nx];
        double T[][] = new double[nt + 1][nx + 1];
        double dt = 1.0 / mt, tc = 1, T0 = 1, kappa = 2e7;
        double len = 100, low = 35, upp = 65;
        double ra = 35, rb = 100, h = rb / nx, h2 = h * h;
        double s0 = dt * kappa * T0 / (ra * ra), g = dt * kappa / h2;

        // Assign the elements in the matrix 2-H_i
        for (int i = 0; i < nx; ++i) {
            d[i] = 2 * (1 + g);
            e[i] = -(1 + 0.5 / (i + 1)) * g;
            c[i] = -(1 - 0.5 / (i + 2)) * g;
        }

        // Modify the first equation from T"=0 at r=0
        d[0] -= 2 * g / 3;
        e[0] += g / 6;

        // Assign the source of the radiation heat
        int na = (int) (len / h);
        int nl = (int) (low / h);
        int nu = (int) (upp / h);
        for (int i = 0; i <= nt; ++i) {
            double t = -dt * i / tc;
            for (int j = nl; j < nu - 1; ++j) {
                s[i][j] = (T0 / Math.pow(ra,2)) * Math.exp(-t/100);
            }
        }

        // Find the temperature field recursively
        for (int i = 1; i <= nt; ++i) {

            // Assign the elements in the matrix 2+H_0
            double d0 = 2 * (1 - g);
            double e0 = (1 + 0.5) * g;
            double c0 = (1 - 0.5) * g;

            // Evaluate b[0] under the condition T"=0 at r=0
            int center = (int) (50 / h); 
            b[0] = d0 * T[i - 1][0] + e0 * T[i - 1][1]
                    + c0 * (4 * T[i - 1][0] - T[i - 1][1]) / 3
                    + s[i - 1][center] + s[i][center];

            // Find the elements in the array b[i]
            for (int j = 1; j < nx; ++j) {

                // Assign the elements in the matrix 2+H_0
                d0 = 2 * (1 - g);
                e0 = (1 + 0.5 / (j + 1)) * g;
                c0 = (1 - 0.5 / (j + 1)) * g;

                // Obtain the elements from the last recursion
                b[j] = d0 * T[i - 1][j] + e0 * T[i - 1][j + 1]
                        + c0 * T[i - 1][j - 1] + s[i - 1][j] + s[i][j];
            }

            // Obtain the solution of the temperature field
            p = tridiagonalLinearEq(d, e, c, b);
            for (int j = 0; j < nx; ++j)
                T[i][j] = p[j];
        }

        // Output the result at every n spatial data points
        System.out.println("x t");
        for (int j = 0; j < nx; j += n) {
            double x = h * (j + 1);
            System.out.println(x + " " + T[nt][j]);
        }
    }

    // Method to solve the tridiagonal linear equation set.

    public static double[] tridiagonalLinearEq(double d[],
            double e[], double c[], double b[]) {
        int m = b.length;
        double w[] = new double[m];
        double y[] = new double[m];
        double z[] = new double[m];
        double v[] = new double[m - 1];
        double t[] = new double[m - 1];

        // Evaluate the elements in the LU decomposition
        w[0] = d[0];
        v[0] = c[0];
        t[0] = e[0] / w[0];
        for (int i = 1; i < m - 1; ++i) {
            w[i] = d[i] - v[i - 1] * t[i - 1];
            v[i] = c[i];
            t[i] = e[i] / w[i];
        }
        w[m - 1] = d[m - 1] - v[m - 2] * t[m - 2];

        // Forward substitution to obtain y
        y[0] = b[0] / w[0];
        for (int i = 1; i < m; ++i)
            y[i] = (b[i] - v[i - 1] * y[i - 1]) / w[i];

        // Backward substitution to obtain z
        z[m - 1] = y[m - 1];
        for (int i = m - 2; i >= 0; --i) {
            z[i] = y[i] - t[i] * z[i + 1];
        }
        return z;
    }
}
