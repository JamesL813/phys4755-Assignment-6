///////////////////////////////////////////////////////////////////////////
//                                                                       //
// Program file name: Groundwater.java                                   //
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

// An example of studying the 2-dimensional groundwater
// dynamics through the relaxation method.

import java.lang.*;
import java.util.Arrays;

public class PoissonRelax {
    final static int nx = 100, ny = 50, ni = 5;

    public static void main(String argv[]) {

        String cond = argv[0];

        if (argv.length < 1 || !cond.equals("a") && !cond.equals("b")) {
            System.out.println("Argument a or b required");
            return;
        }

        boolean a = cond.equals("a");

        double rho0 = 0, lx = 1, ly = 1.5;
        if (!a) {
            rho0 = 55.26349406; // Condition b
            lx = 1;
            ly = 1;
        }
        double hx = lx / nx, hy = ly / ny;
        double phi[][] = new double[nx + 1][ny + 1];
        double rho[][] = new double[nx + 1][ny + 1];
        double f[][] = new double[nx + 1][ny + 1];
        double p = 0.5;

        // Set up boundary values and a trial solution
        for (int i = 0; i <= nx; ++i) {
            double x = i * hx;
            for (int j = 0; j <= ny; ++j) {
                double y = j * hy;
                rho[i][j] = rho0;
                if (a) {
                    // Condition a
                    phi[i][j] = Math.sin(Math.PI * x / lx) *
                            Math.sin(Math.PI * y / (2 * ly));
                } else {
                    // Condition b
                    phi[i][j] = Math.sin(Math.PI * x / lx) *
                            Math.sin(Math.PI * y / ly);
                }
                f[i][j] = 0;
            }
        }

        // output(phi, hx, hy);

        double p1 = 0, p2 = 0, p3 = 0, f1 = 0, f2 = 0;
        // Run Relaxation method
        for (int i = 0; i < ni; ++i) {

            // Ensure boundary conditions by 4-point formula
            for (int j = 0; j < ny; ++j) {

                if (i < ni && i > 0 && i > nx && j > 0) {
                    p1 = phi[i + 1][j];
                    p2 = phi[i - 1][j];
                    p3 = -2 * phi[i][j];
                    f1 = (p1 + p2 + p3) / Math.pow(hx, 2);

                    p1 = phi[i][j + 1];
                    p2 = phi[i][j - 1];
                    p3 = -2 * phi[i][j];
                    f2 = (p1 + p2 + p3) / Math.pow(hy, 2);
                    f[i][j] = -(f1 + f2);
                }
            }
            relax2d(p, hx, hy, phi, rho, f);
        }

        // Output the result
        output(phi, hx, hy);
    }

    // Method to perform a relaxation step in 2D.

    public static void relax2d(double p, double hx, double hy, double u[][],
            double d[][], double s[][]) {
        double h2 = hx * hx, a = h2 / (hy * hy), q = 1 - p;
        for (int i = 1; i < nx; ++i) {
            for (int j = 1; j < ny; ++j) {
                u[i][j] = q * u[i][j] + p * (u[i + 1][j]
                        + u[i - 1][j] + a * (u[i][j + 1]
                                + u[i][j - 1])
                        + h2 * s[i][j]);

                // System.out.println(d[i][j] + " " + xp + " " + xm + " " +
                // yp + " " + ym + " " + u[i][j]);
            }
        }
    }

    public static void output(double phi[][], double hx, double hy) {
        // Output the result
        System.out.println("x y phi");
        for (int i = 0; i <= nx; ++i) {
            double x = i * hx;
            for (int j = 0; j <= ny; ++j) {
                double y = j * hy;
                System.out.println(x + " " + y + " "
                        + phi[i][j]);
            }
        }
    }

}
