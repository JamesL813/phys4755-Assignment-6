///////////////////////////////////////////////////////////////////////////
//                                                                       //
// Program file name: Carlo.java                                         //
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

// An example of Monte Carlo simulation with the
// Metropolis scheme with integrand f(x) = x*x.

import java.lang.*;
import java.util.Random;

public class Carlo {
    static final int nsize = 10000;
    static final int nskip = 15;
    static final int ntotal = nsize * nskip;
    static final int neq = 10000;
    static int iaccept = 0;
    static double x, w, h = 0.4, z = 0.46265167;
    static Random r = new Random();

    public static void main(String argv[]) {
  
        System.out.println("S error acc steps");
        for (int steps = 10000; steps > 0; steps -= 100) {

            x = r.nextDouble();
            w = weight();
            for (int i = 0; i < neq; ++i)
                metropolis();

            double s0 = 0;
            double ds = 0;
            iaccept = 0;
            for (int i = 0; i < ntotal; ++i) {
                metropolis();
                if (i % nskip == 0) {
                    double f = g(x);
                    s0 += f;
                    ds += f * f;
                }
            }
            s0 /= steps;
            ds /= steps;
            ds = Math.sqrt(Math.abs(ds - s0 * s0) / steps);
            s0 *= z;
            ds *= z;
            double accept = 100.0 * iaccept / ntotal;
            System.out.println(s0 + " " + ds + " " + accept + " " + steps);
        }
    }

    public static void metropolis() {
        double xold = x;
        x = x + 2 * h * (r.nextDouble() - 0.5);
        if ((x < 0) || (x > 1))
            x = xold;
        else {
            double wnew = weight();
            if (wnew > w * r.nextDouble()) {
                w = wnew;
                ++iaccept;
            } else
                x = xold;
        }
    }

    public static double weight() {
        return Math.exp(x * x) - 1;
    }

    public static double g(double y) {
        return (Math.exp(-(y * y) / 2));
    }
}
