# include <stdlib.h>
# include <stdio.h>
# include <math.h>

int solve_quartic(double c[5], double s[4]);
int solve_quadric(double c[3], double s[2]);

#define EQN_EPS 1e-9
#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

#ifdef cbrt
#undef cbrt
#endif

double cbrt(double x)
{
    if (x > 0.0)
        return pow(x, 1.0/3.0);
    else
        return -pow(-x, 1.0/3.0);
}
static int is_zero(double x)
{
    return x > -EQN_EPS && x < EQN_EPS;
}
int solve_linear(double c[2], double s[1])
{
    if (is_zero(c[1]))
        return 0;
    s[0] = - c[0] / c[1];
    return 1;
}
int solve_quadric(double c[3], double s[2])
{
    double p, q, D;


// make sure we have a d2 equation

    if (is_zero(c[2]))
        return solve_linear(c, s);


// normal for: x^2 + px + q
    p = c[1] / (2.0 * c[2]);
    q = c[0] / c[2];
    D = p * p - q;

    if (is_zero(D)) {
        // one double root
        s[0] = s[1] = -p;
        return 1;
    }

    if (D < 0.0)
        // no real root
        return 0;

    else {
        // two real roots
        double sqrt_D = sqrt(D);
        s[0] = sqrt_D - p;
        s[1] = -sqrt_D - p;
        return 2;
    }
}
int solve_cubic(double c[4], double s[3])
{
    int	i, num;
    double	sub,
            A, B, C,
            sq_A, p, q,
            cb_p, D;

// normalize the equation:x ^ 3 + Ax ^ 2 + Bx  + C = 0
    A = c[2] / c[3];
    B = c[1] / c[3];
    C = c[0] / c[3];

// substitute x = y - A / 3 to eliminate the quadric term: x^3 + px + q = 0

    sq_A = A * A;
    p = 1.0/3.0 * (-1.0/3.0 * sq_A + B);
    q = 1.0/2.0 * (2.0/27.0 * A *sq_A - 1.0/3.0 * A * B + C);

// use Cardano's formula

    cb_p = p * p * p;
    D = q * q + cb_p;

    if (is_zero(D)) {
        if (is_zero(q)) {
            // one triple solution
            s[0] = 0.0;
            num = 1;
        } else {
            // one single and one double solution
            double u = cbrt(-q);
            s[0] = 2.0 * u;
            s[1] = - u;
            num = 2;
        }
    } else if (D < 0.0) {
        // casus irreductibilis: three real solutions
        double phi = 1.0/3.0 * acos(-q / sqrt(-cb_p));
        double t = 2.0 * sqrt(-p);
        s[0] = t * cos(phi);
        s[1] = -t * cos(phi + M_PI / 3.0);
        s[2] = -t * cos(phi - M_PI / 3.0);
        num = 3;
    } else {
        // one real solution
        double sqrt_D = sqrt(D);
        double u = cbrt(sqrt_D + fabs(q));
        if (q > 0.0)
            s[0] = - u + p / u ;
        else
            s[0] = u - p / u;
        num = 1;
    }

// resubstitute
    sub = 1.0 / 3.0 * A;
    for (i = 0; i < num; i++)
        s[i] -= sub;
    return num;
}
int solve_quartic(double c[5], double s[4])
{
    double	    coeffs[4],
                z, u, v, sub,
                A, B, C, D,
                sq_A, p, q, r;
    int	    i, num;


// normalize the equation:x ^ 4 + Ax ^ 3 + Bx ^ 2 + Cx + D = 0

    A = c[3] / c[4];
    B = c[2] / c[4];
    C = c[1] / c[4];
    D = c[0] / c[4];

// subsitute x = y - A / 4 to eliminate the cubic term: x^4 + px^2 + qx + r = 0

    sq_A = A * A;
    p = -3.0 / 8.0 * sq_A + B;
    q = 1.0 / 8.0 * sq_A * A - 1.0 / 2.0 * A * B + C;
    r = -3.0 / 256.0 * sq_A * sq_A + 1.0 / 16.0 * sq_A * B - 1.0 / 4.0 * A * C + D;

    if (is_zero(r)) {
        // no absolute term:y(y ^ 3 + py + q) = 0
        coeffs[0] = q;
        coeffs[1] = p;
        coeffs[2] = 0.0;
        coeffs[3] = 1.0;

        num = solve_cubic(coeffs, s);
        s[num++] = 0;
    } else {
        // solve the resolvent cubic...
        coeffs[0] = 1.0 / 2.0 * r * p - 1.0 / 8.0 * q * q;
        coeffs[1] = -r;
        coeffs[2] = -1.0 / 2.0 * p;
        coeffs[3] = 1.0;
        (void) solve_cubic(coeffs, s);

        // ...and take the one real solution...
        z = s[0];

        // ...to build two quadratic equations
        u = z * z - r;
        v = 2.0 * z - p;

        if (is_zero(u))
            u = 0.0;
        else if (u > 0.0)
            u = sqrt(u);
        else
            return 0;

        if (is_zero(v))
            v = 0;
        else if (v > 0.0)
            v = sqrt(v);
        else
            return 0;

        coeffs[0] = z - u;
        coeffs[1] = q < 0 ? -v : v;
        coeffs[2] = 1.0;

        num = solve_quadric(coeffs, s);

        coeffs[0] = z + u;
        coeffs[1] = q < 0 ? v : -v;
        coeffs[2] = 1.0;

        num += solve_quadric(coeffs, s + num);
    }

// resubstitute
    sub = 1.0 / 4 * A;
    for (i = 0; i < num; i++)
        s[i] -= sub;

    return num;

}
