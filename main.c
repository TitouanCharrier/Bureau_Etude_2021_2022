#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

int odefunc (double x, const double y[], double f[], void *params)
{
    f[0] = x+2*y[0];   
    return GSL_SUCCESS;
}


int * jac;

int main ()
{
    int dim = 1;
    gsl_odeiv2_system sys = {odefunc, NULL, dim, NULL};

    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
    int i;
    double x0 = 0.0,  xf = 100.0; /* start and end of integration interval */
    double x = x0;
    double y[1] = { 1.0  };  /* initial value */

    for (i = 1; i <= 100; i++)
    {
        double xi = x0 + i * (xf-x0) / 100.0;
        int status = gsl_odeiv2_driver_apply (d, &x, xi, y);

        if (status != GSL_SUCCESS)
        {
            printf ("error, return value=%d\n", status);
            break;
        }

        printf ("%.8e %.8e\n", x, y[0]);
    }

    gsl_odeiv2_driver_free (d);
    return 0;
}

