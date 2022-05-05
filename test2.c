#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <math.h>

#define hbar 1.05e-10
#define m 1e-3
#define V 0
#define pi 3.14
#define L tmax

int f (double x, const double y[], double dy[], void *params_ptr) {
    double E = *(double *) params_ptr;

    dy[0] = y[1];
    dy[1] = -2*m/(hbar*hbar)*(E-V)*y[0];
    dy[2] = y[0]*y[0];

    return GSL_SUCCESS; 
}

int main () {

    int dimension = 3;	
  
    double error_H = 1.e-8;	
    double error_L = 1.e-10;

    const gsl_odeiv_step_type *type_ptr = gsl_odeiv_step_rkf45;

    gsl_odeiv_step *step_ptr = gsl_odeiv_step_alloc (type_ptr, dimension);
    gsl_odeiv_control *control_ptr = gsl_odeiv_control_y_new (error_H, error_L);
    gsl_odeiv_evolve *evolve_ptr = gsl_odeiv_evolve_alloc (dimension);

    gsl_odeiv_system my_system;				

    
    int n = 3;

    //ligne temporelle
    double t, t_next;       
    double tmin, tmax, delta_t;
    tmin = 0.;      
    tmax = 1e-19;    
    delta_t = 1e-22;

    double E = n*n*(hbar*hbar*pi*pi)/(2*m*L*L); 
    double z = sqrt(2/L)*sqrt(2*m*E/hbar);    
    double y[2];
    y[0] = 0.;  
    y[1] = z;  

   //systemes
    my_system.function = f;	
    my_system.jacobian = NULL;	
    my_system.dimension = dimension;	
    my_system.params = &E;	

    t = tmin;   

    double h = 1e-15;

    FILE *fpt; fpt = fopen("MyFile.csv", "w");
    fprintf(fpt,"%lf, %lf,\n", t, y[0]);
    printf ("%.5e %.5e\n", t, y[0]);

    for (t_next = tmin + delta_t; t_next <= tmax; t_next += delta_t) {
        while (t < t_next) {
            gsl_odeiv_evolve_apply (evolve_ptr, control_ptr, step_ptr, &my_system, &t, t_next, &h, y);
        }
        fprintf(fpt,"%.5e, %.5e,\n", t, y[0]);
        printf ("%.5e %.5e %.5e %.5e\n", t, y[0], y[1], y[2]);
    }

    fclose(fpt);   

    //clean
    gsl_odeiv_evolve_free (evolve_ptr);
    gsl_odeiv_control_free (control_ptr);
    gsl_odeiv_step_free (step_ptr);

    return 0;
}
