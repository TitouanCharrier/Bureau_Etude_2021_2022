#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <math.h>

#include "Fonction_de_potentiel.h"

#define hbar 1.05e-1 // yg*nm*nm/ps
#define m 1e-3  // yg
#define pi 3.14

extern int n;
extern double L; // nm

int ShrodingerFunction (double x, const double y[], double dy[], void *params_ptr) {
    double E = *(double *) params_ptr;

    dy[0] = y[1];
    dy[1] = -2*m/(hbar*hbar)*(E-V(x))*y[0];
    dy[2] = y[0]*y[0];

    return GSL_SUCCESS;
}

struct rparams {};

int Solve (const gsl_vector * x, void *params, gsl_vector * f_vector) {

  double New_E = gsl_vector_get (x, 0);
  double New_z = gsl_vector_get (x, 1);

  double error_H = 1.e-8;
  double error_L = 1.e-10;

  gsl_odeiv2_system sys = {ShrodingerFunction, NULL, 3, &New_E};

  gsl_odeiv2_driver * d = 
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45,error_H, error_L, 0.0);

  double y[3];
  y[0] = 0.;
  y[1] = New_z;
  y[2] = 0.;

  //printf("New_E=%f \t New_z=%f\n", New_E, New_z );

  //ligne temporelle
  double t;
  double tmin;
  tmin = 0.;

  t = tmin;

      gsl_odeiv2_driver_apply (d, &t, L, y);

  //printf ("t=%.5e \t y[0]=%.5e \t y[2]=%.5e \t y[2]=%.5e [inside Solve]\n", t, y[0], y[1], y[2]);

  const double y0 = y[0];
  const double y1 = y[2]-1;

  gsl_vector_set (f_vector, 0, y0);
  gsl_vector_set (f_vector, 1, y1);

  gsl_odeiv2_driver_free (d);

  return GSL_SUCCESS;
}

int print_state (size_t iter, gsl_multiroot_fsolver * s) {
  printf ("iter = %3lu \t E = %.3f \t z = %.3f \t "
          "f(x) = %.3e \t Normalisation = %.3e\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1));
  return 0;
}

int main () {

  double E = n*n*(hbar*hbar*pi*pi)/(2*m*L*L); // yg*nm*nm/(ps*ps)
  double z = sqrt(2/L)*sqrt(2*m*E/hbar); // yg^1/2 *nm^1/2 /ps^1/2


  double error_H = 1.e-8;
  double error_L = 1.e-10;

  gsl_odeiv2_system sys = {ShrodingerFunction, NULL, 3, &E};

  gsl_odeiv2_driver * d = 
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45,error_H, error_L, 0.0);

  //ligne temporelle
  double t, t_next;
  double tmin, delta_t;
  tmin = 0.;
  delta_t = 0.01;


  /*----------------------------------------------------------------------------*/

  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;

  int status;
  size_t iter = 0;

  const size_t N = 2;

  struct rparams p;

  gsl_multiroot_function f = {&Solve, N, &p};

  gsl_vector *x = gsl_vector_alloc (N);

  gsl_vector_set (x, 0, E);
  gsl_vector_set (x, 1, z);

  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc (T, 2);
  gsl_multiroot_fsolver_set (s, &f, x);

  //print_state (iter, s);

  do
    {
      iter++;
      gsl_multiroot_fsolver_iterate (s);

      //print_state (iter, s);

      status = gsl_multiroot_test_residual (s->f, 1e-7);
    }
  while (status == GSL_CONTINUE && iter < 50);

  //printf ("status = %s\n", gsl_strerror (status));


  E = gsl_vector_get (s->x, 0);
  z = gsl_vector_get (s->x, 1);

  double E_ev=E/(1.602*100);  // transformation en Ev avec 1ev = 1.602.e-19 J
                              // Z en SI : Z = Z*10^(-12) kg^1/2 * m^1/2 * s^(-1/2)
                              // L en SI = L*10^(-9) m 

  printf("\nResultats :\n\n pour n=%d et L=%f nm : \n\n E = %f Ev, \t z = %f.10^(-12) unités SI \n E = %f J \n", n, L, E_ev, z, E);

  double y[3];
  y[0] = 0.;
  y[1] = z;
  y[2] = 0.;

  t = tmin;

  FILE *fpt; fpt = fopen("MyFile.csv", "w");
  fprintf(fpt,"%lf, %lf,\n", t, y[0]);
  //printf ("%.5e %.5e\n", t, y[0]);

  for (t_next = tmin + delta_t; t_next <= L; t_next += delta_t) {
      while (t < t_next) {
          gsl_odeiv2_driver_apply (d, &t, t_next, y);
      }
      fprintf(fpt,"%lf, %lf,\n", t, y[0]);
      //printf ("t=%.5e \t y[0]=%.5e \t y[2]=%.5e \t y[2]=%.5e\n", t, y[0], y[1], y[2]);
  }

  fclose(fpt);

  //clean
  gsl_odeiv2_driver_free (d);
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);

    return 0;
}
