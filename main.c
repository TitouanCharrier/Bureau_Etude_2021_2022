#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <math.h>

#define hbar 1.05e-1
#define m 1e-3
#define pi 3.14
#define n 1

double V (double x) {
  if (x<=0.5) {
    return 0;
  } else {
    return 50;
  }
}

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

  int dimension = 3;

  double error_H = 1.e-8;
  double error_L = 1.e-10;

  const gsl_odeiv_step_type *type_ptr = gsl_odeiv_step_rkf45;

  gsl_odeiv_step *step_ptr = gsl_odeiv_step_alloc (type_ptr, dimension);
  gsl_odeiv_control *control_ptr = gsl_odeiv_control_y_new (error_H, error_L);
  gsl_odeiv_evolve *evolve_ptr = gsl_odeiv_evolve_alloc (dimension);

  gsl_odeiv_system my_system;

  double L = 1;
  double y[3];
  y[0] = 0.;
  y[1] = New_z;
  y[2] = 0.;

  //printf("New_E=%f \t New_z=%f\n", New_E, New_z );

  //ligne temporelle
  double t, t_next;
  double tmin, delta_t;
  tmin = 0.;
  delta_t = 0.01;


 //systemes
  my_system.function = ShrodingerFunction;
  my_system.jacobian = NULL;
  my_system.dimension = dimension;
  my_system.params = &New_E;

  t = tmin;

  double h = 6.626e-1;

  for (t_next = tmin + delta_t; t_next <= L; t_next += delta_t) {
      while (t < t_next) {
          gsl_odeiv_evolve_apply (evolve_ptr, control_ptr, step_ptr, &my_system, &t, t_next, &h, y);
      }
  }


  //gsl_odeiv_evolve_apply (evolve_ptr, control_ptr, step_ptr, &my_system, &t, L, &h, y);
  printf ("t=%.5e \t y[0]=%.5e \t y[2]=%.5e \t y[2]=%.5e [inside Solve]\n", t, y[0], y[1], y[2]);

  const double y0 = y[0];
  const double y1 = y[2]-1;

  gsl_vector_set (f_vector, 0, y0);
  gsl_vector_set (f_vector, 1, y1);

  gsl_odeiv_evolve_free (evolve_ptr);
  gsl_odeiv_control_free (control_ptr);
  gsl_odeiv_step_free (step_ptr);

  return GSL_SUCCESS;
}

int
print_state (size_t iter, gsl_multiroot_fsolver * s)
{
  printf ("iter = %3u \t E = %.3f \t z = %.3f \t "
          "f(x) = %.3e \t Normalisation = %.3e\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1));
}

int main () {

double L = 1;

double E = n*n*(hbar*hbar*pi*pi)/(2*m*L*L);
double z = sqrt(2/L)*sqrt(2*m*E/hbar);

int dimension = 3;

double error_H = 1.e-8;
double error_L = 1.e-10;


const gsl_odeiv_step_type *type_ptr = gsl_odeiv_step_rkf45;

gsl_odeiv_step *step_ptr = gsl_odeiv_step_alloc (type_ptr, dimension);
gsl_odeiv_control *control_ptr = gsl_odeiv_control_y_new (error_H, error_L);
gsl_odeiv_evolve *evolve_ptr = gsl_odeiv_evolve_alloc (dimension);

gsl_odeiv_system my_system;

//ligne temporelle
double t, t_next;
double tmin, delta_t;
tmin = 0.;
delta_t = 0.01;


/*----------------------------------------------------------------------------*/

const gsl_multiroot_fsolver_type *T;
gsl_multiroot_fsolver *s;

int status;
size_t i, iter = 0;

const size_t N = 2;

struct rparams p;

gsl_multiroot_function f = {&Solve, N, &p};

gsl_vector *x = gsl_vector_alloc (N);

gsl_vector_set (x, 0, E);
gsl_vector_set (x, 1, z);

T = gsl_multiroot_fsolver_hybrids;
s = gsl_multiroot_fsolver_alloc (T, 2);
gsl_multiroot_fsolver_set (s, &f, x);

print_state (iter, s);

do
  {
    iter++;
    gsl_multiroot_fsolver_iterate (s);

    print_state (iter, s);

    status = gsl_multiroot_test_residual (s->f, 1e-7);
  }
while (status == GSL_CONTINUE && iter < 50);

printf ("status = %s\n", gsl_strerror (status));


E = gsl_vector_get (s->x, 0);
z = gsl_vector_get (s->x, 1);

double y[3];
y[0] = 0.;
y[1] = z;
y[2] = 0.;

//systemes
my_system.function = ShrodingerFunction;
my_system.jacobian = NULL;
my_system.dimension = dimension;
my_system.params = &E;

t = tmin;

double h = 6.626e-1;

FILE *fpt; fpt = fopen("MyFile.csv", "w");
fprintf(fpt,"%lf, %lf,\n", t, y[0]);
printf ("%.5e %.5e\n", t, y[0]);

for (t_next = tmin + delta_t; t_next <= L; t_next += delta_t) {
    while (t < t_next) {
        gsl_odeiv_evolve_apply (evolve_ptr, control_ptr, step_ptr, &my_system, &t, t_next, &h, y);
    }
    fprintf(fpt,"%lf, %lf,\n", t, y[0]);
    printf ("t=%.5e \t y[0]=%.5e \t y[2]=%.5e \t y[2]=%.5e\n", t, y[0], y[1], y[2]);
}

fclose(fpt);


    //clean
    gsl_odeiv_evolve_free (evolve_ptr);
    gsl_odeiv_control_free (control_ptr);
    gsl_odeiv_step_free (step_ptr);
    gsl_multiroot_fsolver_free (s);
    gsl_vector_free (x);

    return 0;
}
