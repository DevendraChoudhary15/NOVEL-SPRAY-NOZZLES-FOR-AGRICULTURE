/**
# Inputfile to test the first Version of Lagrange solver together with DNS.
# We initialize a liquid jet (f0 Variable) and a separate droplet (f1) as a VOF Variable,
# as well as severall particles for the Lagrange Solver.
# The droplet f1 will be removed from DNS in the first time step.

We solve first the two-phase Navier--Stokes equations with surface
tension and a momentum-conserving advection of velocity,
than the differential equations for the particles in Lagrange formalism.*/

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "tag.h"
#include "view.h"
#include "save_data.h"

/**
We define the radius of the jet, the initial jet length, the surface tension coefficient,
the jet velocity and Crossflow velocity rougly corresponding to the LPP setup. */

#define radius 0.000225 // 0.225 mm
#define length 0.1e-3 
#define SIGMA 3e-2
#define V_JET 20.7
#define V_CROS 100.0 // 100 m/s

/**
The maximum level of refinement is 10 and the error threshold
on velocity is 10. */

int maxlevel = 10;
double uemax = 10.0;

/**
We set the boundary conditions*/

scalar f0[];
// Jet injection
u.n[back]  = dirichlet(f0[]*V_JET);
u.t[back]  = dirichlet(0);
#if dimension > 2
u.r[back]  = dirichlet(0);
#endif
p[back]    = neumann(0);
f[back]    = f0[];

// Front: free slip wall, symetry boundary (default = nothing to set)
//u.n[front] = neumann(0);
//p[front]   = dirichlet(0);

// Right: open boundary conditions
u.n[right]   = neumann(0);
p[right]     = dirichlet(5.8); // on the right side, impose ambiant air pressure

// From the left comes the crossflow
u.n[left] = dirichlet(V_CROS*(1.0 - (1.0 - z/0.005)*(1.0 - z/0.005)));
u.t[left] = dirichlet(0);
#if dimension >2
u.r[left] = dirichlet(0);
#endif
p[left] = neumann(0); 

/**
The program can take two optional command-line arguments: the maximum
level and the error threshold on velocity. */

int main (int argc, char * argv[])
{
  if (argc > 1)
    maxlevel = atoi (argv[1]);
  if (argc > 2)
    uemax = atof (argv[2]);

  /**
  The initial domain is discretised with $32^3$ grid points. We set
  the origin and domain size. */
  
  init_grid (64);
  origin (-0.001, -0.005, 0);
  size (0.01);

  /**
  We set the density and viscosity of each phase as well as the
  surface tension coefficient and start the simulation. */
  
  rho1 = 997.00, 
  rho2 = 7.19;
  mu1 = 8.94e-4 , mu2 = 1.86e-5;
  f.sigma = SIGMA;

  run();
}

/** Prepare output */
event init_output (t = 0) {

  struct stat st_mov = {0};

  /** If it does not exist, create output folder for movies 
  with permissions for all*/
  if (stat("movies", &st_mov) == -1) {
      mkdir("movies", S_IRWXU | S_IRWXG | S_IRWXO);
  }
}


/**
## Initial conditions */

event init (t = 0) {

  if (!restore (file = "restart")) {

    /**
    We define the VOF Variable for the jet and the separate droplet.
    We use a static refinement down to *maxlevel* in a cylinder 1.2
    times longer than the initial jet and twice the radius. */
    
    refine (z < 1.2*length && sq(x) + sq(y) < 2.*sq(radius) && level < maxlevel);
     
    /** initialize the embedded boundary */

    vertex scalar phi[];
    foreach_vertex() {
      // hole 
      //phi[] = union ( z - 0.5*length , sq(radius) - sq(x) - sq(y));
      // button 
      //phi[] = union ( z - 0.5*length , sq(x) + sq(y) - sq(radius));
      // Wall
      phi[] = 0.0075 - z;
    }
    boundary ({phi});
    fractions (phi, cs, fs);

    /**
    We initialise the auxilliary volume fraction field for a cylinder of
    constant radius. */
    
    fraction (f0, sq(radius) - sq(x) - sq(y));
   
    foreach() {
      f[] = f0[]*(z < length);
    }
    f.refine = f.prolongation = fraction_refine;
    restriction ({f}); // for boundary conditions on levels

    foreach() {
      u.z[] = f[]*V_JET; 
      u.x[] = (1-f[])*V_CROS*(1.0 - (1.0 - z/0.005)*(1.0 - z/0.005));
    }
    boundary ({f,u.z});
  }
  else {
    fprintf(ferr,"restart loaded!!\n");    
  }
}

/**
## Outputs

We log some statistics on the solver. */

event logfile (i++) {
  //if (i == 0)
    fprintf (ferr,
	     "t dt mgp.i mgpf.i mgu.i grid->tn perf.t perf.speed\n");
  fprintf (ferr, "%g %g %d %d %d %ld %g %g\n", 
	   t, dt, mgp.i, mgpf.i, mgu.i,
	   grid->tn, perf.t, perf.speed);
}


/**
We save snapshots of the simulation at regular intervals to
restart or to post-process with [bview](/src/bview). */

event snapshot (t += 1.e-05; t <= 2.e-04) {

  char name[80];
  sprintf (name, "snapshot-%g", t);
  scalar pid[];
  foreach()
    pid[] = fmod(pid()*(npe() + 37), npe());
  boundary ({pid});
  dump (name);
}

/** video output */
event video_output(t += 1.e-05; t <= 2.e-04)
{
  /** Create some videos which can be animated */
  timer timer_movies;
  timer_movies = timer_start();

  view(width = 1920, height = 1080, camera = "front");
  view(fov = 11.0, tx = -0.35, ty = -0.125);
  clear();
  box();
  draw_vof("f");
  draw_vof("f", filled = 1, fc = {0.7, 0.7, 0.7});

  save("movies/movz.mp4");

  view(width = 1920, height = 1080, camera = "left");
  view(fov = 11.0, tx = -0.35, ty = -0.125);
  clear();
  box();
  //draw_vof("f");
  draw_vof("f", filled = 1, fc = {0.7, 0.7, 0.7});

  save("movies/movx.mp4");

  view(width = 1920, height = 1080, camera = "bottom");
  view(fov = 11.0, tx = -0.35, ty = -0.125);
  clear();
  box();
  draw_vof("f");
  draw_vof("f", filled = 1, fc = {0.7, 0.7, 0.7});

  save("movies/movy.mp4");

  double t_elapsed = timer_elapsed(timer_movies);
  fprintf(ferr, "\nOutput VOF + Mesh video!\n");
  fprintf(ferr, "  - time needed: %g s \n", t_elapsed);
}


/**
## Counting droplets

The number and sizes of droplets generated by the atomising jet is a
useful statistics for atomisation problems. This is not a quantity
which is trivial to compute. The *tag()* function is designed to solve
this problem. Any connected region for which *f[] > 1e-3* (i.e. a
droplet) will be identified by a unique "tag" value between 0 and
*n-1*. */

event droplets (i++)
{
  scalar m[];
  foreach(){
    m[] = f[] > 1e-3;
  }
  int n = tag (m);

  /**
  Once each cell is tagged with a unique droplet index, we can easily
  compute the volume *v* and position *b* of each droplet. Note that
  we use *foreach_leaf()* rather than *foreach()* to avoid doing a
  parallel traversal when using OpenMP. This is because we don't have
  reduction operations for the *v* and *b* arrays (yet). */

  double v[n],v1[n];
  coord b[n],k[n],k1[n];
  for (int j = 0; j < n; j++){
    v[j] = v1[j] = b[j].x = b[j].y = b[j].z = k[j].x = k[j].y = k[j].z = k1[j].x = k1[j].y = k1[j].z = 0.;
  }
  foreach_leaf(){
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*f[];
      v1[j] += dv()*(f[] > 1.-1e-3);
      coord p = {x,y,z};
      foreach_dimension(){
	      b[j].x += dv()*f[]*p.x;
	      k[j].x += dv()*f[]*u.x[];
	      k1[j].x += dv()*(f[] > 1.-1e-3)*u.x[];
      }        
    }
  }


  /**
  When using MPI we need to perform a global reduction to get the
  volumes and positions of droplets which span multiple processes. */

#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, v1, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, k, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, k1, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  /**
  Finally we output the volume and position of each droplet to
  standard output. */

#if dimension == 2
  for (int j = 0; j < n; j++)
    fprintf (fout, "%d %g %d %g %g %g\n", i, t,
	     j, v[j], b[j].x/v[j], b[j].y/v[j]);
#else
  for (int j = 0; j < n; j++)
    fprintf (fout, "in event droplets: %d, %g, %d, %g, %g, %g, %g, %g, %g, %g\n", i, t,
	     j, v[j], b[j].x/v[j], b[j].y/v[j], b[j].z/v[j], 
       (v1[j] > 0. ? k1[j].x/v1[j] : k[j].x/v[j]), (v1[j] > 0. ? k1[j].y/v1[j] : k[j].y/v[j]), (v1[j] > 0. ? k1[j].z/v1[j] : k[j].z/v[j]));
#endif
  fflush (fout);
}
  
/**
## Mesh adaptation

We adapt the mesh according to the error on the volume fraction field
and the velocity. */

event adapt(i++) {
  adapt_wavelet ({f,u}, (double[]){0.01,uemax,uemax,uemax}, maxlevel);
}
