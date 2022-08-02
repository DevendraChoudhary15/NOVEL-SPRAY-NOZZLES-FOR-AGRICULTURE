#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "view.h"

/*The actual parameters(SI) are as follows:
radius 5e-5
jet velocity 50
liquid density 848
gas density 34.5
liquid viscosity 2.87e-3
gas viscosity 1.97e-5
surface tension 0.03*/

/**
We define the radius of the jet, the initial jet length, the Reynolds
number and the surface tension coefficient. 
And dimentionless parameters ares used in this simulation*/


#define radius 0.05
#define length 0.025
#define Re 1477.35
#define SIGMA 1.42e-5
#define BG 0.7 // light gray for background
#define DG 0. // dark gray

/**
The default maximum level of refinement is 12 and the error threshold
on velocity is 0.1. */

int maxlevel = 12;
int minlevel = 8;
double uemax = 0.01;


scalar f0[];
u.n[left]  = dirichlet(f0[]*1.);
u.t[left]  = dirichlet(0);
#if dimension > 2
u.r[left]  = dirichlet(0);
#endif
p[left]    = neumann(0);
f[left]    = (y < 0.05);

u.n[right] = neumann(0);
p[right]   = dirichlet(0);

u.n[top] = dirichlet(0);
u.t[top] = neumann(0);


/**
The program can take two optional command-line arguments: the maximum
level and the error threshold on velocity. */

int main (int argc, char * argv[])
{
  if (argc > 1)
    maxlevel = atoi (argv[1]);
  if (argc > 2)
    uemax = atof (argv[2]);
  
  init_grid (256);
  origin (0, 0., 0.);
  size (7.);

  /**
  We set the density and viscosity of each phase as well as the
  surface tension coefficient and start the simulation. */
  
  rho1 = 1., 
	rho2 = 1./24.58;
  mu1 = 2.*radius/Re*rho1, 
  mu2 = mu1/145.7;  
  f.sigma = SIGMA;
	//f[nos] = neumann(0);
  run();
}

/**
## nitial conditions */

event init (t = 0) {
  if (!restore (file = "restart")) {
		
    /**
    We use a static refinement down to *maxlevel* in a cylinder 1.2
    times longer than the initial jet and twice the radius. */
    refine (x < 1. && sq(y) + sq(z) < 2.*sq(radius) && level < maxlevel);
    
    /**
    We initialise the auxilliary volume fraction field for a cylinder of
    constant radius. */
    
    fraction (f0, sq(radius) - sq(y) - sq(z));
    f0.refine = f0.prolongation = fraction_refine;
    restriction ({f0}); // for boundary conditions on levels
    
    /**
    We then use this to define the initial jet and its velocity. */

    foreach() {
      f[] = f0[]*(x < 2.*length);
			//f[] = 0.;	
      u.x[] = f[];
    }
    boundary ({f, u.x});
  }
}

/**
## Outputs

We log some statistics on the solver. */

event logfile (i++) {
  if (i == 0)
    fprintf (ferr,
	     "t dt mgp.i mgpf.i mgu.i grid->tn perf.t perf.speed\n");
  fprintf (ferr, "%g %g %d %d %d %ld %g %g\n", 
	   t, dt, mgp.i, mgpf.i, mgu.i,
	   grid->tn, perf.t, perf.speed);
}


/**
We generate an animation using Basilisk View. */

scalar omega[];

event movie (t += 1e-2)
{
#if dimension == 2
  //scalar omega[];
  vorticity (u, omega);
  //view (fov = 12.1038, quat = {0,0,0,1}, tx = -0.48613, ty = -0.0206686, bg = {1,1,1}, width = 1839, height = 1050, samples = 1);
	view (tx = -0.5);
  clear();
  draw_vof ("f");
  squares ("f", max = 1., min = 0., linear = true, map = cool_warm);
#else // 3D
  scalar pid[];
  foreach()
    pid[] = fmod(pid()*(npe() + 37), npe());
  boundary ({pid}); // not used for the moment
  view (camera = "iso",
	fov = 14.5, tx = -0.418, ty = 0.288,
	width = 1600, height = 1200);
  clear();
  draw_vof ("f");
#endif // 3D
  mirror (n = {0, 1, 0}, alpha = 0.) {
    //draw_vof("f");
    draw_vof("f", edges = true, lw = 1.5, lc = {DG, DG, DG}, filled = 1,
             fc = {BG, BG, BG});
    squares ("omega", linear = true, spread = 10, map = cool_warm);
  }
  save ("movie.mp4");
}

/**
We save snapshots of the simulation at regular intervals to
restart or to post-process with [bview](/src/bview). */

event snapshot (t = 0; t += 0.5; t <= 10) {
  char name[80];
  sprintf (name, "snapshot-%g", t);
  scalar pid[];
  foreach()
    pid[] = fmod(pid()*(npe() + 37), npe());
  boundary ({pid});
  dump (name);
}

/**
## Counting droplets*/


event droplets (t += 0.1)
{
  scalar m[];
  foreach()
    m[] = f[] > 1e-3;
  int n = tag (m);


  double v[n];
  coord b[n];
  for (int j = 0; j < n; j++)
    v[j] = b[j].x = b[j].y = b[j].z = 0.;
  foreach_leaf()
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*f[];
      coord p = {x,y,z};
      foreach_dimension()
	b[j].x += dv()*f[]*p.x;
    }

 /**
 When using MPI we need to perform a global reduction to get the
 volumes and positions of droplets which span multiple processes. */

#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  /**
  Finally we output the volume and position of each droplet to
  standard output. */

  for (int j = 0; j < n; j++)
    fprintf (fout, "%d %g %d %g %g %g\n", i, t,
	     j, v[j], b[j].x/v[j], b[j].y/v[j]);
  fflush (fout);
}
  
/**
## Mesh adaptation*/

event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){0.001,uemax,uemax,uemax}, maxlevel);
}
