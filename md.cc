#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <omp.h>
using namespace std;

//parameters

#define L 50 //lattice length
#define dx 10 //lattice spacing
#define Lx (L/dx) //renormalised length
#define N (Lx*Lx*Lx) // 125 particles (Lx*Lx*Lx) //total number of particles
#define NSTEP 5000000 //total number of time steps

//constants

double const vmax     = 2.7; //T_LJ = 2.4
double const dt       = 0.001;
double const rcut     = 3.0;
double const rcut2    = rcut*rcut;
double const halfdtdt = 0.5*dt*dt;
double const halfdt   = 0.5*dt;

//variables

double x[N];
double y[N];
double z[N];
double vx[N];
double vy[N];
double vz[N];
double fx[N];
double fy[N];
double fz[N];
double fx0[N];
double fy0[N];
double fz0[N];
double xi[N];
double yi[N];
double zi[N];
double T; //reduced temperature
double D; //diffusion constant
int i, j;
double r6, r12;

// file output

FILE *file[10];
const char *dirs[]    = { "./data/x.dat" , "./data/y.dat","./data/z.dat" ,"./data/vx.dat","./data/vy.dat","./data/vz.dat","./data/fx.dat","./data/fy.dat","./data/fz.dat", "./data/temperature_diffusion-coefficient.dat"   };
const char *headers[] = { "#X(t) Coordinates\n#Time(t)\t","#Y(t) Coordinates\n#Time(t)\t","#Z(t) Coordinates\n#Time(t)\t","#X(t) Velocities\n#Time(t)\t","#Y(t) Velocities\n#Time(t)\t","#Z(t) Velocities\n#Time(t)\t","#X(t) Forces\n#Time(t)\t","#Y(t) Forces\n#Time(t)\t","#Z(t) Forces\n#Time(t)\t", "#Time(t)\tTemperature (t)\tDiffusion coefficient (t)" };

//functions

void initialize_force_positions_velocities();
void update_force_positions_velocities(int *);
void lennard_jones_force();
void temperature(double *);
void diffusion_coefficient(double *, double *, double *, int *);
void file_write(int *, double *, double *);
void file_start();
void file_end();

// main program

int main(int argc, char *argv[])
{

  double v2, rms, time_spent, r20, T, D=0; int k=0; // double r2, v2, r20, rms, time_spent;

  srand48(time(0));  // time-based seed for random initial particle velocities

  // Open files
  printf("Preparing output files...\n");
  for(int fi=0; fi<10; fi++ ){
    file[fi] = fopen(dirs[fi],"w");   
  }
  file_start();

  //
  // INITIATE: particle positions & velocities & forces //
  //
  printf("Initiating particle's coordinates and velocities...\n");
  initialize_force_positions_velocities();

  file_write(&k, &T, &D);
  
  ////////////////////////
  // Main Loop (t + dt) //
  ////////////////////////
  printf("Simulating particle's dynamics...\n");

  for(int k=1; k<NSTEP; k++){

    if(k%2000==0){
      printf("\t %f percent\n",((double)k/NSTEP)*100.0); // progress 
    }
    
    //
    // Update: particle positions, velocities & forces
    //
    update_force_positions_velocities( &k );

    //
    // Macroscopic Properties: (i) temperature, (ii) diffusion
    //
    temperature( &T ); 

    diffusion_coefficient( &D, &rms, &r20, &k );

    if (k%1000==0){
    file_write( &k, &T, &D ); // write every 1000 time-steps
    }
  }

  //
  // Finish 
  //
  file_end();

  k=NSTEP;

  srand48(time(0));  // time-based seed for random initial particle velocities

  printf("Temperature: %f\tr2=%f\trms=%f\tDiffusion: %f\n", T, r20/N, sqrt(rms), D);

  return 0;
}


///////////////
// FUNCTIONS //
///////////////

void initialize_force_positions_velocities() 
{
 /* (i) Position initialisation*/
 int l;
 for(int i = 0; i < Lx; i++){
    for(int j = 0; j < Lx; j++){
      for(int k = 0; k < Lx; k++){
        l = (i*Lx*Lx + j*Lx + k);
        x[l] = i*dx;
        y[l] = j*dx;
        z[l] = k*dx;
      }
    }  
  }
 /* (ii) Velocity initialisation*/
 double vxcm, vycm, vzcm = 0.;
 for(int i = 0; i < N; i++){
    vx[i] = vmax*(2.*drand48()-1.); // <V2a>=Kb*T/m  
    vy[i] = vmax*(2.*drand48()-1.);
    vz[i] = vmax*(2.*drand48()-1.);
    vxcm += vx[i]; // summing velocities (for subtracting center of mass V)
    vycm += vy[i];
    vzcm += vz[i];
  }
 /* (iii) Subtract center of mass velocities*/
 vxcm /= N; // averaging velocities (for subtracting center of mass V)
 vycm /= N;
 vzcm /= N;
 for(int i = 0; i < N; i++){
    vx[i] -= vxcm;
    vy[i] -= vycm;
    vz[i] -= vzcm;
  }
  /* (v) Force initialization */
  lennard_jones_force();
  /* (vi) Temperature initialization */
  temperature(&T);
}

void lennard_jones_force() //total force
{ 
  double r2, xmin, ymin, zmin, diff, f;
  
  for(int i = 0; i<N; i++){
    fx[i] = 0;
    fy[i] = 0;
    fz[i] = 0;
  }

  # pragma omp parallel \
  shared ( fx,fy,fz, x,y,z, vx,vy,vz ) \
  private ( i, j, xmin,ymin,zmin, r2, r6, r12, f )

  # pragma omp for
  //Pairwise forces, particles i vs. j 
  for(int i = 0; i<N; i++){
    for(int j = i+1; j<N; j++){

      // Periodic Boundary Conditions
      xmin = ((x[i]-x[j]) - L*round((x[i]-x[j])/L));
      ymin = ((y[i]-y[j]) - L*round((y[i]-y[j])/L));
      zmin = ((z[i]-z[j]) - L*round((z[i]-z[j])/L));

      // Sq. distance, i vs j
      r2 = xmin*xmin + ymin*ymin + zmin*zmin;

      // Ignore i,j particles's forces > rcut2 distance
      if(r2 < rcut2){
        // r6  = r2*r2*r2;
        // r12 = r6*r6;
        f = (48./pow(r2,8) - 24./pow(r2,4)); //LJ-potential

        fx[i] += f*xmin;
        fy[i] += f*ymin;
        fz[i] += f*zmin;
        fx[j] -= f*xmin;
        fy[j] -= f*ymin;
        fz[j] -= f*zmin;
      }
    }
  }
}

void update_force_positions_velocities(int *k) //one time step move dt
{ 
  // UPDATE PARTICLES': Positions (x,y,z), then velocities (vx,vy,vz) according to forces

  // MOVE_R: update particle positions
  for(int i = 0; i < N; i++){
    x[i] += dt*vx[i] + halfdtdt*fx[i]; //updated positions
    y[i] += dt*vy[i] + halfdtdt*fy[i];
    z[i] += dt*vz[i] + halfdtdt*fz[i];
    
    fx0[i] = fx[i]; // fx0, fy0, fz0 store forces at f(t)
    fy0[i] = fy[i];
    fz0[i] = fz[i];
  }
  
  // FORCE: update forces
  lennard_jones_force(); // fx, fy and fz stores forces for f(t+1)
  
  // MOVE_V: update particle velocities
  for(int i = 0; i < N; i++){
    vx[i] += halfdt*(fx0[i] + fx[i]); //updated velocities
    vy[i] += halfdt*(fy0[i] + fy[i]);
    vz[i] += halfdt*(fz0[i] + fz[i]);
  }

  // positions at t0 = k = 1000 stored for diffusion coefficient
  if(*k==1000){ 
    for(int i = 0; i < N; i++){ 
        xi[i]=x[i]; yi[i]=y[i]; zi[i]=z[i];
    }
  } 
}

void temperature(double *T)
{
  // Temperatures calculated over all molecules
  double v2 = 0.;
  for(int i=0; i<N; i++){
    v2 += (vx[i]*vx[i]) + (vy[i]*vy[i]) + (vz[i]*vz[i]);
  }
  *T = v2/(3*N);
}

void diffusion_coefficient(double *D, double *rms, double *r20, int *k)
{
  // Diffusion coefficients calculated over all molecules
  double r2 = 0;
  *r20 = 0;
  for(int i = 0; i < N; i++){
    r2   += x[i]*x[i]   + y[i]*y[i]   + z[i]*z[i];  
    *r20 += xi[i]*xi[i] + yi[i]*yi[i] + zi[i]*zi[i];
  }
  *rms = (r2 - *r20)/N;  
  *D   = *rms/(6.*(*k)*dt); 
}

void file_start(){
  // Files to store datarows (t), 1 particle per ith column
  for(int fi=0; fi<9; fi++ )
  {
    fprintf(file[fi],"%s",headers[fi]);
    for (int i = 0; i < N; i++){
      fprintf(file[fi],  "Molecule %i\t", i); 
    }
    fprintf(file[fi],  "\n");
  }
  // temperature
  int fi = 9;
  fprintf(file[fi], "%s\n", headers[fi]); // time
}

void file_write(int *k, double *T, double *D){

  // @test: to save time
  // Writing datarows (t), 1 particle per ith column
  double *data[] = {x,y,z,vx,vy,vz,fx,fy,fz};

  for(int fi=0; fi<9; fi++ ){
    fprintf(file[fi], "%f\t", *k*dt); // time 
    for (int i = 0; i<N; i++){ // for each molecule
      fprintf(file[fi], "%f\t", data[fi][i]); 
    }
    fprintf(file[fi], "\n");
  }  
  // temperature, diffusion coefficient  
  int fi = 9;
  fprintf(file[fi], "%f\t%f\t%f\n", *k*dt, *T, *D); // time  
}

void file_end(){
  // close output files 
  for(int fi=0; fi<10; fi++ ){
    fclose(file[fi]);
  }  
}
