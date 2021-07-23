#include <omp.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>

using namespace std;

int main()
{
	//example of LBM implementation for channel flow (2D) in C - as introduced in LBM lecture
	//initialization parameters
	//grid parameters
	double height=1.0; //half channel height in [m]
	double length=0.5; //channel length
	double deltaX = 0.025; //grid spacing
	//fluid parameters
	double rho=900; //fluid density [kg/m^3] 
	double viscosity=3e-3; //kinematic fluid viscosity
	//flow parameters
	double Re=100; //Reynolds number
	double meanVelocity=Re*viscosity/(2*height); //mean velocity
	double inflowVelocity[2]={meanVelocity, 0};//inflow velocity --> const
	//time step control
	double machNumber=0.01; //Mach number
	//calculation duration
	long timeSteps=200000; //time steps to go
	//output intervall
	int writeInterval=20000;
	//cpus for shared memory parallelization
	int cpus = 8;
	omp_set_num_threads(cpus); //set number of threads

	//calc lattice dependent variables
	double maxExpectedVelocity=1.5*meanVelocity; //max. expected velocity
	double deltaT=machNumber*deltaX/(sqrt(3.)*maxExpectedVelocity); //time step from Ma condition
	double xsi0=deltaX/deltaT; //molecular velocity
	double speedOfSound=xsi0/sqrt(3.); //isothermal speed of sound (in lattice)
	double omega=pow(speedOfSound,2)*deltaT/(viscosity+0.5*pow(speedOfSound,2)*deltaT); //relaxation parameter for collision
	 //check omega for stability reasons
	  cout << "Omega is: " << omega << endl;
	double epsilon=1e-8; //geometrical tolerance
	long nx=ceil(length/deltaX+epsilon)+1; //number of nodes in x direction
	long ny=ceil(2*height/deltaX+epsilon)+1; //number of nodes in y direction

	//allocate and init matrices: 
	//first array index is x-direction; second array index is y-direction
	double** fluidDensity=new double*[nx]; //matrix for fluid density at every node
	for(int k=0;k<nx;k++)
	{
		fluidDensity[k]=new double[ny];
		for(int l=0;l<ny;l++)
		{
			fluidDensity[k][l]=rho; //init values with given density
		}
	}

	double*** fluidVelocity=new double**[nx]; //matrix for fluid velocity at every node; third array index is velocity direction (0:x-dir;1:y-dir)
	for(int k=0;k<nx;k++)
	{
		fluidVelocity[k]=new double*[ny];
		for(int l=0;l<ny;l++)
		{
			fluidVelocity[k][l]=new double[2];
			fluidVelocity[k][l][0]=0.; //init values with 0. --> flow is at rest at t=0
			fluidVelocity[k][l][1]=0.;
		}
	}

	double weigths[9]={1./9.,1./9.,1./9.,1./9.,1./36.,1./36.,1./36.,1./36.,4./9.}; //weigths of lattice directions for eq

	double**** distribution=new double***[nx]; //distribution at every node; third array index: transport or collision distribution; fourth index: lattice direction
	for(int k=0;k<nx;k++)
	{
		distribution[k]=new double**[ny];
		for(int l=0;l<ny;l++)
		{
			distribution[k][l]=new double*[2];
			distribution[k][l][0]=new double[9];
			distribution[k][l][1]=new double[9];
			//initialize with eq
			for (int i=0;i<9;i++) //all lattice directions
			{
			 //init with rho and u=0;
				distribution[k][l][0][i]=weigths[i]*rho; //transport array
				distribution[k][l][1][i]=weigths[i]*rho; //collision array
			}//end i
		}
	}
	double xsi[9][2];//lattice velocities --> first array index: lattice direction; second index velocity direction (0:x-dir;1:y-dir)
	xsi[0][0]=xsi0; xsi[1][0]=0; xsi[2][0]=-xsi0; xsi[3][0]=0; xsi[4][0]=xsi0; xsi[5][0]=-xsi0; xsi[6][0]=-xsi0; xsi[7][0]=xsi0; xsi[8][0]=0;
	xsi[0][1]=0; xsi[1][1]=xsi0; xsi[2][1]=0; xsi[3][1]=-xsi0; xsi[4][1]=xsi0; xsi[5][1]=xsi0; xsi[6][1]=-xsi0; xsi[7][1]=-xsi0; xsi[8][1]=0;

	//result visualization --> use vtk file here
	//bad style; used lambda here only for simplicity!
	auto writeResults = [](long t, double deltaT, long nx, long ny, double deltaX, double** fluidDensity, double speedOfSound, double*** fluidVelocity)
	{
		string filename("results_unstable" + to_string(static_cast<long>(t)) + ".vtk");

		ofstream results(filename.c_str());

    results << setiosflags(ios::left | ios::showpoint | ios::fixed) << "# vtk DataFile Version 2.0" << endl;
    results << "Lattice Boltzmann solution at timestep " << t * deltaT << endl;
    results << "ASCII" << endl;
    results << "DATASET STRUCTURED_POINTS" << endl;
    results << "DIMENSIONS " << nx << " " << ny << " " << 1 << endl;
    results << "SPACING " << deltaX << " " << deltaX << " " << 0 << endl;
    results << "ORIGIN 0 0 0" << endl;
    results << "POINT_DATA " << nx * ny << endl;
    results << "SCALARS pressure float 1" << endl;
    results << "LOOKUP_TABLE default" << endl;

    for (int l = 0; l < ny; l++)
      for (int k = 0; k < nx; k++)
        results << setprecision(8) << fluidDensity[k][l] * pow(speedOfSound, 2) << endl;

    results << "SCALARS density float 1" << endl;
    results << "LOOKUP_TABLE default" << endl;

    for (int l = 0; l < ny; l++)
      for (int k = 0; k < nx; k++)
        results << setprecision(8) << fluidDensity[k][l] << endl;

    results << "VECTORS velocity float" << endl;

    for (int l = 0; l < ny; l++)
      for (int k = 0; k < nx; k++)
        results << setprecision(8) << fluidVelocity[k][l][0] << " " << setprecision(8) << fluidVelocity[k][l][1] << " " << 0.0 << endl;
		results.close();
	};//end writeResults

	//init time step counter for use outside of loop
	long t = 0;
	//outer loop over all time steps
	for (t=0;t<timeSteps;t++)
	{
	  //transport step for all fluid nodes
	  #pragma omp parallel for
	  for (int k = 1; k < (nx - 1); k++)  //all fluid nodes in x direction
	  {
	  	for (int l = 1; l < (ny - 1); l++) //all fluid nodes in y direction
	  	{
	  	  for (int i = 0; i < 9; i++) //all lattice directions; expected zero velocity
	  	  {
	  	  	int kNN = 0; //index of next neighbor in lattice
	  	  	int lNN = 0; //index of next neighbor in lattice
	  	  
	  	  	switch (i)
	  	  	{
	  	  	case 0:
	  	  	{
	  	  		kNN = k + 1;
	  	  		lNN = l;
	  	  		break;
	  	  	}
	  	  	case 1:
	  	  	{
	  	  		kNN = k;
	  	  		lNN = l + 1;
	  	  		break;
	  	  	}
	  	  	case 2:
	  	  	{
	  	  		kNN = k - 1;
	  	  		lNN = l;
	  	  		break;
	  	  	}
	  	  	case 3:
	  	  	{
	  	  		kNN = k;
	  	  		lNN = l - 1;
	  	  		break;
	  	  	}
	  	  	case 4:
	  	  	{
	  	  		kNN = k + 1;
	  	  		lNN = l + 1;
	  	  		break;
	  	  	}
	  	  	case 5:
	  	  	{
	  	  		kNN = k - 1;
	  	  		lNN = l + 1;
	  	  		break;
	  	  	}
	  	  	case 6:
	  	  	{
	  	  		kNN = k - 1;
	  	  		lNN = l - 1;
	  	  		break;
	  	  	}
	  	  	case 7:
	  	  	{
	  	  		kNN = k + 1;
	  	  		lNN = l - 1;
	  	  		break;
	  	  	}
	  	  	case 8:
	  	  	{
	  	  		kNN = k;
	  	  		lNN = l;
	  	  		break;
	  	  	}
	  	  	}//end switch
	  	  
	  	  	distribution[kNN][lNN][0][i] = distribution[k][l][1][i]; //transport of distribution
	  	  } //end i
	  	} //end l
	  } //end k
  
	  //bottom wall --> half way bounce back
	  #pragma omp parallel for
	  for (int k=1;k<(nx-1);k++) //all nodes in x direction in this line; except corner nodes
	  {
		int l=0; //only bottom wall nodes (first line in y direction)
		//only lattice directions in +y direction --> bounce back
		//i=1
		distribution[k][l+1][0][1]=distribution[k][l][0][3];
		//i=4
		distribution[k+1][l+1][0][4]=distribution[k][l][0][6];
		//i=6
		distribution[k-1][l+1][0][5]=distribution[k][l][0][7];
	  } //end k
	  //do corner nodes
	  distribution[1][1][0][4]=distribution[0][0][0][6]; //k=0; i=4
	  distribution[nx-2][1][0][5]=distribution[nx-1][0][0][7]; //k=nx-1; i=5
  
	  //top wall --> half way bounce back
	  #pragma omp parallel for
	  for(int k=1;k<(nx-1);k++) //all nodes in x direction in this line; except corner nodes
	  {
		int l=ny-1; //only top wall nodes (last line in y direction)
		//only lattice directions in -y direction --> bounce back
		//i=3
		distribution[k][l-1][0][3]=distribution[k][l][0][1];
		//i=6
		distribution[k-1][l-1][0][6]=distribution[k][l][0][4];
		//i=7
		distribution[k+1][l-1][0][7]=distribution[k][l][0][5];
	  } //end k
	  //do corner nodes
	  distribution[1][ny-2][0][7]=distribution[0][ny-1][0][5]; //k=0; i=7
	  distribution[nx-2][ny-2][0][6]=distribution[nx-1][ny-1][0][4]; //k=nx; i=4


	  //inlet nodes --> eq distribution
	  #pragma omp parallel for
	  for (int l=1;l<(ny-1);l++) //all nodes in y direction in this column (expect bottom and top node)
	  {
		int  k=0; //only inlet nodes (first column in x direction)
		for (int i=0;i<9;i++) //all lattice directions; 
		{
		  int kNN=0; //index of next neighbor in lattice
		  int lNN=0; //index of next neighbor in lattice
		  //calc eq
		  double dotProductXsiVelo=xsi[i][0]*inflowVelocity[0]+xsi[i][1]*inflowVelocity[1];
		  double dotProductVelo=pow(inflowVelocity[0],2)+pow(inflowVelocity[1],2);
		  double density=fluidDensity[k+1][l]; //extrapolate density from normal neighbor
		  double eqDistribution=weigths[i]*density*(1.+dotProductXsiVelo/pow(speedOfSound,2)+pow(dotProductXsiVelo,2)/(2.*pow(speedOfSound,4)) - dotProductVelo/(2*pow(speedOfSound,2))); //eq with fixed density and extrapolated velocity
		  //set both distributions at inlet node to eq
		  distribution[k][l][0][i]=eqDistribution;
		  distribution[k][l][1][i]=eqDistribution;
		  switch (i) //do transport for +x directions
		  {
			case 0:
			{
				kNN=k+1; 
				lNN=l;
				break;
			}
			case 4:
			{
				kNN=k+1; 
				lNN=l+1;
				break;
			}
			case 7:
			{
				kNN=k+1; 
				lNN=l-1;   
				break;
			}
			default: continue; //continue in loop for other directions
		  }//end switch
		  //do transport if necessary
		  distribution[kNN][lNN][0][i]=distribution[k][l][1][i]; //transport of distribution
		} //end i
	  } //end l

	  //outlet nodes --> eq distribution
	  #pragma omp parallel for
	  for (int l=1;l<(ny-1);l++) //all nodes in y direction in this column (expect bottom and top node)
	  {
		int  k=nx-1; //only inlet nodes (first column in x direction)
		for (int i=0;i<9;i++) //all lattice directions; 
		{
		  int kNN=0; //index of next neighbor in lattice
		  int lNN=0; //index of next neighbor in lattice
		  //calc eq
		  double veloX=fluidVelocity[k-1][l][0];	//extrapolate velocity
		  double veloY=fluidVelocity[k-1][l][1];
		  double dotProductXsiVelo=xsi[i][0]*veloX+xsi[i][1]*veloY;
		  double dotProductVelo=pow(veloX,2)+pow(veloY,2);
		  double eqDistribution=weigths[i]*rho*(1.+dotProductXsiVelo/pow(speedOfSound,2)+pow(dotProductXsiVelo,2)/(2.*pow(speedOfSound,4)) - dotProductVelo/(2*pow(speedOfSound,2))); //eq with fixed density and extrapolated velocity
		  //set both distributions at inlet node to eq
		  distribution[k][l][0][i]=eqDistribution;
		  distribution[k][l][1][i]=eqDistribution;
		  switch (i) //do transport for -x directions
		  {
			case 2:
			{
				kNN=k-1; 
				lNN=l;
				break;
			}
			case 5:
			{
			  kNN=k-1; 
			  lNN=l+1;
			  break;
			}
			case 6:
			{
			  kNN=k-1; 
			  lNN=l-1;   
			  break;
			}
			default: continue; //continue in loop for other directions
			}//end switch
		  //do transport if necessary
		  distribution[kNN][lNN][0][i]=distribution[k][l][1][i]; //transport of distribution
		} //end i
	  } //end l
  
	  //calc moments
	  for(int k=0;k<nx;k++)  //all nodes in x direction
	  {
		for (int l=1;l<(ny-1);l++) //all nodes in y direction; except of wall nodes
		{
		  fluidDensity[k][l]=0.; //reset value
		  fluidVelocity[k][l][0]=0.; //reset value
		  fluidVelocity[k][l][1]=0.; //reset value
		  for (int i=0;i<9;i++) //all lattice directions
		  {
		   fluidDensity[k][l]    =fluidDensity[k][l]+distribution[k][l][0][i]; //density
		   fluidVelocity[k][l][0]=fluidVelocity[k][l][0]+xsi[i][0]*distribution[k][l][0][i]; //momentum x direction
		   fluidVelocity[k][l][1]=fluidVelocity[k][l][1]+xsi[i][1]*distribution[k][l][0][i]; //momentum x direction
		  } //end i

		  fluidVelocity[k][l][0]=fluidVelocity[k][l][0]/fluidDensity[k][l]; //calc velocity from momentum
		  fluidVelocity[k][l][1]=fluidVelocity[k][l][1]/fluidDensity[k][l]; //calc velocity from momentum
		} //end l
	  } //end k
  
	  //collision step
	  #pragma omp parallel for
	  for (int k=1;k<(nx-1);k++)  //all fluid nodes in x direction
	  {
		for (int l=1;l<(ny-1);l++) //all fluid nodes in y direction
		{
		  for (int i=0;i<9;i++) //all lattice directions;  
		  {
			//calc eq for this direction
			double dotProductXsiVelo=xsi[i][0]*fluidVelocity[k][l][0]+xsi[i][1]*fluidVelocity[k][l][1];
			double dotProductVelo=pow(fluidVelocity[k][l][0],2)+pow(fluidVelocity[k][l][1],2);
			double eqDistribution=weigths[i]*fluidDensity[k][l]*(1.+dotProductXsiVelo/pow(speedOfSound,2)+pow(dotProductXsiVelo,2)/(2.*pow(speedOfSound,4)) - dotProductVelo/(2*pow(speedOfSound,2)));
			distribution[k][l][1][i]=distribution[k][l][0][i]+omega*(eqDistribution-distribution[k][l][0][i]); //collision step
		  } //end i
		} //end l
	  } //end k
  
	  //do convergence measure --> 
	  if (t % 1000 == 0)
	  {
		//use node at half channel height after 2/3 of channel length
		int nxRes = 2 * nx / 3;
		int nyRes = ny / 2;
		double yPosition = (nyRes - 0.5) * deltaX; //y position of nyRes in channel
		double analyticalSolution = 1.5 * meanVelocity * (2 * yPosition/height - pow(yPosition,2) / pow(height,2)); //analytical solution for x-velocity
		double err = sqrt(pow(fluidVelocity[nxRes][nyRes][0] - analyticalSolution,2) + pow(fluidVelocity[nxRes][nyRes][1], 2)) / analyticalSolution * 100; //calc error as RMS value
		cout << "time step: "<< t << " error [%]: " << err << endl;
		cout << "numerical solution " << fluidVelocity[nxRes][nyRes][0] << endl << endl;
	  }
	  //write results for appropiate time steps
	  if (t % writeInterval == 0) writeResults(t, deltaT, nx, ny, deltaX, fluidDensity, speedOfSound, fluidVelocity);
	} //end time loop

	//write final result
	writeResults(t, deltaT, nx, ny, deltaX, fluidDensity, speedOfSound, fluidVelocity);
}//end main



