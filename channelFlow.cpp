#include <omp.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>

#include <vector>

using namespace std;

int main(int argc, char *argv[])
{
  //----------------------------------------------------------------
  //      Parameter Definitions
  //----------------------------------------------------------------
  double u_0 = 0.5 / 10.0;            //max speed of the wall
  double u_w = 0.0;             //current speed of the wall 

  double n = M_PI / 10.0;
  //-----------------------------Geometry---------------------------------
	double height=1.0;            //half channel height in [m]
	double length=0.5;            //channel length
	double deltaX = 0.01 * atof(argv[1]);        //grid spacing
	double epsilon=1e-8;          //geometrical tolerance

	long nx=ceil(length/deltaX+epsilon)+1;   //number of nodes in x direction
	long ny=ceil(2*height/deltaX+epsilon)+1; //number of nodes in y direction

  //-----------------------------Time---------------------------------
  double deltaT = 5e-04 * pow(atof(argv[1]),2);
	long timeSteps=25.0 / deltaT;        //time steps to go
	int writeInterval=10000;      // time interval to write (nt = pi)
  int minTimeStep= 20000000;      // minimal timestep to start writing

  cout << "MaxInterval: " << timeSteps << endl;

	int cpus = 16;                 //number of parallel threads
	omp_set_num_threads(cpus);    //set number of threads

  //-----------------------------Fluid---------------------------------
	double maxExpectedVelocity=u_0; //max. expected velocity
  
	double rho=900;               //fluid density [kg/m^3] 
	double viscosity=3e-3;        //kinematic fluid viscosity

  //-----------------------------Calc---------------------------------
  double xsi0=deltaX/deltaT; //molecular velocity
  double machNumber = sqrt(3.) * deltaT * maxExpectedVelocity / deltaX;
  double speedOfSound = xsi0 / sqrt(3.);
  double omega = pow(speedOfSound, 2) * deltaT / (viscosity + 0.5 * pow(speedOfSound, 2) * deltaT);


  //----------------------------------------------------------------
  //      Write Parameters to StdOut
  //----------------------------------------------------------------

	cout << "DOMAIN" << endl;
	cout << "height " << height << endl;
	cout << "length " << length << endl;
	cout << "deltaX " << deltaX << endl;
	cout << "nx " << nx << endl;
	cout << "ny " << ny << endl << endl;
  
	cout << "FLUID" << endl;
	cout << "rho " << rho << endl;
	cout << "viscosity " << viscosity << endl << endl;

	cout << "TIME" << endl;
	cout << "timeSteps " << timeSteps << endl;
	cout << "deltaT " << deltaT << endl << endl;

	//check omega for stability reasons
	cout << "STABILITY" << endl;
	cout << "omega " << omega << endl;
	cout << "speedOfSound " << speedOfSound << endl;
	cout << "machNumber " << machNumber << endl << endl;

  //----------------------------------------------------------------
  //      Data Structure Definition
  //----------------------------------------------------------------

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

  //----------------------------------------------------------------
  //      VTK Writer
  //----------------------------------------------------------------
  
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

  //----------------------------------------------------------------
  //      LBM Algorithm
  //----------------------------------------------------------------

	//init time step counter for use outside of loop
	long t = 0;
  vector<double> verror;
  vector<vector<double>> vvelo_write;

	//outer loop over all time steps
	for (t=0;t<timeSteps;t++)
	{
    //calculate time dependent wall speed
    u_w = u_0 * cos(n * t * deltaT);

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

    //-----------------------------Periodic BC---------------------------------
    #pragma omp parallel for 
    for (int l = 1; l < ny - 1; ++l) { //loop over all nodes in y direction
      // - Left Boundary Transport Step
      for (int i = 0; i < 9; i++) //all lattice directions; expected zero velocity
      {
        int kNN = 0; //index of next neighbor in lattice
        int lNN = 0; //index of next neighbor in lattice
        int k = 0;
      
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
          kNN = nx - 1;
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
          kNN = nx - 1;
          lNN = l + 1;
          break;
        }
        case 6:
        {
          kNN = nx - 1;
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
      }//end i 
      for (int i = 0; i < 9; i++) //all lattice directions; expected zero velocity
      {
        int kNN = 0; //index of next neighbor in lattice
        int lNN = 0; //index of next neighbor in lattice
        int k = nx - 1; 
        switch (i)
        {
        case 0:
        {
          kNN = 0;
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
          kNN = 0;
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
          kNN = 0;
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
    }
  
    //-----------------------------Bottom: Moving Wall (Bounce Back)---------------------------------
    double veloFactor = 2 * u_w / pow(speedOfSound,2);

    #pragma omp parallel for
    for (int k=1;k<(nx-1);k++) //all nodes in x direction in this line; except corner nodes
    {
	    int l=0; //only bottom wall nodes (first line in y direction)
    	//only lattice directions in +y direction --> bounce back
    	//i=1
    	distribution[k][l+1][0][1] = distribution[k][l][0][3] + weigths[3] * xsi[3][0] * fluidDensity[k][l+1] * veloFactor;
    	//i=4
    	distribution[k+1][l+1][0][4]=distribution[k][l][0][6] + weigths[6] * xsi[6][0] * fluidDensity[k+1][l+1] * veloFactor;
    	//i=6
    	distribution[k-1][l+1][0][5]=distribution[k][l][0][7] + weigths[7] * xsi[7][0] * fluidDensity[k-1][l+1] * veloFactor;
    } //end k
    
    //left corner
    int l = 0;
    int k = 0;
    double factor = 2 * (fluidDensity[0][0] / pow(speedOfSound,2)) * u_w;
    distribution[k][l+1][0][1] = distribution[k][l][0][3] + weigths[3] * xsi[3][0] * fluidDensity[k][l+1] * veloFactor;
    distribution[k+1][l+1][0][4] = distribution[k][l][0][6] + weigths[6] * xsi[6][0] * fluidDensity[k+1][l+1] * veloFactor;
    distribution[nx-1][l+1][0][5] = distribution[k][l][0][7] + weigths[7] * xsi[7][0] * fluidDensity[nx-1][l+1] * veloFactor;

    //right corner
    k = nx - 1;
    factor = 2 * (fluidDensity[nx-1][0] / pow(speedOfSound,2)) * u_w;
    distribution[k][l+1][0][1] = distribution[k][l][0][3] + weigths[3] * xsi[3][0] * fluidDensity[k][l+1] * veloFactor;
    distribution[0][l+1][0][4]= distribution[k][l][0][6] + weigths[6] * xsi[6][0] * fluidDensity[0][l+1] * veloFactor;
    distribution[k-1][l+1][0][5]= distribution[k][l][0][7] + weigths[7] * xsi[7][0] * fluidDensity[k-1][l+1] * veloFactor;
  
    //-----------------------------Top (Specular Reflection)---------------------------------
    #pragma omp parallel for
    for (int k = 0; k < nx; ++k) {
      //fix y to top boundary
      int l = ny - 1; 
      // specular reflection of wall distributions --> positive y-direction
      distribution[k][l - 1][0][3] = distribution[k][l][0][1];  //i = 1
      distribution[k][l - 1][0][6] = distribution[k][l][0][5];  //i = 6
      distribution[k][l - 1][0][7] = distribution[k][l][0][4];  //i = 7
    }
    //----------------------------Moments----------------------------------
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
  
    //-----------------------------Collision Step---------------------------------
	  #pragma omp parallel for
	  for (int k=0; k<nx; k++)  //all fluid nodes in x direction
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

    //convergence measure
    double error = 0.0;
    double nus;
    int testk = nx /2 + 1;
    for (int i = 1; i < ny-1; ++i) {
      nus = (deltaX * (float(i) - 0.5)) * sqrt(n / (2 * viscosity));
      error += pow(fabs(u_0 * exp(-nus) * cos(n * t * deltaT - nus) - fluidVelocity[testk][i][0]),2);
    }
    error = sqrt(error / (ny - 2));
    verror.push_back(error);
	  //do convergence measure --> 
	  if (t % 1000 == 0)
	  {
		//use node at half channel height after 2/3 of channel length
		cout << "time step: "<< t << " error [%]: "<< error << endl;
	  }
	  //write results for appropiate time steps
	  if ((t % writeInterval == 0) && (t > minTimeStep)){
      writeResults(t, deltaT, nx, ny, deltaX, fluidDensity, speedOfSound, fluidVelocity);
      vector<double> vvelo;
      int testk = nx /2 + 1;
      vvelo.push_back(t);
      for (int i = 0; i < ny; ++i) {
        vvelo.push_back(fluidVelocity[testk][i][0]);         
      }
      vvelo_write.push_back(vvelo);
    }
	} //end time loop

	//write final result
	//writeResults(t, deltaT, nx, ny, deltaX, fluidDensity, speedOfSound, fluidVelocity);

  //write velocity vector
  ofstream ofile;
  ofile.open("velo.log");
  for (int i = 0; i < vvelo_write.size(); ++i) {
    for (int j = 0; j < vvelo_write[0].size(); ++j) {
      ofile << vvelo_write[i][j] << ',';
    }  
    ofile << endl;
  }
  ofile.close();

  //write error vector
  ofile.open("error.log");
  for (int i = 0; i < verror.size(); ++i) {
    ofile << i * deltaT << "," << verror[i] << endl;
  }
  ofile.close();
}//end main


