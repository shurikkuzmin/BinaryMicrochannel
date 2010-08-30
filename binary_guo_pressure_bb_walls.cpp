//This file is intended to try to use Zou-He pressure BC for both fields.

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>

//Domain size
const int NY=51;
const int NX=751;

//Time steps
const int N=40000;
const int NOUTPUT=1000;

//Fields and populations
double f[NX][NY][9], f2[NX][NY][9], feq[NX][NY][9], g[NX][NY][9], g2[NX][NY][9],geq[NX][NY][9];
double rho[NX][NY],ux[NX][NY],uy[NX][NY],phase[NX][NY];

//Laplacians and gradients
double laplace[NX][NY],gradx[NX][NY],grady[NX][NY];

//Pressure boundary conditions
double rho_inlet=1.0135;
double rho_outlet=1.0;
double phase_inlet=1.0;
double phase_outlet=1.0;

double force_x=0.000;
double force_y=0.000;

//Binary-liquid parameters
double aconst=0.04;
double kconst=0.04;
double gammaconst=1.0;
double tau_liq=2.5;
double tau_gas=0.7;

//Wall wettability parameter
double wall_gradient=0.0;

//BGK relaxation parameter
double omega=1.0;
int width=10;

//Magic Irina's parameters
double omegaginzburg=8.0*(2.0-omega)/(8.0-omega);
double omegamat[]={1.0,1.0,1.0,omega,omega,omega,1.0,omegaginzburg,omegaginzburg};

//Underlying lattice parameters
double weights[]={4.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};
int cx[]={0,1,0,-1,0,1,-1,-1,1};
int cy[]={0,0,1,0,-1,1,1,-1,-1};
double gmat[]={1.0,-2.0,-2.0,-2.0,-2.0,4.0,4.0,4.0,4.0};
int compliment[]={0,3,4,1,2,7,8,5,6};
float wxx[] = {0.0, 1.0/3.0, -1.0/6.0, 1.0/3.0, -1.0/6.0, -1.0/24.0, -1.0/24.0, -1.0/24.0, -1.0/24.0};
float wyy[] = {0.0, -1.0/6.0, 1.0/3.0, -1.0/6.0, 1.0/3.0, -1.0/24.0, -1.0/24.0, -1.0/24.0, -1.0/24.0};
float wxy[] = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0/4.0, -1.0/4.0, 1.0/4.0, -1.0/4.0};

float gradstencilx[9]={0.0,4.0/12.0,0.0,-4.0/12.0,0.0,
			          1.0/12.0,-1.0/12.0,-1.0/12.0,1.0/12.0};

float gradstencily[9]={0.0,0.0,4.0/12.0,0.0,-4.0/12.0,
			          1.0/12.0,1.0/12.0,-1.0/12.0,-1.0/12.0};

float laplacestencil[9]={-20.0/6.0,4.0/6.0,4.0/6.0,4.0/6.0,4.0/6.0,
					   1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0};

//Matrix M
double M[9][9];


void matrix_init()
{
	for (int iCoor=0;iCoor<9;iCoor++)
	{
		M[0][iCoor]=1.0;
		M[1][iCoor]=cx[iCoor]*sqrt(3.0);
		M[2][iCoor]=cy[iCoor]*sqrt(3.0);
		M[3][iCoor]=(cx[iCoor]*cx[iCoor]-1.0/3.0)*3.0/sqrt(2.0);
		M[4][iCoor]=cx[iCoor]*cy[iCoor]*3.0;
		M[5][iCoor]=(cy[iCoor]*cy[iCoor]-1.0/3.0)*3.0/sqrt(2.0);
		M[6][iCoor]=gmat[iCoor]/2.0;
		M[7][iCoor]=gmat[iCoor]*cx[iCoor]*sqrt(1.5)/2.0;
		M[8][iCoor]=gmat[iCoor]*cy[iCoor]*sqrt(1.5)/2.0;
	}
}

void writedensity(std::string const & fName)
{
	std::string fullName = "./tmp/" + fName+ ".dat";
	std::ofstream fout(fullName.c_str());
	fout.precision(10);

	for (int iY=NY-1; iY>=0; --iY)
	{
		for (int iX=0; iX<NX; ++iX)
			fout<<rho[iX][iY]<<" ";
		fout<<"\n";
	}

}

void writephase(std::string const & fName)
{
	std::string fullName = "./tmp/" + fName+ ".dat";
	std::ofstream fout(fullName.c_str());
	fout.precision(10);

	for (int iY=NY-1; iY>=0; --iY)
	{
		for (int iX=0; iX<NX; ++iX)
			fout<<phase[iX][iY]<<" ";
		fout<<"\n";
	}

}


void writevelocity(std::string const & fName)
{
	std::string fullName = "./tmp/" + fName+ ".dat";
	std::ofstream fout(fullName.c_str());
	fout.precision(10);

	for (int iY=NY-1; iY>=0; --iY)
	{
		for (int iX=0; iX<NX; ++iX)
			fout<<ux[iX][iY]<<" ";
		fout<<"\n";
	}

}

void writepopulationslice(std::string const & fName,int iPop)
{
	std::stringstream fullNamef;
    std::stringstream fullNameg;

	fullNamef<<"./tmp/"<<fName<<"f"<<iPop<<".dat";
	std::ofstream fout(fullNamef.str().c_str());
	fullNameg<<"./tmp/"<<fName<<"g"<<iPop<<".dat";
	std::ofstream gout(fullNameg.str().c_str());
	fout.precision(10);
	gout.precision(10);

	for (int iY=NY-1; iY>=0; --iY)
	{
		for (int iX=0; iX<NX; ++iX)
		{
			fout<<f[iX][iY][iPop]<<" ";
            gout<<g[iX][iY][iPop]<<" ";
		}
		fout<<"\n";
		gout<<"\n";
	}

}



void init()
{

	//Phase initialization prior any equilibrium functions calculations
    for(int iX=0;iX<NX;iX++)
		for(int iY=0; iY<NY; iY++)
		{
			if ( (iX>=(NX-1)/3) && (iX<=2*(NX-1)/3) && (iY>=width) && (iY<=NY-width-1) )
            {
                phase[iX][iY]=-1.0;
            }
            else
            {
                phase[iX][iY]=1.0;
            }

			if (iX==0)
			{
				phase[iX][iY]=phase_inlet;
			}

			if (iX==NX-1)
			{
				phase[iX][iY]=phase_outlet;
			}

		}

    //Wall gradients stuff needs to be calculated outside the loop
    for(int iX=0;iX<NX;iX++)
    {
		//First-order accuracy
		phase[iX][0]=phase[iX][1]-wall_gradient;
        //Second-order accuracy
        //phase[iX][iY]=(4.0*phase[iX][iY+1]-phase[iX][iY+2]-2.0*wall_gradient)/3.0;

		//First-order accuracy
		phase[iX][NY-1]=phase[iX][NY-2]-wall_gradient;
		//Second-order accuracy
        //phase[iX][iY]=(4.0*phase[iX][iY-1]-phase[iX][iY-2]-2.0*wall_gradient)/3.0;
    }


	//Bulk nodes initialization
	for(int iX=1;iX<NX-1;iX++)
		for(int iY=1;iY<NY-1;iY++)
		{


			double laplace_temp=0.0;
			double gradx_temp=0.0;
			double grady_temp=0.0;
			for(int k=0;k<9;k++)
			{
				int iX2=(iX+cx[k]+NX) % NX;
				int iY2=(iY+cy[k]+NY) % NY;
				laplace_temp+=laplacestencil[k]*phase[iX2][iY2];
				gradx_temp+=gradstencilx[k]*phase[iX2][iY2];
				grady_temp+=gradstencily[k]*phase[iX2][iY2];
			}

            laplace[iX][iY]=laplace_temp;
            gradx[iX][iY]=gradx_temp;
            grady[iX][iY]=grady_temp;

			//Initialization of the macroscopic fields
			rho[iX][iY]=1.0;
			ux[iX][iY]=0.0;
			uy[iX][iY]=0.0;


			double phase_temp=phase[iX][iY];
			double dense_temp=rho[iX][iY];
			double ux_temp=ux[iX][iY];
			double uy_temp=uy[iX][iY];

			double feq;
			double geq;
			double sum=0.0;
			double sum_phase=0.0;
			double phase_square=phase_temp*phase_temp;
			double pressure_bulk=dense_temp/3.0+aconst*(-0.5*phase_square+3.0/4.0*phase_square*phase_square)-kconst*phase_temp*laplace_temp;
			double chemical_pot=gammaconst*(aconst*(-phase_temp+phase_temp*phase_temp*phase_temp)-kconst*laplace_temp);

			for (int k=1; k<9; k++)
			{
				feq=weights[k]*(3.0*pressure_bulk+3.0*dense_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
								+4.5*dense_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+2.0*ux_temp*uy_temp*cx[k]*cy[k]))
					+kconst*(wxx[k]*gradx_temp*gradx_temp+wyy[k]*grady_temp*grady_temp+wxy[k]*gradx_temp*grady_temp);
				geq=weights[k]*(3.0*chemical_pot+3.0*phase_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
								+4.5*phase_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp
												 +2.0*ux_temp*uy_temp*cx[k]*cy[k]));
				sum+=feq;
				sum_phase+=geq;

				f[iX][iY][k]=feq;
				g[iX][iY][k]=geq;
			}
			f[iX][iY][0]=dense_temp-sum;
			g[iX][iY][0]=phase_temp-sum_phase;

		}


	//Wall nodes initialization
    for(int iX=0;iX<NX;iX++)
	{
		rho[iX][0]=1.0;

		for(int iPop=0;iPop<9;iPop++)
		{
			f[iX][0][iPop]=weights[iPop]*rho[iX][0];
			g[iX][0][iPop]=weights[iPop]*phase[iX][0];

		}

		rho[iX][NY-1]=1.0;

		for(int iPop=0;iPop<9;iPop++)
		{
			f[iX][NY-1][iPop]=weights[iPop]*rho[iX][NY-1];
			g[iX][NY-1][iPop]=weights[iPop]*phase[iX][NY-1];
		}
	}
    //Inlet and outlet nodes initialization
    //TODO:: Do it according BC conditions
    for(int iY=1;iY<NY-1;iY++)
	{
		rho[0][iY]=rho_inlet;
        rho[NX-1][iY]=rho_outlet;

		for(int iPop=0;iPop<9;iPop++)
		{
			f[0][iY][iPop]=weights[iPop]*rho[0][iY];
			g[0][iY][iPop]=weights[iPop]*phase[0][iY];

			f[NX-1][iY][iPop]=weights[iPop]*rho[NX-1][iY];
			g[NX-1][iY][iPop]=weights[iPop]*phase[NX-1][iY];
		}

	}



}


void collide_bulk()
{
    //The phase field should be calculated prior the laplacians
    //Calculation of phase for bulk nodes
    for(int iX=1;iX<NX-1;iX++)
        for(int iY=1;iY<NY-1;iY++)
		{
            phase[iX][iY]=0.0;
            for(int iPop=0;iPop<9;iPop++)
   				phase[iX][iY]+=g[iX][iY][iPop];
		}

    //Update phase field for inlet and outlet
    for(int iY=1;iY<NY-1;iY++)
    {
        phase[0][iY]=phase_inlet;
        phase[NX-1][iY]=phase_outlet;
    }


    //Calculation of the wall phases because of the wall gradients
    for(int iX=0;iX<NX;iX++)
    {
        //First-order accuracy
		phase[iX][0]=phase[iX][1]-wall_gradient;
        //Second-order accuracy
        //phase[iX][0]=(4.0*phase[iX][1]-phase[iX][2]-2.0*wall_gradient)/3.0;

		//First-order accuracy
		phase[iX][NY-1]=phase[iX][NY-2]-wall_gradient;
		//Second-order accuracy
        //phase[iX][NY-1]=(4.0*phase[iX][NY-2]-phase[iX][NY-3]-2.0*wall_gradient)/3.0;
    }



    for(int iX=1;iX<NX-1;iX++)
        for(int iY=1;iY<NY-1;iY++)
		{

			//Construction equilibrium
			rho[iX][iY]=0.0;
			ux[iX][iY]=0.0;
			uy[iX][iY]=0.0;

			for(int iPop=0;iPop<9;iPop++)
			{
				rho[iX][iY]+=f[iX][iY][iPop];
				ux[iX][iY]+=f[iX][iY][iPop]*cx[iPop];
				uy[iX][iY]+=f[iX][iY][iPop]*cy[iPop];
			}

			ux[iX][iY]=(ux[iX][iY]+0.5*force_x)/rho[iX][iY];
			uy[iX][iY]=(uy[iX][iY]+0.5*force_y)/rho[iX][iY];


			double laplace_temp=0.0;
			double gradx_temp=0.0;
			double grady_temp=0.0;
			for(int k=0;k<9;k++)
			{
				int iX2=(iX+cx[k]+NX) % NX;
				int iY2=(iY+cy[k]+NY) % NY;
				laplace_temp+=laplacestencil[k]*phase[iX2][iY2];
				gradx_temp+=gradstencilx[k]*phase[iX2][iY2];
				grady_temp+=gradstencily[k]*phase[iX2][iY2];
			}

            laplace[iX][iY]=laplace_temp;
            gradx[iX][iY]=gradx_temp;
            grady[iX][iY]=grady_temp;

			double phase_temp=phase[iX][iY];
			double dense_temp=rho[iX][iY];
			double ux_temp=ux[iX][iY];
			double uy_temp=uy[iX][iY];

			double sum=0.0;
			double sum_phase=0.0;
			double phase_square=phase_temp*phase_temp;
			double pressure_bulk=dense_temp/3.0+aconst*(-0.5*phase_square+3.0/4.0*phase_square*phase_square)-kconst*phase_temp*laplace_temp;
			double chemical_pot=gammaconst*(aconst*(-phase_temp+phase_temp*phase_temp*phase_temp)-kconst*laplace_temp);

			double force[9];

			//Obtain force population
            for (int k=0;k<9;k++)
            {
				force[k]=weights[k]*(1.0-0.5*omega)*(3.0*force_x*cx[k]+3.0*force_y*cy[k]+
                            9.0*((cx[k]*cx[k]-1.0/3.0)*force_x*ux_temp+cx[k]*cy[k]*(force_x*uy_temp+force_y*ux_temp)+
 							(cy[k]*cy[k]-1.0/3.0)*force_y*uy_temp));
            }
			for (int k=1; k<9; k++)
			{
				feq[iX][iY][k]=weights[k]*(3.0*pressure_bulk+3.0*dense_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
								+4.5*dense_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+2.0*ux_temp*uy_temp*cx[k]*cy[k]))
				+kconst*(wxx[k]*gradx_temp*gradx_temp+wyy[k]*grady_temp*grady_temp+wxy[k]*gradx_temp*grady_temp);
				geq[iX][iY][k]=weights[k]*(3.0*chemical_pot+3.0*phase_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
								+4.5*phase_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+2.0*ux_temp*uy_temp*cx[k]*cy[k]));
				sum+=feq[iX][iY][k];
				sum_phase+=geq[iX][iY][k];
			}

			feq[iX][iY][0]=dense_temp-sum;
			geq[iX][iY][0]=phase_temp-sum_phase;

            double tau_temp=tau_gas+(phase_temp+1.0)/2.0*(tau_liq-tau_gas);
            double omega_temp=1.0/tau_temp;

			for(int k=0; k < 9; k++)
			{
				f2[iX][iY][k]=f[iX][iY][k]*(1.0-omega_temp)+omega_temp*feq[iX][iY][k]+force[k];
				g2[iX][iY][k]=g[iX][iY][k]*(1.0-omega_temp)+omega_temp*geq[iX][iY][k];
			}


		}

    //Specify the densities and velocities for inlet and outlet
    for(int iY=1;iY<NY-1;iY++)
    {
        rho[0][iY]=rho_inlet;
        rho[NX-1][iY]=rho_outlet;
        ux[0][iY]=ux[1][iY];
        uy[0][iY]=uy[1][iY];
        ux[NX-1][iY]=ux[NX-2][iY];
        uy[NX-1][iY]=uy[NX-2][iY];
    }




}


void guo_binary_construction(int iY)
{
    //Construct equilibrium by using previously defined laplacians and gradients
    //Inlet part
    double phase_temp=phase[0][iY];
	double dense_temp=rho[0][iY];
	double ux_temp=ux[1][iY];
	double uy_temp=uy[1][iY];


	double phase_square=phase_temp*phase_temp;
	double laplace_temp=laplace[1][iY];
	double gradx_temp=gradx[1][iY];
	double grady_temp=grady[1][iY];

//    if ((iY==1)||(iY==NY-2))
//    {
//        std::cout.precision(10);
//        std::cout<<gradx_temp<<" "<<grady_temp<<" "<<laplace_temp<<"\n";
//        std::cout<<"Velocities"<<ux_temp<<" "<<uy_temp<<"\n";
//    }


	double pressure_bulk=dense_temp/3.0+aconst*(-0.5*phase_square+3.0/4.0*phase_square*phase_square)-kconst*phase_temp*laplace_temp;
	double chemical_pot=gammaconst*(aconst*(-phase_temp+phase_temp*phase_temp*phase_temp)-kconst*laplace_temp);

    double long_rho_term[9],long_phase_term[9];
    for(int k=1;k<9;k++)
    {
        long_rho_term[k]=3.0*pressure_bulk+3.0*dense_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
								+4.5*dense_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+2.0*ux_temp*uy_temp*cx[k]*cy[k]);
        long_phase_term[k]=3.0*chemical_pot+3.0*phase_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
								+4.5*phase_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp
												 +2.0*ux_temp*uy_temp*cx[k]*cy[k]);
    }
	feq[0][iY][1]=weights[1]*long_rho_term[1]
			+kconst*(wxx[1]*gradx_temp*gradx_temp+wyy[1]*grady_temp*grady_temp+wxy[1]*gradx_temp*grady_temp);
    feq[0][iY][5]=weights[5]*long_rho_term[5]
			+kconst*(wxx[5]*gradx_temp*gradx_temp+wyy[5]*grady_temp*grady_temp+wxy[5]*gradx_temp*grady_temp);
	feq[0][iY][8]=weights[8]*long_rho_term[8]
			+kconst*(wxx[8]*gradx_temp*gradx_temp+wyy[8]*grady_temp*grady_temp+wxy[8]*gradx_temp*grady_temp);

    geq[0][iY][1]=weights[1]*long_phase_term[1];
    geq[0][iY][5]=weights[5]*long_phase_term[5];
    geq[0][iY][8]=weights[8]*long_phase_term[8];

	double tau_temp=tau_gas+(phase_temp+1.0)/2.0*(tau_liq-tau_gas);
	double omega_temp=1.0/tau_temp;


	f2[0][iY][1]=feq[0][iY][1]+(1.0-omega_temp)*(f[1][iY][1]-feq[1][iY][1]);
	f2[0][iY][5]=feq[0][iY][5]+(1.0-omega_temp)*(f[1][iY][5]-feq[1][iY][5]);
	f2[0][iY][8]=feq[0][iY][8]+(1.0-omega_temp)*(f[1][iY][8]-feq[1][iY][8]);

	g2[0][iY][1]=geq[0][iY][1]+(1.0-omega_temp)*(g[1][iY][1]-geq[1][iY][1]);
	g2[0][iY][5]=geq[0][iY][5]+(1.0-omega_temp)*(g[1][iY][5]-geq[1][iY][5]);
	g2[0][iY][8]=geq[0][iY][8]+(1.0-omega_temp)*(g[1][iY][8]-geq[1][iY][8]);

    //Outlet part
    phase_temp=phase[NX-1][iY];
	dense_temp=rho[NX-1][iY];
	ux_temp=ux[NX-2][iY];
	uy_temp=uy[NX-2][iY];


	phase_square=phase_temp*phase_temp;
	laplace_temp=laplace[NX-2][iY];
	gradx_temp=gradx[NX-2][iY];
	grady_temp=grady[NX-2][iY];

	pressure_bulk=dense_temp/3.0+aconst*(-0.5*phase_square+3.0/4.0*phase_square*phase_square)-kconst*phase_temp*laplace_temp;
	chemical_pot=gammaconst*(aconst*(-phase_temp+phase_temp*phase_temp*phase_temp)-kconst*laplace_temp);

    for(int k=1;k<9;k++)
    {
        long_rho_term[k]=3.0*pressure_bulk+3.0*dense_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
								+4.5*dense_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp+2.0*ux_temp*uy_temp*cx[k]*cy[k]);
        long_phase_term[k]=3.0*chemical_pot+3.0*phase_temp*(cx[k]*ux_temp+cy[k]*uy_temp)
								+4.5*phase_temp*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp+(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp
												 +2.0*ux_temp*uy_temp*cx[k]*cy[k]);
    }

	feq[NX-1][iY][3]=weights[3]*long_rho_term[3]
			+kconst*(wxx[3]*gradx_temp*gradx_temp+wyy[3]*grady_temp*grady_temp+wxy[3]*gradx_temp*grady_temp);
    feq[NX-1][iY][6]=weights[6]*long_rho_term[6]
			+kconst*(wxx[6]*gradx_temp*gradx_temp+wyy[6]*grady_temp*grady_temp+wxy[6]*gradx_temp*grady_temp);
	feq[NX-1][iY][7]=weights[7]*long_rho_term[7]
			+kconst*(wxx[7]*gradx_temp*gradx_temp+wyy[7]*grady_temp*grady_temp+wxy[7]*gradx_temp*grady_temp);

    geq[NX-1][iY][3]=weights[3]*long_phase_term[3];
    geq[NX-1][iY][6]=weights[6]*long_phase_term[6];
    geq[NX-1][iY][7]=weights[7]*long_phase_term[7];

	tau_temp=tau_gas+(phase_temp+1.0)/2.0*(tau_liq-tau_gas);
	omega_temp=1.0/tau_temp;


	f2[NX-1][iY][3]=feq[NX-1][iY][3]+(1.0-omega_temp)*(f[NX-2][iY][3]-feq[NX-2][iY][3]);
	f2[NX-1][iY][6]=feq[NX-1][iY][6]+(1.0-omega_temp)*(f[NX-2][iY][6]-feq[NX-2][iY][6]);
	f2[NX-1][iY][7]=feq[NX-1][iY][7]+(1.0-omega_temp)*(f[NX-2][iY][7]-feq[NX-2][iY][7]);

	g2[NX-1][iY][3]=geq[NX-1][iY][3]+(1.0-omega_temp)*(g[NX-2][iY][3]-geq[NX-2][iY][3]);
	g2[NX-1][iY][6]=geq[NX-1][iY][6]+(1.0-omega_temp)*(g[NX-2][iY][6]-geq[NX-2][iY][6]);
	g2[NX-1][iY][7]=geq[NX-1][iY][7]+(1.0-omega_temp)*(g[NX-2][iY][7]-geq[NX-2][iY][7]);

}

void guo_pressure()
{
    for(int iY=1;iY<NY-1;iY++)
    {
        guo_binary_construction(iY);
    }


}



void update_bounce_back()
{
	//BB nodes density and velocity specification
	for(int iX=0;iX<NX;iX++)
	{
		int iXtop=(iX+1+NX)%NX;
		int iXbottom=(iX-1+NX)%NX;

		f2[iX][0][2]=f2[iX][1][4];
		f2[iX][0][5]=f2[iXtop][1][7];
		f2[iX][0][6]=f2[iXbottom][1][8];

		f2[iX][NY-1][4]=f2[iX][NY-2][2];
		f2[iX][NY-1][7]=f2[iXbottom][NY-2][5];
		f2[iX][NY-1][8]=f2[iXtop][NY-2][6];

		//BB for the scalar phase??? Ask Irina about it
		g2[iX][0][2]=g2[iX][1][4];
		g2[iX][0][5]=g2[iXtop][1][7];
		g2[iX][0][6]=g2[iXbottom][1][8];

		g2[iX][NY-1][4]=g2[iX][NY-2][2];
		g2[iX][NY-1][7]=g2[iXbottom][NY-2][5];
		g2[iX][NY-1][8]=g2[iXtop][NY-2][6];


		rho[iX][0]=1.0;
		rho[iX][NY-1]=1.0;
		ux[iX][0]=0.0;ux[iX][NY-1]=0.0;
		uy[iX][0]=0.0;uy[iX][NY-1]=0.0;
	}
}


int main(int argc, char* argv[])
{

    if (argc!=1)
        width=atoi(argv[1]);

    std::cout<<"Width="<<width<<"\n";

    matrix_init();
    init();

	for(int counter=0;counter<=N;counter++)
	{

        collide_bulk();
        guo_pressure();

        {
  			std::stringstream filewritepopulations;
            filewritepopulations<<std::fixed;
 			std::stringstream counterconvert;
 			counterconvert<<counter;

            filewritepopulations<<"pop"<<std::string(6-counterconvert.str().size(),'0')<<counter;

            //for (int iPop=0;iPop<9;iPop++)
            //    writepopulationslice(filewritepopulations.str(),iPop);
        }

        update_bounce_back();

		//Streaming
		for(int iX=1;iX<NX-1;iX++)
			for(int iY=1;iY<NY-1;iY++)
				for(int iPop=0;iPop<9;iPop++)
				{
					int iX2=(iX-cx[iPop]+NX)%NX;
					int iY2=(iY-cy[iPop]+NY)%NY;
					f[iX][iY][iPop]=f2[iX2][iY2][iPop];
                    g[iX][iY][iPop]=g2[iX2][iY2][iPop];
				}

		//Writing files
		if (counter%NOUTPUT==0)
		{
			std::cout<<counter<<"\n";

			std::stringstream filewritedensity;
 			std::stringstream filewritevelocity;
 			std::stringstream filewritephase;
 			std::stringstream counterconvert;
 			counterconvert<<counter;

 			filewritedensity<<std::fixed;
			filewritevelocity<<std::fixed;
			filewritephase<<std::fixed;

			filewritedensity<<"density"<<std::string(6-counterconvert.str().size(),'0')<<counter;
			filewritevelocity<<"velocity"<<std::string(6-counterconvert.str().size(),'0')<<counter;
            filewritephase<<"phase"<<std::string(6-counterconvert.str().size(),'0')<<counter;

 			writedensity(filewritedensity.str());
			writevelocity(filewritevelocity.str());
            writephase(filewritephase.str());
		}


	}

	return 0;
}
