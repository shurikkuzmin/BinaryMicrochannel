//This file is intended to try to use Zou-He pressure BC for both fields.


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>

//Domain size
const int NY=31;
const int NX=201;

//Time steps
const int N=1000;
const int NOUTPUT=100;

//Fields and populations
double f[NX][NY][9], f2[NX][NY][9], feq[NX][NY][9], g[NX][NY][9], g2[NX][NY][9],geq[NX][NY][9];
double rho[NX][NY],ux[NX][NY],uy[NX][NY],phase[NX][NY];

//Laplacians and gradients
double laplace[NX][NY],gradx[NX][NY],grady[NX][NY];

//Pressure boundary conditions
double rho_inlet=1.003;
double rho_outlet=1.0;
double phase_inlet=1.0;
double phase_outlet=1.0;

double force_x=0.0001;
double force_y=0.000;

//Binary-liquid parameters
double aconst=0.04;
double kconst=0.04;
double gammaconst=1.0;

//Wall wettability parameter
double wall_gradient=0.5;

//BGK relaxation parameter
double omega=1.0;

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

void init()
{

	//Phase initialization prior any equilibrium functions calculations
    for(int iX=0;iX<NX;iX++)
		for(int iY=0; iY<NY; iY++)
		{
			if ( (iX>=(NX-1)/4) && (iX<=3*(NX-1)/4) && (iY>=5) && (iY<=NY-6) )
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

			if (iY==0)
			{
				//First-order accuracy
				phase[iX][iY]=phase[iX][iY+1]-wall_gradient;
                //Second-order accuracy
                //phase[iX][iY]=(4.0*phase[iX][iY+1]-phase[iX][iY+2]-2.0*wall_gradient)/3.0;
			}

			if (iY==NY-1)
			{
			    //First-order accuracy
			    phase[iX][iY]=phase[iX][iY-1]-wall_gradient;
			    //Second-order accuracy
                //phase[iX][iY]=(4.0*phase[iX][iY-1]-phase[iX][iY-2]-2.0*wall_gradient)/3.0;
            }
		}

	//Bulk nodes initialization
	for(int iX=0;iX<NX;iX++)
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

}


void collide_bulk()
{
    //The phase field should be calculated prior the laplacians
    for(int iX=0;iX<NX;iX++)
        for(int iY=1;iY<NY-1;iY++)
		{
            phase[iX][iY]=0.0;
            for(int iPop=0;iPop<9;iPop++)
   				phase[iX][iY]+=g[iX][iY][iPop];
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



    for(int iX=0;iX<NX;iX++)
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


			for(int k=0; k < 9; k++)
			{
				f2[iX][iY][k]=f[iX][iY][k]*(1.0-omega)+omega*feq[iX][iY][k]+force[k];
				g2[iX][iY][k]=g[iX][iY][k]*(1.0-omega)+omega*geq[iX][iY][k];
			}


		}

}

void guo_binary_construction(int iX)
{
    //Construct equilibrium by using previously defined laplacians and gradients
    //Bottom part
    double phase_temp=phase[iX][0];
	double dense_temp=rho[iX][0];
	double phase_square=phase_temp*phase_temp;
	double laplace_temp=laplace[iX][1];
	double gradx_temp=gradx[iX][1];
	double grady_temp=grady[iX][1];

	double pressure_bulk=dense_temp/3.0+aconst*(-0.5*phase_square+3.0/4.0*phase_square*phase_square)-kconst*phase_temp*laplace_temp;
	double chemical_pot=gammaconst*(aconst*(-phase_temp+phase_temp*phase_temp*phase_temp)-kconst*laplace_temp);

	feq[iX][0][2]=3.0*weights[2]*pressure_bulk
			+kconst*(wxx[2]*gradx_temp*gradx_temp+wyy[2]*grady_temp*grady_temp+wxy[2]*gradx_temp*grady_temp);
    feq[iX][0][5]=3.0*weights[5]*pressure_bulk
			+kconst*(wxx[5]*gradx_temp*gradx_temp+wyy[5]*grady_temp*grady_temp+wxy[5]*gradx_temp*grady_temp);
	feq[iX][0][6]=3.0*weights[6]*pressure_bulk
			+kconst*(wxx[6]*gradx_temp*gradx_temp+wyy[6]*grady_temp*grady_temp+wxy[6]*gradx_temp*grady_temp);

    geq[iX][0][2]=3.0*weights[2]*chemical_pot;
    geq[iX][0][5]=3.0*weights[5]*chemical_pot;
    geq[iX][0][6]=3.0*weights[6]*chemical_pot;

	f2[iX][0][2]=feq[iX][0][2]+(1.0-omega)*(f[iX][1][2]-feq[iX][1][2]);
	f2[iX][0][5]=feq[iX][0][5]+(1.0-omega)*(f[iX][1][5]-feq[iX][1][5]);
	f2[iX][0][6]=feq[iX][0][6]+(1.0-omega)*(f[iX][1][6]-feq[iX][1][6]);

	g2[iX][0][2]=geq[iX][0][2]+(1.0-omega)*(g[iX][1][2]-geq[iX][1][2]);
	g2[iX][0][5]=geq[iX][0][5]+(1.0-omega)*(g[iX][1][5]-geq[iX][1][5]);
	g2[iX][0][6]=geq[iX][0][6]+(1.0-omega)*(g[iX][1][6]-geq[iX][1][6]);

    //Top part
    phase_temp=phase[iX][NY-1];
    dense_temp=rho[iX][NY-1];
	phase_square=phase_temp*phase_temp;
	laplace_temp=laplace[iX][NY-2];
	gradx_temp=gradx[iX][NY-2];
	grady_temp=grady[iX][NY-2];

	pressure_bulk=dense_temp/3.0+aconst*(-0.5*phase_square+3.0/4.0*phase_square*phase_square)-kconst*phase_temp*laplace_temp;
	chemical_pot=gammaconst*(aconst*(-phase_temp+phase_temp*phase_temp*phase_temp)-kconst*laplace_temp);

	feq[iX][NY-1][4]=3.0*weights[4]*pressure_bulk
            +kconst*(wxx[4]*gradx_temp*gradx_temp+wyy[4]*grady_temp*grady_temp+wxy[4]*gradx_temp*grady_temp);
    feq[iX][NY-1][7]=3.0*weights[7]*pressure_bulk
			+kconst*(wxx[7]*gradx_temp*gradx_temp+wyy[7]*grady_temp*grady_temp+wxy[7]*gradx_temp*grady_temp);
	feq[iX][NY-1][8]=3.0*weights[8]*pressure_bulk
			+kconst*(wxx[8]*gradx_temp*gradx_temp+wyy[8]*grady_temp*grady_temp+wxy[8]*gradx_temp*grady_temp);

    geq[iX][NY-1][4]=3.0*weights[4]*chemical_pot;
    geq[iX][NY-1][7]=3.0*weights[7]*chemical_pot;
    geq[iX][NY-1][8]=3.0*weights[8]*chemical_pot;

	f2[iX][NY-1][4]=feq[iX][NY-1][4]+(1.0-omega)*(f[iX][NY-2][4]-feq[iX][NY-2][4]);
	f2[iX][NY-1][7]=feq[iX][NY-1][7]+(1.0-omega)*(f[iX][NY-2][7]-feq[iX][NY-2][7]);
	f2[iX][NY-1][8]=feq[iX][NY-1][8]+(1.0-omega)*(f[iX][NY-2][8]-feq[iX][NY-2][8]);

	g2[iX][NY-1][4]=geq[iX][NY-1][4]+(1.0-omega)*(g[iX][NY-2][4]-geq[iX][NY-2][4]);
	g2[iX][NY-1][7]=geq[iX][NY-1][7]+(1.0-omega)*(g[iX][NY-2][7]-geq[iX][NY-2][7]);
	g2[iX][NY-1][8]=geq[iX][NY-1][8]+(1.0-omega)*(g[iX][NY-2][8]-geq[iX][NY-2][8]);


}

void guo_hydrodynamics_mixed_construction(int iX)
{
    //Construct equilibrium by using previously defined laplacians and gradients
    //Bottom part
    double phase_temp=phase[iX][0];
	double dense_temp=rho[iX][0];

	double pressure_bulk=dense_temp/3.0;
	double chemical_pot=phase_temp/3.0;

	feq[iX][0][2]=3.0*weights[2]*pressure_bulk;
    feq[iX][0][5]=3.0*weights[5]*pressure_bulk;
	feq[iX][0][6]=3.0*weights[6]*pressure_bulk;

    geq[iX][0][2]=3.0*weights[2]*chemical_pot;
    geq[iX][0][5]=3.0*weights[5]*chemical_pot;
    geq[iX][0][6]=3.0*weights[6]*chemical_pot;

	f2[iX][0][2]=feq[iX][0][2]+(1.0-omega)*(f[iX][1][2]-feq[iX][1][2]);
	f2[iX][0][5]=feq[iX][0][5]+(1.0-omega)*(f[iX][1][5]-feq[iX][1][5]);
	f2[iX][0][6]=feq[iX][0][6]+(1.0-omega)*(f[iX][1][6]-feq[iX][1][6]);

	g2[iX][0][2]=geq[iX][0][2]+(1.0-omega)*(g[iX][1][2]-geq[iX][1][2]);
	g2[iX][0][5]=geq[iX][0][5]+(1.0-omega)*(g[iX][1][5]-geq[iX][1][5]);
	g2[iX][0][6]=geq[iX][0][6]+(1.0-omega)*(g[iX][1][6]-geq[iX][1][6]);

    //Top part
    phase_temp=phase[iX][NY-1];
    dense_temp=rho[iX][NY-1];

	pressure_bulk=dense_temp/3.0;
	chemical_pot=phase_temp/3.0;

	feq[iX][NY-1][4]=3.0*weights[4]*pressure_bulk;
    feq[iX][NY-1][7]=3.0*weights[7]*pressure_bulk;
	feq[iX][NY-1][8]=3.0*weights[8]*pressure_bulk;

    geq[iX][NY-1][4]=3.0*weights[4]*chemical_pot;
    geq[iX][NY-1][7]=3.0*weights[7]*chemical_pot;
    geq[iX][NY-1][8]=3.0*weights[8]*chemical_pot;

	f2[iX][NY-1][4]=feq[iX][NY-1][4]+(1.0-omega)*(f[iX][NY-2][4]-feq[iX][NY-2][4]);
	f2[iX][NY-1][7]=feq[iX][NY-1][7]+(1.0-omega)*(f[iX][NY-2][7]-feq[iX][NY-2][7]);
	f2[iX][NY-1][8]=feq[iX][NY-1][8]+(1.0-omega)*(f[iX][NY-2][8]-feq[iX][NY-2][8]);

	g2[iX][NY-1][4]=geq[iX][NY-1][4]+(1.0-omega)*(g[iX][NY-2][4]-geq[iX][NY-2][4]);
	g2[iX][NY-1][7]=geq[iX][NY-1][7]+(1.0-omega)*(g[iX][NY-2][7]-geq[iX][NY-2][7]);
	g2[iX][NY-1][8]=geq[iX][NY-1][8]+(1.0-omega)*(g[iX][NY-2][8]-geq[iX][NY-2][8]);


}

void guo_hydrodynamics_construction(int iX)
{
    //Construct equilibrium by using previously defined laplacians and gradients
    //Bottom part
    double phase_temp=phase[iX][0];
	double dense_temp=rho[iX][0];
    double ux_temp=ux[iX][1];
    double uy_temp=uy[iX][1];


	double pressure_bulk=dense_temp/3.0;
	double chemical_pot=phase_temp/3.0;
    double feqeq[9],geqeq[9];

    for(int k=0;k<9;k++)
    {
        feqeq[k]=weights[k]*dense_temp*(1.0+3.0*(cx[k]*ux_temp+cy[k]*uy_temp)
                                        +4.5*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp
                                               +2.0*cx[k]*cy[k]*ux_temp*uy_temp
                                               +(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp));
        geqeq[k]=weights[k]*phase_temp*(1.0+3.0*(cx[k]*ux_temp+cy[k]*uy_temp)
                                        +4.5*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp
                                               +2.0*cx[k]*cy[k]*ux_temp*uy_temp
                                               +(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp));
    }

	feq[iX][0][2]=3.0*weights[2]*pressure_bulk;
    feq[iX][0][5]=3.0*weights[5]*pressure_bulk;
	feq[iX][0][6]=3.0*weights[6]*pressure_bulk;

    geq[iX][0][2]=3.0*weights[2]*chemical_pot;
    geq[iX][0][5]=3.0*weights[5]*chemical_pot;
    geq[iX][0][6]=3.0*weights[6]*chemical_pot;

	f2[iX][0][2]=feq[iX][0][2]+(1.0-omega)*(f[iX][1][2]-feqeq[2]);
	f2[iX][0][5]=feq[iX][0][5]+(1.0-omega)*(f[iX][1][5]-feqeq[5]);
	f2[iX][0][6]=feq[iX][0][6]+(1.0-omega)*(f[iX][1][6]-feqeq[6]);

	g2[iX][0][2]=geq[iX][0][2]+(1.0-omega)*(g[iX][1][2]-geqeq[2]);
	g2[iX][0][5]=geq[iX][0][5]+(1.0-omega)*(g[iX][1][5]-geqeq[5]);
	g2[iX][0][6]=geq[iX][0][6]+(1.0-omega)*(g[iX][1][6]-geqeq[6]);

    //Top part
    phase_temp=phase[iX][NY-1];
    dense_temp=rho[iX][NY-1];
    ux_temp=ux[iX][NY-2];
    uy_temp=uy[iX][NY-2];

	pressure_bulk=dense_temp/3.0;
	chemical_pot=phase_temp/3.0;

    for(int k=0;k<9;k++)
    {
        feqeq[k]=weights[k]*dense_temp*(1.0+3.0*(cx[k]*ux_temp+cy[k]*uy_temp)
                                        +4.5*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp
                                               +2.0*cx[k]*cy[k]*ux_temp*uy_temp
                                               +(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp));
        geqeq[k]=weights[k]*phase_temp*(1.0+3.0*(cx[k]*ux_temp+cy[k]*uy_temp)
                                        +4.5*((cx[k]*cx[k]-1.0/3.0)*ux_temp*ux_temp
                                               +2.0*cx[k]*cy[k]*ux_temp*uy_temp
                                               +(cy[k]*cy[k]-1.0/3.0)*uy_temp*uy_temp));
    }


	feq[iX][NY-1][4]=3.0*weights[4]*pressure_bulk;
    feq[iX][NY-1][7]=3.0*weights[7]*pressure_bulk;
	feq[iX][NY-1][8]=3.0*weights[8]*pressure_bulk;

    geq[iX][NY-1][4]=3.0*weights[4]*chemical_pot;
    geq[iX][NY-1][7]=3.0*weights[7]*chemical_pot;
    geq[iX][NY-1][8]=3.0*weights[8]*chemical_pot;

	f2[iX][NY-1][4]=feq[iX][NY-1][4]+(1.0-omega)*(f[iX][NY-2][4]-feqeq[4]);
	f2[iX][NY-1][7]=feq[iX][NY-1][7]+(1.0-omega)*(f[iX][NY-2][7]-feqeq[7]);
	f2[iX][NY-1][8]=feq[iX][NY-1][8]+(1.0-omega)*(f[iX][NY-2][8]-feqeq[8]);

	g2[iX][NY-1][4]=geq[iX][NY-1][4]+(1.0-omega)*(g[iX][NY-2][4]-geqeq[4]);
	g2[iX][NY-1][7]=geq[iX][NY-1][7]+(1.0-omega)*(g[iX][NY-2][7]-geqeq[7]);
	g2[iX][NY-1][8]=geq[iX][NY-1][8]+(1.0-omega)*(g[iX][NY-2][8]-geqeq[8]);

}


void update_guo_walls()
{
   	for(int iX=0;iX<NX;iX++)
	{
		rho[iX][0]=rho[iX][1];
		rho[iX][NY-1]=rho[iX][NY-2];

		ux[iX][0]=0.0;
		uy[iX][0]=0.0;
		ux[iX][NY-1]=0.0;
		uy[iX][NY-1]=0.0;

        guo_binary_construction(iX);
        //guo_hydrodynamics_mixed_construction(iX);
        //guo_hydrodynamics_construction(iX);
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
		//phase[iX][0]=0.5;
		//phase[iX][NY-1]=0.5;
		ux[iX][0]=0.0;ux[iX][NY-1]=0.0;
		uy[iX][0]=0.0;uy[iX][NY-1]=0.0;
	}
}

void update_half_bounce_back()
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

		//Equilibrium values for the phase field
		g2[iX][0][2]=weights[2]*phase[iX][0];
		g2[iX][0][5]=weights[5]*phase[iX][0];
		g2[iX][0][6]=weights[6]*phase[iX][0];

		g2[iX][NY-1][4]=weights[4]*phase[iX][NY-1];
		g2[iX][NY-1][7]=weights[7]*phase[iX][NY-1];
		g2[iX][NY-1][8]=weights[8]*phase[iX][NY-1];


		rho[iX][0]=1.0;
		rho[iX][NY-1]=1.0;
		ux[iX][0]=0.0;ux[iX][NY-1]=0.0;
		uy[iX][0]=0.0;uy[iX][NY-1]=0.0;
	}
}



int main(int argc, char* argv[])
{

    matrix_init();
    init();


	for(int counter=0;counter<=N;counter++)
	{

        collide_bulk();
        //update_half_bounce_back();
        update_bounce_back();
        //update_guo_walls();

		//Streaming
		for(int iX=0;iX<NX;iX++)
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
