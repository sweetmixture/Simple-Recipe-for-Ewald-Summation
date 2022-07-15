#include <iostream>
#include <iomanip>
#include <cmath>

using std::cout;
using std::endl;

#define TO_EV 14.39964390675221758120

double dp( double* r1, double* r2 )
{
	return r1[0]*r2[0] + r1[1]*r2[1] + r1[2]*r2[2];
}

double norm( double* r )
{
	return sqrt(dp(r,r));
}
double norm( double x, double y, double z )
{
	return sqrt(x*x+y*y+z*z);
}

void GetCart( const double Lx, const double Ly, const double Lz,
	const double* frac_x, const double* frac_y, const double* frac_z,
	const int index, double* out )
{
	out[0] = Lx*frac_x[index];
	out[1] = Ly*frac_y[index];
	out[2] = Lz*frac_z[index];
}

int main()
{
	// UnitCell Info
	double Lx,Ly,Lz;
	double V;
	Lx = Ly = Lz = 4.2112;		// In Angs		// Lz = 6;
	Lx = 3.89;
	Ly = 4.3;
	Lz = 3.2;
	V = Lx*Ly*Lz;			// In Angs**3

	double Real_E, Reci_E;
	Real_E = Reci_E = 0.;		// In eV

	// Ion Info
	const int NumberOfIons = 8;

	double frac_x[NumberOfIons], frac_y[NumberOfIons], frac_z[NumberOfIons];
	double cart_x[NumberOfIons], cart_y[NumberOfIons], cart_z[NumberOfIons];
	double charge[NumberOfIons];

	double force_x[NumberOfIons] = {0.};
	double force_y[NumberOfIons] = {0.};
	double force_z[NumberOfIons] = {0.};

	charge[0] = +2;	
	charge[1] = +2; 
	charge[2] = +2; 
	charge[3] = +2;

	charge[4] = -2; 
	charge[5] = -2; 
	charge[6] = -2; 
	charge[7] = -2; 

	frac_x[0] = frac_y[0] = frac_z[0] = 0.;
	frac_x[0] = 0.05;

	frac_x[1] = 0.5;
	frac_y[1] = 0.52;
	frac_z[1] = 0.;
	
	frac_x[2] = frac_z[2] = 0.5;
	frac_y[2] = 0.;

	frac_x[3] = 0.;
	frac_y[3] = 0.52;
	frac_z[3] = 0.5;

	frac_x[4] = frac_y[4] = frac_z[4] = 0.5;

	frac_x[5] = 0.5;
	frac_y[5] = frac_z[5] = 0.;
	frac_z[5] = 0.04;
	
	frac_x[6] = frac_z[6] = 0.;
	frac_y[6] = 0.5;
	frac_x[6] = -0.2;
	frac_z[6] = 0.32;

	frac_x[7] = frac_y[7] = 0.;
	frac_z[7] = 0.5;
/*
Mg core 0.00 0.00 0.00
Mg core 0.50 0.52 0.00
Mg core 0.50 0.00 0.50
Mg core 0.00 0.52 0.50

O  core 0.50 0.50 0.50
O  core 0.50 0.00 0.00
O  core 0.00 0.50 0.00
O  core 0.00 0.00 0.50


0.00000    0.00000    0.00000
0.00000    0.50000    0.50000
0.50000    0.00000    0.50000
0.50000    0.50000    0.00000

0.50000    0.00000    0.00000
0.50000    0.50000    0.50000
0.00000    0.00000    0.50000
0.00000    0.50000    0.00000
*/



	// EwaldSum Info
	double rcut, gcut;	// cut offs
	rcut = gcut = 40.;
	// ---------------------------> Set to 50 Angs.. not sure how to optimise efficient Truncate Radii

	double Sigma = (Lx+Ly+Lz)/5.6;

	//// //// //// MAIN
	
	// Working Variables
	double dr[3];		// working r1 - r2 vectors in RealSpace
	double gr[3];
	double nx,ny,nz;	// variables for translational cell images;
	nx = ny = nz = 10;

	double g_2;

	int iter = 1;
	for(int i=0;i<NumberOfIons;i++)
	{	for(int j=0;j<NumberOfIons;j++)
		{
			// Get Cartesian Coord
			cart_x[i] = Lx*frac_x[i];
			cart_y[i] = Ly*frac_y[i];
			cart_z[i] = Lz*frac_z[i];

			cart_x[j] = Lx*frac_x[j];
			cart_y[j] = Ly*frac_y[j];
			cart_z[j] = Lz*frac_z[j];
	
			dr[0] = cart_x[i] - cart_x[j];
			dr[1] = cart_y[i] - cart_y[j];
			dr[2] = cart_z[i] - cart_z[j];	// r_i - r_j;

			// Space Sum
			for(int rnx = -nx; rnx <= +nx; rnx++)
			{
				for(int rny = -ny; rny <= +ny; rny++)
				{
					for(int rnz = -nz; rnz <= +nz; rnz++)
					{
						// =================================================================================================  REAL SUM
						if( norm(rnx*Lx,rny*Ly,rnz*Ly) < rcut )		// Within the spherical radius (under rcut)
						{
							if( rnx == 0 && rny == 0 && rnz == 0 )	// In Central Unit Cell, Choosing the Same Ion
							{	
								if ( i != j )
								{
									dr[0] = cart_x[i] - cart_x[j];
									dr[1] = cart_y[i] - cart_y[j];
									dr[2] = cart_z[i] - cart_z[j];	// r_i - r_j;

									dr[0] = dr[0] - rnx*Lx;
									dr[1] = dr[1] - rny*Ly;
									dr[2] = dr[2] - rnz*Lz;

								Real_E += (0.5*charge[i]*charge[j])/norm(dr) * (erfc(norm(dr)/Sigma)) * TO_EV;
							
								}
							}
							else
							{
								dr[0] = cart_x[i] - cart_x[j];
								dr[1] = cart_y[i] - cart_y[j];
								dr[2] = cart_z[i] - cart_z[j];	// r_i - r_j;
							
								dr[0] = dr[0] - rnx*Lx;
								dr[1] = dr[1] - rny*Ly;
								dr[2] = dr[2] - rnz*Lz;

								Real_E += (0.5*charge[i]*charge[j])/norm(dr) * (erfc(norm(dr)/Sigma)) * TO_EV;

								if ( std::isinf(Real_E) )
								{
									cout << " SUM IS INF !!!\n";
									printf("Txyz  : %d\t%d\t%d\n",rnx,rny,rnz);
									printf("dr    : %12.6lf%12.6lf%12.6lf\n",dr[0],dr[1],dr[2]);
									printf("normdr / normdr-1 : %12.6lf%12.5lf\n",norm(dr),1/norm(dr));
									
									return 0;
								}
								//printf("Real E: %20.6lf\n",Real_E);
							}
						}
						// =================================================================================================  REAL SUM

						/* ***													*** */

						// =================================================================================================  RECIPROCAL SUM

						gr[0] = 2.*M_PI/Lx * rnx;
						gr[1] = 2.*M_PI/Ly * rny;
						gr[2] = 2.*M_PI/Lz * rnz;

						if( norm(gr) < gcut )
						{

							dr[0] = cart_x[i] - cart_x[j];
							dr[1] = cart_y[i] - cart_y[j];
							dr[2] = cart_z[i] - cart_z[j];	// r_i - r_j;
						/*
							if ( std::isinf(1./(norm(gr)*norm(gr))) )
							{
								cout << "1./g**2 is inf : " << norm(gr) << endl;
								printf("lattice translation indices : %d\t%d\t%d\n",rnx,rny,rnz);
								//return 0;
							}
						*/
							if( rnx == 0 && rny == 0 && rnz == 0 ) {}	// in Self at the end
							else
							{

							g_2 = dp(gr,gr);
							Reci_E += TO_EV * ((2.*M_PI)/V) * (charge[i]*charge[j]) * exp(-0.25*Sigma*Sigma*g_2)/g_2 * cos( dp(gr,dr) );
							//cout << std::setprecision(8) << Reci_E << endl;
							}
						}
						// =================================================================================================  RECIPROCAL SUM
			
						//printf("iter/real/reci : %12.d%12.6lf%12.6lf\n",iter,Real_E,Reci_E);
						iter++;
					}//rnz
				}//rny
			}//rnz
		}// Ion 'i'
	}// Ion 'j'

	printf("Real sum : %12.8lf (eV) \n", Real_E);
	printf("Reci sum : %12.8lf (eV) \n", Reci_E);

	double Reci_self = 0.;
	for(int i=0;i<NumberOfIons;i++)
	{	Reci_self += -(charge[i]*charge[i])/Sigma/sqrt(M_PI)*TO_EV;
	}

	cout << "Reci_self (eV) : " << std::setprecision(16) << Reci_self << endl;
	cout << endl;
	cout << "Recl      (eV) : " << std::setprecision(16) << Real_E << endl;
	cout << "Reci      (eV) : " << std::setprecision(16) << Reci_self + Reci_E << endl;
	cout << "total     (eV) : " << std::setprecision(16) << Real_E + Reci_E + Reci_self << endl;
	
	printf("Cell Info (Angs): %12.4lf\t%12.4lf\t%12.4lf\n",Lx,Ly,Lz);
	for(auto i=0;i<NumberOfIons;i++)
	{
		printf("Ion %d           : %12.6lf\t%12.6lf\t%12.6lf\t%12.6lf\n",i+1,charge[i],frac_x[i],frac_y[i],frac_z[i]);
	}	
		


	cout << "\n================================================ calculating forces!!\n";

	double dr_2;
	double fx_tmp, fy_tmp, fz_tmp;

	double sd[6] = {0.};

	for(int i=0;i<NumberOfIons;i++)
	{	for(int j=0;j<NumberOfIons;j++)
		{
			// Get Cartesian Coord
			cart_x[i] = Lx*frac_x[i];
			cart_y[i] = Ly*frac_y[i];
			cart_z[i] = Lz*frac_z[i];

			cart_x[j] = Lx*frac_x[j];
			cart_y[j] = Ly*frac_y[j];
			cart_z[j] = Lz*frac_z[j];

			// Space Sum
			for(int rnx = -nx; rnx <= +nx; rnx++)
			{
				for(int rny = -ny; rny <= +ny; rny++)
				{
					for(int rnz = -nz; rnz <= +nz; rnz++)
					{
						// =================================================================================================  REAL SUM
						if( norm(rnx*Lx,rny*Ly,rnz*Ly) < rcut )		// Within the spherical radius (under rcut)
						{
							if( rnx == 0 && rny == 0 && rnz == 0 )	// In Central Unit Cell, Choosing the Same Ion
							{	
								if ( i != j )
								{
									dr[0] = cart_x[i] - cart_x[j];
									dr[1] = cart_y[i] - cart_y[j];
									dr[2] = cart_z[i] - cart_z[j];	// r_i - r_j;

									dr[0] = dr[0] - rnx*Lx;
									dr[1] = dr[1] - rny*Ly;
									dr[2] = dr[2] - rnz*Lz;

									dr_2 = dp(dr,dr);

									/*
										equation by Grad (E_short) 
									*/

									force_x[i] += TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[0]/norm(dr);
									force_y[i] += TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[1]/norm(dr);
									force_z[i] += TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[2]/norm(dr);

									force_x[j] -= TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[0]/norm(dr);
									force_y[j] -= TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[1]/norm(dr);
									force_z[j] -= TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[2]/norm(dr);
									
							// Strain derivative - https://doi.org/10.1016/j.cpc.2015.01.003, Knuth et al
							sd[0] += TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[0]/norm(dr) * dr[0];
							sd[1] += TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[1]/norm(dr) * dr[1];
							sd[2] += TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[2]/norm(dr) * dr[2];

							sd[3] += TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[1]/norm(dr) * dr[2];
							sd[4] += TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[0]/norm(dr) * dr[2];
							sd[5] += TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[0]/norm(dr) * dr[1];
								}
							}
							else
							{
								dr[0] = cart_x[i] - cart_x[j];
								dr[1] = cart_y[i] - cart_y[j];
								dr[2] = cart_z[i] - cart_z[j];	// r_i - r_j;
							
								dr[0] = dr[0] - rnx*Lx;
								dr[1] = dr[1] - rny*Ly;
								dr[2] = dr[2] - rnz*Lz;

								dr_2 = dp(dr,dr);												
								force_x[i] += TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[0]/norm(dr);
								force_y[i] += TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[1]/norm(dr);
								force_z[i] += TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[2]/norm(dr);

								force_x[j] -= TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[0]/norm(dr);
								force_y[j] -= TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[1]/norm(dr);
								force_z[j] -= TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[2]/norm(dr);
							

							// Strain derivative - https://doi.org/10.1016/j.cpc.2015.01.003, Knuth et al
							sd[0] += TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[0]/norm(dr) * dr[0];
							sd[1] += TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[1]/norm(dr) * dr[1];
							sd[2] += TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[2]/norm(dr) * dr[2];

							sd[3] += TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[1]/norm(dr) * dr[2];
							sd[4] += TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[0]/norm(dr) * dr[2];
							sd[5] += TO_EV*(0.5*charge[i]*charge[j])*((-2./Sigma/sqrt(M_PI))*(exp(-dr_2/Sigma/Sigma)/norm(dr))-(erfc(norm(dr)/Sigma)/dr_2)) * dr[0]/norm(dr) * dr[1];

							}
						}
						// =================================================================================================  REAL SUM

						/* ***													*** */

						// =================================================================================================  RECIPROCAL SUM

						gr[0] = 2.*M_PI/Lx * rnx;
						gr[1] = 2.*M_PI/Ly * rny;
						gr[2] = 2.*M_PI/Lz * rnz;	// this has to be modified to deal with a,b,c != M_PI/2.

						if( norm(gr) < gcut )
						{

							dr[0] = cart_x[i] - cart_x[j];
							dr[1] = cart_y[i] - cart_y[j];
							dr[2] = cart_z[i] - cart_z[j];	// r_i - r_j;

							if( rnx == 0 && rny == 0 && rnz == 0 ) {}	// in Self at the end
							else
							{

								g_2 = dp(gr,gr);
								dr_2= dp(dr,dr);
								/*
									equation by Grad (E_long)
								*/
								force_x[i] += TO_EV*((2.*M_PI)/V)*(charge[i]*charge[j])*exp(-0.25*Sigma*Sigma*g_2)/g_2*-sin(dp(gr,dr)) * gr[0];
								force_y[i] += TO_EV*((2.*M_PI)/V)*(charge[i]*charge[j])*exp(-0.25*Sigma*Sigma*g_2)/g_2*-sin(dp(gr,dr)) * gr[1];
								force_z[i] += TO_EV*((2.*M_PI)/V)*(charge[i]*charge[j])*exp(-0.25*Sigma*Sigma*g_2)/g_2*-sin(dp(gr,dr)) * gr[2];

								force_x[j] -= TO_EV*((2.*M_PI)/V)*(charge[i]*charge[j])*exp(-0.25*Sigma*Sigma*g_2)/g_2*-sin(dp(gr,dr)) * gr[0];
								force_y[j] -= TO_EV*((2.*M_PI)/V)*(charge[i]*charge[j])*exp(-0.25*Sigma*Sigma*g_2)/g_2*-sin(dp(gr,dr)) * gr[1];
								force_z[j] -= TO_EV*((2.*M_PI)/V)*(charge[i]*charge[j])*exp(-0.25*Sigma*Sigma*g_2)/g_2*-sin(dp(gr,dr)) * gr[2];


								/// Derivatives chain w.r.t r vector (position vectors)

							        sd[0]      += TO_EV*((2.*M_PI)/V)*(charge[i]*charge[j])*exp(-0.25*Sigma*Sigma*g_2)/g_2*-sin(dp(gr,dr)) * gr[0] * dr[0];
							        sd[1]      += TO_EV*((2.*M_PI)/V)*(charge[i]*charge[j])*exp(-0.25*Sigma*Sigma*g_2)/g_2*-sin(dp(gr,dr)) * gr[1] * dr[1];
							        sd[2]      += TO_EV*((2.*M_PI)/V)*(charge[i]*charge[j])*exp(-0.25*Sigma*Sigma*g_2)/g_2*-sin(dp(gr,dr)) * gr[2] * dr[2];

								sd[3]	   += TO_EV*((2.*M_PI)/V)*(charge[i]*charge[j])*exp(-0.25*Sigma*Sigma*g_2)/g_2*-sin(dp(gr,dr)) * gr[1] * dr[2];
								sd[4]	   += TO_EV*((2.*M_PI)/V)*(charge[i]*charge[j])*exp(-0.25*Sigma*Sigma*g_2)/g_2*-sin(dp(gr,dr)) * gr[0] * dr[2];
								sd[5]	   += TO_EV*((2.*M_PI)/V)*(charge[i]*charge[j])*exp(-0.25*Sigma*Sigma*g_2)/g_2*-sin(dp(gr,dr)) * gr[0] * dr[1];	// derivative w.r.t real (r_i) !

								/// Derivative chain w.r.t. G vector (reciprocal lattice vectors)


								sd[0] += TO_EV*((2.*M_PI)/V)*(charge[i]*charge[j])
										* ( -2.*exp(-0.25*Sigma*Sigma*g_2)*gr[0]/g_2/g_2*cos(dp(gr,dr))
										   -0.5*exp(-0.25*Sigma*Sigma*g_2)*gr[0]/g_2*Sigma*Sigma*cos(dp(gr,dr))
										     -  exp(-0.25*Sigma*Sigma*g_2)/g_2*dr[0]*sin(dp(gr,dr)) ) * -gr[0];

								sd[1] += TO_EV*((2.*M_PI)/V)*(charge[i]*charge[j])
										* ( -2.*exp(-0.25*Sigma*Sigma*g_2)*gr[1]/g_2/g_2*cos(dp(gr,dr))
										   -0.5*exp(-0.25*Sigma*Sigma*g_2)*gr[1]/g_2*Sigma*Sigma*cos(dp(gr,dr))
										     -  exp(-0.25*Sigma*Sigma*g_2)/g_2*dr[1]*sin(dp(gr,dr)) ) * -gr[1];

								sd[2] += TO_EV*((2.*M_PI)/V)*(charge[i]*charge[j])
										* ( -2.*exp(-0.25*Sigma*Sigma*g_2)*gr[2]/g_2/g_2*cos(dp(gr,dr))
										   -0.5*exp(-0.25*Sigma*Sigma*g_2)*gr[2]/g_2*Sigma*Sigma*cos(dp(gr,dr))
										     -  exp(-0.25*Sigma*Sigma*g_2)/g_2*dr[2]*sin(dp(gr,dr)) ) * -gr[2];
// swap indices - w.r.t. real strain derivative sum
								sd[3] += TO_EV*((2.*M_PI)/V)*(charge[i]*charge[j])
										* ( -2.*exp(-0.25*Sigma*Sigma*g_2)*gr[2]/g_2/g_2*cos(dp(gr,dr))
										   -0.5*exp(-0.25*Sigma*Sigma*g_2)*gr[2]/g_2*Sigma*Sigma*cos(dp(gr,dr))
										     -  exp(-0.25*Sigma*Sigma*g_2)/g_2*dr[2]*sin(dp(gr,dr)) ) * -gr[1];

								sd[4] += TO_EV*((2.*M_PI)/V)*(charge[i]*charge[j])
										* ( -2.*exp(-0.25*Sigma*Sigma*g_2)*gr[2]/g_2/g_2*cos(dp(gr,dr))
										   -0.5*exp(-0.25*Sigma*Sigma*g_2)*gr[2]/g_2*Sigma*Sigma*cos(dp(gr,dr))
										     -  exp(-0.25*Sigma*Sigma*g_2)/g_2*dr[2]*sin(dp(gr,dr)) ) * -gr[0];

								sd[5] += TO_EV*((2.*M_PI)/V)*(charge[i]*charge[j])
										* ( -2.*exp(-0.25*Sigma*Sigma*g_2)*gr[1]/g_2/g_2*cos(dp(gr,dr))
										   -0.5*exp(-0.25*Sigma*Sigma*g_2)*gr[1]/g_2*Sigma*Sigma*cos(dp(gr,dr))
										     -  exp(-0.25*Sigma*Sigma*g_2)/g_2*dr[1]*sin(dp(gr,dr)) ) * -gr[0];

								/// Derivatives chain w.r.t Cell Volume
							        sd[0]      += -1.* TO_EV * ((2.*M_PI)/V/V) * V * (charge[i]*charge[j]) * exp(-0.25*Sigma*Sigma*g_2)/g_2 * cos( dp(gr,dr) );
							        sd[1]      += -1.* TO_EV * ((2.*M_PI)/V/V) * V * (charge[i]*charge[j]) * exp(-0.25*Sigma*Sigma*g_2)/g_2 * cos( dp(gr,dr) );
							        sd[2]      += -1.* TO_EV * ((2.*M_PI)/V/V) * V * (charge[i]*charge[j]) * exp(-0.25*Sigma*Sigma*g_2)/g_2 * cos( dp(gr,dr) );
							}
						}
						// =================================================================================================  RECIPROCAL SUM

					}//rnz
				}//rny
			}//rnz
		}
	}

	cout << "Final raw derivatives\n";
	for(auto i=0;i<NumberOfIons;i++)
	{
		printf("Ion %d           : %12.6lf\t%12.6lf\t%12.6lf\t%12.6lf\n",i+1,charge[i],force_x[i],force_y[i],force_z[i]);
	}

	cout << " END ================================================================================================= Finalise\n";



	cout << "dE/de1(xx) : " <<  std::setprecision(16)  << sd[0] << endl;
	cout << "dE/de2(yy) : " <<  std::setprecision(16)  << sd[1] << endl;
	cout << "dE/de3(zz) : " <<  std::setprecision(16)  << sd[2] << endl;
	cout << "dE/de4(yz) : " <<  std::setprecision(16)  << sd[3] << endl;
	cout << "dE/de5(xz) : " <<  std::setprecision(16)  << sd[4] << endl;
	cout << "dE/de6(xy) : " <<  std::setprecision(16)  << sd[5] << endl;



	return 0;
}









/*
cout << "i/j  : " << i << '/' << j << endl;		
printf("ion %d %12.4lf%12.4lf%12.4lf\n",i,cart_x[i],cart_y[i],cart_z[i]);
printf("ion %d %12.4lf%12.4lf%12.4lf\n",j,cart_x[j],cart_y[j],cart_z[j]);
*/
