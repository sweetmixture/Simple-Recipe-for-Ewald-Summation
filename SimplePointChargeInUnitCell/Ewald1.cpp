#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

#define TO_EV 14.39964390675221758120  

double dp( double rx, double ry, double rz, double gx, double gy, double gz )
{
	return rx*gx + ry*gy + rz*gz;
}

int main()
{
	// Define Unit Cell ... Single Point Charge at 0.5 0.5 0.5

	double Lx, Ly, Lz;
	Lx = 4.;//4.2112;
	Ly = 4.;//4.2112;
	Lz = 4.;//4.2112;

	double x_frac, y_frac, z_frac;
	x_frac = 0.5;
	y_frac = 0.5;
	z_frac = 0.5;

	double x, y ,z;
	x = Lx*x_frac;
	y = Ly*y_frac;
	z = Lz*z_frac;

	double V = Lx*Ly*Lz;
	double q = -2.;

	double Sigma = (Lx + Ly + Lz)/3.;


	double Real_E = 0.;
	int nx, ny, nz;
	nx = ny = nz = 50;

	double rx, ry ,rz;

	int IterReal = 0;
	for(int rnx = -nx ; rnx <= +nx ; rnx++)
	{
		for(int rny = -ny ; rny <= +ny ; rny++)
		{
			for(int rnz = -nz ; rnz <= +nz ; rnz++)
			{
				if( rnx == 0 && rny == 0 && rnz == 0 )
				{	continue;	}
				else
				{
					IterReal++;
					rx = x - x - rnx * Lx;
					ry = y - y - rny * Ly;
					rz = z - z - rnz * Lz;

					Real_E += 0.5 * (q*q) / pow(rx*rx + ry*ry + rz*rz,0.5) * erfc(pow(rx*rx + ry*ry + rz*rz,0.5)/Sigma);
					//Real_E += 0.5 * (q*q) / pow(rx*rx + ry*ry + rz*rz,0.5);// * erfc(pow(rx*rx + ry*ry + rz*rz,0.5)/Sigma);
					//cout << "Iter / REAL_E : " << IterReal << "\t" << Real_E*TO_EV << endl;
				}
			}
		}
	} 

	double Reci_E = 0.;
	double gx, gy, gz, gg;

	nx = ny = nz = 50;
	for(int rnx = -nx ; rnx <= +nx ; rnx++)
	{
		for(int rny = -ny ; rny <= +ny ; rny++)
		{
			for(int rnz = -nz ; rnz <= +nz ; rnz++)
			{
				if( rnx == 0 && rny == 0 && rnz == 0 )
				{	continue;	}
				else
				{
					gx = 2.*M_PI*rnx/Lx;
					gy = 2.*M_PI*rny/Ly;
					gz = 2.*M_PI*rnz/Lz;

					gg = dp(gx,gy,gz,gx,gy,gz);
					//cout << gg << endl;
					Reci_E += 2.*M_PI/V * (q*q) * exp(-0.25*Sigma*Sigma*gg) / gg * cos(dp(x,y,z,gx,gy,gz));
					//Reci_E += 2.*M_PI/V * (q*q) * exp(-0.25*Sigma*Sigma*gg) / gg * cos(dp(x,y,z,gx,gy,gz));
					//cout << "nx ny nz / RECI_E : " << rnx << " " << rny << " " << rnz << " " <<Reci_E*TO_EV << endl;
				}
			}
		}
	}

	// Self
	//double Self_E = q*q / (Sigma/sqrt(M_PI));

	cout << "Real/Reci: " << Real_E*TO_EV << "\t" << Reci_E*TO_EV << endl;
	cout << "SUM : " << TO_EV*(Real_E + Reci_E) << endl;

	return 0;
}
