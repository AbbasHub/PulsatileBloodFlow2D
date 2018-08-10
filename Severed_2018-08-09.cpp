
#include <iostream>
#include <cmath>

#include <cstdio>	//for I/O (either this or stdio.h)
#include <stdio.h>	//for I/O (either this or cstdio)
#include <stdlib.h>

#include <sstream>	//to convert number to string and vice versa
#include <string>

#include <ctime>		//for CPU time

#include <complex>
#include "TECIO.h"	//Include "tecio.f90"	//#include <io.h>		//for I/O

using namespace std;

/*!**********************************************************
!
!	The lattice_Boltzmann method is implemented in this code
	to solve single-phase fluid flow through a compliant vessel.

!	This code was very crudely converted from Fortran 90 to C++
	As such, it strives for readability and ease in adding other
	functionalities rather than optimization. The structures of the
	arrays are kept as 2D or 3D matrices to enhance the readability.

!	written/updated by Abbas Fakhai (09/15/2017)-(08/09/2018)

!	(2018-02-20): converted from FORTRAN 90 to C++
!	(2018-05-19): changed Fobj[Nx+2][Ny] to Fobj[Nx+2][Ny+2]
!	(2018-05-29): modified Bouzidi and other BCs
!	(2018-08-08): many other modifications and improvements
!***********************************************************/

int t, tf, step, t_start, t_sever, t_beat;

const int hsize = 1;

const int L0 = hsize * 32;

//----------------------------------------------------------------------------
const int Ny = L0 - 2 * (hsize - 1);

const int Yw1 = (0)      + 1 * hsize;	// stationary wall index at the bottom
const int Yw2 = (Ny - 1) - 1 * hsize;	// stationary wall index at the top

const int L_D = 10;	// L_D = length/diameter ratio

const int Nx = 1 + L_D * (Yw2 - Yw1 - 1);	//updated (2018-06-09)

const int X0 = (Nx - 1) / 2;	//updated (2018-06-09)
const int Y0 = (Ny - 1) / 2;	//updated (2018-06-09)
//----------------------------------------------------------------------------

const int t_propagation = int((Nx - 1.)*sqrt(3.) - 1)*hsize;

const int ex[9] = { 0, 1, 0, -1, 0, 1, -1, -1, 1 };
const int ey[9] = { 0, 0, 1, 0, -1, 1, 1, -1, -1 };
const double Wa[9] = { 16 / 36.L, 4 / 36.L, 4 / 36.L, 4 / 36.L, 4 / 36.L, 1 / 36.L, 1 / 36.L, 1 / 36.L, 1 / 36.L };

const double pi = 3.14159265358979L;

const double Rho0 = 1.L / pow(hsize, 3);

// you may change the following parameters:
//----------------------------------------------------------
const bool deformable = 1;	//0 = rigid; 1 = deformable
const bool is_severed = 1;	//0 = normal vessel; 1 = severed vessel

const double alpha = 0.01L;	// wall compliance factor

const double tau = 0.75L;	// 0.75L;	// 0.55L;
const double mu = Rho0 * (tau - 0.5L) / 3.L;
//----------------------------------------------------------
//you may change the parameters above

const double s8 = 1.L / tau;
const double s5 = 1.L;	// (2.L - s8) / (8.L - s8) * 8.L;
const double S[9] = { 1L, 1.0L, 1.L, 1L, s5, 1L, s5, s8, s8 };

double Umax, Re, Wo;
double omega;
double p0_in;	//mean pressure at the inlet (inlet pressure without the oscillatory part)
double p0_out;	//mean pressure at the outlet
double Delta_p;	//pressure difference along the vessel [Delta_p = -grad_p*(Nx-1)]
double p_tissue, p_oscillatory;	//2018-05-31: obtain from physical units (via conversion)

double g[9][Nx + 2][Ny + 2], geq[9];	// modified the 2nd and 3rd dimensions of g (02/18/18)
double P[Nx][Ny], Ux[Nx][Ny], Uy[Nx][Ny];

double Fobj[Nx + 2][Ny + 2];// modified the 1st dimension (02/19/18) and 2nd dimension (05/17/18)

double Vw1[Nx], Vw2[Nx];	//   velocity of the bottom and top walls
double yr1[Nx], yr2[Nx];	// y-location of the bottom and top walls
double y1new[Nx], y2new[Nx];

int FreshNodes, KilledNodes;

//**********************************************************************

int Nb1, Nb2;	//Number of boundary nodes (might change as the vessel walls move)

int jB[9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };	//Bounce - Back directions (ADJUSTED FOR C++)

struct Border_Node_f
{
	int X, Y;			//indices of the boundary nodes (inside the solid object)
	double Delta[8];	//Delta for interpolation
};
Border_Node_f *Borders1, *Borders2;	//pointers to be allocated later

/* --------------------------------------------------------- */
/* === function declarations === */
/* === Moving_BC: */
void convert_to_physical_units_2D(double &Pphy_Plbm, double &Lphy_Llbm);	//2018-06-09
void Initialize_Fobj_for_Vessel_Walls();
void Find_or_Update_Boundary_Nodes();
void Update_Boundary_Nodes_Bottom();
void Update_Boundary_Nodes_Top();
void Find_Delta(int mA, double mB, double Y1, double &Delta);
void Update_Fobj_for_Vessel_Walls();
void Fill_Fluid_Node(int X, int Y, int Ffrac[][3]);
void Fill_Inlet_node(int X, int Y);		//2018-05-17
void Fill_Outlet_node(int X, int Y);	//2018-05-17
void Fill_Fluid_Node_g_new(int X, int Y, int Ffrac[][3]);
void Fresh_Macroscopic_Values(int X, int Y);

void Link_Bounce_Back();
void Bouzidi_linear();
void Bouzidi_quadratic();	//2018-05-30
void Multireflection_BC();	//2018-06-18

/* === Vessel_Deformable: */
void Boundary_Conditions();
void Streaming();

void Initialize_Yr_and_Vw_and_p();
void Get_Theoretical_Wall_Locaion_and_pressure();	//2018-05-31
void Setup_Simulation_Parameters();		//2018-05-31
void Initialize_P_U_g();
void Equilibrium_g(double P, double U, double V);
void Incompressible_LBM_g();
void Collision_g();
void Inlet_ZouHe();
void Outlet_ZouHe();
void Macroscopic_Properties_g();
void DataPoints_ASCII();	//2018-08-09
void CONVERT(double IN[9], double(*OUT)[9]);
void RECONVERT(double IN[9], double(*OUT)[9]);
void BGK_Collision(int X, int Y);	//2018-06-04
void MRT_Collision(int X, int Y);	//2018-06-04
void Calculate_Pressure_and_Move_Walls();
void Calculate_Pressure_and_Move_Wall_1_Bottom();
void Calculate_Pressure_and_Move_Wall_2_Top();
void Poiseuille_Deformable();			//2018-04-20
void Poiseuille_Deformable_pressure();	//2018-04-20
void Calculate_Shear_Stress_LBM(double Pphy_Plbm, double Lphy_Llbm);	//2018-06-09
double Stress_Tensor(int Xf, int Yf);		//2018-06-10
double Stress_Tensor_BGK(int Xf, int Yf);	//2018-06-06
double Stress_Tensor_MRT(int Xf, int Yf);	//2018-06-13

/* === Utility functions: */
double MinVal1D(double arr[Nx]);
double MaxVal1D(double arr[Nx]);
double MinVal2D(double arr[Nx][Ny]);
double MaxVal2D(double arr[Nx][Ny]);
/* --------------------------------- */

/*!**********************************************************
!
!			updated by Abbas Fakhai 08/09/2018
!	(Postdoctoral fellow at the University of Pennsylvania)
!
!***********************************************************/

int main()
{

	clock_t begin_time = clock();

	//==========================================
	double Pphy_Plbm;	//to convert pressure to physical units (2018-05-29)
	double Lphy_Llbm;	//to convert length   to physical units (2018-06-09)
	convert_to_physical_units_2D(Pphy_Plbm, Lphy_Llbm);	//2018-06-09

	Setup_Simulation_Parameters();

	Initialize_Yr_and_Vw_and_p();
	//==========================================
	//= these routines are for the moving walls
	Initialize_Fobj_for_Vessel_Walls();

	Find_or_Update_Boundary_Nodes();	//fixed a bug (03/22/2018)
//= these routines are for the moving walls
//==========================================

	Initialize_P_U_g();

	//= == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
	printf("\n    tf    Lx   Ly   Yw1  Y0  Yw2  alpha   tau   Umax     Re     Wo    p_tissue \n\n");
	printf("%7i%5i%5i%5i%5i%5i%7.4f%7.3f%7.3f%8.2f%7.3f%9.3f \n\n", tf, Nx - 1, Yw2 - Yw1 - 1, Yw1, Y0, Yw2, alpha, tau, Umax, Re, Wo, p_tissue);

	printf("%6s%6s%4s%8s%8s%11s%11s%11s%11s\n\n", "t", "nF", "nS", "P_min", "P_max", "Ux_max", "Fobj_1-2", "Vw1_min", "Vw2_max");
	//= == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==

	for (t = 0; t <= tf; t++) {

		if (isnan(P[X0][Y0])) {
			printf("!!! THE PROGRAM DIVERGED AT t = %i \n", t);
			abort();	// exit(1);
		}
		else if (t%step == 0) {
			DataPoints_ASCII();
			if (t % (tf / 10) == 0 & deformable) {
				Poiseuille_Deformable(); //for steady flow in a deformable vessel (2018-04-20)
			}

			printf("%7i%4i%4i%9.4f%9.4f%11.3e%11.3e%11.3e%11.3e\n", t, FreshNodes, KilledNodes, MinVal2D(P), MaxVal2D(P), MaxVal2D(Ux), Fobj[X0][Y0] - Fobj[X0][Y0 + 1], MinVal1D(Vw1), MaxVal1D(Vw2));
		}

		Incompressible_LBM_g();

		Calculate_Shear_Stress_LBM(Pphy_Plbm, Lphy_Llbm);	//(2018-06-09)

		if (deformable) {
			Calculate_Pressure_and_Move_Walls();
		}

	}

	if (deformable) {
		Poiseuille_Deformable_pressure(); //for steady flow in a deformable vessel (2018-04-20)
	}

	clock_t end_time = clock();
	printf("\n*****************************************\n");
	printf("Time Elapsed: %12.1f", (end_time - begin_time) / double(CLOCKS_PER_SEC));
	printf("\n*****************************************\n");

	return (0);
}

//**********************************************************************

void convert_to_physical_units_2D(double &Pphy_Plbm, double &Lphy_Llbm)	//2018-06-09
{
	double D_lbm = Yw2 - Yw1 - 1.;
	double D_phy = 0.2e-3;		//Diameter = 0.2 mm
	double L_phy = L_D * D_phy;	//Length of the vessel

	double dx_lbm = 1;
	double dx_phy = dx_lbm * D_phy / D_lbm;

	double rho_lbm = Rho0;
	double rho_phy = 1060;	//blood density [kg/m3]

	double mu_lbm = mu;
	double mu_phy = 1.043e-3;	//1 cP = 1e-3 Pa.s = 1 kg/(m.s)

	double dt_lbm = 1;
	double dt_phy = dt_lbm * mu_lbm / mu_phy * (rho_phy*D_phy*D_phy) / (rho_lbm*D_lbm*D_lbm);

	double p_in_mmHg  = 20.8L;	//mean inlet  pressure (mmHg)
	double p_out_mmHg = 20.6L;	//mean outlet pressure (mmHg)
	if (is_severed) {
		p_in_mmHg  = 2.0L;
		p_out_mmHg = 0.0L;
	}

	double p_in_phy  = p_in_mmHg  * 101325.L / 760.L;	// 760 mmHg = 101325 Pa
	double p_out_phy = p_out_mmHg * 101325.L / 760.L;	// 760 mmHg = 101325 Pa
	double p_in_lbm = p_in_phy * dt_phy / dt_lbm * mu_lbm / mu_phy;

	double gradp_phy = (p_out_phy - p_in_phy) / L_phy;
	double gradp_lbm = gradp_phy * dx_phy / dx_lbm * p_in_lbm / p_in_phy;	//pressure gradient
	double Dp_lbm = -gradp_lbm * (Nx - 1);									//pressure drop

	double heart_beat_phy = 1.L;	//beat pulase duration is 1 sec
	double heart_beat_lbm = heart_beat_phy * dt_lbm / dt_phy;

	printf("=============================================================================\n");
	printf("Each iteration in LBM corresponds to %e seconds in physical units.\n", dt_phy);
	printf("Therefore the total physical time in this simulation is %7.1f seconds.\n\n", heart_beat_phy);
	printf("A blood pressure of %2.0f mmHg is %f in LB units (p0_in = %f lu)\n\n", p_in_mmHg, p_in_lbm, p_in_lbm);
	printf("The pressure grad of %6.1f Pa/m is equal to %f in LB units.\n", gradp_phy, gradp_lbm);
	printf("Therefore the pressure difference in LBM should be [pin - pout = %f] \n", Dp_lbm);
	printf("=============================================================================\n");

	Pphy_Plbm = p_in_phy / p_in_lbm;	//conversion factor for pressure
	Lphy_Llbm = D_phy / D_lbm;			//conversion factor for pressure

	p0_in = p_in_lbm;
	p0_out = p0_in - Dp_lbm;

//	double U_phy_2d = -gradp_phy * pow(D_phy, 2) / (8. * mu_phy);
	double U_phy_3d = -gradp_phy * pow(D_phy, 2) / (16 * mu_phy);

//	double Q_phy_2d = -gradp_phy * (pow(D_phy, 3) * 1.) / (12. * mu_phy);
//	double Q_phy_3d = -gradp_phy * (pow(D_phy, 4) * pi) / (144 * mu_phy);

	printf("heart_beat_lbm = %f, Ux_physical = %f \n", heart_beat_lbm, U_phy_3d);

	t_beat = (int)heart_beat_lbm;

}

//**********************************************************************

void Setup_Simulation_Parameters()
{

	p_tissue = p0_in;	//changed from p0_out to p0_in to correctly measure diameter (2018-06-01)

	p_oscillatory = (p0_in - p0_out);
	if (is_severed) {
		p_oscillatory = (p0_in - p0_out) / 10;
	}

	Delta_p = p0_out - p0_in;
//*******************************************************

	t_start = t_propagation * 2;

	tf = (t_beat + t_start) / 100 * 100;
	step = tf / 100;
//	tf = 20000;
//	step = tf / 10;

	t_sever = tf / 4 * 4;
	if (is_severed) {
		t_sever = 0;	//tf/2
	}

//*******************************************************
	omega = 2.L*pi / t_beat;

	double U_pulse = -p_oscillatory * sin(omega) / (2 * mu) * ((Y0 + 0.5L) - (Yw1 + 0.5L))*((Y0 + 0.5L) - (Yw2 - 0.5L));
	
	printf("t_beat = %i, U_pulse = %f \n", t_beat, U_pulse);
	printf("p0_in  = %f, p_oscillatory = %f \n", p0_in, p_oscillatory);
	printf("p0_out = %f, p_tissue = %f \n", p0_out, p_tissue);
	printf("dy_inlet = %f, dy_outlet = %f \n", (p0_in + p_oscillatory - p_tissue) / alpha, (p0_out - p_oscillatory - p_tissue) / alpha);
	printf("=============================================================================\n");

	if (U_pulse > 0.1L) {
		printf("Oscillatory acceleration is too large. U_pulse = %f \n", U_pulse);
		printf("Increase temporal resolution (or increase t_beat) \n");
		exit(2);
	}

	Umax = -Delta_p / (Nx - 1) * pow(Yw2 - Yw1 - 1, 2.L) / (8 * mu);

	Re = Rho0*Umax / mu / 2.L*(Yw2 - Yw1 - 1);					// Reynolds  number
	Wo = (Yw2 - Yw1 - 1) / 2.L * pow(omega * Rho0 / mu, 0.5L);	// Womersley number

}

//**********************************************************************

void Initialize_Yr_and_Vw_and_p()	//modified 2018-05-31
{

	int X;

	if (deformable) {	//the vessel walls are compliant
		Get_Theoretical_Wall_Locaion_and_pressure();
	}
	else {	//the vessel walls are rigid
		for (X = 0; X < Nx; X++) {
			yr1[X] = Yw1 + 0.5L;	//initial vessel wall location at the bottom
			yr2[X] = Yw2 - 0.5L;	//initial vessel wall location at the top
			for (int Y = 0; Y < Ny; Y++) {
				P[X][Y] = p0_in + Delta_p * X / (Nx - 1.L);
			}
		}
	}

	for (X = 0; X < Nx; X++) {
		y1new[X] = yr1[X];
		y2new[X] = yr2[X];
	}

	for (X = 0; X < Nx; X++) {
		Vw1[X] = 0;
		Vw2[X] = 0;
	}

}

//**********************************************************************

void Get_Theoretical_Wall_Locaion_and_pressure()	//2018-05-31
{

	double yr1_in = Yw1 + 0.5L - (p0_in - p_tissue) / alpha;
	double yr2_in = Yw2 - 0.5L + (p0_in - p_tissue) / alpha;

	double yr1_out = Yw1 + 0.5L - (p0_out - p_tissue) / alpha;
	double yr2_out = Yw2 - 0.5L + (p0_out - p_tissue) / alpha;

	if ((p0_in + p_oscillatory - p_tissue) / alpha > Yw1) {	//inlet expands too much
		printf("displacement (Dp/alpha)=%f is too large. \n Increase alpha, decrease Delta_p, or expand wall location \n", (p0_in + p_oscillatory - p_tissue) / alpha);
		exit(10);
	}
	if ((p0_out - p_oscillatory - p_tissue) / alpha > Ny/4) { //outlet collapses too much
		printf("displacement (Dp/alpha)=%f is too large. \n Increase alpha, decrease Delta_p, or expand wall location \n", (p0_out - p_oscillatory - p_tissue) / alpha);
		exit(10);
	}

	if (yr1_in < 1) {
		printf("yr1_in is out of bounds: yr1[0]=%f \n", yr1_in);
		exit(1);
	}
	if (yr2_in > Ny - 2) {
		printf("yr2_in is out of bounds: yr2[0]=%f \n", yr2_in);
		exit(2);
	}
	if (yr1_out < 1) {
		printf("yr1_out is out of bounds: yr1[Nx-1]=%f \n", yr1_out);
		exit(3);
	}
	if (yr2_out > Ny - 2) {
		printf("yr2_out is out of bounds: yr2[Nx-1]=%f \n", yr2_out);
		exit(4);
	}

	double R0 = (yr2_in - yr1_in) / 2;
	double RL = (yr2_out - yr1_out) / 2;
	for (int X = 0; X < Nx; X++) {
		double Rx = (pow(RL, 4) - pow(R0, 4)) * X / (Nx - 1.L) + pow(R0, 4);
		Rx = pow(Rx, 0.25L);

		yr1[X] = (Y0 + 0.5L) - Rx;	//updated (2018-06-09)
		yr2[X] = (Y0 + 0.5L) + Rx;	//updated (2018-06-09)

		for (int Y = 0; Y < Ny; Y++) {
			P[X][Y] = (yr2[X] - (Yw2 - 0.5))*alpha + p_tissue;
		}
	}

}

//**********************************************************************

void Initialize_P_U_g()
{

	int X, Y, i;

	for (X = 0; X < Nx; X++) {
		for (Y = 0; Y < Ny; Y++) {
			Ux[X][Y] = 0;
			Uy[X][Y] = 0;
		}
	}
//==Initializing the velocity profile (due to a pressure gradient)
	for (X = 0; X < Nx; X++) {
		for (Y = (int)ceil(yr1[X] - 0.01); Y <= (int)floor(yr2[X] + 0.01); Y++) {
			double dp_dx = 0;
			if (X == 0) {
				dp_dx = P[X + 1][Y] - P[X][Y];
			}
			else if (X == Nx-1) {
				dp_dx = P[X][Y] - P[X - 1][Y];
			}
			else {
				dp_dx = (P[X + 1][Y] - P[X - 1][Y])/2;
			}
			Ux[X][Y] = dp_dx / (2 * mu) * (Y - yr1[X])*(Y - yr2[X]);
		}
	}

	for (X = 0; X < Nx+2; X++) {
		for (Y = 0; Y < Ny+2; Y++) {
			for (i = 0; i < 9; i++) {
				g[i][X][Y] = 0;
			}
		}
	}

	for (X = 0; X < Nx; X++) {
		for (Y = 0; Y < Ny; Y++) {

			if (Fobj[X + 1][Y + 1] < 1) {	//solid node
				P [X][Y] = 0;
				Ux[X][Y] = 0;
				Uy[X][Y] = 0;
			}
			else {
				Equilibrium_g(P[X][Y], Ux[X][Y], Uy[X][Y]);

				for (i = 0; i < 9; i++) {
					g[i][X + 1][Y + 1] = geq[i];
				}
			}

		}
	}

}

//**********************************************************************

void Equilibrium_g(double P, double U, double V)
{

	double U2, eU;

	U2 = U*U + V*V;

	for (int i = 0; i < 9; i++) {
		eU = ex[i] * U + ey[i] * V;
		geq[i] = Wa[i] * (P + Rho0 / 3 * (eU * (3.L + 4.5L*eU) - 1.5L*U2));
	}

}

//**********************************************************************

void Incompressible_LBM_g()
{

	Collision_g();

	Boundary_Conditions();

	Streaming();

//==BC_after_streaming
	Inlet_ZouHe();
	Outlet_ZouHe();
//=================

	Macroscopic_Properties_g();

}

//**********************************************************************

void Collision_g()
{

	for (int X = 0; X < Nx; X++) {
		for (int Y = 0; Y < Ny; Y++) {

			if (Fobj[X + 1][Y + 1] < 1) continue;	//solid node

			Equilibrium_g(P[X][Y], Ux[X][Y], Uy[X][Y]);

			//*******************	COLLISION(g)	*************
		//	BGK_Collision(X, Y);
			MRT_Collision(X, Y);
			//*******************	COLLISION(g)	*************

		}
	}

}

//**********************************************************************

void BGK_Collision(int X, int Y)
{

	for (int i = 0; i < 9; i++) {
		g[i][X + 1][Y + 1] = g[i][X + 1][Y + 1] * (1.L - s8) + geq[i] * s8;
	}

}

//**********************************************************************

void MRT_Collision(int X, int Y)
{

	int i;
	double tmp[9];

	for (i = 0; i < 9; i++) tmp[i] = g[i][X + 1][Y + 1] - geq[i];
	CONVERT(tmp, &geq);
	for (i = 0; i < 9; i++) tmp[i] = geq[i] * S[i];
	RECONVERT(tmp, &geq);
	for (i = 0; i < 9; i++) g[i][X + 1][Y + 1] = g[i][X + 1][Y + 1] - geq[i];

}

//**********************************************************************

void Boundary_Conditions()
{

//	Link_Bounce_Back();
//	Bouzidi_linear();
	Bouzidi_quadratic();
//	Multireflection_BC();

}

//**********************************************************************

void Streaming()
{

	int X, Y, i;
	double fnew[8][Nx][Ny];

	for (X = 0; X < Nx; X++) {
		for (Y = 0; Y < Ny; Y++) {
			for (i = 1; i < 9; i++) {
				fnew[i - 1][X][Y] = g[i][X + 1 - ex[i]][Y + 1 - ey[i]];
			}
		}
	}

	for (X = 0; X < Nx; X++) {
		for (Y = 0; Y < Ny; Y++) {
			for (i = 1; i < 9; i++) {
				g[i][X + 1][Y + 1] = fnew[i - 1][X][Y];
			}
		}
	}

}

//**********************************************************************

void Inlet_ZouHe()
{

	int X, Y;
	double Uin, Pin;

	Pin = p0_in; 
	if (t >= t_start) Pin = p0_in + p_oscillatory * sin(omega*(t + 1 - t_start));

	X = 1;
	for (Y = (int)ceil(yr1[0] - 0.01) + 1; Y <= (int)floor(yr2[0] + 0.01) + 1; Y++) {
//	for (Y = 1; Y <= Ny; Y++) {

		Uin = Pin - g[0][X][Y] - g[2][X][Y] - 2 * g[3][X][Y] - g[4][X][Y] - 2 * g[6][X][Y] - 2 * g[7][X][Y];
		Uin = Uin * 3.L / Rho0;

		g[1][X][Y] = g[3][X][Y] + 2.L * Rho0 / 9.L * Uin;
		g[5][X][Y] = Rho0 / 18.L * Uin - 0.5L * (g[2][X][Y] - g[4][X][Y]) + g[7][X][Y];
		g[8][X][Y] = Rho0 / 18.L * Uin + 0.5L * (g[2][X][Y] - g[4][X][Y]) + g[6][X][Y];

	}

}

//**********************************************************************

void Outlet_ZouHe()
{

	int X, Y;
	double Uout, Pout;

	Pout = p0_out;
	if (t >= t_start + t_propagation) Pout = p0_out + p_oscillatory * sin(omega*(t + 1 - t_start - t_propagation));
	if (t > t_sever) Pout = 0;

	X = Nx;
	for (Y = (int)ceil(yr1[Nx - 1] - 0.01) + 1; Y <= (int)floor(yr2[Nx - 1] + 0.01) + 1; Y++) {
//	for (Y = 1; Y <= Ny; Y++) {

		Uout = g[0][X][Y] + 2 * g[1][X][Y] + g[2][X][Y] + g[4][X][Y] + 2 * g[5][X][Y] + 2 * g[8][X][Y] - Pout;
		Uout = Uout * 3.L / Rho0;

		g[3][X][Y] = g[1][X][Y] - 2.L * Rho0 / 9.L * Uout;
		g[6][X][Y] = -Rho0 / 18.L * Uout - 0.5L * (g[2][X][Y] - g[4][X][Y]) + g[8][X][Y];
		g[7][X][Y] = -Rho0 / 18.L * Uout + 0.5L * (g[2][X][Y] - g[4][X][Y]) + g[5][X][Y];

	}

}

//**********************************************************************

void Macroscopic_Properties_g()
{

	int X, Y, i;

	for (X = 0; X < Nx; X++) {
		for (Y = 0; Y < Ny; Y++) {

			if (Fobj[X + 1][Y + 1] < 1) {	//solid node
				P [X][Y] = 0;
				Ux[X][Y] = 0;
				Uy[X][Y] = 0;
			}
			else {
				P[X][Y] = 0;
				for (i = 0; i < 9; i++) {
					P[X][Y] += g[i][X + 1][Y + 1];
				}
				Ux[X][Y] = 0;
				Uy[X][Y] = 0;
				for (i = 1; i < 9; i++) {
					Ux[X][Y] += g[i][X + 1][Y + 1] * ex[i];
					Uy[X][Y] += g[i][X + 1][Y + 1] * ey[i];
				}
				Ux[X][Y] = 3 * Ux[X][Y] / Rho0;
				Uy[X][Y] = 3 * Uy[X][Y] / Rho0;
			}

		}
	}

}

//**********************************************************************

void DataPoints_ASCII()
{

	string FileName, Id;
	int static counter = -1;

	counter = counter + 1;
	Id = to_string(counter);

	FileName = "deformable_vessel_" + Id;

	FILE *fp11;
	fp11 = fopen(FileName.c_str(), "w");
	if (fp11 == NULL) {
		printf("cannot open file %s \n", FileName.c_str());
		exit(11);
	}
	fprintf(fp11, "Variables = X, Y, Ux, Uy, P, Fobj \n");
	fprintf(fp11, "Zone t =\" %6.6i \", F=Point, I= %i, J= %i \n", t, Ny, Nx);

	//--Write out the flow field data
		for (int X = 0; X < Nx; X++) {
			for (int Y = 0; Y < Ny; Y++) {
			fprintf(fp11, "%f %f %e %e %e %e \n", float(X) / L0, float(Y) / L0, Ux[X][Y], Uy[X][Y], P[X][Y], Fobj[X + 1][Y + 1]);
		}
	}

	fclose(fp11);

}

//**********************************************************************

void CONVERT(double IN[9], double(*OUT)[9])
{

	(*OUT)[0] = IN[0] + IN[1] + IN[2] + IN[3] + IN[4] + IN[5] + IN[6] + IN[7] + IN[8];
	(*OUT)[1] = -IN[1] - IN[2] - IN[3] - IN[4] + (IN[5] + IN[6] + IN[7] + IN[8]) * 2 - IN[0] * 4;
	(*OUT)[2] = IN[5] + IN[6] + IN[7] + IN[8] - (IN[1] + IN[2] + IN[3] + IN[4]) * 2 + IN[0] * 4;
	(*OUT)[3] = IN[1] - IN[3] + IN[5] - IN[6] - IN[7] + IN[8];
	(*OUT)[4] = IN[5] - IN[6] - IN[7] + IN[8] - (IN[1] - IN[3]) * 2;
	(*OUT)[5] = IN[2] - IN[4] + IN[5] + IN[6] - IN[7] - IN[8];
	(*OUT)[6] = IN[5] + IN[6] - IN[7] - IN[8] - (IN[2] - IN[4]) * 2;
	(*OUT)[7] = IN[1] - IN[2] + IN[3] - IN[4];
	(*OUT)[8] = IN[5] - IN[6] + IN[7] - IN[8];

}

//**********************************************************************

void RECONVERT(double IN[9], double(*OUT)[9])
{

	double C0, C7, C8;

	C0 = IN[0] / 9;
	C7 = IN[7] / 4;
	C8 = IN[8] / 4;

	(*OUT)[0] = C0 - (IN[1] - IN[2]) / 9;
	(*OUT)[1] = C0 - (IN[1] + 2 * IN[2]) / 36 + (IN[3] - IN[4]) / 6 + C7;
	(*OUT)[2] = C0 - (IN[1] + 2 * IN[2]) / 36 + (IN[5] - IN[6]) / 6 - C7;
	(*OUT)[3] = C0 - (IN[1] + 2 * IN[2]) / 36 - (IN[3] - IN[4]) / 6 + C7;
	(*OUT)[4] = C0 - (IN[1] + 2 * IN[2]) / 36 - (IN[5] - IN[6]) / 6 - C7;
	(*OUT)[5] = C0 + (IN[2] + 2 * IN[1]) / 36 + (IN[3] + IN[5]) / 6 + (IN[4] + IN[6]) / 12 + C8;
	(*OUT)[6] = C0 + (IN[2] + 2 * IN[1]) / 36 - (IN[3] - IN[5]) / 6 - (IN[4] - IN[6]) / 12 - C8;
	(*OUT)[7] = C0 + (IN[2] + 2 * IN[1]) / 36 - (IN[3] + IN[5]) / 6 - (IN[4] + IN[6]) / 12 + C8;
	(*OUT)[8] = C0 + (IN[2] + 2 * IN[1]) / 36 + (IN[3] - IN[5]) / 6 + (IN[4] - IN[6]) / 12 - C8;

}

//**********************************************************************

void Calculate_Pressure_and_Move_Walls()
{

	Calculate_Pressure_and_Move_Wall_1_Bottom();
	Calculate_Pressure_and_Move_Wall_2_Top();

//==now update wall locations, etc.
	Update_Fobj_for_Vessel_Walls();

	Find_or_Update_Boundary_Nodes();

}

//**********************************************************************

void Calculate_Pressure_and_Move_Wall_1_Bottom()
{

	int Xf, Yf;
	double Ps;

	for (Xf = 0; Xf < Nx; Xf++) {

		Yf = (int)ceil(yr1[Xf]);	//adjacent fluid node

		Ps =  P[Xf][Y0] - p_tissue;
	//	Ps = (P[Xf][Yf] + P[Xf][Yf + 1]) / 2 - p_tissue;	//pressure change at node(Xf, Yf)
	//	Ps = (P[Xf][Yf + 1] + P[Xf][Yf + 2]) / 2 - p_tissue;	//pressure change at node(Xf, Yf)
	//	Ps = (P[Xf][Yf + 1] + P[Xf][Yf + 2] + P[Xf][Yf + 3]) / 3 - p_tissue;	//pressure change at node(Xf, Yf)

		y1new[Xf] = (Yw1 + 0.5L) - Ps / alpha;

		if (abs(y1new[Xf] - yr1[Xf]) > 1) {
			printf("\n y1new - yr1 > 1, i.e. wall velocity (Vw1) is larger than 1.0");
			printf("\n This is because Ps (Delta_p) is too large. Increase alpha and rerun.");
			printf("\n t, Xf, Yf, y1, y1new: %i %i %i %f %f \n", t, Xf, Yf, yr1[Xf], y1new[Xf]);
			printf("FreshNodes, KilledNodes: %i %i \n", FreshNodes, KilledNodes);
			printf("P[Xf][Yf], P[Xf][Yf + 1], P[Xf][Yf + 2]: %f %f %f \n", P[Xf][Yf], P[Xf][Yf + 1], P[Xf][Yf + 2]);
			exit(7);
		}

	}
	for (int X = 0; X < Nx; X++) {
		Vw1[X] = y1new[X] - yr1[X];
		yr1[X] = y1new[X];
	}

}

//**********************************************************************

void Calculate_Pressure_and_Move_Wall_2_Top()
{

	int Xf, Yf;
	double Ps;

	for (Xf = 0; Xf < Nx; Xf++) {

		Yf = (int)floor(yr2[Xf]);	//adjacent fluid node

		Ps =  P[Xf][Y0 + 1] - p_tissue;//Ps =  P[Xf][Y0] - p_tissue;
	//	Ps = (P[Xf][Yf - 1] + P[Xf][Yf]) / 2 - p_tissue;	//pressure change at node(Xf, Yf)
	//	Ps = (P[Xf][Yf - 2] + P[Xf][Yf - 1]) / 2.L - p_tissue;	//pressure change at node(Xf, Yf)
	//	Ps = (P[Xf][Yf - 3] + P[Xf][Yf - 2] + P[Xf][Yf - 1]) / 3 - p_tissue;	//pressure change at node(Xf, Yf)

		y2new[Xf] = (Yw2 - 0.5L) + Ps / alpha;

		if (abs(y2new[Xf] - yr2[Xf]) > 1) {
			printf("\n y2new - yr2 > 1, i.e. wall velocity (Vw2) is larger than 1.0");
			printf("\n This is because Ps (Delta_p) is too large. Increase alpha and rerun.");
			printf("\n t, Xf, Yf, y2, y2new: %i %i %i %f %f \n", t, Xf, Yf, yr2[Xf], y2new[Xf]);
			exit(8);
		}

	}
	for (int X = 0; X < Nx; X++) {
		Vw2[X] = y2new[X] - yr2[X];
		yr2[X] = y2new[X];
	}

}

//**********************************************************************

void Poiseuille_Deformable()	//(2018 - 04 - 20)
{

	string FileName;
	string id;
	int static counter = -1;

	counter = counter + 1;
	id = to_string(counter);

	FileName = "Poiseuille_deformable_" + id;

	FILE *fp4;
	fp4 = fopen(FileName.c_str(), "w");
	if (fp4 == NULL) {
		printf("cannot open file %s \n", FileName.c_str());
		exit(4);
	}
	fprintf(fp4, "Variables = X, Y1_ana, y1r, Y2_ana, y2r \n");
	fprintf(fp4, "Zone t =\" %6.6i \" \n", t);

	double R0 = (yr2[0] - yr1[0]) / 2;
	double RL = (yr2[Nx - 1] - yr1[Nx - 1]) / 2;
	for (int X = 0; X < Nx; X++) {
		double Rx = (pow(RL, 4) - pow(R0, 4)) * X / (Nx - 1.L) + pow(R0, 4);
		Rx = pow(Rx, 0.25L);

		fprintf(fp4, "%i %e %e %e %e \n", X, (Y0 + 0.5L) - Rx, yr1[X], (Y0 + 0.5L) + Rx, yr2[X]);
	}

	fclose(fp4);

}

//**********************************************************************

void Poiseuille_Deformable_pressure()	//(2018 - 04 - 20)
{

	FILE *fp5;
	fp5 = fopen("Poiseuille_deformable_p.txt", "w");
	if (fp5 == NULL) {
		printf("cannot open file Poiseuille_deformable_p.txt \n");
		exit(5);
	}
	fprintf(fp5, "Variables = X, P_ana, P_LBM \n");

	double R0 = (yr2[0] - yr1[0]) / 2;
	double RL = (yr2[Nx - 1] - yr1[Nx - 1]) / 2;
	for (int X = 0; X < Nx; X++) {
		double Rx = (pow(RL, 4) - pow(R0, 4)) * X / (Nx - 1.L) + pow(R0, 4);
		Rx = pow(Rx, 0.25L);
		double Pana = ((Y0 + 0.5L) + Rx - (Yw2 - 0.5)) * alpha + p_tissue;	//updated (2018-06-09)
	//	double Pana = (yr2[X] - (Yw2 - 0.5))*alpha + p_tissue;	another way to calculate P_LBM

		fprintf(fp5, "%i %e %e \n", X, Pana, P[X][Y0]);
	}

	fclose(fp5);

}

//**********************************************************************

void Calculate_Shear_Stress_LBM(double Pphy_Plbm, double Lphy_Llbm)	//2018-06-09
{

	if (tf > 100) {
		if (t % (tf / 100) != 0) return;
	}

	string FileName;
	string id;
	int static counter = -1;

	counter = counter + 1;
	id = to_string(counter);

	FileName = "Wall_Shear_Stress_" + id;

	FILE *fp2;
	fp2 = fopen(FileName.c_str(), "w");
	if (fp2 == NULL) {
		printf("cannot open file %s \n", FileName.c_str());
		exit(2);
	}
	fprintf(fp2, "Variables = X, WSS_LBM \n");
	fprintf(fp2, "Zone t =\" %6.6i \" \n", t);

	float static wss_avg[Nx];	//to calculate average wall-shear stress over time
	int static nW[Nx];			//to calculate average wall-shear stress over time
	for (int Xf = 0; Xf < Nx; Xf++) {
		wss_avg[Xf] = 0.;
		nW[Xf] = 0;
	}

	for (int Xf = 0; Xf < Nx; Xf++) {

		int Yf = (int)ceil(yr1[Xf]);	//bottom wall for now (maybe add upper wall later)

		double h = Yf - yr1[Xf];

		if (h < 0.05 || h > 0.95) continue;	//added 2018-06-09
	//	if (h < 1.e-6) Yf = Yf + 1;	//added 04/30/2018 to avoid devision by 0
	//	if (Fobj[Xf + 1][Yf + 1] < 1) {	//solid node
	//		printf("Calculate_Shear_Stress_LBM: Fobj(Xf,Yf)<1, %9.3f \n", Fobj[Xf + 1][Yf + 1]);
	//	}
		//===========================================
		double WSS_1 = Stress_Tensor(Xf, Yf);
		double WSS_2 = Stress_Tensor(Xf, Yf + 1);
	//	double WSS_3 = Stress_Tensor(Xf, Yf + 2);

		double WSS = (1 + h) * WSS_1 - h * WSS_2;	//linear extrapolation
	//	double WSS = (1 + h)*(2 + h) / 2 * WSS_1 - h*(2 + h)*WSS_2 + h*(1 + h) / 2 * WSS_3;	//quadratic extrapolation
		//=================================================================================
		wss_avg[Xf] = wss_avg[Xf] + (float)WSS;
		nW[Xf] ++;
		//=====================================

		fprintf(fp2, "%e %e \n", Xf*Lphy_Llbm*1000, WSS*Pphy_Plbm*10);	// 1 N/m2 = 10 Dyne/cm2
	}
	if (t == tf) {	//semi-analytical shear stress at the bottom wall
		fprintf(fp2, "Zone t = \"analytical\" \n");
		for (int Xf = 1; Xf < Nx - 1; Xf++) {
			double Txy_ana = 0.5 * (yr1[Xf] - yr2[Xf]) * (P[Xf + 1][Y0] - P[Xf - 1][Y0]) / 2;
			fprintf(fp2, "%e %e \n", Xf*Lphy_Llbm*1000, Txy_ana*Pphy_Plbm*10);
		}

	//= write the average wall shear stress
		FILE *fp3;
		fp3 = fopen("WSS_avg.txt", "w");
		if (fp3 == NULL) {
			printf("cannot open file WSS_avg.txt \n");
			exit(3);
		}
		fprintf(fp3, "Variables = x(mm), WSS_avg \n");
		fprintf(fp3, "Zone t = \"Avg WSS\" \n");
		for (int Xf = 0; Xf < Nx; Xf++) {
			if (nW[Xf] == 0) continue;
			fprintf(fp3, "%e %e \n", Xf*Lphy_Llbm*1000, wss_avg[Xf]/nW[Xf]*Pphy_Plbm*10);
		}
		fclose(fp3);

	}

	fclose(fp2);

}

//**********************************************************************

double Stress_Tensor(int Xf, int Yf)	//2018-06-06
{

	double WSS = Stress_Tensor_BGK(Xf, Yf);
//	double WSS = Stress_Tensor_MRT(Xf, Yf);	//not much different than BGK

	return(WSS);

}

//**********************************************************************

double Stress_Tensor_BGK(int Xf, int Yf)	//2018-06-06
{
	double mx, my;

	//= calculating the slope of the wall segment
	if (Xf == 0) {
	//	mx = 1.;
	//	my = yr1[Xf + 1]- yr1[Xf];
		mx = 2.;
		my = -yr1[Xf + 2] + 4 * yr1[Xf + 1] - 3 * yr1[Xf];
	}
	else if (Xf == Nx - 1) {
	//	mx = 1.;
	//	my = yr1[Xf] - yr1[Xf - 1];
		mx = 2.;
		my = -yr1[Xf] + 4 * yr1[Xf - 1] - 3 * yr1[Xf - 2];
	}
	else {
		mx = 2.;	//slope of the wall
		my = yr1[Xf + 1] - yr1[Xf - 1];
	}

	//= calculating the vector normal to the wall
	double nx = -my / (mx*mx + my*my);
	double ny =  mx / (mx*mx + my*my);

	//= calculating the stress tensor
	Equilibrium_g(P[Xf][Yf], Ux[Xf][Yf], Uy[Xf][Yf]);

	double T_xx = 0;
	double T_xy = 0;
	double T_yy = 0;
	for (int i = 1; i < 9; i++) {
		T_xx += (g[i][Xf + 1][Yf + 1] - geq[i]) * ex[i] * ex[i];
		T_xy += (g[i][Xf + 1][Yf + 1] - geq[i]) * ex[i] * ey[i];
		T_yy += (g[i][Xf + 1][Yf + 1] - geq[i]) * ey[i] * ey[i];
	}

	//= calculating Wall Shear Stress (WSS)
	double WSS = T_xx * nx*(1 - nx*nx) + T_xy * ny*(1 - 2 * nx*nx) - T_yy * nx*ny*ny;
	WSS = WSS * (0.5L / tau - 1) * 3;	//2018-05-08 (*3 added to correct for g=f/3)

	//= calculating normal stress (2018-06-13)
//	double WSS = T_yy * ny*(1 - ny*ny) + T_xy * nx*(1 - 2 * ny*ny) - T_xx * ny*nx*nx;
//	WSS = WSS * (0.5L / tau - 1) * 3;	//2018-05-08 (*3 added to correct for g=f/3)

	return(WSS);

}

//**********************************************************************

double Stress_Tensor_MRT(int Xf, int Yf)	//2018-06-13
{
	double mx, my;

	//= calculating the slope of the wall segment
	if (Xf == 0) {
		mx = 2.;
		my = -yr1[Xf + 2] + 4 * yr1[Xf + 1] - 3 * yr1[Xf];
	}
	else if (Xf == Nx - 1) {
		mx = 2.;
		my = -yr1[Xf] + 4 * yr1[Xf - 1] - 3 * yr1[Xf - 2];
	}
	else {
		mx = 2.;
		my = yr1[Xf + 1] - yr1[Xf - 1];
	}

	//= calculating the vector normal to the wall
	double nx = -my / (mx*mx + my*my);
	double ny =  mx / (mx*mx + my*my);

	//= calculating the stress tensor
	Equilibrium_g(P[Xf][Yf], Ux[Xf][Yf], Uy[Xf][Yf]);

	double tmp[9];
	for (int i = 0; i < 9; i++) tmp[i] = g[i][Xf + 1][Yf + 1] - geq[i];
	CONVERT(tmp, &geq);
	for (int i = 0; i < 9; i++) tmp[i] = geq[i] * S[i];
	RECONVERT(tmp, &geq);

	double T_xx = 0;
	double T_xy = 0;
	double T_yy = 0;
	for (int i = 1; i < 9; i++) {
		T_xx += geq[i] * ex[i] * ex[i];
		T_xy += geq[i] * ex[i] * ey[i];
		T_yy += geq[i] * ey[i] * ey[i];
	}

	//= calculating Wall Shear Stress (WSS)
	double WSS = T_xx * nx*(1 - nx*nx) + T_xy * ny*(1 - 2 * nx*nx) - T_yy * nx*ny*ny;
	WSS = WSS * (0.5L - tau) * 3;	//2018-05-08 (*3 added to correct for g=f/3)

	//= calculating normal stress (2018-06-13)
//	double WSS = T_yy * ny*(1 - ny*ny) + T_xy * nx*(1 - 2 * ny*ny) - T_xx * ny*nx*nx;
//	WSS = WSS * (0.5L - tau) * 3;	//2018-05-08 (*3 added to correct for g=f/3)

	return(WSS);

}

//**********************************************************************

double MinVal1D(double arr[Nx])
{
	double minval = arr[0];

	for (int i = 1; i < Nx; i++) {
		if (arr[i] < minval) minval = arr[i];
	}

	return(minval);

}

//**********************************************************************

double MaxVal1D(double arr[Nx])
{
	double maxval = arr[0];

	for (int i = 1; i < Nx; i++) {
		if (arr[i] > maxval) maxval = arr[i];
	}

	return(maxval);

}

//**********************************************************************

double MinVal2D(double arr[Nx][Ny])
{
	double minval = arr[0][Y0];

	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			if (Fobj[i + 1][j + 1] < 1) { continue; } //avoid solid nodes (2018-06-01)
			if (arr[i][j] < minval) minval = arr[i][j];
		}
	}

	return(minval);

}

//**********************************************************************

double MaxVal2D(double arr[Nx][Ny])
{
	double maxval = arr[0][Y0];

	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			if (Fobj[i + 1][j + 1] < 1) { continue; } //avoid solid nodes (2018-06-01)
			if (arr[i][j] > maxval) maxval = arr[i][j];
		}
	}

	return(maxval);

}

/*!************************************************
!
!	subroutines for curved boundary treatment
!
!*************************************************/

void Initialize_Fobj_for_Vessel_Walls()	//updated 2018-06-09
{

	int X, Y;

//= if Fobj < 1 then it is treated as a solid node (Fobj = 1 is fluid)
	for (X = 0; X < Nx; X++) {
		for (Y = -1; Y <= Y0; Y++) {
			Fobj[X + 1][Y + 1] = (yr1[X] - (Y0 + 0.5L)) / (Y - (Y0 + 0.5L));
		}
		for (Y = Y0 + 1; Y < Ny + 1; Y++) {
			Fobj[X + 1][Y + 1] = (yr2[X] - (Y0 + 0.5L)) / (Y - (Y0 + 0.5L));
		}
	}
//= now take care of the end points
	for (Y = 0; Y < Ny + 2; Y++) {
		Fobj[0][Y] = Fobj[1][Y] * 2 - Fobj[2][Y];
		Fobj[Nx + 1][Y] = Fobj[Nx][Y] * 2 - Fobj[Nx - 1][Y];
	}

}

//**********************************************************************

void Find_or_Update_Boundary_Nodes()
{
	Update_Boundary_Nodes_Bottom();

	Update_Boundary_Nodes_Top();

}

//**********************************************************************

void Update_Boundary_Nodes_Bottom()	//BUG FIXED: change X to X+1 in Fobj (04/05/2018)
{

	int i, X, Y, Yx, n;
	double Delta[8];

	if (Borders1 != 0) delete[] Borders1;

//--first find the number of boundary nodes at the bottom wall
	X = 0;
	Y = (int)floor(yr1[X]);
	if (Fobj[X + 1][Y + 1] >= 1) { //BUG IN C++: This is a fluid node (2018-06-09)
		Y = Y - 1;
	}
	n = 1;	//bug: must be 1 not 0 (03/22/2018)
	for (X = 1; X < Nx; X++) {
		Yx = (int)floor(yr1[X]);
		if (Fobj[X + 1][Yx + 1] >= 1) { //BUG IN C++: This is a fluid node (2018-06-09)
			Yx = Yx - 1;
		}
		if (Yx != Y) {
			n = n + 1;
		}
		n = n + 1;
		Y = Yx;
	}
	Nb1 = n;
	Borders1 = new Border_Node_f[n];

//--now store the x and y values of the border nodes in Borders1
	X = 0;
	Y = (int)floor(yr1[X]);
	if (Fobj[X + 1][Y + 1] >= 1) { //BUG IN C++: This is a fluid node (2018-06-09)
		Y = Y - 1;
	}
	n = 0;
	Borders1[n].X = X;
	Borders1[n].Y = Y;
	//*********** find Delta ***********
	for (i = 0; i < 8; i++) Delta[i] = 2;	//making sure only Boundary Nodes with Delta < 1 are considered
	if (Fobj[X + 2][Y + 1] >= 1) {	//the neighboring cell is a fluid cell
		Find_Delta(0, yr1[X + 1] - yr1[X], yr1[X] - Y, Delta[0]);	//for alpha = 1
	}
	Delta[1] = 1 - (yr1[X] - Y);									//for alpha = 2
	if (Fobj[X + 2][Y + 2] >= 1) {	//the neighboring cell is a fluid cell
		Find_Delta(1, yr1[X + 1] - yr1[X], yr1[X] - Y, Delta[4]);	//for alpha = 5
	}
	//**********************************
	for (i = 0; i < 8; i++) Borders1[n].Delta[i] = Delta[i];
	for (X = 1; X < Nx - 1; X++) {
		Yx = (int)floor(yr1[X]);
		if (Fobj[X + 1][Yx + 1] >= 1) { //BUG IN C++: This is a fluid node (2018-06-09)
			Yx = Yx - 1;
		}
		if (Yx != Y) {
			n = n + 1;
			for (i = 0; i < 8; i++) Delta[i] = 2;	//making sure only Boundary Nodes with Delta < 1 are considered
			if (Yx > Y) {
				Borders1[n].X = X;
				Borders1[n].Y = Y;
				Find_Delta(-1, yr1[X] - yr1[X - 1], yr1[X] - Y, Delta[5]);		//for alpha = 6
			}
			else {
				Borders1[n].X = X - 1;
				Borders1[n].Y = Yx;
				Find_Delta(1, yr1[X] - yr1[X - 1], yr1[X - 1] - Yx, Delta[4]);	//for alpha = 5
			}
			for (i = 0; i < 8; i++) Borders1[n].Delta[i] = Delta[i];
		}
		n = n + 1;
		Borders1[n].X = X;
		Borders1[n].Y = Yx;
		//*********** find Delta ***********
		for (i = 0; i < 8; i++) Delta[i] = 2;	//making sure only Boundary Nodes with Delta < 1 are considered
		if (Fobj[X + 2][Yx + 1] >= 1) {	//the neighboring cell is a fluid cell
			Find_Delta(0, yr1[X + 1] - yr1[X], yr1[X] - Yx, Delta[0]);	//for alpha = 1
		}
		Delta[1] = 1 - (yr1[X] - Yx);									//for alpha = 2
		if (Fobj[X][Yx + 1] >= 1) {	//the neighboring cell is a fluid cell
			Find_Delta(0, yr1[X] - yr1[X - 1], yr1[X] - Yx, Delta[2]);	//for alpha = 3
		}
		if (Fobj[X + 2][Yx + 2] >= 1) {	//the neighboring cell is a fluid cell
			Find_Delta(1, yr1[X + 1] - yr1[X], yr1[X] - Yx, Delta[4]);	//for alpha = 5
		}
		if (Fobj[X][Yx + 2] >= 1) {	//the neighboring cell is a fluid cell
			Find_Delta(-1, yr1[X] - yr1[X - 1], yr1[X] - Yx, Delta[5]);	//for alpha = 6
		}
		//**********************************
		for (i = 0; i < 8; i++) Borders1[n].Delta[i] = Delta[i];

		Y = Yx;
	}
	X = Nx - 1;
	Yx = (int)floor(yr1[X]);
	if (Fobj[X + 1][Yx + 1] >= 1) { //BUG IN C++: This is a fluid node (2018-06-09)
		Yx = Yx - 1;
	}
	if (Yx != Y) {
		n = n + 1;
		for (i = 0; i < 8; i++) Delta[i] = 2;	//making sure only Boundary Nodes with Delta < 1 are considered
		if (Yx > Y) {
			Borders1[n].X = X;
			Borders1[n].Y = Y;
			Find_Delta(-1, yr1[X] - yr1[X - 1], yr1[X] - Y, Delta[5]);		//for alpha = 6
		}
		else {
			Borders1[n].X = X - 1;
			Borders1[n].Y = Yx;
			Find_Delta(1, yr1[X] - yr1[X - 1], yr1[X - 1] - Yx, Delta[4]);	//for alpha = 5
		}
		for (i = 0; i < 8; i++) Borders1[n].Delta[i] = Delta[i];
	}
	n = n + 1;
	Borders1[n].X = X;
	Borders1[n].Y = Yx;
	//*********** find Delta ***********
	for (i = 0; i < 8; i++) Delta[i] = 2;	//making sure only Boundary Nodes with Delta < 1 are considered
	Delta[1] = 1 - (yr1[X] - Yx);									//for alpha = 2
	if (Fobj[X][Yx + 1] >= 1) {	//the neighboring cell is a fluid cell
		Find_Delta(0, yr1[X] - yr1[X - 1], yr1[X] - Yx, Delta[2]);	//for alpha = 3
	}
	if (Fobj[X][Yx + 2] >= 1) {	//the neighboring cell is a fluid cell
		Find_Delta(-1, yr1[X] - yr1[X - 1], yr1[X] - Yx, Delta[5]);	//for alpha = 6
	}
	//**********************************
	for (i = 0; i < 8; i++) Borders1[n].Delta[i] = Delta[i];

}

//**********************************************************************

void Update_Boundary_Nodes_Top()	//BUG FIXED: change X to X+1 in Fobj (04/05/2018)
{

	int i, X, Y, Yx, n;
	double Delta[8];

	if (Borders2 != 0) delete[] Borders2;

//--first find the number of boundary nodes at the top wall
	X = 0;
	Y = (int)ceil(yr2[X]);
	if (Fobj[X + 1][Y + 1] >= 1) { //BUG IN C++: This is a fluid node (2018-06-09)
		Y = Y + 1;
	}
	n = 1;	//bug: must be 1 not 0 (03/22/2018)
	for (X = 1; X < Nx; X++) {
		Yx = (int)ceil(yr2[X]);
		if (Fobj[X + 1][Yx + 1] >= 1) { //BUG IN C++: This is a fluid node (2018-06-09)
			Yx = Yx + 1;
		}
		if (Yx != Y) {
			n = n + 1;
		}
		n = n + 1;
		Y = Yx;
	}
	Nb2 = n;
	Borders2 = new Border_Node_f[n];

//--now store the x and y values of the border nodes in Borders2
	X = 0;
	Y = (int)ceil(yr2[X]);
	if (Fobj[X + 1][Y + 1] >= 1) { //BUG IN C++: This is a fluid node (2018-06-09)
		Y = Y + 1;
	}
	n = 0;
	Borders2[n].X = X;
	Borders2[n].Y = Y;
	//*********** find Delta ***********
	for (i = 0; i < 8; i++) Delta[i] = 2;	//making sure only Boundary Nodes with Delta < 1 are considered
	if (Fobj[X + 2][Y + 1] >= 1) {	//the neighboring cell is a fluid cell
		Find_Delta(0, yr2[X + 1] - yr2[X], yr2[X] - Y, Delta[0]);	//for alpha = 1
	}
	Delta[3] = 1 - (Y - yr2[X]);									//for alpha = 4
	if (Fobj[X + 2][Y] >= 1) {	//the neighboring cell is a fluid cell
		Find_Delta(-1, yr2[X + 1] - yr2[X], yr2[X] - Y, Delta[7]);	//for alpha = 8
	}
	//**********************************
	for (i = 0; i < 8; i++) Borders2[n].Delta[i] = Delta[i];
	for (X = 1; X < Nx - 1; X++) {
		Yx = (int)ceil(yr2[X]);
		if (Fobj[X + 1][Yx + 1] >= 1) { //BUG IN C++: This is a fluid node (2018-06-09)
			Yx = Yx + 1;
		}
		if (Yx != Y) {
			n = n + 1;
			for (i = 0; i < 8; i++) Delta[i] = 2;	//making sure only Boundary Nodes with Delta < 1 are considered
			if (Yx > Y) {
				Borders2[n].X = X - 1;
				Borders2[n].Y = Yx;
				Find_Delta(-1, yr2[X] - yr2[X - 1], yr2[X - 1] - Yx, Delta[7]);	//for alpha = 8
			}
			else {
				Borders2[n].X = X;
				Borders2[n].Y = Y;
				Find_Delta(1, yr2[X] - yr2[X - 1], yr2[X] - Y, Delta[6]);		//for alpha = 7
			}
			for (i = 0; i < 8; i++) Borders2[n].Delta[i] = Delta[i];
		}
		n = n + 1;
		Borders2[n].X = X;
		Borders2[n].Y = Yx;
		//*********** find Delta ***********
		for (i = 0; i < 8; i++) Delta[i] = 2;	//making sure only Boundary Nodes with Delta < 1 are considered
		if (Fobj[X + 2][Yx + 1] >= 1) {	//the neighboring cell is a fluid cell
			Find_Delta(0, yr2[X + 1] - yr2[X], yr2[X] - Yx, Delta[0]);	//for alpha = 1
		}
		if (Fobj[X][Yx + 1] >= 1) {	//the neighboring cell is a fluid cell
			Find_Delta(0, yr2[X] - yr2[X - 1], yr2[X] - Yx, Delta[2]);	//for alpha = 3
		}
		Delta[3] = 1 - (Yx - yr2[X]);									//for alpha = 4	(BUG fixed 1/8/18)
		if (Fobj[X][Yx] >= 1) {	//the neighboring cell is a fluid cell
			Find_Delta(1, yr2[X] - yr2[X - 1], yr2[X] - Yx, Delta[6]);	//for alpha = 7
		}
		if (Fobj[X + 2][Yx] >= 1) {	//the neighboring cell is a fluid cell
			Find_Delta(-1, yr2[X + 1] - yr2[X], yr2[X] - Yx, Delta[7]);	//for alpha = 8
		}
		//**********************************
		for (i = 0; i < 8; i++) Borders2[n].Delta[i] = Delta[i];

		Y = Yx;
	}
	X = Nx - 1;
	Yx = (int)ceil(yr2[X]);
	if (Fobj[X + 1][Yx + 1] >= 1) { //BUG IN C++: This is a fluid node (2018-06-09)
		Yx = Yx + 1;
	}
	if (Yx != Y) {
		n = n + 1;
		for (i = 0; i < 8; i++) Delta[i] = 2;	//making sure only Boundary Nodes with Delta < 1 are considered
		if (Yx > Y) {
			Borders2[n].X = X - 1;
			Borders2[n].Y = Yx;
			Find_Delta(-1, yr2[X] - yr2[X - 1], yr2[X] - Yx, Delta[7]);	//for alpha = 8
		}
		else {
			Borders2[n].X = X;
			Borders2[n].Y = Y;
			Find_Delta(1, yr2[X] - yr2[X - 1], yr2[X] - Y, Delta[6]);	//for alpha = 7
		}
		for (i = 0; i < 8; i++) Borders2[n].Delta[i] = Delta[i];
	}
	n = n + 1;
	Borders2[n].X = X;
	Borders2[n].Y = Yx;
	//*********** find Delta ***********
	for (i = 0; i < 8; i++) Delta[i] = 2;	//making sure only Boundary Nodes with Delta < 1 are considered
	if (Fobj[X][Yx + 1] >= 1) {	//the neighboring cell is a fluid cell
		Find_Delta(0, yr2[X] - yr2[X - 1], yr2[X] - Yx, Delta[2]);	//for alpha = 3
	}
	Delta[3] = 1 - (Yx - yr2[X]);									//for alpha = 4	(BUG fixed 1 / 8 / 18)
	if (Fobj[X][Yx] >= 1) {	//the neighboring cell is a fluid cell
		Find_Delta(1, yr2[X] - yr2[X - 1], yr2[X] - Yx, Delta[6]);	//for alpha = 7
	}
	//**********************************
	for (i = 0; i < 8; i++) Borders2[n].Delta[i] = Delta[i];

}

//**********************************************************************

void Find_Delta(int mA, double mB, double Y1, double &Delta)
{

	Delta = 1 - abs(Y1 / (mA - mB));

	if (Delta < 0) { //because of stupid C++ (round-off error leads to a small negative number)
	//	printf("Delta=%e \n", Delta);
		Delta = 0;
	}

}

//**********************************************************************

void Update_Fobj_for_Vessel_Walls()
{

	int X, Y, c1, c2;
	int Ffrac[3][3];
	double Fobj_old[Nx + 2][Ny + 2];

	for (X = 0; X < Nx + 2; X++) {
		for (Y = 0; Y < Ny + 2; Y++) {
			Fobj_old[X][Y] = Fobj[X][Y];
		}
	}

	Initialize_Fobj_for_Vessel_Walls();	//updated 2018-06-09

	//= = check and see if any solid node has become a fluid node
	c1 = 0;
	c2 = 0;
	for (X = 1; X <= Nx; X++) {
		for (Y = 1; Y <= Ny; Y++) {		//changed from Y<Ny-1 to Y<=Ny (2018-05-17)
			if (Fobj_old[X][Y] < 1) {	//it was a solid node
				if (Fobj[X][Y] >= 1) {	//it has become a fluid node

					c1 = c1 + 1;

					for (int i = -1; i <= 1; i++) {
						for (int j = -1; j <= 1; j++) {
							Ffrac[i + 1][j + 1] = int(Fobj_old[X + i][Y + j]);
						}
					}

					Fill_Fluid_Node(X - 1, Y - 1, Ffrac);

				}
			}
			if (Fobj_old[X][Y] >= 1) {	//it was a fluid node
				if (Fobj[X][Y] < 1) {	//it has become a solid node
					c2 = c2 + 1;
				}
			}
		}
	}
	FreshNodes = c1; //+ FreshNodes
	KilledNodes = c2; //+ KilledNodes

}

//**********************************************************************

void Fill_Fluid_Node(int X, int Y, int Ffrac[][3])	//modified 2018 - 05 - 17
{

	if (X == 0) {		//bug fixed(2018 - 05 - 17)
		Fill_Inlet_node(X, Y);
	}
	else if(X == Nx-1) {	//bug fixed(2018 - 05 - 17)
		Fill_Outlet_node(X, Y);
	}
	else {
		Fill_Fluid_Node_g_new(X, Y, Ffrac);	//(01 / 16 / 2018)
	}

	Fresh_Macroscopic_Values(X, Y);	//(01 / 29 / 2018)

}

//**********************************************************************

void Fill_Inlet_node(int X, int Y)	//written 2018-05-17
{

	if (Y < Y0) {	//new inlet - fluid node at the bottom - left
		g[0][X + 1][Y + 1] = g[0][X + 1][Y + 2];
		g[1][X + 1][Y + 1] = g[1][X + 1][Y + 2];
		g[2][X + 1][Y + 1] = g[2][X + 1][Y + 2];
		g[3][X + 1][Y + 1] = g[3][X + 1][Y + 2];
	//	g[4][X + 1][Y + 1] = g[4][X + 1][Y + 2];
		g[5][X + 1][Y + 1] = g[5][X + 1][Y + 2];
		g[6][X + 1][Y + 1] = g[6][X + 1][Y + 2];
	//	g[7][X + 1][Y + 1] = g[7][X + 1][Y + 2];
		g[8][X + 1][Y + 1] = g[8][X + 1][Y + 2];
	}
	else{	//new inlet - fluid node at the top - left
		g[0][X + 1][Y + 1] = g[0][X + 1][Y];
		g[1][X + 1][Y + 1] = g[1][X + 1][Y];
	//	g[2][X + 1][Y + 1] = g[2][X + 1][Y];
		g[3][X + 1][Y + 1] = g[3][X + 1][Y];
		g[4][X + 1][Y + 1] = g[4][X + 1][Y];
		g[5][X + 1][Y + 1] = g[5][X + 1][Y];
	//	g[6][X + 1][Y + 1] = g[6][X + 1][Y];
		g[7][X + 1][Y + 1] = g[7][X + 1][Y];
		g[8][X + 1][Y + 1] = g[8][X + 1][Y];
	}

}

//**********************************************************************

void Fill_Outlet_node(int X, int Y)	//written 2018-05-17
{

	if (Y < Y0) {	//new outlet - fluid node at the bottom - right
		g[0][X + 1][Y + 1] = g[0][X + 1][Y + 2];
		g[1][X + 1][Y + 1] = g[1][X + 1][Y + 2];
		g[2][X + 1][Y + 1] = g[2][X + 1][Y + 2];
		g[3][X + 1][Y + 1] = g[3][X + 1][Y + 2];
	//	g[4][X + 1][Y + 1] = g[4][X + 1][Y + 2];
		g[5][X + 1][Y + 1] = g[5][X + 1][Y + 2];
		g[6][X + 1][Y + 1] = g[6][X + 1][Y + 2];
		g[7][X + 1][Y + 1] = g[7][X + 1][Y + 2];
	//	g[8][X + 1][Y + 1] = g[8][X + 1][Y + 2];
	}
	else {			//new outlet - fluid node at the top - right
		g[0][X + 1][Y + 1] = g[0][X + 1][Y];
		g[1][X + 1][Y + 1] = g[1][X + 1][Y];
	//	g[2][X + 1][Y + 1] = g[2][X + 1][Y];
		g[3][X + 1][Y + 1] = g[3][X + 1][Y];
		g[4][X + 1][Y + 1] = g[4][X + 1][Y];
	//	g[5][X + 1][Y + 1] = g[5][X + 1][Y];
		g[6][X + 1][Y + 1] = g[6][X + 1][Y];
		g[7][X + 1][Y + 1] = g[7][X + 1][Y];
		g[8][X + 1][Y + 1] = g[8][X + 1][Y];
	}

}

//**********************************************************************

void Fill_Fluid_Node_g_new(int X, int Y, int Ffrac[][3])	//bug here (4/4/18) changed Y -> Y+1
{

	int SumFrac = 0;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			SumFrac += Ffrac[i][j];
		}
	}

	if (SumFrac == 0) {
		printf("Fill_Fluid_Node_g_new: SumFrac=0!!! program terminated...");
		exit(3);
	}

	for (int i = 0; i < 9; i++) { //BUG HERE!!! (02/26/2018): MUST BE (-ex(i),-ey(i)) INSTEAD OF (-I,-I)
		if (Ffrac[1 - ex[i]][1 - ey[i]] != 1) {	//the neighboring node was a solid node
			g[i][X + 1][Y + 1] = (g[i][X][Y] * Ffrac[0][0] + g[i][X + 1][Y] * Ffrac[1][0] + g[i][X + 2][Y] * Ffrac[2][0]
				+ g[i][X][Y + 1] * Ffrac[0][1] + g[i][X + 2][Y + 1] * Ffrac[2][1]
				+ g[i][X][Y + 2] * Ffrac[0][2] + g[i][X + 1][Y + 2] * Ffrac[1][2] + g[i][X + 2][Y + 2] * Ffrac[2][2]) / SumFrac;
		}
		else {	//the neihgboring node was a fluid node, therefore its PDF has been streamed into this node
				//DO NOTHING(make sure in the "streaming" PDFs are streamed into the solid nodes)
		}
	}

}

//**********************************************************************

void Fresh_Macroscopic_Values(int X, int Y)	//bug here (4/4/18) changed Y -> Y+1
{

	P[X][Y] = 0;
	for (int i = 0; i < 9; i++) {
		P[X][Y] += g[i][X + 1][Y + 1];
	}
	Ux[X][Y] = 0;
	Uy[X][Y] = 0;
	for (int i = 1; i < 9; i++) {
		Ux[X][Y] += g[i][X + 1][Y + 1] * ex[i];
		Uy[X][Y] += g[i][X + 1][Y + 1] * ey[i];
	}
	Ux[X][Y] = 3 * Ux[X][Y] / Rho0;
	Uy[X][Y] = 3 * Uy[X][Y] / Rho0;

}

//**********************************************************************

void Link_Bounce_Back()
{

	int n, I, j, X, Y, Xf, Yf;

	for (n = 0; n < Nb1; n++) {
		X = Borders1[n].X + 1; //bug fixed (04/03/2018): +1 for g
		Y = Borders1[n].Y + 1; //bug fixed (04/03/2018): +1 for g
		for (I = 1; I <= 8; I++) {
			j = jB[I];
			Xf = X + ex[I];
			Yf = Y + ey[I];
			g[I][X][Y] = g[j][Xf][Yf]; //+ 2 * Rho0 * Wa[I] * (ey[I] * Vw1[X-1]);	//for incompressible LBM(2017 / 11 / 19)
		}
	}

	for (n = 0; n < Nb2; n++) {
		X = Borders2[n].X + 1; //bug fixed (04/03/2018): +1 for g
		Y = Borders2[n].Y + 1; //bug fixed (04/03/2018): +1 for g
		for (I = 1; I <= 8; I++) {
			j = jB[I];
			Xf = X + ex[I];
			Yf = Y + ey[I];
			g[I][X][Y] = g[j][Xf][Yf]; //+ 2 * Rho0 * Wa[I] * (ey[I] * Vw2[X - 1]);	//for incompressible LBM(2017 / 11 / 19)
		}
	}

}

//**********************************************************************

void Bouzidi_linear()
{

	int n, I, j, X, Y, X1, Y1, X2, Y2;
	double D;

	for (n = 0; n < Nb1; n++) {
		X = Borders1[n].X + 1; //bug fixed (04/03/2018): +1 for g
		Y = Borders1[n].Y + 1; //bug fixed (04/03/2018): +1 for g
		for (I = 1; I <= 8; I++) {
			D = Borders1[n].Delta[I - 1];
			if (D < 1) {
				j = jB[I];
				X1 = X + ex[I];
				Y1 = Y + ey[I];

				//= = Bouzidi(2001) :
				if (D < 0.5) {
					X2 = X + ex[I] * 2;
					Y2 = Y + ey[I] * 2;
				//== NEW BUG FIXED (2018-05-08): X2==0 or X2==Nx+1
					if (X2 == 0 || X2 == Nx + 1) X2 = X1;
				//-- New bug found (2018-05-16): (X2,Y2) might not be a fluid node!
					if (Fobj[X2][Y2] < 1) {	//avoid accessing unknown data
						X2 = X1;
						Y2 = Y1;
					}

					g[I][X][Y] = 2 * D * g[j][X1][Y1] + (1 - 2 * D) * g[j][X2][Y2];
				}
				else {
					g[I][X][Y] = 0.5L / D * (g[j][X1][Y1] + (2 * D - 1) * g[I][X1][Y1]);
				}

			}
		}
	}

	for (n = 0; n < Nb2; n++) {
		X = Borders2[n].X + 1; //bug fixed (04/03/2018): +1 for g
		Y = Borders2[n].Y + 1; //bug fixed (04/03/2018): +1 for g
		for (I = 1; I <= 8; I++) {
			D = Borders2[n].Delta[I - 1];
			if (D < 1) {
				j = jB[I];
				X1 = X + ex[I];
				Y1 = Y + ey[I];

				//= = Bouzidi(2001) :
				if (D < 0.5) {
					X2 = X + ex[I] * 2;
					Y2 = Y + ey[I] * 2;
				//== NEW BUG FIXED (2018-05-08): X2==0 or X2==Nx+1
					if (X2 == 0 || X2 == Nx + 1) X2 = X1;
				//-- New bug found (2018-05-16): (X2,Y2) might not be a fluid node!
					if (Fobj[X2][Y2] < 1) {	//avoid accessing unknown data
						X2 = X1;
						Y2 = Y1;
					}

					g[I][X][Y] = 2 * D * g[j][X1][Y1] + (1 - 2 * D) * g[j][X2][Y2];
				}
				else {
					g[I][X][Y] = 0.5L / D * (g[j][X1][Y1] + (2 * D - 1) * g[I][X1][Y1]);
				}

			}
		}

	}

}

//**********************************************************************

void Bouzidi_quadratic()
{

	int n, I, j, X, Y, X1, Y1, X2, Y2, X3, Y3;
	double D;

	for (n = 0; n < Nb1; n++) {
		X = Borders1[n].X + 1; //bug fixed (04/03/2018): +1 for g
		Y = Borders1[n].Y + 1; //bug fixed (04/03/2018): +1 for g
		for (I = 1; I <= 8; I++) {
			D = Borders1[n].Delta[I - 1];
			if (D < 1) {
				j = jB[I];
				X1 = X + ex[I];
				Y1 = Y + ey[I];
				X2 = X + ex[I] * 2;
				Y2 = Y + ey[I] * 2;

				//= = Bouzidi(2001):
				if (D < 0.5) {
					X3 = X + ex[I] * 3;
					Y3 = Y + ey[I] * 3;
					if (X2 == 0 || X2 == Nx + 1) X2 = X1;
					if (X3 == -1 || X3 == Nx + 2) X3 = X1;
					if (X3 == 0 || X3 == Nx + 1) X3 = X2;
					if (Fobj[X2][Y2] < 1) {
						X2 = X1;
						Y2 = Y1;
					}
					if (Fobj[X3][Y3] < 1) {
						X3 = X2;
						Y3 = Y2;
					}

					g[I][X][Y] = g[j][X1][Y1] * (1 + 2 * D) * D
						+ g[j][X2][Y2] * (1 - 2 * D) * (1 + 2 * D)
						- g[j][X3][Y3] * (1 - 2 * D) * D;
				}
				else {
					if (X2 == 0 || X2 == Nx + 1) X2 = X1;
					if (Fobj[X2][Y2] < 1) {
						X2 = X1;
						Y2 = Y1;
					}

					g[I][X][Y] = (g[j][X1][Y1]
						- g[I][X1][Y1] * (1 - 2 * D) * (1 + 2 * D)
						+ g[I][X2][Y2] * (1 - 2 * D) * D) / (D*(1 + 2 * D));
				}

			}
		}
	}

	for (n = 0; n < Nb2; n++) {
		X = Borders2[n].X + 1;
		Y = Borders2[n].Y + 1;
		for (I = 1; I <= 8; I++) {
			D = Borders2[n].Delta[I - 1];
			if (D < 1) {
				j = jB[I];
				X1 = X + ex[I];
				Y1 = Y + ey[I];
				X2 = X + ex[I] * 2;
				Y2 = Y + ey[I] * 2;

				//= = Bouzidi(2001):
				if (D < 0.5) {
					X3 = X + ex[I] * 3;
					Y3 = Y + ey[I] * 3;
					if (X2 == 0 || X2 == Nx + 1) X2 = X1;
					if (X3 == -1 || X3 == Nx + 2) X3 = X1;
					if (X3 == 0 || X3 == Nx + 1) X3 = X2;
					if (Fobj[X2][Y2] < 1) {
						X2 = X1;
						Y2 = Y1;
					}
					if (Fobj[X3][Y3] < 1) {
						X3 = X2;
						Y3 = Y2;
					}

					g[I][X][Y] = g[j][X1][Y1] * (1 + 2 * D) * D
						+ g[j][X2][Y2] * (1 - 2 * D) * (1 + 2 * D)
						- g[j][X3][Y3] * (1 - 2 * D) * D;
				}
				else {
					if (X2 == 0 || X2 == Nx + 1) X2 = X1;
					if (Fobj[X2][Y2] < 1) {
						X2 = X1;
						Y2 = Y1;
					}

					g[I][X][Y] = (g[j][X1][Y1]
						- g[I][X1][Y1] * (1 - 2 * D) * (1 + 2 * D)
						+ g[I][X2][Y2] * (1 - 2 * D) * D) / (D*(1 + 2 * D));
				}

			}
		}

	}

}

//**********************************************************************

void Multireflection_BC()	//2018-06-18
{

	int n, I, j, X, Y, X1, Y1, X2, Y2, X3, Y3;
	double D;
	double kt, ks, kf, kstar;
	double k1, k2, k3, k4, k5;

	for (n = 0; n < Nb1; n++) {
		X = Borders1[n].X + 1;
		Y = Borders1[n].Y + 1;
		for (I = 1; I <= 8; I++) {
			D = Borders1[n].Delta[I - 1];
			if (D < 1) {
				j = jB[I];
				X1 = X  + ex[I];
				Y1 = Y  + ey[I];
				X2 = X1 + ex[I];
				Y2 = Y1 + ey[I];
				X3 = X2 + ex[I];
				Y3 = Y2 + ey[I];

				//= = Ginzburg and d'Humieres (2003):
				kt = 3 + 2 * D;
				ks = 1 + D*(6 + 4 * D);
				kf = 2 * (1 + D)*(1 + D);

			//	kstar = 2 * ks / kf;
				kstar = 1 + D / 15 * (3 * (5 + D) + 4 * sqrt(3.)*(1 - D));

				k1 = kstar*kf / ks - 1;
				k2 = 2 - kstar;
				k3 = k1 - kstar*(kt - 2) / ks;
				k4 = 2 - kstar*kt / ks;
				k5 = k1 - k3 - 1;

				if (X2 == 0 || X2 == Nx + 1) {
					X2 = X1;
				//	Y2 = Y1;
				}
				if (X3 == -1 || X3 == Nx + 2) {
					X3 = X1;
				//	Y3 = Y1;
				}
				if (X3 == 0 || X3 == Nx + 1) {
					X3 = X2;
				//	Y3 = Y2;
				}

				if (Fobj[X2][Y2] < 1) {
					X2 = X1;
					Y2 = Y1;
				}
				if (Fobj[X3][Y3] < 1) {
					X3 = X2;
					Y3 = Y2;
				}

				g[I][X][Y] = k1*g[j][X1][Y1] + k2*g[j][X2][Y2] + k3*g[j][X3][Y3]
					+ k4*g[I][X1][Y1] + k5*g[I][X2][Y2];
			}
		}
	}

	for (n = 0; n < Nb2; n++) {
		X = Borders2[n].X + 1;
		Y = Borders2[n].Y + 1;
		for (I = 1; I <= 8; I++) {
			D = Borders2[n].Delta[I - 1];
			if (D < 1) {
				j = jB[I];
				X1 = X + ex[I];
				Y1 = Y + ey[I];
				X2 = X1 + ex[I];
				Y2 = Y1 + ey[I];
				X3 = X2 + ex[I];
				Y3 = Y2 + ey[I];

				//= = Ginzburg and d'Humieres (2003):
				kt = 3 + 2 * D;
				ks = 1 + D*(6 + 4 * D);
				kf = 2 * (1 + D)*(1 + D);

			//	kstar = 2 * ks / kf;
				kstar = 1 + D / 15 * (3 * (5 + D) + 4 * sqrt(3.)*(1 - D));

				k1 = kstar*kf / ks - 1;
				k2 = 2 - kstar;
				k3 = k1 - kstar*(kt - 2) / ks;
				k4 = 2 - kstar*kt / ks;
				k5 = k1 - k3 - 1;

				if (X2 == 0 || X2 == Nx + 1) {
					X2 = X1;
				//	Y2 = Y1;
				}
				if (X3 == -1 || X3 == Nx + 2) {
					X3 = X1;
				//	Y3 = Y1;
				}
				if (X3 == 0 || X3 == Nx + 1) {
					X3 = X2;
				//	Y3 = Y2;
				}

				if (Fobj[X2][Y2] < 1) {
					X2 = X1;
					Y2 = Y1;
				}
				if (Fobj[X3][Y3] < 1) {
					X3 = X2;
					Y3 = Y2;
				}

				g[I][X][Y] = k1*g[j][X1][Y1] + k2*g[j][X2][Y2] + k3*g[j][X3][Y3]
					+ k4*g[I][X1][Y1] + k5*g[I][X2][Y2];

			}
		}

	}

}
