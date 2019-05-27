//****************************************************
//
// Nonlinear Pendulum Solver
//
// Zachary Cotman
// April 27, 2018
//
// To do:
//   make the solver adaptive
//
// See Project Journal April 23rd for notes on the 4th order Runge-Kutta
//   that is being implemented here
//
// Compare to the version from april 25th!
//
//****************************************************

#include <stdio.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;

//_______________________________________________________
// FUNCTIONS

// damped driven pendulum first-order diff eqs
void ddp_y_dot(double *y_dot, double *y, double t, void *params_ptr);


void RK4_step( double t0, double d_t, 
               double *y, int n, void *params_ptr,
			   void (*y_dot)(double *, double *, double, void *) );		

//_______________________________________________________
// STRUCTURES

// damped driven pendulum parameters
typedef struct
{
	double omega_drive;
	double omega_0;
	double beta;
	double gamma;
	double delta;
}
dpp_params;

// some constants associated with the dpp
// const int n_dpp_params = 5;  // don't need
const int n_dpp_eqs    = 2;

const double pi = 4.*atan(1);

// So... I have parameters, functions, constants,
//   associated with the dpp.
// I would prefer to wrap these all in a class...
// But making that class usable for ANY diff eq...
// And need to be able to pass the class properly without
//   the solver needing to know anything about the class...
// Hrmmmm

//_______________________________________________________
// MAIN
int main()
{
	//********************************************************************
	// The Function:  
	//   phi_ddot = -2*beta*phi_dot - omega_0^2*sin(phi) 
	//              + gamma*omega_0^2*cos(omega_drive*t + delta)
	//
	//   Straight from Taylor's Classical Mechanics, Chapter 12
	//
	// To test the solver play with these parameters!
	//********************************************************************
	dpp_params dpp_params_1;
		dpp_params_1.omega_drive = 2.*pi; // 2 pi
		dpp_params_1.omega_0     = 1.5 * dpp_params_1.omega_drive;
		dpp_params_1.beta        = dpp_params_1.omega_0/4;
		dpp_params_1.gamma       = 1.105; // <- determines if CHAOS ENSUES
		dpp_params_1.delta       = 0;
		
	void *params_ptr = &dpp_params_1;
	
	double ics[n_dpp_eqs] = {0.};
	ics[0] = 0;   // phi_0
	ics[1] = 0;   // phi_dot_0
	
	double d_t   = 0.0001;  // initial time step in seconds
	
	// time equal to a number of periods of the driving force
	double tmax;  
		double driving_period = 2.*pi/dpp_params_1.omega_drive;
		double periods = 100;
		tmax = periods * driving_period;
		
	//*********************************************************************
	// file output and header
	ofstream out;
	out.open("diff_eq_apr_27.dat");
	out << "# omega_drive = " << dpp_params_1.omega_drive << endl
	    << "# omega_0  = "    << dpp_params_1.omega_0     << endl
		<< "# beta  = "       << dpp_params_1.beta        << endl
		<< "# gamma  = "      << dpp_params_1.gamma       << endl
		<< "# delta = "       << dpp_params_1.delta       << endl
		<< "#" << endl;
	out << "phi_0 = "     << ics[0] << endl
	    << "phi_dot_0 = " << ics[1] << endl
		<< "#" << endl;
	out << "# t  theta  theta_dot" << endl;
	
	
	// setup an array for the solver to take some steps before checking
	//   to see if it has moved forward far enough in time to be finished
		
		int n_steps = 1000; 
		int n_eqs = n_dpp_eqs; // CHANGE IF USING OTHER DIFF EQS
	
		double t[n_steps]  = {0.};
		double t0; // the current time
	
		// 2-D array, oh boy! haha
		double x[n_eqs][n_steps] = {0.};
	
		// save the initial conditions to the ends of the arrays
		for(int i = 0; i < n_eqs; ++i)
		{
			x[i][n_steps-1] = ics[i];
		}
	
	//********************************************************************
	// SOLVE!
	// 
	// int counter = 0;
	
	while(t[n_steps-1] < tmax)
	{
		// Get the initial time and "positions" from the end values
		//   of the previous loop of the solver
		//
        // Note: I saved the initial conditions to the ends of the 
		//   arrays so that it works right even for the first loop!
		for(int i = 0; i < n_eqs; ++i)
		{
			x[i][0] = x[i][n_steps-1];
		}
		t[0] = t[n_steps-1];
		
		//*************************************************************
		// To do: Add the "ic" dimensionality
		for(int i = 1; i < n_steps; ++i)
		{
			// current positions and time 
			for(int k = 0; k < n_eqs; ++k)
				{ ics[k] = x[k][i-1]; }
		    t0 = t[i-1];
		  
		    // Take a step forward!
		    RK4_step( t0, d_t, ics, n_eqs, params_ptr, ddp_y_dot );
			
			// assign step to array
			for(int k = 0; k < n_eqs; ++k)
				{ x[k][i] = ics[k]; }
			
			t[i]  = t[n_steps-1] + i*d_t;
		}
		
		//**************************************************************
		// Printing every 100 steps
		for(int i = 0; i < n_steps; ++i)
	    {
			if(i % 100 == 0)
			{
				out << scientific << setprecision(8)
					<< setw(10) << t[i]  << "  ";
				for(int k = 0; k < n_eqs; ++k)
				{
					out << setw(10) << x[k][i] << " ";
				}
				out << endl;
			}
		}
		
		// I was just interested in watching the solver loop through things
		// ++counter;
		// cout << "counter = " << counter 
		//      << "   time = " << t[n_steps-1]    << endl;
	}

	out.close();
	
	return 0;
}

//_______________________________________________________
// FUNCTION DEFINITIONS


//***************************************
// Runge-Kutta 4
//***************************************
// RK4_step( t0, d_t, y, n, params_ptr, dpp_y_dot );
void RK4_step( double t0, double d_t, 
               double *y, int n, void *params_ptr,
			   void (*y_dot)(double *, double *, double, void *) )
{
		// Note: n is the number of first-order diff eqs
		double k1[n] = {0.};
		double k2[n] = {0.};
		double k3[n] = {0.};
		double k4[n] = {0.};
		
		// save the initial y values for reuse
		double y0[n] = {0.};
		for(int i = 0; i < n; ++i)
		  { y0[i] = y[i]; }
		
		// k1
		  for(int i = 0; i < n; ++i)
			{ y[i] = y0[i]; }
		  double t = t0;
		y_dot(k1, y, t, params_ptr);
		
		
		// k2
		  for(int i = 0; i < n; ++i)
			{ y[i] = y0[i] + d_t/2. * k1[i]; }
		  t = t0 + d_t/2.;
	    y_dot(k2, y, t, params_ptr);
		
		
		// k3
		  for(int i = 0; i < n; ++i)
			{ y[i] = y0[i] + d_t/2. * k2[i]; }
		  t = t0 + d_t/2.;
		y_dot(k3, y, t, params_ptr);
		
		
		// k4
		  for(int i = 0; i < n; ++i)
			{ y[i] = y0[i] + d_t * k3[i]; }
		  t = t0 + d_t;
		y_dot(k4, y, t, params_ptr);
		
		// final values for y1 and y2
		for(int i = 0; i < n; ++i)
		{
			y[i] = y0[i] + d_t/6. * ( k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
		}
	
	
    return;
}


//***************************************
// damped driven pendulum
//***************************************
void ddp_y_dot(double *y_dot, double *y, double t, void *params_ptr)
{
	
	double omega_drive = ((dpp_params *) params_ptr)->omega_drive;
	double omega_0     = ((dpp_params *) params_ptr)->omega_0;
	double beta        = ((dpp_params *) params_ptr)->beta;
	double gamma       = ((dpp_params *) params_ptr)->gamma;
	double delta       = ((dpp_params *) params_ptr)->delta;	
		
	y_dot[0] = y[1];
	y_dot[1] = -2*beta*y[1] - omega_0*omega_0*sin(y[0]) + gamma*omega_0*omega_0*cos(omega_drive*t + delta);
		
	return;
}