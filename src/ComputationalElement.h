/*
This is a class that describes a basic computational element for a numerical integration of the Helmholtz equation. 
More precisely, given the wave equation \ddot{p} - c^2 \Delta p = 0, subject to the boundary condition that 
in a plane, p is everywhere 0 except in a small area. In this area it oscillates as p = p_0 exp(-i omega t), 
with p_0 constant over the area. Only outgoing waves are allowed. Then we can easily solve this problem using 
Green's functions. The solution at a point in space X depends on the location of the little nonzero area, the 
normal to the plane in which everything vanishes, the area of the element, omega/c and p_0.  Now, if we have 
loads of these elements, each with their own normal, etc. but all oscillating exactly in phase, then it makes 
no sense to store the value of p_0 for each of these. SO:
     
The class consists of 
Data:
  Position
  Area
  Normal
Methods:
  Usual constructors and assignations
  Pressure output: given a position in space and omega/c, then this will give the pressure at that position, 
  assuming that the input pressure was 1.0. You are free to multiply result by p_0 to get the true pressure.

Notes on notation:
All variables must be prefixed. 
First prefix indicates scope: c = class, g = global, a = argument, l = loop. 
Second prefix indicates type: d = double, i = integer.  

EDITS:
Written on 14.05.2009. 

*/

#ifndef COMPUTATIONAL_ELEMENT
#define COMPUTATIONAL_ELEMENT

#include <cmath>

const double PI2=6.28318530718;

class Computational_Element
{
	private:
		double cd_position[3];
		double cd_normal[3];
		double cd_area;
		
	public:
		Computational_Element(){} // default constructor doesn't do anything.
		
		Computational_Element(double ad_posx,double ad_posy,double ad_posz,double ad_normx,double ad_normy,double ad_normz,double ad_area)
		{
			cd_position[0]=ad_posx;cd_position[1]=ad_posy;cd_position[2]=ad_posz;
			double ld_Norm = sqrt(ad_normx*ad_normx +ad_normy*ad_normy +ad_normz*ad_normz);  
			cd_normal[0]=ad_normx/ld_Norm;cd_normal[1]=ad_normy/ld_Norm;cd_normal[2]=ad_normz/ld_Norm; 
			// Make sure the normal is a unit vector.
			cd_area=ad_area; 
		}

		void Set_Position(double ad_posx,double ad_posy,double ad_posz)
		{
			cd_position[0]=ad_posx;cd_position[1]=ad_posy;cd_position[2]=ad_posz;
		}

		void Set_Normal(double ad_normx,double ad_normy,double ad_normz)
		{
			double ld_Norm = sqrt(ad_normx*ad_normx +ad_normy*ad_normy +ad_normz*ad_normz);  
			cd_normal[0]=ad_normx/ld_Norm;cd_normal[1]=ad_normy/ld_Norm;cd_normal[2]=ad_normz/ld_Norm; 
			// Make sure the normal is a unit vector.
		}

		void Set_Area(double ad_area)
		{
			cd_area = ad_area;
		}

		void Output_Position(double &ad_posx,double &ad_posy,double &ad_posz)
		{
			ad_posx=cd_position[0];ad_posy=cd_position[1];ad_posz=cd_position[2];
		}

		void Output_Normal(double &ad_normx,double &ad_normy,double &ad_normz)
		{
			ad_normx=cd_normal[0];ad_normy=cd_normal[1];ad_normz=cd_normal[2];
		}

		void Output_Area(double &ad_area)
		{
			ad_area=cd_area;
		}

		void Calculate_Pressure(double ad_X, double ad_Y, double ad_Z, double ad_omega_div_c, double &ad_real_pressure, double &ad_imag_pressure);
			
};


#endif
