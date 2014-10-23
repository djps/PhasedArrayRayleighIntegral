/*
This is the base class for the transducer element. Each independent element in the true physical system 
corresponds to one object of this class. Each object will be broken into a load of computational_elements 
for the purpose of calculating the pressure.   This class is just a base class, because I want to make it 
possible to have unfocused, cylindrically focused and spherically focused elements, each of which will be 
a different class with its own way of defining computational_elements. 
However, to limit things a little, I assume that each element is a part of an annulus (possibly with 0 
internal radius). So each can be defined by the internal and external radius and 2 angles. 
NOTE Actually, that's not quite true. Although I name the variables internal_ and external_radius, really 
they're just doubles, as are the angles. They are not used within this class - they are only input or 
output. So you can actually just 

Note that I use the coordinate system,
    x
    | 
    |
    |_______Z
   /
  /
 /  
y

and angles theta are measured from the +x direction around towards the +y direction.  In general, I'll have 
the focus at the origin and the transducer in -z region, but this is not necessary at the moment. What matters 
here is just the coordinate system.

Notes on notation:
All variables must be prefixed. 
First prefix indicates scope: c = class, g = global, a = argument, l = loop. 
Second prefix indicates type: d = double, i = integer.  

EDITS:
Written on 14.05.2009. 

*/

#ifndef TRANSDUCER_ELEMENT
#define TRANSDUCER_ELEMENT

#include "ComputationalElement.h"
#include <new>
#include <exception>
#include <iostream>
#include <cstdlib>

class Transducer_Element
{
	protected:
	
		double cd_internal_radius, cd_external_radius;
		double cd_theta1, cd_theta2;
		int ci_number_of_radials, ci_number_of_angles; // number of bits in which each direction is divided.
		double cd_real_pressure, cd_imag_pressure;
		Computational_Element *ccomp_array_of_computational_elements;

	public:      
   
		// Have a pure virtual generator to stop people instantiating this, and to force them to write a Generate_Elements() function.  
		virtual void Generate_Elements() =0;

		// Calculate the pressure at a given point, for a given value of omega/c:
		void Calculate_Pressure(double ad_X, double ad_Y, double ad_Z, double ad_omega_div_c, double &ad_real_pressure, double &ad_imag_pressure)
		{
			// Note: We add up the pressure contributions from each of the computational elements assumning an input pressure 
			// of 1.0. Then we multiply by the true pressure. 
			double ld_real_pressure = 0.0;
			double ld_imag_pressure = 0.0;
			double ld_temp_real = 0.0;
			double ld_temp_imag = 0.0;
			int li_count = 0;
			char lc_tmp; 
			for (int li_r = 0; li_r < ci_number_of_radials; li_r++)
			{
				for (int li_theta = 0; li_theta < ci_number_of_angles; li_theta++)
				{
					ccomp_array_of_computational_elements[li_count].Calculate_Pressure(ad_X, ad_Y, ad_Z, ad_omega_div_c,ld_temp_real,ld_temp_imag); // calculate pressure for each of the comp. elements.
					ld_real_pressure += ld_temp_real; 
					ld_imag_pressure += ld_temp_imag;                
					li_count++;
				}
			}
			ad_real_pressure = cd_real_pressure*ld_real_pressure - cd_imag_pressure*ld_imag_pressure ;   
			ad_imag_pressure = cd_real_pressure*ld_imag_pressure + cd_imag_pressure*ld_real_pressure ;
			/*
			 std::cout << "TransdElement  " <<   ad_real_pressure <<"  " << ad_imag_pressure << std::endl;   
			 std::cout << "TransdElement  " <<   cd_real_pressure <<"  " << cd_imag_pressure << std::endl;   
			 std::cout << "TransdElement  " <<   ld_real_pressure <<"  " << ld_imag_pressure << std::endl;
			 std::cout << "Number of elements " << ci_number_of_radials   <<"  " <<ci_number_of_angles << std::endl;
			 std::cin >> lc_tmp;
			*/
		}      

		// Constructors;
		Transducer_Element(double ad_internal_radius, double ad_external_radius, double ad_theta1,  double ad_theta2, int ai_number_of_radials, int ai_number_of_angles, double ad_real_pressure, double ad_imag_pressure)
		{ // No real need for a constructor like this, considering this class is pure virtual.
			cd_internal_radius = ad_internal_radius; cd_external_radius = ad_external_radius;
			cd_theta1= ad_theta1  ;  cd_theta2 = ad_theta2;
			ci_number_of_radials =  ai_number_of_radials ;  ci_number_of_angles = ai_number_of_angles; 
			cd_real_pressure = ad_real_pressure ;  cd_imag_pressure = ad_imag_pressure ;          
			ccomp_array_of_computational_elements = new Computational_Element[ci_number_of_angles*ci_number_of_radials];
		}

		virtual ~Transducer_Element() 
		// virtual destructor, ensures that if I delete objects of classes derived from Transducer_Element, it deletes correctly.
		{
			delete[] ccomp_array_of_computational_elements;
		}

		// pretty stupid, but just because otherwise the destructor may delete something that is NULL.
		Transducer_Element() 
		{
			ccomp_array_of_computational_elements = new Computational_Element[1];
		}

		// functions to read and write the pressures.
		void Input_Pressure(double ad_real_pressure, double ad_imag_pressure)
		{
			cd_real_pressure = ad_real_pressure ;  cd_imag_pressure = ad_imag_pressure ;
		}

		void Input_Radii(double ad_internal_radius, double ad_external_radius )
		{     
			cd_internal_radius = ad_internal_radius;
			cd_external_radius = ad_external_radius;
		}

		void Input_Angles(double ad_theta1, double ad_theta2)
		{     
			cd_theta1 = ad_theta1;
			cd_theta2 = ad_theta2;
		}
		
		// This resizes the mesh. But it doesn't re-generate the mesh. Basically, 
		// I want to allow the user to re-generate the mesh when he or she feels like it, 
		// not only by changing the inputs. And if the user has to explicitly call the 
		// re-generator, it's silly to include it in other routines, because then it will be used too often. 
		void Input_Element_Numbers(int ai_number_of_radials, int ai_number_of_angles)
		{ 
			ci_number_of_radials = ai_number_of_radials;
			ci_number_of_angles = ai_number_of_angles;
		 
		#ifdef DEBUG
		 //std::cout << ci_number_of_radials << "  " << ai_number_of_radials << std::endl; 
		 //std::cout << ci_number_of_angles << "  " << ai_number_of_angles << std::endl;
		#endif 

			delete[] ccomp_array_of_computational_elements;
			try
			{
		//std::cout << "About to try " << std::endl;
				ccomp_array_of_computational_elements = new  Computational_Element[ci_number_of_angles*ci_number_of_radials];		
			}
			catch (std::exception& e)
			{
				std::cout << e.what() << std::endl;
				char tmp;
				std::cin >> tmp;
				exit(1);
			}
			catch (...){
				std::cout << "Failed to allocate memory. Not sure why. " << std::endl;
				char tmp;
				std::cin >> tmp;
				exit(1);
			}
		}

		// NB: This outputs the value of the pressure on the transducer. Not the more useful pressure, i.e. 
		// the one at a given point in space. 
		void Output_Pressure(double &ad_real_pressure, double &ad_imag_pressure)
		{          
			ad_real_pressure = cd_real_pressure ;  ad_imag_pressure = cd_imag_pressure ;
			//std::cout << "Pressure on the transducer element "<< cd_real_pressure << "  " << cd_imag_pressure << std::endl;
		}

		void Output_Radii(double &ad_internal_radius, double &ad_external_radius )
		{     
			ad_internal_radius = cd_internal_radius;
			ad_external_radius = cd_external_radius;
		}

		void Output_Angles(double &ad_theta1, double &ad_theta2)
		{     
			ad_theta1 = cd_theta1;
			ad_theta2 = cd_theta2;
		}

		void Output_Element_Numbers(int &ai_number_of_radials, int &ai_number_of_angles)
		{  
			ai_number_of_radials = ci_number_of_radials;
			ai_number_of_angles = ci_number_of_angles;         
		}

		// Calculates the centroid position of this element. This may seem a little dumb, but I like it. 
		void Output_Centroid(double &ad_x, double &ad_y, double &ad_z)
		{ 
			double ld_x = 0.0;
			double ld_y = 0.0;
			double ld_z = 0.0;
			ad_x = 0.0; ad_y = 0.0; ad_z = 0.0;
			for (int li_ind = 0; li_ind < (ci_number_of_radials*ci_number_of_angles); li_ind++)
			{
				ccomp_array_of_computational_elements[li_ind].Output_Position(ld_x,ld_y,ld_z); 
				ad_x += ld_x;
				ad_y += ld_y;
				ad_z += ld_z;
			}
			ad_x /= double(ci_number_of_radials*ci_number_of_angles);
			ad_y /= double(ci_number_of_radials*ci_number_of_angles);
			ad_z /= double(ci_number_of_radials*ci_number_of_angles);
		} 

		// takes in an index and outputs the position of the associated computational element.
		void Output_Position_Of_Computational_Element(double &ad_x, double &ad_y, double &ad_z, int ai_index_of_comp_element)
		{ 
			if (ai_index_of_comp_element < ci_number_of_angles*ci_number_of_radials) ccomp_array_of_computational_elements[ai_index_of_comp_element].Output_Position(ad_x,ad_y,ad_z) ; 
		} 

		// Calculates the area element. This may seem a little dumb, but I like it. 
		double Output_Area(void)
		{ 
			double ld_area = 0.0;
			double rd_area = 0.0;
			for (int li_ind = 0; li_ind < (ci_number_of_radials*ci_number_of_angles); li_ind++)
			{
				ccomp_array_of_computational_elements[li_ind].Output_Area(ld_area); 
				rd_area += ld_area;
			}
			return rd_area;
		} 
      
};

#endif

