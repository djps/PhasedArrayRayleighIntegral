/*
This class specifies the base class TransducerElement to the case where the curvature of the element face is 
spherical. I assume that the focus of this sphere is at the origin, and that the transducer is entirely 
contained in the z<0 plane. 

Note that I use the coordinate system,
    x
    | 
    |
    |_______Z
   /
  /
 /  
y

and angles theta are measured from the +x direction around towards the +y direction.  

Notes on notation:
All variables must be prefixed. 
First prefix indicates scope: c = class, g = global, a = argument, l = loop. 
Second prefix indicates type: d = double, i = integer.  
 
EDITS:
Written on 28.04.2010. 

*/


#ifndef CIRCULAR_PISTON_TRANSDUCER_ELEMENT
#define CIRCULAR_PISTON_TRANSDUCER_ELEMENT

#include "TransducerElement.h"
#include "ComputationalElement.h"
#include <iostream>

class Circular_Piston_Element : public Transducer_Element
{ 
/* Basically just make a flat piston element. Defined by position of centre, normal direction and radius. 
  radius is part of Transducer_Element, while position and normal are specialized.   
 */
   
	public: 
		double position[3]; 
		double normal[3]; 

	public:
		// generate the computational mesh.
		void Generate_Elements(); 

		// Constructors;
		Circular_Piston_Element(double ad_external_radius, int ai_number_of_radials, int ai_number_of_angles, double ad_real_pressure, double ad_imag_pressure, double posX, double posY, double posZ, double normX, double normY, double normZ)
		{
			cd_internal_radius = 0.0; cd_external_radius = ad_external_radius;
			cd_theta1= 0.0  ;  cd_theta2 = PI2;
			ci_number_of_radials =  ai_number_of_radials ;  ci_number_of_angles = ai_number_of_angles; 
			cd_real_pressure = ad_real_pressure ;  cd_imag_pressure = ad_imag_pressure ;
			position[0]=posX;
			position[1]=posY;
			position[2]=posZ;
			normal[0]=normX;
			normal[1]=normY;
			normal[2]=normZ;
			ccomp_array_of_computational_elements = new Computational_Element[ci_number_of_angles*ci_number_of_radials];
			this->Generate_Elements(); // Once everything is set up, call the function to make the elements. 
		}

		// stupid default constructor. Probably will get used though. 
		Circular_Piston_Element()
		{ 
			ccomp_array_of_computational_elements = new Computational_Element[1];
		}

		// stupid default destructor. Doesn't do anything. Reason is that all non-automatic variables are deleted by 
		// the Transducer_Element class. 
		~Circular_Piston_Element()
		{ 
#ifdef DEBUG2
		std::cout << "Deleting Circular_Piston_Element" << std::endl;
#endif  
		}

		void Input_Position(double posX, double posY, double posZ)
		{  
#ifdef DEBUG2
			  std::cout << posX << "  " << posY << "  " << posZ << std::endl;
			  std::cout << position << std::endl;
			  std::cout << sizeof(position[0]) << ",  " << sizeof(posX) << std::endl;
#endif
		  position[0]=posX;
#ifdef DEBUG2
			std::cout << posY << "  " << std::endl;
#endif
		  position[1]=posY;
#ifdef DEBUG2
		  std::cout << posY << "  " << std::endl;
#endif
		  position[2]=posZ;		 
#ifdef DEBUG2
		  std::cout << posZ << "  " << std::endl;
#endif
		}
	
		void Output_Position(double &posX, double &posY, double &posZ)
		{
			posX=position[0];
			posY=position[1];
			posZ=position[2];		 
		}
		
		void Input_Normal(double normX, double normY, double normZ)
		{
			normal[0]=normX;
			normal[1]=normY;
			normal[2]=normZ;		 
		}
		
		void Output_Normal(double &normX, double &normY, double &normZ)
		{
			normX=normal[0];
			normY=normal[1];
			normZ=normal[2];		 
		}

};



#endif


