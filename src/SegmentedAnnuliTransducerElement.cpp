/*
Implement the generate elements function in here.  


This requires a bit of geometry. We have a vector that is normal to the piston. Therefore the vector a formed by cross-producting normal with the x-direction unit vector gives us a vector  in the plane of the piston (unless this vanishes, in which case, normal = x direction, so a can be set to y-direction.) Then b = a X n will give a second vector, orthogonal to both. Then a,b give us coordinate vectors in the plane of the piston. Let theta = 0 correspond to direction a, and increment theta in the positive b direction. We should normalize both a and b to the radius of the piston, because then we can just add cosines/sines of them to get the required location.  

EDITS:
Written on 28.04.2010. 


*/

#include "CircularPistonTransducerElement.h"

void Circular_Piston_Element::Generate_Elements()
{// OK, this generates the computational elements. We need the right number of them, they must be in the right place and have the right normal. Go for it!
   // Note that "position" of the element will be the central point. Just seems easiest.
	double a[3], b[3];
	if ((fabs(normal[1])>1.e-9)||(fabs(normal[2])>1.e-9)){ // i.e. check that it doesn't point in the x-direction. 
		a[0] = 0.0; a[1] = normal[2]; a[2] = -normal[1];
		b[0] = -(normal[1]*normal[1]+normal[2]*normal[2]); b[1] = normal[0]*normal[1]; b[2] = normal[0]*normal[2];
	}
	else{ 
		std::cout << normal[0] << "  "<< normal[1] << "  "<< normal[2] << std::endl; 
		
		a[0] = 0.0; a[1] =1.0; a[2] = 0.0;
		b[0] = 0.0; b[1] =0.0; b[2] = 1.0;
	}
	// normalise to the length of the radius. 
	double tmp;
	tmp = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
	a[0]*=cd_external_radius/tmp;a[1]*=cd_external_radius/tmp;a[2]*=cd_external_radius/tmp;
	tmp = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
	b[0]*=cd_external_radius/tmp;b[1]*=cd_external_radius/tmp;b[2]*=cd_external_radius/tmp;
#ifdef DEBUG
	//std::cout<<"Hello" << std::endl;
	if ((a[0]==0.0)&&(b[0]==0)){
		std::cout << normal[0] << "  "<< normal[1] << "  "<< normal[2] << std::endl; 
		std::cout << a[0] << "  "<< a[1] << "  "<< a[2] << std::endl; 
		std::cout << b[0] << "  "<< b[1] << "  "<< b[2] << std::endl;
	}
#endif 
	double ld_dtheta = (cd_theta2-cd_theta1)/double(ci_number_of_angles);
   	char lc_tmp;
   	double ld_x, ld_y, ld_z; // useful variables. 
   	int li_count = 0; // counter variable. 
   	for (int li_r = 0; li_r < ci_number_of_radials; li_r++){
      for (int li_theta = 0; li_theta < ci_number_of_angles; li_theta++)
      {
		  ld_x = (li_r+0.5)/double(ci_number_of_radials) * (a[0]*cos(cd_theta1+(li_theta+0.5)*ld_dtheta) + b[0]*sin(cd_theta1+(li_theta+0.5)*ld_dtheta) );
		  ld_y = (li_r+0.5)/double(ci_number_of_radials) * (a[1]*cos(cd_theta1+(li_theta+0.5)*ld_dtheta) + b[1]*sin(cd_theta1+(li_theta+0.5)*ld_dtheta) );
		  ld_z = (li_r+0.5)/double(ci_number_of_radials) * (a[2]*cos(cd_theta1+(li_theta+0.5)*ld_dtheta) + b[2]*sin(cd_theta1+(li_theta+0.5)*ld_dtheta) );
		  ccomp_array_of_computational_elements[li_count].Set_Position(position[0]+ld_x,position[1]+ld_y,position[2]+ld_z);
          ccomp_array_of_computational_elements[li_count].Set_Normal(normal[0],normal[1],normal[2]); // Easy - it's a piston, so it's got a simple normal.
		  ccomp_array_of_computational_elements[li_count].Set_Area(ld_dtheta*(cd_external_radius/double(ci_number_of_radials))*(li_r+0.5)*(cd_external_radius/double(ci_number_of_radials)));  // i.e. dtheta * dr * r. Note there is no curvature.    
         li_count++;
         //std::cout << "Area " << ld_dtheta*ld_dr*(cd_internal_radius+(li_r+0.5)*ld_dr) << std::endl;
         //std::cout << "Position " <<  ld_x << ", "<<  ld_y << ", "<<  ld_z << ", " << std::endl;
      }      
   }
   //std::cin >> lc_tmp;
  
   
}




