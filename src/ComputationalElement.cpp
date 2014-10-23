/*
Here we implement any complicated methods associated with ComputationalElement.h

EDITS: Written on 14.05.2009. 

*/

#include "ComputationalElement.h"
#include <iostream>

#define PI 3.141592653589793238

void Computational_Element::Calculate_Pressure(double ad_X, double ad_Y, double ad_Z, double ad_omega_div_c, double &ad_real_pressure, double &ad_imag_pressure){
/* 
 Calculate the pressure at the point X,Y,Z, assuming the input pressure is 1.0  Scale result accordingly. 
 */

   double ld_normal_dot_displacement = cd_normal[0]*(ad_X - cd_position[0]) + cd_normal[1]*(ad_Y - cd_position[1]) + cd_normal[2]*(ad_Z - cd_position[2]); 
   // calculate Norm dot (X - position).
   
   double ld_distance = sqrt((ad_X - cd_position[0])*(ad_X - cd_position[0]) + (ad_Y - cd_position[1])*(ad_Y - cd_position[1]) + (ad_Z - cd_position[2])*(ad_Z - cd_position[2]) ); 
   // calculate distance between source and measuring point.
   
   double ld_cos_phase_shift  = cos(ad_omega_div_c*ld_distance); 
   
   double ld_sin_phase_shift  = sin(ad_omega_div_c*ld_distance);
   
   ad_real_pressure = -(cd_area/(2.0*PI))*(ld_normal_dot_displacement/(ld_distance*ld_distance))*(3.0*ld_cos_phase_shift/ld_distance + ad_omega_div_c*ld_sin_phase_shift);   
   
   ad_imag_pressure = (cd_area/(2.0*PI))*(ld_normal_dot_displacement/(ld_distance*ld_distance))*(ad_omega_div_c*ld_cos_phase_shift - 3.0*ld_sin_phase_shift/ld_distance);        

   char lc_tmp; 
   //std::cout << (cd_area/(2.0*PI)) << " "<< ld_normal_dot_displacement << " "<< ld_distance << " "<< 3.0*ld_cos_phase_shift/ld_distance << " " << ad_omega_div_c*ld_sin_phase_shift << std::endl;
   //std::cout << "press " << ad_real_pressure << " " << ad_imag_pressure << std::endl;
   //std::cin >>  lc_tmp;
}

