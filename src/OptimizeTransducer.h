/*
This class is the base class for all transducers and contains all the stuff about optimization that we need. 
We have a lot of  transducer elements that will be specialized in derived classes. We also have information 
about the region over which optimization occurs, and the manner of optimization. All these things are indep of the 
actual transducer, so can be contained in here.

Note that I use the coordinate system,
    x
    | 
    |
    |_______Z
   /
  /
 /  
y

and angles theta are measured from the +x direction around towards the +y direction. This header goes in every file, 
even if it isn't relevant.   

Notes on notation:
All variables must be prefixed. 
First prefix indicates scope: c = class, g = global, a = argument, l = loop. 
Second prefix indicates type: d = double, i = integer.  

EDITS: 
Written 28.04.2010.

*/

//#define DEBUG_PSEUDO_INVERSE
//#define DEBUG_OPTIMIZATION

#ifndef OPTIMIZE_TRANSDUCER
#define OPTIMIZE_TRANSDUCER

extern "C"{ 
// solve a general least-squares
   void dgelss_(const int *M,const int *N, const int *nrhs, double *A, const int *lda, double *B, const int *ldb, double *S, const double *rcond, int *rank, double *work, const int *lwork, int *info);
// generate a square matrix C = B^T B.
   void dsyrk_(const char *UPLO, const char *TRANS, const int *N,const int *K,const double *ALPHA,const double *A,const int *LDA,const double *BETA,double *C,const int *LDC);
// multiply a vector by a matrix.
   void dgemv_(const char *TRANS, const int *M, const int *N, const double *alpha, const double *cpd_transfer_matrix, const int *K,const double *cpd_desired_result_of_optimization, const int *INC, const double *beta,  double *x, const int *INC2);
   // Factorize a matrix into Q T Q^T.
   void dsytrd_(const char *UPLO, const int *N, double *C, const int *LDA, double *D, double *E, double *TAU, double *WORK, const int *LWORK,  int *INFO);
   void dsysv_(const char *UPLO, const int *N, const int *NRHS,  double *C, const int *LDA, int* IPIV,  double *x, const int *LDX , double *WORK, const int *LWORK,  int *INFO);
   void dcopy_(const int *N,const double *D,const int *INC,double *D2,const int *INC2);
   void dsterf_(const int *N, double *D2, double *E2,  int *INFO);
   void dormtr_(const char *SIDE, const char *UPLO,const char *TRANS, const int *M, const int *N, const  double *A, const int *LDA,  const double *TAU,  double *C, const int *LDC, double *WORK, const int *LWORK, int *INFO ) ;
   void dptsv_(const int *N, const int *NRHS, double *D, double *E, double *x, const int *INC, int *INFO );
   double dnrm2_(const int *N, const double *D2, const int *INC);
   double ddot_(const int *N, const double *D2, const int *INC, const double *E2, const int *INC2);
   void dscal_(const int *N, const double *minus_lambda, double *lpd_optimized_pressure, const int *INC);
   void daxpy_(const int *N, const double *alpha, const double *x, const int *INC, double *y, const int *INC2);
   void dposv_(const char *UPLO, const int *N, const int *NRHS, double *C, const int *LDA, double *x, const int *LDB, int *INFO ); 
}

#include "TransducerElement.h"
#include "ComputationalElement.h"
#include "CircularPistonTransducerElement.h"
#include <iostream>
#include "Fitter.h"
#include "array3d.cpp"
#include <complex>
#include <cstring>

using namespace std;

class Optimize_Transducer : public Function
{

	public:
		Transducer_Element **cte_array_of_elements; 
		// NOTE this must point to an array of pointers. It's the best way to allow access to arrays of derived elements. 
		int ci_number_of_elements;
		int ci_number_of_angles,ci_number_of_radial_bits; 
		// Each Transducer_Element will probably need two numbers of this type, either referring to r and theta or to x and y. 
		double cd_omega_div_c; 
		// This needs to be somewhere, and will be constant over the transducer so this seems like the best place for it. 
      
		Array3D<double> *ca_real_array_for_optimization; // Must assign this before it is useful. Can only do this when the optimization data are in place.
		Array3D<double> *ca_imag_array_for_optimization; // Must assign this before it is useful. Can only do this when the optimization data are in place.
		Array3D<double> *ca_energy_in_optimization_region; // This is used in the conjugate gradient method. I don't want to have to keep assigning and deleting large arrays. Or even small arrays.
      
      
		// Optimization data.
		double cd_position_of_focus[3]; // I'm assuming that we want a single focus. 
		double cd_value_at_focus; // This is the value of the pressure that we want at the focus. If we make this just a magnitude, then there will be a phase invariance to the whole thing and we cannot find a minimum. So we force this to be real-valued.
		double cd_xrange_of_optimized_region[2]; // The minimum and maximum values of x over which the optimization function looks.      
		double cd_yrange_of_optimized_region[2]; // The minimum and maximum values of y over which the optimization function looks.      
		double cd_zrange_of_optimized_region[2]; // The minimum and maximum values of z over which the optimization function looks.
		int ci_x_resolution_optimization; // The number of steps in the x-direction of optimization.
		int ci_y_resolution_optimization;
		int ci_z_resolution_optimization;
      
		double *cpd_optimization_parameters; // Things like the width of exponential premultipliers, the weighting of the maximization of pressure, etc. 
		int ci_number_of_optimization_parameters; 

		// These are quite useful and aren't really optimization parameters. 
		double cd_desired_width_of_focus,cd_scale_factor_for_wide_focus,cd_bias_factor_for_wide_focus; 
      
      
	public:
	
		Optimize_Transducer()
		{ 
			// default. Not very interesting. 
			cte_array_of_elements = new Transducer_Element*[1];          
		}
      
		~Optimize_Transducer()
		{  
#ifdef DEBUG
            std::cout << "Deleting Optimize_Transducer" << std::endl;
#endif 
			for (int p = 0; p < ci_number_of_elements; p++) delete cte_array_of_elements[p];
			delete[] cte_array_of_elements;
		}
      
		virtual void Setup_Transducer(string &filename) = 0; // Pure virtual. Stops people initializing this class.  

		void Setup_Optimization(string &filename); // Grab information from the file and use it to set up the optimization stuff.

		void Setup_Pressures(); // Once the transducer and optimization data are in, generate an array of  3-d arrays of pressures over the region of interest. One 3-d array for each transducer.  

		void Generate_Initial_Pressures(double *parameters); // Uses phase matching to try generate the best initial parameters. Otherwise you have to start at dumb ones that the user has input. 

		void Generate_Pressures_Using_Pseudoinverse(double *lpd_pressures) ; // Uses pseudoinverse to generate the initial pressures. Needs Setup_Optimization and Setup_Pressures to have completed. 

		double function(int N, double *parameters); // This is the function that is called during optimization. The parameters will be the real and imaginary parts of all the pressures. N should be twice the number of transducer elements. 
			
		// Functions related to the fancy1 method of optimization, although not exclusively.
		void Generate_Desired_Output_Vector_OLD(double *aa_desired_output_vector);
		void Generate_Desired_Output_Vector(double *aa_desired_output_vector);
		void Generate_Spatial_Vector_For_Optimization(double *apd_spatial_vector_x, double *apd_spatial_vector_y, double *apd_spatial_vector_z, double *apd_position_of_focus);
		void Generate_Transfer_Matrix(double *apd_transfer_matrix);
		void Generate_Pressures_Using_Fancy_Method1(double *apd_pressures); 
		void Generate_Pressures_Using_Fancy_Method2(double *apd_pressures);
		void Generate_Pressures_Using_Fancy_Method3(double *apd_pressures);
		void Generate_Pressures_Using_Piston_Time_Reversal(double *apd_pressures);   
            
		// New type of input/output.
		virtual void Input_Number_Of_Elements(int ai_number_of_elements){
			ci_number_of_elements = ai_number_of_elements ;
			//std::cout << "I am here " << std::endl;
			//std::cout << cspher_array_of_elements << std::endl;
			delete[]  cte_array_of_elements;
			//std::cout << "I am here  now" << std::endl;
			cte_array_of_elements = new Transducer_Element*[ci_number_of_elements];		
		}
      
		void Output_Number_Of_Elements(int &ai_number_of_elements){
			ai_number_of_elements = ci_number_of_elements ;         
		}

		void Input_Omega_Div_C(double ad_omega_div_c)  {
			cd_omega_div_c= ad_omega_div_c;
		}

		void Output_Omega_Div_C(double &ad_omega_div_c) {
			ad_omega_div_c= cd_omega_div_c;
		}

		 
		// Area is only dependent on the computational elements, so Output_Area need not be virtual 
		// in Transducer_Element - the same function applies regardless of element shape.
		void Output_Area(void){	
			double area = 0.0;
			for (int ind = 0; ind < ci_number_of_elements; ind++){
				area += cte_array_of_elements[ind]->Output_Area();
			}
			std::cout << "Area " << area << std::endl;
		}
       
      
		// Inputs that just call routines implemented elsewhere.
		void Input_Pressure(double ad_real_pressure, double ad_imag_pressure, int ai_index_of_element) { 
			if (ai_index_of_element < ci_number_of_elements) 
			{
				cte_array_of_elements[ai_index_of_element]->Input_Pressure(ad_real_pressure, ad_imag_pressure);
			}
		}

		void Input_Radii(double ad_internal_radius, double ad_external_radius, int ai_index_of_element) { 
			// NOTE that even though these are labelled radii, they are just numbers and are used as x,y 
			// coords in the rectangular elements.  
			if (ai_index_of_element < ci_number_of_elements){
				cte_array_of_elements[ai_index_of_element]->Input_Radii(ad_internal_radius, ad_external_radius) ;
			}
		}

		void Input_Angles(double ad_theta1, double ad_theta2, int ai_index_of_element)
		{ 
			// NOTE that even though these are labelled angles, they are just numbers and are used as x,y 
			// coords in the rectangular elements.  
			if (ai_index_of_element < ci_number_of_elements){
				cte_array_of_elements[ai_index_of_element]->Input_Angles(ad_theta1, ad_theta2);
			}
		}

		void Input_Element_Numbers(int ai_number_of_radials, int ai_number_of_angles, int ai_index_of_element)
		{ 
			// NOTE that even though these are labelled radius and angle, they are just numbers and are used as 
			// numbers of x and y coords in rect case.  
			if (ai_index_of_element < ci_number_of_elements){
			cte_array_of_elements[ai_index_of_element]->Input_Element_Numbers(ai_number_of_radials, ai_number_of_angles);
			}
		}

		// Outputs that just call routines implemented elsewhere.
		void Output_Pressure(double &ad_real_pressure, double &ad_imag_pressure, int ai_index_of_element)
		{ 
			if (ai_index_of_element < ci_number_of_elements) 
			{
				cte_array_of_elements[ai_index_of_element]->Output_Pressure(ad_real_pressure, ad_imag_pressure);
			}
		}

		void Output_Radii(double &ad_internal_radius, double &ad_external_radius, int ai_index_of_element)
		{ 
			if (ai_index_of_element < ci_number_of_elements) 
			{
				cte_array_of_elements[ai_index_of_element]->Output_Radii(ad_internal_radius, ad_external_radius) ;
			}
		}

		void Output_Angles(double &ad_theta1, double &ad_theta2, int ai_index_of_element)
		{ 
			if (ai_index_of_element < ci_number_of_elements) 
			{
				cte_array_of_elements[ai_index_of_element]->Output_Angles(ad_theta1, ad_theta2);
			}
		}

		void Output_Element_Numbers(int &ai_number_of_radials, int &ai_number_of_angles, int ai_index_of_element)
		{ 
			if (ai_index_of_element < ci_number_of_elements) 
			{
				cte_array_of_elements[ai_index_of_element]->Output_Element_Numbers(ai_number_of_radials, ai_number_of_angles);
			}
		}

		// Generate computational elements.
		void Generate_Single_Element(int ai_index_of_element)
		{
			if (ai_index_of_element < ci_number_of_elements) 
			{
				cte_array_of_elements[ai_index_of_element]->Generate_Elements();
			}
		}

		void Generate_All_Elements()
		{
			for (int ai_index_of_element = 0; ai_index_of_element < ci_number_of_elements; ai_index_of_element++) 
			{
				cte_array_of_elements[ai_index_of_element]->Generate_Elements();
			}
		}

		void Calculate_Pressure(double ad_X, double ad_Y, double ad_Z,  double &ad_real_pressure, double &ad_imag_pressure)
		{ 
			// Calculate the pressure in space.
			ad_real_pressure = 0.0;
			ad_imag_pressure = 0.0;
			double ld_real_pressure = 0.0;
			double ld_imag_pressure = 0.0;                  
			for (int ai_index_of_element = 0; ai_index_of_element < ci_number_of_elements; ai_index_of_element++) 
			{
				cte_array_of_elements[ai_index_of_element]->Calculate_Pressure(ad_X, ad_Y, ad_Z, cd_omega_div_c, ld_real_pressure, ld_imag_pressure);
				ad_real_pressure +=ld_real_pressure ; 
				ad_imag_pressure +=ld_imag_pressure ; 
				//std::cout << ld_real_pressure  << " "  << ld_imag_pressure  << std::endl;
				//cspher_array_of_elements[ai_index_of_element].Output_Pressure(ld_real_pressure, ld_imag_pressure);
			} 
			//std::cout << "Pressure calculated at point " << ad_X << ", "<< ad_Y << ", "<< ad_Z << " was  " << ad_real_pressure << " + i" <<       ad_imag_pressure << std::endl;    
		}

		void Output_All_Pressures()
		{ 
			// Output the pressure at the transducer face. 
			double ld_real_pressure, ld_imag_pressure;
			std::cout << "*****************PRESSURES FOR EACH ELEMENT ****************" << std::endl;
			for (int li_index_of_element = 0 ;li_index_of_element   < ci_number_of_elements; li_index_of_element++ ) 
			{
				this->Output_Pressure(ld_real_pressure, ld_imag_pressure, li_index_of_element);
				std::cout <<  li_index_of_element << "  " <<    ld_real_pressure << "  " <<  ld_imag_pressure << std::endl;
			}
		}
        
};

#endif



