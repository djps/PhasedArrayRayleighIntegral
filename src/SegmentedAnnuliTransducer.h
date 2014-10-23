/*
This class defines a  spherical bowl that is divided into annuli and segments.  I assume that the centre of this 
sphere is at the origin, and that the transducer is entirely contained in the z<0 plane.  The elements are forced 
to be equal area, so that to generate them, all we need to do is input the number of annuli and number of segments

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
Written on 29.04.2010

NOTES:
This class inherits almost everything from Optimize_Transducer. The only things that are new are setting up the 
layout of the elements and initializing the pointers to the elements. We also need a radius of curvature, and we
 shall need several methods of laying out the transducers. These go into the Setup_Transducer routine. 
*/

#ifndef CIRCULAR_PISTON_TRANSDUCER
#define CIRCULAR_PISTON_TRANSDUCER

extern "C"{ 
// solve a general least-squares
   void dgelss_(const int *M,const int *N, const int *nrhs, double *A, const int *lda, double *B, const int *ldb, double *S, const double *rcond, int *rank, double *work, const int *lwork, int *info);
// generate a square matrix C = B^T B.
   void dsyrk_(const char *UPLO, const char *TRANS, const int *N,const int *K,const double *ALPHA,const double *A,const int *LDA,const double *BETA,double *C,const int *LDC);
// multiply a vector by a matrix.
   void dgemv_(const char *TRANS, const int *M, const int *N, const double *alpha, const double *cpd_transfer_matrix, const int *K,const double *cpd_desired_result_of_optimization, const int *INC, const double *beta,  double *x, const int *INC2);
   // Factorize a matrix into Q T Q^T.
   void dsytrd_(const char *UPLO, const int *N, double *C, const int *LDA, double *D, double *E, double *TAU, double *WORK, const int *LWORK,  int *INFO);
   void dcopy_(const int *N,const double *D,const int *INC,double *D2,const int *INC2);
   void dsterf_(const int *N, double *D2, double *E2,  int *INFO);
   void dormtr_(const char *SIDE, const char *UPLO,const char *TRANS, const int *M, const int *N, const  double *A, const int *LDA,  const double *TAU,  double *C, const int *LDC, double *WORK, const int *LWORK, int *INFO ) ;
   void dptsv_(const int *N, const int *NRHS, double *D, double *E, double *x, const int *INC, int *INFO );
   double dnrm2_(const int *N, const double *D2, const int *INC);
   double ddot_(const int *N, const double *D2, const int *INC, const double *E2, const int *INC2);
   void dposv_(const char *UPLO, const int *N, const int *NRHS, double *C, const int *LDA, double *x, const int *LDB, int *INFO ); 
}

#include "CircularPistonTransducerElement.h"
#include "TransducerElement.h"
#include "OptimizeTransducer.h"
#include "ComputationalElement.h"
#include <iostream>
#include "array3d.cpp"
#include <complex>

using namespace std;

class Circular_Piston_Transducer : public Optimize_Transducer
{

	private: 
		double cd_radius_of_sphere;
		double cd_radius_of_piston;
		double cd_minimum_separation;
		double cd_radius_of_transducer;
		double cd_radius_of_imaging_aperture;
	
	
	public:
	  
		void Setup_Transducer(string &filename); // Grab information from the file and use it to set up the transducer. 

		void Input_Number_Of_Elements(int ai_number_of_elements)
		{
			ci_number_of_elements = ai_number_of_elements ;
			//std::cout << "I am here " << std::endl;
			//std::cout << cspher_array_of_elements << std::endl;
			delete[]  cte_array_of_elements;
			//std::cout << "I am here  now" << std::endl;
			cte_array_of_elements = new Transducer_Element*[ci_number_of_elements];		
			for (int p = 0; p < ci_number_of_elements;p++) cte_array_of_elements[p] = new Circular_Piston_Element;
		}  
      
};

#endif



