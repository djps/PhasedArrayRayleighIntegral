#ifndef FUNCTION
#define FUNCTION

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdio>
#include <cmath>
using namespace std;


//#include "runner.hpp"
/*
HPP file. 
1) Defines the class "Function", which contains a function that gives a positive real number. This function is what we are minimizing; doesn't really need 
to be positive, as long as it is what we're minimizing.
2) Defines the routine which does the minimization. Sorry these comments aren't very good. I'm not paying much attention.
*/

class Function
{
   public: 
     virtual double function(int N, double *parameters) = 0; // pure virtual.
};

class LineSearch : public Function
{ // A meta-function which is used in moving a multidimensional function along its gradient. Essentially, we map the multidimensional function to a 1-dimensional function.
  public:
   static double *InitialPars, *Gradient;
   static int NumberOfElements;
   double *tempvec; //useful for calculating.
   Function *MultiDimensionalFunction; // i.e. a pointer to the multidimensional function.
   
   virtual double function(int N, double *parameters); // the LineSearch's function.
   
   void Initialize(int N, double *init, double *grad); // set up the function
   
   void Delete(); // delete stuff from memory.
};


class  Tester : public Function
{ // a straightforward class which is just an implementation of a Function
   virtual double function(int N, double *parameters);
};



// Useful constants.
const double SMALL = 1.e-10;
const double GOLD = 0.38197; // the golden mean.
const int MaxCount = 1000;
const double TOL1d = 1.e-16; // For 1-d minimization  -- some day later, should allow the Function itself to have its own Tolerance. 
const double TOLMultid = 1.e-16; // For multi-d minimization

double Minimizer(int NumberOfFunctions, Function **FunctionWeWantToMinimize, int NumberOfParameters, double *Parameters);

double MinimizerConjGrad(int NumberOfFunctions, Function **FunctionWeWantToMinimize, int NumberOfParameters, double *Parameters);

double Evaluator(int NumberOfFunctions, Function **FunctionWeWantToMinimize, int NumberOfParameters, double *Parameters);



#endif
