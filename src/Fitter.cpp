/*
CPP file.
This file contains the basic functions that fit data to a given function. More specifically, we have:
1) A class, which contains a function that returns a positive real value, measuring the deviation between the data point and the fitted function. I guess
that a derivative function may also be useful, although for now 
*/

/* Minimizing a function is pretty tricky. I'm going to use a steepest descent method coupled with a line-search algorithm.
That is, if we have a multidimensional function, we shoot in the direction of its (numerically calculated) gradient, with the shooting length
determined by a 1-dimensional minimization. This means that we need to define a new type of Function, which contains the multi-dimensional parameters, and the gradient,
and is a function of a single variable.
1-dimensional minimization works as follows:
 Bracket an area: figure out which direction to go from the initial point to move downhill. Then expand the region (by doubling) until it encloses a minimum.
 Do a Golden-mean bracketting with this initial region to find the minimum. 
*/

#include "Fitter.h"

// Class Tester things

double Tester::function(int N, double *pars)
{ 
   if (N==2)
   {
       double x = pars[0], y=pars[1];
       return (sin(x)*sin(y));
   }
   else return 0.0;
}



// Class linesearch things;

double *LineSearch::InitialPars, *LineSearch::Gradient;
int LineSearch::NumberOfElements = 0;


double LineSearch::function(int N, double *parameters)
{  // This uses the function that is already inside the Linesearch class.
   // LineSearch is a 1-d function, so N = 1. If not, just quit happily.
   if (N == 1)
   { 
//     cout << "Inside LineSearch " << N << " " << parameters[0] << endl;
     if (NumberOfElements ==0) return 0.0; // i.e. quit if it's not been initialized.
//     double *tempvec = new double [NumberOfElements];
     double lamb = parameters[0];
     for (int p = 0; p < NumberOfElements; p++)
      {  
         tempvec[p] = InitialPars[p] + lamb*Gradient[p];
//	 cout << tempvec[p] << ", ";
      }
//      cout << "Temp " << endl;
     double output =  (MultiDimensionalFunction->function(NumberOfElements, tempvec)); // i.e. use the internal Function's function to determine the result.
//     delete tempvec;
//     cout << "out " <<  output << endl;
     return output;  
   }
   else return 0.0;
}
   
void LineSearch::Initialize(int N, double *init, double *grad)
{ // set up the function
   NumberOfElements = N;
   InitialPars = new double [N];
   Gradient = new double [N];
   for (int p = 0; p<N; p++)
   {
     InitialPars[p] = init[p];
     Gradient[p] = grad[p];
   }
}


void LineSearch::Delete()
{ // delete stuff from memory.
   if (NumberOfElements != 0) 
   {  
     delete InitialPars;
     delete Gradient;
   }
}












double Minimizer(int NumberOfFunctions, Function **FunctionWeWantToMinimize, int NumberOfParameters, double *Parameters)
{ // Minimizes the sum of the Functions, each of which takes in an array of parameters. 
  // Steps  down the gradient; apparently this is not awfully efficient.
  int count =0; // keep track of how many times we do a loop.  
   cout << scientific; 
   cout << setprecision(8);
  if (NumberOfParameters == 1) // then we use the 1-d minimization.
  {
 //   cout << "Here" << endl;
    double fa=0.0,fb=0.0,fc=0.0,ftmp;   // these are the values of the function at the start, middle and end of the region.
    double xa = Parameters[0], lamb = SMALL;  // x0 is the initial position. lamb is the step-size.
    double xb , xc,xtmp;

      // First, we want to determine the starting region. 
      // Evaluate the direction in which the function decreases
      xb = xa+lamb;
      xc = xa-lamb;  // i.e. look to the left and the right of the starting point.
      for (int p = 0; p<NumberOfFunctions; p++)
       {
         fa += FunctionWeWantToMinimize[p]->function(1,&xa);
         fb += FunctionWeWantToMinimize[p]->function(1,&xb);
         fc += FunctionWeWantToMinimize[p]->function(1,&xc);
       }
       // Now we choose a direction. We go to decreasing x if fc < fa. Otherwise, x increases.
       if (fc<fa) lamb = -lamb; // i.e. lamb is now negative.
       if ((fc > fa) && (fb > fa)) 
         {  // then we are at a minimum.
//	   cout << endl;
	   cout << "at min  " << fa << ", " <<  fb << ", " << fc << endl;
	   cout << "at min  " << xa << ", " <<  xb << ", " << xc << endl;
	    Parameters[0] = xa;
 //         cout << "  asaf  " << fa << endl;
	    return fa;
 //         cout << "  asaf  " << fa << endl;
	 }
       xb = xa + lamb; 
       xc = xa + 2.0*lamb;
       // reevaluate the functions
       fb = 0.0;
       fc = 0.0;
       for (int p = 0; p<NumberOfFunctions; p++)
              {
                 fb += FunctionWeWantToMinimize[p]->function(1,&xb);
                 fc += FunctionWeWantToMinimize[p]->function(1,&xc);
              }
       // We stop looping when we have a minimum; that is, when fb < fa and fb < fc. 
       // Since we use the old value of xc for the new value of xb, it is not possible to get to the stage that fb > fa, so we needn't test this. Just stop when fb < fc.
        while ((fb>fc)&&(fabs(xc-xa)<MaxCount))  // i.e. don't let things get too large...
	{
	    lamb *= 2.0;
	    // Use old values of xc, fc as new values of xb, fb. i.e. double the size of the region that you are encompassing.
	    xb = xc;
	    fb = fc; 
	    xc = xa + 2.0*lamb;         
	    fc = 0.0;
            for (int p = 0; p<NumberOfFunctions; p++)
                {
                   fc += FunctionWeWantToMinimize[p]->function(1,&xc);
                }
           if (fabs(xc-xa)>MaxCount)  
	          {
		    cout << "ungracious exit from domain doubling " << xa << " " << xb << " " << xc << endl;
		  }
//          cout << "Determining bracket " << xa << " " << fa << " "  << xb << " " << fb << " " << xc << " " << fc << endl;

 	    
	 }
// OK, now we have the starting region, and the value of the function at the two ends and a central point. 
// Of course, I favour xa slightly, since that was my original input. I favour xa by starting with it. However, to guarantee convergence, you need to look at golden-mean intervals
// on both sides of xc. 
// Start by looking at the point GOLD*(xb-xa) away from xa. If the function at this point is less than fb, then
// we change our region so that the old value of xb is the new boundary point xc, and the new point is the middle point. If the function is greater than fb, then this new point
// replaces xa (i.e. also a boundary point), and we continue. Note that it may happen (very unlikely) that the new point is larger than fa, and there may be a minimum between
// this point and xa. However, I'm going to ignore this case, since there will certainly be a minimum in the new region that I've defined above, and it's probably a better
// minimum.
// Then do the same thing for the point GOLD*(xc-xb) away from xb. (Clearly we are trying to move as little as possible away from xa, to find the closest minimum...)
        while (fabs(xc-xa) > TOL1d)  // this is the condition for halting.
	{ 
	   xtmp = xa + (xb-xa)*GOLD;
	   ftmp = 0.0;
           for (int p = 0; p<NumberOfFunctions; p++)
                {
                   ftmp += FunctionWeWantToMinimize[p]->function(1,&xtmp);
                }
	   if (ftmp >= fb) // then xtmp replaces xa.
	     {
	       xa = xtmp; 
	       fa = ftmp;
	     }
	   else // then xtmp replaces xb, and xb replaces xc.
	    { 
	       xc = xb;
	       fc = fb;
	       xb = xtmp;
	       fb = ftmp;
	    }
	   // Now approach minimum from other side.
	   xtmp = xb + (xc-xb)*GOLD;
	   ftmp = 0.0;
           for (int p = 0; p<NumberOfFunctions; p++)
                {
                   ftmp += FunctionWeWantToMinimize[p]->function(1,&xtmp);
                }
	   if (ftmp >= fb) // then xtmp replaces xc.
	     {
	       xc = xtmp; 
	       fc = ftmp;
	     }
	   else // then xtmp replaces xb, and xb replaces xa.
	    { 
	       xa = xb;
	       fa = fb;
	       xb = xtmp;
	       fb = ftmp;
	    }
	   count++ ;
//	   cout << "X's" << ", " << xa  << ", " <<  xb << ", " <<  xc << std::endl;
//	   cout << count << ", " << fa  << ", " <<  fb << ", " <<  fc << std::endl;
	   if (count > MaxCount) 
	      {
	         cout << "Didnt converge here. " << std::endl;
		 xa = xb;
	      }
  	}
// Finally, we have the minimum. The new value of xa must replace the parameter. The return of this routine is the value of the funtion at xa.
   Parameters[0] = xa;
   return fa;
   }  // End of 1-dimensional minimization.

  else if (NumberOfParameters >1) // i.e. multidimensional minimization.
   {
       cout << "In minimizer" << endl;  
     // Set up some things which are useful right at the start.
       double Change = 100.0;  // Change measures how different the new value is from the old one. i.e. setting it to 100 forces us to go through the loop at least once.
       double fnew;
       double fa = 0.0;
       double *fgrad = new double [NumberOfParameters];
       double *delPars = new double [NumberOfParameters]; // the parameters, shifted slightly to evaluate the gradient.
       LineSearch *LocalLineSearch = new LineSearch [NumberOfFunctions]; // i.e. make a meta-function out of each function.
       Function **PointToFunctions = new Function *[NumberOfFunctions]; // Something to take note of the addresses of the local line search functions
       
       // Initialise the LineSearch variables.
       LineSearch::InitialPars = new double[NumberOfParameters];
       LineSearch::Gradient = new double[NumberOfParameters];
       LineSearch::NumberOfElements = NumberOfParameters;
       for (int p = 0; p < NumberOfFunctions; p++)
       {
         LocalLineSearch[p].MultiDimensionalFunction = FunctionWeWantToMinimize[p];  // now local line search[p] contains a pointer the the p'th function that we want to minimize.
	 PointToFunctions[p] = &(LocalLineSearch[p]);
         LocalLineSearch[p].tempvec = new double[NumberOfParameters];  // declare this memory once rather than each time we call it.
       }
       count = 0;
       while (Change > TOLMultid)
       {   
         // first thing to do is to determine the value of the function at the given point, and its gradient. 
	 cout << count << " Parameters  ";
         for (int p = 0; p < NumberOfParameters; p++)
         {
           delPars[p] = Parameters[p];
	   cout << Parameters[p] << ", " ;
         }
	 cout << endl;
	 fa = 0.0;
         for (int p = 0; p < NumberOfFunctions; p++)  // Value of function
          {
	    fa += FunctionWeWantToMinimize[p]->function(NumberOfParameters,Parameters);
          }   
         for (int q = 0; q < NumberOfParameters; q++) // value of function slightly displaced.
	   {
	      fgrad[q] = 0.0;
	      delPars[q] += SMALL; // shift it a little
	      for (int p = 0; p < NumberOfFunctions; p++)
	      { 
	          fgrad[q] += FunctionWeWantToMinimize[p]->function(NumberOfParameters,delPars);
	      }
	      delPars[q] = Parameters[q]; // shift it back. 
	      fgrad[q] = (fgrad[q]-fa)/SMALL; // now make it into a gradient.
	      // Finally, insert the information into LineSearch class.
	      LineSearch::InitialPars[q] = Parameters[q];
	      LineSearch::Gradient[q] = fgrad[q];	      
	      cout << "Pars " << q << ", " << LineSearch::InitialPars[q] << ", " << Parameters[q] << endl; 
	      cout << "Gradient " << q << ", " << LineSearch::Gradient[q] << ", " << fgrad[q] << endl; 
	   }
	   
	  // Ok, we've determined the gradient, now we want to minimize the function by moving along the gradient. Need a variable that tells us how far along the gradient:
	  double del = 0.0;
	  cout << "Call the 1-d minimizer  " << endl;
	  fnew = Minimizer(NumberOfFunctions, PointToFunctions, 1, &del); 
	  cout << "Before and after " << fa << "  " << fnew << ",  " << del << endl;
	  // fnew is the new value of the function, lamb is the quantity that we've moved along the gradient.
	  for (int q = 0; q < NumberOfParameters; q++) 
	   {
	     Parameters[q] += del*fgrad[q]; // Update the parameters.
	   }
	  Change  = fabs(fnew - fa); // how much has the function changed?
	  count++;
	  cout << count << " Change  " << Change <<  endl;
	  if (count > MaxCount) Change = -1000;
	 }
	 std::cout << "Change is " << Change << std::endl;
// Reclaim some memory.
       for (int p = 0; p < NumberOfFunctions; p++)
         {
            delete LocalLineSearch[p].tempvec;
         }
	 delete PointToFunctions;
	 delete LocalLineSearch;
	 delete delPars;
	 delete fgrad;
	 delete LineSearch::InitialPars;
	 delete LineSearch::Gradient;
	 return  fa;
    }
       
       


  return 0.0;

}



double Evaluator(int NumberOfFunctions, Function **FunctionWeWantToMinimize, int NumberOfParameters, double *Parameters)
{ // Evaluates the sum of the input functions at the parameters.
 double output = 0.0; 
 for (int p  = 0; p< NumberOfFunctions; p++)
  {
     output += FunctionWeWantToMinimize[p] -> function(NumberOfParameters, Parameters);
  }
  return output;
}






























double MinimizerConjGrad(int NumberOfFunctions, Function **FunctionWeWantToMinimize, int NumberOfParameters, double *Parameters)
{ // Minimizes the sum of the Functions, each of which takes in an array of parameters. 
// Uses conjugate gradient method.
  
  int count =0; // keep track of how many times we do a loop.  
   cout << scientific; 
   cout << setprecision(8);
  if (NumberOfParameters == 1) // then we use the 1-d minimization.
  {
    cout << "Here" << endl;
    double fa=0.0,fb=0.0,fc=0.0,ftmp;   // these are the values of the function at the start, middle and end of the region.
    double xa = Parameters[0], lamb = SMALL;  // x0 is the initial position. lamb is the step-size.
    double xb , xc,xtmp;

      // First, we want to determine the starting region. 
      // Evaluate the direction in which the function decreases
      xb = xa+lamb;
      xc = xa-lamb;  // i.e. look to the left and the right of the starting point.
      for (int p = 0; p<NumberOfFunctions; p++)
       {
         fa += FunctionWeWantToMinimize[p]->function(1,&xa);
         fb += FunctionWeWantToMinimize[p]->function(1,&xb);
         fc += FunctionWeWantToMinimize[p]->function(1,&xc);
       }
       // Now we choose a direction. We go to decreasing x if fc < fa. Otherwise, x increases.
       if (fc<fa) lamb = -lamb; // i.e. lamb is now negative.
       if ((fc > fa) && (fb > fa)) 
         {  // then we are at a minimum.
//	   cout << endl;
	       cout << "at min  " << fa << ", " <<  fb << ", " << fc << endl;
	       cout << "at min  " << xa << ", " <<  xb << ", " << xc << endl;
	        Parameters[0] = xa;
 //         cout << "  asaf  " << fa << endl;
	        return fa;
 //         cout << "  asaf  " << fa << endl;
	       }
          //cout << "Values " << fa << " " << fb << " " << fc << endl;
       xb = xa + lamb; 
       xc = xa + 2.0*lamb;
       cout << "in 1-d fitter " << xa << " " << lamb << endl;  
       
       // reevaluate the functions
       fb = 0.0;
       fc = 0.0;
       for (int p = 0; p<NumberOfFunctions; p++)
              {
                 fb += FunctionWeWantToMinimize[p]->function(1,&xb);
                 fc += FunctionWeWantToMinimize[p]->function(1,&xc);
              }
              //cout << "Values " << fa << " " << fb << " " << fc << endl;
              //cout << "X's  " << xa << "  " << xb <<  " " << xc << endl;
       // We stop looping when we have a minimum; that is, when fb < fa and fb < fc. 
       // Since we use the old value of xc for the new value of xb, it is not possible to get to the stage that fb > fa, so we needn't test this. Just stop when fb < fc.
        while ((fb>fc)&&(fabs(xc-xa)<MaxCount))  // i.e. don't let things get too large...
	{
	    lamb *= 2.0;
	    // Use old values of xc, fc as new values of xb, fb. i.e. double the size of the region that you are encompassing.
	    xb = xc;
	    fb = fc; 
	    xc = xa + 2.0*lamb;         
	    fc = 0.0;
            for (int p = 0; p<NumberOfFunctions; p++)
                {
                   fc += FunctionWeWantToMinimize[p]->function(1,&xc);
                }
                cout << "Values " << fa << " " << fb << " " << fc << endl;
                //cout << "X's  " << xa << "  " << xb <<  " " << xc << endl;
                if (fabs(xc-xa)>MaxCount)  
	          {
		    cout << "ungracious exit from domain doubling " << xa << " " << xb << " " << xc << endl;
		  }
//          cout << "Determining bracket " << xa << " " << fa << " "  << xb << " " << fb << " " << xc << " " << fc << endl;

 	    
	 }
   // cout << "*******************************Have found starting region" << endl;
   // cout <<  "X's  " << xa << "  " << xb <<  " " << xc << endl;
// OK, now we have the starting region, and the value of the function at the two ends and a central point. 
// Of course, I favour xa slightly, since that was my original input. I favour xa by starting with it. However, to guarantee convergence, you need to look at golden-mean intervals
// on both sides of xc. 
// Start by looking at the point GOLD*(xb-xa) away from xa. If the function at this point is less than fb, then
// we change our region so that the old value of xb is the new boundary point xc, and the new point is the middle point. If the function is greater than fb, then this new point
// replaces xa (i.e. also a boundary point), and we continue. Note that it may happen (very unlikely) that the new point is larger than fa, and there may be a minimum between
// this point and xa. However, I'm going to ignore this case, since there will certainly be a minimum in the new region that I've defined above, and it's probably a better
// minimum.
// Then do the same thing for the point GOLD*(xc-xb) away from xb. (Clearly we are trying to move as little as possible away from xa, to find the closest minimum...)
        while (fabs(xc-xa) > TOL1d)  // this is the condition for halting.
	{ 
     // cout << "********************************************************INSIDE GOLDEN MEAN BIT *******************************************" << endl;
     // cout << "Values " << fa << " " << fb << " " << fc << endl;
     // cout <<  "X's  " << xa << "  " << xb <<  " " << xc << endl;
	   xtmp = xa + (xb-xa)*GOLD;
	   ftmp = 0.0;
           for (int p = 0; p<NumberOfFunctions; p++)
                {
                   ftmp += FunctionWeWantToMinimize[p]->function(1,&xtmp);
                }
	   if (ftmp >= fb) // then xtmp replaces xa.
	     {
	       xa = xtmp; 
	       fa = ftmp;
	     }
	   else // then xtmp replaces xb, and xb replaces xc.
	    { 
	       xc = xb;
	       fc = fb;
	       xb = xtmp;
	       fb = ftmp;
	    }
	   // Now approach minimum from other side.
	   xtmp = xb + (xc-xb)*GOLD;
	   ftmp = 0.0;
           for (int p = 0; p<NumberOfFunctions; p++)
                {
                   ftmp += FunctionWeWantToMinimize[p]->function(1,&xtmp);
                }
	   if (ftmp >= fb) // then xtmp replaces xc.
	     {
	       xc = xtmp; 
	       fc = ftmp;
	     }
	   else // then xtmp replaces xb, and xb replaces xa.
	    { 
	       xa = xb;
	       fa = fb;
	       xb = xtmp;
	       fb = ftmp;
	    }
	   count++ ;
//	   cout << "X's" << ", " << xa  << ", " <<  xb << ", " <<  xc << std::endl;
//	   cout << count << ", " << fa  << ", " <<  fb << ", " <<  fc << std::endl;
	   if (count > MaxCount) 
	      {
	         cout << "Didnt converge here. " << std::endl;
		 xa = xb;
	      }
  	}
// Finally, we have the minimum. The new value of xa must replace the parameter. The return of this routine is the value of the funtion at xa.
   Parameters[0] = xa;
   return fa;
   }  // End of 1-dimensional minimization.

  else if (NumberOfParameters >1) // i.e. multidimensional minimization.
   {  // This is largely copied from numerical recipes. 
  // It may be tricky to read, because they use short names for their vectors
/* Basically, the method works as follows:
1) start with g = -grad (func), h=g, xi = g, fa = f(initial point).
2) Enter loop : Minimize f in direction xi; update p to new point, output the new value of the function.
3) If no significant change, quit. Else, fa = new value of f, xi = grad(func), gg = g\dot g, dgg = (xi + g)\dot xi.
4) If gg = 0, then either quit or check that it's not a saddle point. Otherwise, gamma = dgg/gg.
5) Finally,  g = -xi, h = g + gamma*h, xi = h. Then back to 2).
*/
      //char lc_tmp;
      //cin >> lc_tmp;
       cout << "In minimizer" << endl;  
     // Set up some things which are useful right at the start.
       double Change = 100.0;  // Change measures how different the new value is from the old one. i.e. setting it to 100 forces us to go through the loop at least once.
       double fnew, gg, dgg, gamma;
       double fa = 0.0;
       double *g = new double [NumberOfParameters]; // This is the "g" vector used by conjugate gradients.
       double *h = new double [NumberOfParameters]; // This is the "h" vector used in conj grad.
       double *xi = new double [NumberOfParameters];
       double *delPars = new double [NumberOfParameters]; // the parameters, shifted slightly to evaluate the gradient.
       LineSearch *LocalLineSearch = new LineSearch [NumberOfFunctions]; // i.e. make a meta-function out of each function.
       Function **PointToFunctions = new Function *[NumberOfFunctions]; // Something to take note of the addresses of the local line search functions
       
       // Initialise the LineSearch variables. This is not part of numerical recipes; this is part of the class-based minimization. 
       LineSearch::InitialPars = new double[NumberOfParameters];
       LineSearch::Gradient = new double[NumberOfParameters];
       LineSearch::NumberOfElements = NumberOfParameters;
       for (int p = 0; p < NumberOfFunctions; p++)
       {
         LocalLineSearch[p].MultiDimensionalFunction = FunctionWeWantToMinimize[p];  // now local line search[p] contains a pointer the the p'th function that we want to minimize.
	   PointToFunctions[p] = &(LocalLineSearch[p]);
         LocalLineSearch[p].tempvec = new double[NumberOfParameters];  // declare this memory once rather than each time we call it.
       }
       count = 0;
// Now start doing what the numerical recipes guys do.
      
 // first thing to do is to determine the value of the function at the given point, and its gradient. 
         for (int p = 0; p < NumberOfParameters; p++)
         {
           delPars[p] = Parameters[p];
         }
   	  fa = Evaluator(NumberOfFunctions,FunctionWeWantToMinimize,  NumberOfParameters, delPars); // Evaluate the function.   
        cout << "Value of function " << fa << endl;
        for (int q = 0; q < NumberOfParameters; q++) // Evaluate the function slightly displaced, i.e. the gradient.
	   {
	      g[q] = 0.0;
	      delPars[q] += SMALL; // shift it a little
            g[q] = Evaluator(NumberOfFunctions,FunctionWeWantToMinimize,  NumberOfParameters, delPars); // Evaluate the function.   
	      delPars[q] = Parameters[q]; // shift it back.
	      g[q] = (fa - g[q])/SMALL; // now make it into the negative of the  gradient.
            h[q] = g[q];
            xi[q] = g[q]; 
	      // Finally, insert the information into LineSearch class.
	      LineSearch::InitialPars[q] = Parameters[q];
	      LineSearch::Gradient[q] = xi[q];	      
	      //cout << "Pars " << q << ", " << LineSearch::InitialPars[q] << ", " << Parameters[q] << endl; 
	      //cout << "Gradient " << q << ", " << LineSearch::Gradient[q] << ", " << xi[q] << endl; 
	   }

// Now we enter the loop.
       while (Change > TOLMultid) 
       {   
  	     double del = 0.0; // This measures how far we move along xi to minimize the function. 
	     //cout << "Call the 1-d minimizer  " << endl;
	     fnew = MinimizerConjGrad(NumberOfFunctions, PointToFunctions, 1, &del); //minimize in the direction of xi. 
	     //cout << "Before and after " << fa << "  " << fnew << ",  " << del << endl;
	     // fnew is the new value of the function, del is the quantity that we've moved along the gradient.
	     for (int q = 0; q < NumberOfParameters; q++) 
	      {
	        Parameters[q] += del*xi[q]; // Update the parameters.
              delPars[q] = Parameters[q]; // 
	      }
	     Change  = fabs(fnew - fa); // how much has the function changed?
           gg = 0.0;
           dgg = 0.0;
           fa = fnew;
          // cout << "New value of error function " << fa << endl;
           for (int q = 0; q < NumberOfParameters; q++) // Evaluate the function slightly displaced, i.e. the gradient.
       	   {
	           xi[q] = 0.0;
	           delPars[q] += SMALL; // shift it a little
                 xi[q] = Evaluator(NumberOfFunctions,FunctionWeWantToMinimize,  NumberOfParameters, delPars); // Evaluate the function.   
	           delPars[q] = Parameters[q]; // shift it back.
	           xi[q] = (xi[q]-fa)/SMALL; // now make it into  the  gradient.
                 gg += g[q]*g[q];
                 dgg += (xi[q]+g[q])*xi[q];
	         }
           if (gg <= TOLMultid) 
              {  // If the gradient is very small, then we should either quit, or do a check whether we are at a saddle.
                 cout << "gradient was very small " << gg << endl;
                 Change = 0.0;  // Just quit.
              }
           else
             { // i.e. if the gradient is nonzero
                 gamma = dgg/gg; 
                 for (int q = 0; q < NumberOfParameters; q++)  
                    {
                       g[q] = -xi[q];
                       h[q] *= gamma; // Multiply by gamma
                       h[q] += g[q];   // shift it by g[q].
                       xi[q] = h[q];
	          // Finally, insert the new xi into LineSearch class.
       	           LineSearch::InitialPars[q] = Parameters[q];
	                 LineSearch::Gradient[q] = xi[q];	     
                    cout << xi[q] << ",  " ; 
                    }
                    cout << endl;
              }  // end of the else. 

	     count++;
	     cout << count << " Change  " << Change <<  endl;
	     if (count > MaxCount) Change = -1000;
	 }
	 std::cout << "Change is " << Change << std::endl;
// Reclaim some memory.
       for (int p = 0; p < NumberOfFunctions; p++)
         {
            delete LocalLineSearch[p].tempvec;
         }
	 delete PointToFunctions;
	 delete LocalLineSearch;
	 delete delPars;
	 delete g;
       delete h;
       delete xi;
	 delete LineSearch::InitialPars;
	 delete LineSearch::Gradient;
	 return  fa;
    }
       
       


  return 0.0;

}


