/*
OK, now we start getting into the less generic stuff. This class does the following:

1) Read in data from a file and initialize the transducer. 
2) Read in data from a file and initialize the optimization routine (i.e. figure out the weighting function).
3) Run through the optimization. 
4) Output pressure to a data file at a variety of points.

This is all handled by passing in an argv array from the command line. The arguments are as follows.
argv[0]: Program name.
argv[1]: Transducer type: Type 1 is a spherical bowl with circular elements.
argv[2]: Filename containing info for outputting
argv[3]: Filename containing info for optimizing
argv[4]: Filename containing info for setting up the transducer
argv[5]: Filename to which pressure map in space will be saved.
argv[6]: Filename to which the pressure on each element will be saved.
argv[7]: Optimization option: 1 = pseudoinverse
argv[8]:
argv[9]:
argv[10]:
argv[11]:
argv[12]:
argv[13]:
argv[14]:

*/

#ifndef RUNNER
#define RUNNER

#define SAVE_THRESHOLD 0.001

#include "CircularPistonTransducerElement.h"
#include "TransducerElement.h"
#include "ComputationalElement.h"
#include "CircularPistonTransducer.h"
#include "OptimizeTransducer.h"
#include "Fitter.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <cstdio>
#include <cmath>

using namespace std;

const int EXPECTED_ARGS = 8;

class Runner
{ 
	private: 
		string Outputter;// = "Outputter.dat";
		string Optimizer;// = "Optimizer.dat";
		string Transducer;// = "Transducer.dat";
		Function **lf_pointer_to_array_of_functions ; // Note: only 1 transducer over which to optimize things.
		Optimize_Transducer *my_transducer;  // Pointer to the base class; any actual transducer is derived from Optimize_Transducer.

	public:
		Runner(int argc, char *argv[]);
		void Normalize_Pressures(double *pressures, int N);
}; 

#endif


