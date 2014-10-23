#include "Runner.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <cstdio>
#include <cmath>

using namespace std;

int main(int argc, char *argv[]){ 

/*

main file, main.cpp, a wrapper which executes the runner script.
 
.. cpp:function:: bool namespaced::theclass::method(int arg1, std::string arg2)

   Describes a method with parameters and types
 
 */
   
   // instantiate Runner object with the handle my_transducer.
   Runner my_transducer(argc, argv);
   
   return 0;
}