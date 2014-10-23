#include "CircularPistonTransducer.h"
 
int Consistency_Checker(double *x, double *y, double *z, double r, double separation, int N){
/* Returns 1 for consistent, 0 if inconsistent. Check that the arrays of positions 
  (x,y,z) [0..N-1] are consistent, in the sense that pistons placed at the positions 
  and having radius r are at least separation away from each other. Actually, I cheat. 
  The geometry of measuring the distance between the edges of two pistons on the surface 
  of a transducer is a little tricky, so all I do is ensure that the centres are more 
  than 2r+separation from each other. If this is the case, then they are at least 
  2r+separation from each other when properly placed, i.e. noncoplanar. 
 */
	double dist2 = (2*r+separation)*(2*r+separation);
	int ret = 1;
	for (int i = 0; i<N; i++){
		for (int j = 0; j<N; j++){
			if (i!=j){
				if (((x[i]-x[j])*(x[i]-x[j]) +(y[i]-y[j])*(y[i]-y[j]) +(z[i]-z[j])*(z[i]-z[j])) < dist2) {
					std::cout << i << "  " << j << std::endl;
					std::cout << x[i] << "  " << y[i] << "  " << z[i] << std::endl;
					std::cout << x[j] << "  " << y[j] << "  " << z[j] << std::endl;
					ret=0;
				} 
			}
		} 
	} 
	return ret;
}



void Randomly_Move_Pistons_Normal(double *x, double *y, double *z, double r_sphere, int N,double amount){ 
/* Shift each transducer centre by a random offset in both x and y, normally distributed with a std 
   deviation of 1, then multiplied by "amount" to get a distance. I use the Box-Muller transform 
   to get normally distributed values.
 */
	srand(time(NULL));
	double xShift,yShift;
	double rand1,rand2;
	for (int p = 0; p<N; p++){
		rand1 = double(rand()%10000)/double(10000);
		rand2 = double(rand()%10000)/double(10000);
		xShift = amount*sqrt(-2.0*log(rand1))*cos(PI2*rand2);
		yShift = amount*sqrt(-2.0*log(rand1))*sin(PI2*rand2);
		x[p] += xShift;
		y[p] += yShift;
		z[p] = -sqrt(r_sphere*r_sphere -x[p]*x[p] - y[p]*y[p]);
	}
}


void Setup_Regular_Hexagonal_Straight_Lines(double *&x, double *&y, double *&z, double r_piston, double r_sphere, double separation, double r_aperture,double r_transducer,  int &N){ 
// TOO TRICKY FOR NOW! Fill up the transducer with hexagonally packed pistons by placing a line of pistons across the y=0 axis, making sure the minimum spacing between elements is separation, they don't overlap the external radius which is r_transducer, and they don't overlap the internal radius, which is r_aperture.  Output number of elements to N.  
// Actually, this is slightly false: Instead of being separation apart, their centres are 2r+separation apart in a straight line. This is not quite the same, because of the curvature, but it's good enough.
//int NumberInFirstLine = int((r_transducer));
// NB: If you overwrite N, you better resize x,y,z otherwise you will get a seg fault.

}

void Setup_Regular_Dumb_Packing_Straight_Lines(double *&x, double *&y, double *&z, double r_piston, double r_sphere, double separation, double r_aperture,double r_transducer,  int &N){ 
/* Fill up the transducer with pistons by placing a line of pistons across the y=0 axis, making 
   sure the minimum spacing between elements is separation, they don't overlap the external radius 
   which is r_transducer, and they don't overlap the internal radius, which is r_aperture.  
   Output number of elements to N.  
   ACTUALLY!!! This is not quite right. What I'm doing here is placing points regularly on a square 
   and then projecting back to the sphere. It's not great, but it'll do. 
   So this is how we do it. The separation between centres is s = 2r_piston+separation. 
   Start at ind_x = -rtranscucer/s and ind_y =  -rtranscucer/s. Increment them and test to see 
   whether the relevant point makes it into the sphere.
 */
	int Ntemp;
	double s  = 2*r_piston+separation;
	Ntemp = int(r_transducer/s); 
	N=0; // reset the counter.
	delete[]  x;
	delete[]  y;
	delete[]  z;
	x = new double [(2*Ntemp+1)*(2*Ntemp+1)]; // It certainly isn't bigger than this. Allow the extra space to just be wasted.
	y = new double [(2*Ntemp+1)*(2*Ntemp+1)]; // It certainly isn't bigger than this. Allow the extra space to just be wasted.
	z = new double [(2*Ntemp+1)*(2*Ntemp+1)]; // It certainly isn't bigger than this. Allow the extra space to just be wasted.
	double radial_position,xpos,ypos;
	for (int ind_x = -Ntemp; ind_x <=Ntemp; ind_x++){
		xpos =  double(ind_x)*s;
		for (int ind_y = -Ntemp; ind_y <=Ntemp; ind_y++){
			ypos = double(ind_y)*s;
			radial_position = sqrt(xpos*xpos + ypos*ypos);
			if ( ((radial_position)<(r_transducer-r_piston)) && ((radial_position)>(r_aperture-r_piston))){
				//std::cout << N << " " << (2*Ntemp+1)*(2*Ntemp+1) <<  endl;
				x[N]=xpos;
				y[N]=ypos;
				z[N]=-sqrt(r_sphere*r_sphere - xpos*xpos - ypos*ypos);
				N++;
			}
		}
	}
}

void Circular_Piston_Transducer::Setup_Transducer(string &filename){
/* Grab information from the file and use it to set up the transducer.
   filename is the name of the file with the transducer parameters in it. This should contain, in order:
   1) 5 doubles: Radius of curvature of the transducer, radius of piston, minimum separation, radius of 
      transducer in projection, radius of interior hole in transducer.  
   2) integer: Number of elements that the transducer contains. 
   3) 2 integer: Number of radial and angular bits that each transducer element is cut into
   4) 2 double: f and c. From this, make omega/c = 2pi*f/c. 
   5) integer: option - allows you to choose between either loading positions directly, randomly assigning, etc.
   6) IF option is appropriate: 4 double: one line for each transducer element, x, y, pressure_real,pressure_imag. 
      Of course these pressures might not be used, but load them in anyway.
  */
	ifstream Parameters; 
	Parameters.open(filename.c_str());
#ifdef DEBUG
	std::cout << filename << std::endl;
#endif
	if(!Parameters.is_open()) 
{
	cout << "1 Couldn't open your parameter file. Sorry. No success this time. " << filename << endl;      
}     
	char comma; 
	// unimportant variable, used when reading in comma-separated lists.
	string discardMe; 
	// similar to comma. 
   
	Parameters >> cd_radius_of_sphere >> comma >> cd_radius_of_piston >> comma >> cd_minimum_separation  >> comma  >> cd_radius_of_transducer >> comma >> cd_radius_of_imaging_aperture;  // radius of sphere  making up the transducer
	getline(Parameters,discardMe); 
	// go to next line, discarding all comments on rest of line.
	
#ifdef DEBUG
	std::cout << cd_radius_of_sphere <<  comma << cd_radius_of_piston << comma << cd_minimum_separation  <<  comma << cd_radius_of_transducer << comma << cd_radius_of_imaging_aperture << std::endl;
#endif

	int li_number_of_elements;
	Parameters >> li_number_of_elements; 
	Input_Number_Of_Elements(li_number_of_elements); 
	// NOTE must call this otherwise you never assign memory to the different elements and get SIGSEGVs. 
	
	getline(Parameters,discardMe);
	
#ifdef DEBUG
	std::cout << ci_number_of_elements  << std::endl;
#endif
	
	Parameters >> ci_number_of_radial_bits >> comma >> ci_number_of_angles;
	getline(Parameters,discardMe); 
	
#ifdef DEBUG
	std::cout << ci_number_of_radial_bits << comma << ci_number_of_angles << std::endl;
#endif
	
	double ld_freq,ld_speed;
	Parameters >> ld_freq >> comma >>  ld_speed;
	getline(Parameters,discardMe); 
#ifdef DEBUG
	std::cout << ld_freq << "  " << ld_speed << std::endl;
#endif
	
	cd_omega_div_c = (PI2*ld_freq/ld_speed); 
	int option;
	Parameters >> option;
	getline(Parameters,discardMe); 
	
	// OK, now it's time to break into different cases. 
	// But note: all the different cases do is give different positions of the 
	// transducer elements. The actual generation of things can come later and 
	// will be the same for all cases. So do it like that! Don't use constructors 
	// for this - you might want to setup transducers multiple times and you don't 
	// want to have to quit the program to do that.
	
	// Set up arrays to hold the positions and pressure.
	double *ld_x,*ld_y,*ld_z,*ld_pressure_real,*ld_pressure_imag;
	ld_x = new double[ci_number_of_elements];
	ld_y = new double[ci_number_of_elements];
	ld_z = new double[ci_number_of_elements];
	ld_pressure_real = new double[ci_number_of_elements];
	ld_pressure_imag = new double[ci_number_of_elements];
	
	switch (option){
	case 1: // In this case, we load up the positions from the remainder of the file.
		for (int li_ind = 0; li_ind < ci_number_of_elements ; li_ind++){
			Parameters >> ld_x[li_ind] >> comma >> ld_y[li_ind]>> comma >> ld_z[li_ind] >> comma >> ld_pressure_real[li_ind] >> comma >> ld_pressure_imag[li_ind];
			getline(Parameters,discardMe);
			
#ifdef DEBUG
		//std::cout << li_ind << comma <<  ld_x[li_ind] << comma << ld_y[li_ind] << comma <<  ld_z[li_ind] << comma <<ld_pressure_real[li_ind] << comma << ld_pressure_imag[li_ind] << std::endl;
#endif
			
		}
		break;
	case 2:
		Setup_Regular_Dumb_Packing_Straight_Lines(ld_x, ld_y,ld_z, cd_radius_of_piston, cd_radius_of_sphere, cd_minimum_separation ,cd_radius_of_imaging_aperture  ,cd_radius_of_transducer ,  li_number_of_elements);
		Input_Number_Of_Elements(li_number_of_elements); 
		// NOTE must call this otherwise you never assign memory to the different elements and get SIGSEGVs. 
		cout << " Number of elements " << li_number_of_elements<<endl;
		break;
	}
	
#ifdef DEBUG
	std::cout << "LOADED" << std::endl;
//TEMP stuff.
//	Transducer_Element *tmp = new Circular_Piston_Element;
//	dynamic_cast<Circular_Piston_Element*>(tmp)->Input_Position(ld_x[0],ld_y[0],1.0);
//	std::cout << "Test succeeded " << std::endl;
#endif

    // Check the consistency of the newly generated grid. 
    if (!Consistency_Checker(ld_x, ld_y, ld_z, cd_radius_of_piston, cd_minimum_separation, li_number_of_elements)){
		std::cerr<< "Failed consistency check" << std::endl;
		//exit(1);
	}	
	
	// OK, now we know the positions, etc. Run through all the elements and load this info.
	
	for (int li_ind = 0; li_ind < ci_number_of_elements ; li_ind++){ 
		// OK, I could use function calls for these - I implemented the functions in a previous version 
		// but I don't see what the gain would be.
		//cout << li_ind << " " << ld_x[li_ind] <<  " "<< ld_y[li_ind] <<  " "<< ld_z[li_ind] <<  endl ;
		
		dynamic_cast<Circular_Piston_Element*>( cte_array_of_elements[li_ind])->Input_Position(ld_x[li_ind],ld_y[li_ind],ld_z[li_ind]);
		dynamic_cast<Circular_Piston_Element*>( cte_array_of_elements[li_ind])->Input_Normal(ld_x[li_ind],ld_y[li_ind],ld_z[li_ind]);
		cte_array_of_elements[li_ind]->Input_Radii(0.0, cd_radius_of_piston);
		cte_array_of_elements[li_ind]->Input_Element_Numbers(ci_number_of_radial_bits, ci_number_of_angles);
		cte_array_of_elements[li_ind]->Input_Angles(0.0, PI2);
		cte_array_of_elements[li_ind]->Input_Pressure(ld_pressure_real[li_ind], ld_pressure_imag[li_ind]); 
	}
	this->Generate_All_Elements(); 
	// Now all of the transducer elements set themselves up. i.e. they assign info to their computational elements. 
	Parameters.close();
}
