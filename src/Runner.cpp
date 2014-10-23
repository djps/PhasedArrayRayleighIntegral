/*
 Here, I implement the different commands for setting up and optimizing the transducer. 
*/      

#include "Runner.h"
      
using namespace std;
      
void Output_Data(Optimize_Transducer *my_transducer, string parameterFile, string outputFile){
	// MODIFIED! 18.08.2010: Instead of outputting all data, only output data where the intensity 
	// is at least SAVE_THRESHOLD of the value at the focus. 
	double ld_xrange_of_output[2];
	double ld_yrange_of_output[2];
	double ld_zrange_of_output[2];
	int li_x_resolution_output,li_y_resolution_output,li_z_resolution_output;
	ifstream Parameters; 
	Parameters.open(parameterFile.c_str());
	if (!Parameters.is_open()) 
	{
		cerr << "Couldn't open your outputting-parameter file. Sorry. No success this time. " << endl;      
	}     
	char Comma ; // unimportant variable, used when reading in comma-separated lists.
	string DiscardMe; // similar to Comma.
	Parameters >> ld_xrange_of_output[0] >> Comma >> ld_xrange_of_output[1];     
//   cout << cd_xrange_of_output[0] << Comma << cd_xrange_of_output[1]   << endl;
	getline(Parameters,DiscardMe);
	Parameters >> ld_yrange_of_output[0] >> Comma >> ld_yrange_of_output[1];     
//   cout << cd_yrange_of_output[0] << Comma << cd_yrange_of_output[1]   << endl;
	getline(Parameters,DiscardMe);
	Parameters >> ld_zrange_of_output[0] >> Comma >> ld_zrange_of_output[1];     
//   cout << cd_zrange_of_output[0] << Comma << cd_zrange_of_output[1]   << endl;
	getline(Parameters,DiscardMe);
	Parameters >> li_x_resolution_output; 
//   cout << ci_x_resolution_output  << endl;
	getline(Parameters,DiscardMe);
	Parameters >> li_y_resolution_output; 
//   cout << ci_y_resolution_output  << endl;
	getline(Parameters,DiscardMe);
	Parameters >> li_z_resolution_output;
//   cout << ci_z_resolution_output  << endl;
	Parameters.close();
	
	// Now do the outputting
	ofstream lo_output; 
	lo_output.open(outputFile.c_str());
	if (!lo_output.is_open()) 
	{
		cerr << "Couldn't open your output file. Sorry. No success this time. " << endl;      
	}     
	double ld_dx = ((li_x_resolution_output >0) ? (ld_xrange_of_output[1]-ld_xrange_of_output[0])/double(li_x_resolution_output):0.0);
	double ld_dy = ((li_y_resolution_output >0) ? (ld_yrange_of_output[1]-ld_yrange_of_output[0])/double(li_y_resolution_output):0.0);
	double ld_dz = ((li_z_resolution_output >0) ? (ld_zrange_of_output[1]-ld_zrange_of_output[0])/double(li_z_resolution_output):0.0);
	double ld_x, ld_y, ld_z, ld_press_real, ld_press_imag, ld_press_magn; 
	
	my_transducer->Calculate_Pressure(my_transducer->cd_position_of_focus[0],my_transducer->cd_position_of_focus[1],my_transducer->cd_position_of_focus[2], ld_press_real, ld_press_imag); // Get the value at the focus. 
	double ld_max = ld_press_real*ld_press_real  + ld_press_imag*ld_press_imag; // Get the intensity at the focus
	
	ld_x = ld_xrange_of_output[0];
	for (int li_x = 0; li_x <(li_x_resolution_output+1) ; li_x++)
	{ 
		ld_y = ld_yrange_of_output[0];
		for (int li_y = 0; li_y <(li_y_resolution_output+1) ; li_y++)
		{
			ld_z = ld_zrange_of_output[0];
			for (int li_z = 0; li_z <(li_z_resolution_output+1) ; li_z++)
			{
				my_transducer->Calculate_Pressure(ld_x, ld_y, ld_z, ld_press_real, ld_press_imag);
				ld_press_magn = ld_press_real*ld_press_real  + ld_press_imag*ld_press_imag;
				if (ld_press_magn > SAVE_THRESHOLD*ld_max) {
					lo_output << ld_x << ", "<<   ld_y <<", "<<  ld_z << ", "<< ld_press_real << ", "<< ld_press_imag <<", "<< ld_press_real*ld_press_real  + ld_press_imag*ld_press_imag << "\n"; // only save data if the ratio with the steered intensity exceeds the save threshold
					if (ld_press_magn > ld_max) ld_max = ld_press_magn; // If you've got a small estimate of the max pressure value, you output too much. This way, we  can update it as we go. This should reduce file size for little effort.  
				}
				ld_z +=ld_dz;
			}
			ld_y += ld_dy;
		}
		ld_x +=ld_dx;
	}
	lo_output.close();
}




void Output_Data_Binary(Optimize_Transducer *my_transducer,string parameterFile,string outputFile){ 

/* Dumps the information in binary format. Format is 3 integers: nx, ny, nz. 
   Then nx*ny*nz values representing the real parts, then another set for the imaginary parts. 
   No point in outputting the magnitudes, because we're going to have to process things a lot anyway.  
   I'm hoping this leads to much smaller files and faster running of the program, as there is much less I/O.
 */
 
	double ld_xrange_of_output[2];
	double ld_yrange_of_output[2];
	double ld_zrange_of_output[2];
	int li_x_resolution_output, li_y_resolution_output, li_z_resolution_output;
	ifstream Parameters; 
	Parameters.open(parameterFile.c_str());
	if (!Parameters.is_open()) 
	{
		cerr << "Couldn't open your outputting-parameter file. Sorry. No success this time. " << endl;      
	}     
	char Comma ; // unimportant variable, used when reading in comma-separated lists.
	string DiscardMe; // similar to Comma.
	Parameters >> ld_xrange_of_output[0] >> Comma >> ld_xrange_of_output[1] ;     
//   cout << cd_xrange_of_output[0] << Comma << cd_xrange_of_output[1]   << endl;
	getline(Parameters,DiscardMe);
	Parameters >> ld_yrange_of_output[0] >> Comma >> ld_yrange_of_output[1] ;     
//   cout << cd_yrange_of_output[0] << Comma << cd_yrange_of_output[1]   << endl;
	getline(Parameters,DiscardMe);
	Parameters >> ld_zrange_of_output[0] >> Comma >> ld_zrange_of_output[1] ;     
//   cout << cd_zrange_of_output[0] << Comma << cd_zrange_of_output[1]   << endl;
	getline(Parameters,DiscardMe);
	Parameters >> li_x_resolution_output ; 
//   cout << ci_x_resolution_output  << endl;
	getline(Parameters,DiscardMe);
	Parameters >> li_y_resolution_output ; 
//   cout << ci_y_resolution_output  << endl;
	getline(Parameters,DiscardMe);
	Parameters >> li_z_resolution_output ;
//   cout << ci_z_resolution_output  << endl;
	Parameters.close();
	
	// Now generate the data. 
	Array3D<double> *la_real_pressures = new Array3D<double>;
	Array3D<double> *la_imag_pressures = new Array3D<double>;
	la_real_pressures->Resize_Array(li_x_resolution_output + 1, li_y_resolution_output + 1,li_z_resolution_output + 1);
	la_imag_pressures->Resize_Array(li_x_resolution_output + 1, li_y_resolution_output + 1,li_z_resolution_output + 1);
	double ld_dx = ((li_x_resolution_output >0) ? (ld_xrange_of_output[1]-ld_xrange_of_output[0])/double(li_x_resolution_output):0.0);
	double ld_dy = ((li_y_resolution_output >0) ? (ld_yrange_of_output[1]-ld_yrange_of_output[0])/double(li_y_resolution_output):0.0);
	double ld_dz = ((li_z_resolution_output >0) ? (ld_zrange_of_output[1]-ld_zrange_of_output[0])/double(li_z_resolution_output):0.0);
	double ld_x, ld_y, ld_z, ld_press_real, ld_press_imag, ld_press_magn; 
	ld_x = ld_xrange_of_output[0];
	for (int li_x = 0; li_x <(li_x_resolution_output+1) ; li_x++)
	{ 
		std::cout << li_x << ", " << std::flush;
		ld_y = ld_yrange_of_output[0];
		for (int li_y = 0; li_y <(li_y_resolution_output+1) ; li_y++)
		{
			ld_z = ld_zrange_of_output[0];
			for (int li_z = 0; li_z <(li_z_resolution_output+1) ; li_z++)
			{
				my_transducer->Calculate_Pressure(ld_x, ld_y, ld_z, ld_press_real, ld_press_imag);
				//lo_output << ld_x << ", "<<   ld_y <<", "<<  ld_z << ", "<< ld_press_real << ", "<< ld_press_imag <<", "<< ld_press_real*ld_press_real  + ld_press_imag*ld_press_imag << "\n";
				la_real_pressures->Input_Value(li_x,li_y,li_z,ld_press_real);
				la_imag_pressures->Input_Value(li_x,li_y,li_z,ld_press_imag);
				ld_z +=ld_dz;
			}
			ld_y += ld_dy;
		}
		ld_x +=ld_dx;
	}
	std::cout << std::endl;
	
	// Now do the outputting
	ofstream lo_output; 
	lo_output.open(outputFile.c_str(), ios::out  | ios::binary);
	
	if (!lo_output.is_open()) 
	{
		cerr << "Couldn't open your output file. Sorry. No success this time. " << endl;      
	}    
	int nx =li_x_resolution_output+1;  
	int ny =li_y_resolution_output+1;  
	int nz =li_z_resolution_output+1;  
	lo_output.write((char*) &nx, sizeof(int));
	lo_output.write((char*) &ny, sizeof(int));
	lo_output.write((char*) &nz, sizeof(int));
	lo_output.write((char*) la_real_pressures->array, sizeof(double)*(nx*ny*nz));
	lo_output.write((char*) la_imag_pressures->array, sizeof(double)*(nx*ny*nz));
	
	lo_output.close();
}


void Output_Other_Data(Optimize_Transducer *my_transducer, string outputFile){

/* 
Outputs the pressure configuration on the transducer bowl. 
It makes sense to output just the centroid coordinates, real and imaginary parts 
of the pressure, because then I can use the output as an input to another simulation. 
 */
 
	ofstream lo_output; 
	lo_output.open(outputFile.c_str());
	if(!lo_output.is_open()) 
	{
		cerr << "Couldn't open your output file. Sorry. No success this time. " << endl;      
	}
	double ld_real_pressure, ld_imag_pressure;     
	double ld_x, ld_y, ld_z; 
	for (int li_element_index = 0; li_element_index < my_transducer->ci_number_of_elements; ++li_element_index)
	{ // cycle over all elements. 
		my_transducer->cte_array_of_elements[li_element_index]->Output_Pressure(ld_real_pressure, ld_imag_pressure) ; // Get the pressure on this transducer element.
		my_transducer->cte_array_of_elements[li_element_index]->Output_Centroid(ld_x, ld_y, ld_z ) ; // Get the centroid. 
		lo_output << ld_x << ", "<<   ld_y <<", "<<  ld_z << ", "<< ld_real_pressure << ", "<< ld_imag_pressure << "\n";
	}
	lo_output.close();
}



Runner::Runner(int argc, char *argv[]){

/*
	argv[0]  program name
	argv[1]  Option: 1 means that we use circular pistons.
	argv[2]  Input file containing info about outputting
	argv[3]  Input file containing info about optimizing
	argv[4]  Input file containing info about transducer
	argv[5]  Output file containing pressure
	argv[6]  Output file containing transducer
	argv[7]: Optimization option: 1 = pseudoinverse
 */
     
	if (argc != EXPECTED_ARGS)
	{
		cerr<< "Wrong no of args. You gave " << argc <<  " and should have given " << EXPECTED_ARGS << endl;
		exit(1);
	}
	
	switch (atoi(argv[1]))
	{
	case 1: 
		my_transducer = new Circular_Piston_Transducer;
		break;
	}
   
	// File names
    Outputter = argv[2];
    Optimizer = argv[3];
    Transducer = argv[4];
   
	// Set up the transducer.
	lf_pointer_to_array_of_functions  = new Function*[1] ; 
	lf_pointer_to_array_of_functions[0] = my_transducer; // point it to the transducer.
	cout << "Setting up transducer ................." << flush;
	my_transducer->Setup_Transducer(Transducer);
	cout << " DONE" << endl;
	cout << "Setting up optimization ................." <<flush;
	my_transducer->Setup_Optimization(Optimizer);
	cout << " DONE" << endl;
	//cout << "Setting up output ................."  << flush;
	//my_transducer->Setup_Output(Outputter);
	//cout << " DONE" << endl;
   
    cout << "Output will be to file named " << argv[5] << endl;
    string ls_filename = argv[5]; // Basically I use this to cast argv[5] into a string. I don't like casts inside function calls.
	//  my_transducer->Change_Output_File_Name(ls_filename); 
   
    cout << "Output of transducer stuff  will be to file named " << argv[6] << endl;
    string ls_filename2 = argv[6]; // Basically I use this to cast argv[6] into a string. I don't like casts inside function calls.
	//my_transducer->Change_Other_Output_File_Name(ls_filename2); 
   
    cout << "Outputting initial transducer stuff to file.... " << flush;
    Output_Other_Data(dynamic_cast<Optimize_Transducer*>(my_transducer),ls_filename2);
    cout << "DONE!" << endl;
   
	int li_number_of_elements;
	my_transducer->Output_Number_Of_Elements(li_number_of_elements);
	my_transducer->Output_Area();
   
	double *pressures = new double [2*li_number_of_elements];
	switch (atoi(argv[7]))
	{
		case -1: // Use the pressures that were loaded in the transducer file. i.e. don't do anything.
			break;		   
		case 0: // Use phase-coherence at the desired focus.
			cout << "Optimizing pressures... " << flush;
			my_transducer->Generate_Initial_Pressures(pressures);
			break;
		case 1: // Pseudoinverse. Need to determine the transfer matrix first. 
			std::cout << "Getting information for transfer matrix................." << flush;
			my_transducer->Setup_Pressures();
			std::cout << "DONE ! " <<std::endl;
			cout << "Optimizing pressures... " << flush;
			my_transducer->Generate_Pressures_Using_Pseudoinverse(pressures);  // Get the starting guess using pseudoinverse.
			break;
		case 2:
			std::cout << "Getting information for transfer matrix................." << flush;
			my_transducer->Setup_Pressures();
			std::cout << "DONE ! " <<std::endl;
			cout << "Optimizing pressures... " << flush;
			my_transducer->Generate_Pressures_Using_Fancy_Method2(pressures); // Try it and see!
			break;
		case 3:
			std::cout << "Getting information for transfer matrix................." << flush;
			my_transducer->Setup_Pressures();
			std::cout << "DONE ! " <<std::endl;
			cout << "Optimizing pressures... " << flush;
			my_transducer->Generate_Pressures_Using_Fancy_Method3(pressures); // Try it and see!
			break;
		case 4: // Piston time-reversal. No need for transfer matrix,
			cout << "Optimizing pressures... " << flush;
			my_transducer->Generate_Pressures_Using_Piston_Time_Reversal(pressures); // Try it and see!
			break;
		default:
			std::cout << "Getting information for transfer matrix................." << flush;
			my_transducer->Setup_Pressures();
			std::cout << "DONE ! " <<std::endl;
			cout << "Optimizing pressures... " << flush;
			my_transducer->Generate_Pressures_Using_Pseudoinverse(pressures);  // Get the starting guess using pseudoinverse.				   
	}
	cout << "DONE!" << endl;

	Normalize_Pressures(pressures, li_number_of_elements); // Normalize them so that the maximum amplitude is 1. 
   
	// Now copy the pressures into the transducer. 
	for (int li_ind = 0; li_ind<li_number_of_elements; li_ind++)
	{
		my_transducer->Input_Pressure(pressures[li_ind], pressures[li_ind+li_number_of_elements], li_ind);
	}
	cout << "Outputting data to file ..... " << flush;   
	Output_Data(dynamic_cast<Optimize_Transducer*>(my_transducer),Outputter,ls_filename) ; 	
	cout << "DONE! " << endl;
	cout << "Outputting the transducer data to file ..... " << flush;   
	Output_Other_Data(dynamic_cast<Optimize_Transducer*>(my_transducer),ls_filename2);
	cout << "DONE! " << endl;
   
	delete my_transducer; 
   
 }




void Runner::Normalize_Pressures(double *pressures, int N){ 

/* 
Straightforward: run through array which consists of N real parts and N corresponding 
imag parts, normalize maximum magnitude to 1.
*/
 
	double ld_tmp_max = 0.0;
	double ld_real, ld_imag, ld_magn;
	for(int li_ind = 0; li_ind < N; li_ind++)
	{
		ld_real = pressures[li_ind];
		ld_imag = pressures[li_ind+N];
		ld_magn = ld_real*ld_real + ld_imag*ld_imag;
		if(ld_magn > ld_tmp_max) ld_tmp_max = ld_magn;
	}
	// Now normalize.  
	for (int li_ind = 0; li_ind < N; li_ind++)
	{
		pressures[li_ind] /= sqrt(ld_tmp_max/10.0); // possible divide by zero. Oh well...
		pressures[li_ind+N] /= sqrt(ld_tmp_max/10.0); 
	}    
}