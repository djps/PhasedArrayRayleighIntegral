#include "OptimizeTransducer.h"
 
void Optimize_Transducer::Setup_Optimization(string &filename){  
// Grab information from the file and use it to set up the optimization stuff.
/*
   OK, here we are grabbing info for the optimization. This info is stored in the current class, so it's really quite simple. 
   The expected info, in order, is:
   1) double: 3 values, x,y,z of focus 
   2) double: pressure value at focus. 
   3) double: 2 values, x_min and x_max for optimizing region. 
   4) double: same for y. 
   5) double: same for z.       
   6) int: number of x values for optimization grid. 
   7) int: number of y values for optimization grid.
   8) int: number of z values for optimization grid.   
   9) int: number of parameters. The fancy1 method is designed to be somewhat adaptable and you can pass loads of parameters to it. However, the first of these will be overwritten by lambda, the weighting function for relating the two optimizations, i.e. maximizing pressure output at transducer vs making the nicest possible focus.   
   10) array of comma separated doubles: the parameters for fancy1.
	11) 3 doubles: cd_desired_width_of_focus,cd_scale_factor_for_wide_focus,cd_bias_factor_for_wide_focus
		  
	DEPRECATED
	9) double:  tolerance value. Anything below it is ignored completely. Only used in the deprecated conj-grad optimization
	10) double: tolerance value.  anything above it is penalized heavily.  Only used in the deprecated conj-grad optimization
	11) double: condition number. Maximum allowable condition number for the fancy1 method of optimizing the problem.
  */
	ifstream Parameters; 
	Parameters.open(filename.c_str());
	if(!Parameters.is_open()) 
	{
      cout << "Couldn't open your optimization parameter file. Sorry. No success this time. " << endl;      
	}     
	char Comma ; // unimportant variable, used when reading in comma-separated lists.
	string DiscardMe; // similar to Comma. 
	Parameters >> cd_position_of_focus[0] >> Comma >>cd_position_of_focus[1] >> Comma >>cd_position_of_focus[2] ; 
	cout << "FOCUS " << cd_position_of_focus[0] << " " << cd_position_of_focus[1] << " "  << cd_position_of_focus[2] << endl; 
	getline(Parameters,DiscardMe);
	Parameters >>   cd_value_at_focus;  
//   cout <<cd_value_at_focus << endl;
	getline(Parameters,DiscardMe);
	Parameters >> cd_xrange_of_optimized_region[0] >> Comma >>cd_xrange_of_optimized_region[1] ;       
	//cout << cd_xrange_of_optimized_region[0] << "  " << cd_xrange_of_optimized_region[1]<< endl;
	getline(Parameters,DiscardMe);
	Parameters >>  cd_yrange_of_optimized_region[0] >> Comma >>cd_yrange_of_optimized_region[1] ;       
	//cout << cd_yrange_of_optimized_region[0] << "  " << cd_yrange_of_optimized_region[1]<< endl;
	getline(Parameters,DiscardMe);
	Parameters >> cd_zrange_of_optimized_region[0] >> Comma >> cd_zrange_of_optimized_region[1] ; 
	//cout << cd_zrange_of_optimized_region[0] << "  " << cd_zrange_of_optimized_region[1]<< endl;
	getline(Parameters,DiscardMe);
	Parameters >> ci_x_resolution_optimization ; 
	//cout << ci_x_resolution_optimization<< endl;
	getline(Parameters,DiscardMe);
	Parameters >>  ci_y_resolution_optimization ;
	//cout << ci_y_resolution_optimization<< endl;
	getline(Parameters,DiscardMe);
	Parameters >>  ci_z_resolution_optimization ;
	//cout << ci_z_resolution_optimization<< endl;
	getline(Parameters,DiscardMe);
   
	// This is the only place to read in optimization parameters. 
	Parameters >> ci_number_of_optimization_parameters;
	getline(Parameters,DiscardMe);
	cpd_optimization_parameters = new double [ci_number_of_optimization_parameters];
	for (int ind = 0; ind < ci_number_of_optimization_parameters-1; ind++) Parameters >> cpd_optimization_parameters[ind] >> Comma;  
	Parameters >> cpd_optimization_parameters[ci_number_of_optimization_parameters-1];
	getline(Parameters,DiscardMe);
	Parameters >> cd_desired_width_of_focus >> Comma >> cd_scale_factor_for_wide_focus >> Comma >> cd_bias_factor_for_wide_focus; 

	Parameters.close(); 
   
}


double Optimize_Transducer::function(int N, double *parameters){ 
 
	/* This is the function that is called during optimization. The parameters will be the real and 
	   imaginary parts of all the pressures. N should be twice the number of transducers. 
	 */
	// There will be a bunch of parameters. These are all stored in cpd_optimization_parameters, but it helps to name them.
	int ci_number_of_maxima = int(cpd_optimization_parameters[0]); 
	// When analyzing with conj grad, we need to know how many maxima to consider.
	double cd_soft_tolerance = cpd_optimization_parameters[1]; 
	// If magnitudes away from focus are less than this, don't include them in optimization at all. 
	// If in excess of this, don't get too worried, but gently guide them down.
	double cd_hard_tolerance= cpd_optimization_parameters[2]; 
	// If anything is larger than this, do whatever you can to lower it. 
      
	if(N != (2*ci_number_of_elements)){
		cout << "Trouble in the optimization function"; 
		return 0; 
	}
	else
	{
		char lc_tmp;  
		double ld_dx = (cd_xrange_of_optimized_region[1]-cd_xrange_of_optimized_region[0])/double(ci_x_resolution_optimization);
		double ld_dy = (cd_yrange_of_optimized_region[1]-cd_yrange_of_optimized_region[0])/double(ci_y_resolution_optimization);
		double ld_dz = (cd_zrange_of_optimized_region[1]-cd_zrange_of_optimized_region[0])/double(ci_z_resolution_optimization);
		double ld_error_function = 0.0; // keeps track of the deviation from what we want.
		double ld_x, ld_y, ld_z, ld_press_real, ld_press_imag, ld_press_magn;
		   
		double ld_pos_max_x,ld_pos_max_y,ld_pos_max_z,ld_max; // keep track of the position of the maximum.
		ld_max = 0.0;
		ld_pos_max_x =  cd_xrange_of_optimized_region[0];
		ld_pos_max_y =  cd_yrange_of_optimized_region[0];
		ld_pos_max_z =  cd_zrange_of_optimized_region[0];
		for (int li_x = 0; li_x <(ci_x_resolution_optimization+1) ; li_x++)
		{ 
			//std::cin >> lc_tmp; 
			for (int li_y = 0; li_y <(ci_y_resolution_optimization+1) ; li_y++)
			{
				for (int li_z = 0; li_z <(ci_z_resolution_optimization+1) ; li_z++)
				{
					ld_press_real = 0.0;
					ld_press_imag = 0.0;
					for (int li_ind = 0; li_ind < ci_number_of_elements; li_ind++) 
					{ 
						// Here, we run over the transducers and get the pressures. 
						ld_press_real += (parameters[li_ind]*ca_real_array_for_optimization[li_ind].Output_Value(li_x,li_y,li_z) - parameters[li_ind+ci_number_of_elements] *ca_imag_array_for_optimization[li_ind].Output_Value(li_x,li_y,li_z));
						ld_press_imag += (parameters[li_ind]*ca_imag_array_for_optimization[li_ind].Output_Value(li_x,li_y,li_z) + parameters[li_ind+ci_number_of_elements] *ca_real_array_for_optimization[li_ind].Output_Value(li_x,li_y,li_z));
					}
					ca_energy_in_optimization_region[0].Input_Value(li_x, li_y, li_z,  (ld_press_real*ld_press_real + ld_press_imag*ld_press_imag));
				}
			}
		}
		// OK, so now we have the energy in the optimizing region. Let's get the top 5 maxima
		int *lpi_x_pos = new int[ci_number_of_maxima];         
		int *lpi_y_pos = new int[ci_number_of_maxima];         
		int *lpi_z_pos = new int[ci_number_of_maxima];         
		double *lpd_x_pos = new double[ci_number_of_maxima];         
		double *lpd_y_pos = new double[ci_number_of_maxima];         
		double *lpd_z_pos = new double[ci_number_of_maxima];         
		double *lpd_value = new double[ci_number_of_maxima];
				
		ca_energy_in_optimization_region[0].Find_Local_Maxima(ci_number_of_maxima, lpi_x_pos,lpi_y_pos,lpi_z_pos,lpd_value );
		for (int ind = 0; ind < ci_number_of_maxima; ind++)
		{
			lpd_x_pos[ind] = cd_xrange_of_optimized_region[0] + ld_dx*lpi_x_pos[ind]; 
			lpd_y_pos[ind] = cd_yrange_of_optimized_region[0] + ld_dy*lpi_y_pos[ind]; 
			lpd_z_pos[ind] = cd_zrange_of_optimized_region[0] + ld_dz*lpi_z_pos[ind]; 
		} 
   
		// OK, we want the main maximum to be at the focus, we want it to be as large as possible, and we want all other peaks to be less than soft_tolerance. If any other peak is larger than hard_tolerance, we whack it hard. But not too hard. What has sometimes happened is that we whack the sidelobes too hard, the transducer switches off, then when I force the pressure output to be normalized to 1, suddenly I get basically a random distribution. No to that! 
    
#ifdef DEBUG_OPTIMIZATION
   std::cout << "Error function: " << ld_error_function << ", "; 
#endif   
   
		// First, the largest max MUST be at the focus.
		ld_error_function += 1000.0*sqrt((lpd_x_pos[0] - cd_position_of_focus[0])*(lpd_x_pos[0] - cd_position_of_focus[0]) +(lpd_y_pos[0] - cd_position_of_focus[1])*(lpd_y_pos[0] - cd_position_of_focus[1]) +(lpd_z_pos[0] - cd_position_of_focus[2])*(lpd_z_pos[0] - cd_position_of_focus[2]));
   
#ifdef DEBUG_OPTIMIZATION
   std::cout <<  ld_error_function << ", "; 
#endif   

		// Next, it should be equal to the desired value. But this is less important than having the max at focus. 
		//ld_error_function += 10.0*fabs(cd_value_at_focus - lpd_value[0]);
		 
#ifdef DEBUG_OPTIMIZATION
   //std::cout <<  ld_error_function << ", "; 
#endif   

		// 3) All other peaks must be lower than hard_ and soft_tolerance. Recall that these things are ratios, so it's hard to figure out the best coefficient to add this term to the rest.  
		for (int ind = 1; ind < ci_number_of_maxima ; ind++)
		{
			if (lpd_value[ind]/lpd_value[0] >= cd_hard_tolerance)  ld_error_function += 1000.0*(lpd_value[ind]/lpd_value[0] - cd_hard_tolerance);
			if (lpd_value[ind]/lpd_value[0] >= cd_soft_tolerance)  ld_error_function += 50.0*(lpd_value[ind]/lpd_value[0] - cd_soft_tolerance); // Note, this is added even if the peak is also above hard_tolerance. This is deliberate and enforces continuity, which is always nice. 
		}
		
#ifdef DEBUG_OPTIMIZATION
   std::cout <<  ld_error_function << std::endl; 
   std::cout << cd_value_at_focus  << "  " <<  lpd_value[0]  << std::endl; 
#endif   

	return ld_error_function; 
	}
}




void Optimize_Transducer::Setup_Pressures(){ 
	// Here, we generate a 3-d array for each transducer element. This gives us the pressure in the region of interest, assuming p_0 =  1.0 over the relevant transducer. By linearity, this is useful.
   // NB: I've messed this up in the past. I have previously simply used whatever pressure each element was initialized with, not a pressure of 1.0. This is because I'm not as bright as my mom tells people. Still, it's stupid to force the user to initialize the pressures to 1.0 - it is good that the user can input his/her own choice of pressure to generate an output. So what I'm adding is that this routine: locally stores the value of the pressure on the element, changes the element's value to 1.0, calculates the pressure, then restores the input value. 
   ca_real_array_for_optimization = new Array3D<double> [ci_number_of_elements] ;
   ca_imag_array_for_optimization = new Array3D<double> [ci_number_of_elements] ;
   ca_energy_in_optimization_region = new Array3D<double> [1]; // I know it looks stupid. But this is the easiest way, 'cos I can just copy other code I've already stress-tested.
   
   for (int li_ind = 0; li_ind < ci_number_of_elements; li_ind++)
   {
      ca_real_array_for_optimization[li_ind].Resize_Array(ci_x_resolution_optimization+1,ci_y_resolution_optimization+1,ci_z_resolution_optimization+1);
      ca_imag_array_for_optimization[li_ind].Resize_Array(ci_x_resolution_optimization+1,ci_y_resolution_optimization+1,ci_z_resolution_optimization+1);
   }
   ca_energy_in_optimization_region[0].Resize_Array(ci_x_resolution_optimization+1,ci_y_resolution_optimization+1,ci_z_resolution_optimization+1); 
   // i.e. make it the correct size. 
   
   
   // generate regions on the grid.
   double ld_dx = (cd_xrange_of_optimized_region[1]-cd_xrange_of_optimized_region[0])/double(ci_x_resolution_optimization);
   double ld_dy = (cd_yrange_of_optimized_region[1]-cd_yrange_of_optimized_region[0])/double(ci_y_resolution_optimization);
   double ld_dz = (cd_zrange_of_optimized_region[1]-cd_zrange_of_optimized_region[0])/double(ci_z_resolution_optimization);
   
   double ld_tmp_real_pressure,ld_tmp_imag_pressure;  
   
   double ld_x, ld_y, ld_z, ld_real_pressure, ld_imag_pressure;
   std::cout << ci_number_of_elements <<  ":  " << flush; 
   for (int li_ind = 0; li_ind < ci_number_of_elements; li_ind++)
   {
      std::cout << li_ind << ", "  << flush ;
      // Store the pressure on this element, then set the element pressure to 1.0
      cte_array_of_elements[li_ind]->Output_Pressure(ld_tmp_real_pressure,ld_tmp_imag_pressure);
      cte_array_of_elements[li_ind]->Input_Pressure(1.0,0.0);      
         ld_x = cd_xrange_of_optimized_region[0];
         for (int li_x = 0; li_x <(ci_x_resolution_optimization+1) ; li_x++)
         {
            ld_y = cd_yrange_of_optimized_region[0];
            for (int li_y = 0; li_y <(ci_y_resolution_optimization+1) ; li_y++)
            {
               ld_z = cd_zrange_of_optimized_region[0];
               for (int li_z = 0; li_z <(ci_z_resolution_optimization+1) ; li_z++)
               {
                  cte_array_of_elements[li_ind]->Calculate_Pressure(ld_x,  ld_y, ld_z, cd_omega_div_c, ld_real_pressure, ld_imag_pressure);
                  ca_real_array_for_optimization[li_ind].Input_Value(li_x,li_y,li_z,ld_real_pressure);
                  ca_imag_array_for_optimization[li_ind].Input_Value(li_x,li_y,li_z,ld_imag_pressure);
                  ld_z +=ld_dz; 
                  #ifdef DEBUG_PSEUDO_INVERSE
                 // cout << "Pressures at point " << li_x <<", "<< li_y <<", "<< li_z << ": " << ld_real_pressure << ", "  << ld_imag_pressure<< endl;
                     //cin >> lc_tmp; 
                  #endif   

               }
               ld_y += ld_dy;
            }
            ld_x +=ld_dx;
         }
         // Restore the pressure on the element to what it was originally. 
      cte_array_of_elements[li_ind]->Input_Pressure(ld_tmp_real_pressure,ld_tmp_imag_pressure);      
   }
}

void Optimize_Transducer::Generate_Initial_Pressures(double *parameters)
{ // Generates the initial guess for the pressures based on distances and phases. writes the results into parameters, which is assumed to be long enough.
   double ld_x, ld_y, ld_z; // position of centroid of elements.
   double ld_dist; 
   for (int li_ind = 0; li_ind < ci_number_of_elements; li_ind++)
   {
      cte_array_of_elements[li_ind]->Output_Centroid(ld_x, ld_y, ld_z);
      ld_dist = sqrt((ld_x - cd_position_of_focus[0])*(ld_x - cd_position_of_focus[0]) + (ld_y - cd_position_of_focus[1])*(ld_y - cd_position_of_focus[1])+ (ld_z - cd_position_of_focus[2])*(ld_z - cd_position_of_focus[2]));
      parameters[li_ind] = cos(cd_omega_div_c * ld_dist);
      parameters[li_ind+ci_number_of_elements] = -sin(cd_omega_div_c * ld_dist);
   }
   std::cout <<  "Focus " << cd_position_of_focus[0] << " "<< cd_position_of_focus[1] << " "<< cd_position_of_focus[2] << " " << std::endl;
}







void Optimize_Transducer::Generate_Pressures_Using_Pseudoinverse(double *apd_pressures)
{ // Uses the pseudoinverse to generate the pressures. Writes the corresponding pressures into the array, real first then imaginary. Array must be twice as long as the number of transducer elements.  
   int li_number_of_points_in_space = (ci_x_resolution_optimization+1)*(ci_y_resolution_optimization+1)*(ci_z_resolution_optimization+1) ; 
   // Note the +1 is quite important. Derives from my way of describing a resolution, i.e. resolution is number of divisions, not number of regions.
   double *cpd_desired_result_of_optimization,  *cpd_optimized_pressure, *cpd_array_for_transfer_matrix;  
   try{
   		cpd_desired_result_of_optimization = new double[2*li_number_of_points_in_space]; 
   		cpd_optimized_pressure = new double[2*ci_number_of_elements];
   		cpd_array_for_transfer_matrix = new double[4*ci_number_of_elements*li_number_of_points_in_space];
   }
   catch (bad_alloc&)
   {
	   cerr << "Error allocating memory in Generate_Pressures_Using_Pseudoinverse." << endl;
	   exit(-1);
   }

 
   char lc_tmp; // Just a rubbish variable that is useful if you want to pause until keypressed. I should figure out the right way to do this!
   
   // Now we step through the grid where we want to optimize the solution. A lot of the following apparent garbage is just for the sake of figuring out where we are in space so that we can look for the grid point that corresponds to the focus. 
   double ld_dx = (cd_xrange_of_optimized_region[1]-cd_xrange_of_optimized_region[0])/double(ci_x_resolution_optimization);
   double ld_dy = (cd_yrange_of_optimized_region[1]-cd_yrange_of_optimized_region[0])/double(ci_y_resolution_optimization);
   double ld_dz = (cd_zrange_of_optimized_region[1]-cd_zrange_of_optimized_region[0])/double(ci_z_resolution_optimization);
   double ld_x, ld_y, ld_z, ld_real_pressure, ld_imag_pressure; 
   
   double ld_dist_to_focus ; // just a temporary starting point. 
   double ld_min_dist_to_focus = 10000.0; // just a temporary starting point.
   double ld_actual_position_of_focus[3] ; // a temporary variable just intended for output to screen. Basically, if the user uses a grid on which the focus doesn't appear, then we have to choose a nearby point at which to position the peak. This variable contains that point.  
   
   int li_counter, li_index_of_focus;
   li_counter = 0;
   ld_x = cd_xrange_of_optimized_region[0];
   // NOTE REALLY IMPORTANT: The order of x,y,z in the following loop defines the map from 3-d space to 1-d array. So it MUST be the same as that used below in generating the transfer matrix. Otherwise you get garbage. You don't want garbage. 
   for (int li_x = 0; li_x <(ci_x_resolution_optimization+1) ; li_x++)
   {
      ld_y = cd_yrange_of_optimized_region[0];
      for (int li_y = 0; li_y <(ci_y_resolution_optimization+1) ; li_y++)
      {
         ld_z = cd_zrange_of_optimized_region[0];
         for (int li_z = 0; li_z <(ci_z_resolution_optimization+1) ; li_z++)
         {
            cpd_desired_result_of_optimization[li_counter] = 0.0; 
            cpd_desired_result_of_optimization[li_counter + li_number_of_points_in_space] = 0.0; 
            ld_dist_to_focus = (ld_x-cd_position_of_focus[0])*(ld_x-cd_position_of_focus[0])+(ld_y-cd_position_of_focus[1])*(ld_y-cd_position_of_focus[1])+(ld_z-cd_position_of_focus[2])*(ld_z- cd_position_of_focus[2]); // OK, this is the squared distance, but it's good enough.
            if(ld_dist_to_focus < ld_min_dist_to_focus) 
            {
               ld_min_dist_to_focus = ld_dist_to_focus ; //update the minimum.
               li_index_of_focus = li_counter; // get the appropriate index.
               ld_actual_position_of_focus[0] = ld_x;
               ld_actual_position_of_focus[1] = ld_y;
               ld_actual_position_of_focus[2] = ld_z;
               #ifdef DEBUG_PSEUDO_INVERSE
                   //     cout << "Current position  " <<  ld_x << " " << ld_y << " " << ld_z << " " << endl;
                    //    cout << "Dist of this point to focus = " << ld_min_dist_to_focus<< endl;
                    //    cout << "Current index " << li_counter << endl;
               #endif   
            } 
            
            ++li_counter;
            ld_z +=ld_dz; 
         }
         ld_y += ld_dy;
      }
      ld_x +=ld_dx;
   }
   cpd_desired_result_of_optimization[li_index_of_focus] = cd_value_at_focus;
#ifdef EXTENDED_FOCUS
// This makes the focus a little larger in the z direction 
   cpd_desired_result_of_optimization[li_index_of_focus+1] = 0.8*cd_value_at_focus; 
   cpd_desired_result_of_optimization[li_index_of_focus-1] = 0.8*cd_value_at_focus;
   cpd_desired_result_of_optimization[li_index_of_focus+2] = 0.4*cd_value_at_focus; 
   cpd_desired_result_of_optimization[li_index_of_focus-2] = 0.4*cd_value_at_focus;
#endif 
   /*   think it makes sense to force the "desired focus" to be equal to the actual focus. */
   cd_position_of_focus[0] = ld_actual_position_of_focus[0];  
   cd_position_of_focus[1] = ld_actual_position_of_focus[1];  
   cd_position_of_focus[2] = ld_actual_position_of_focus[2];  
            
#ifdef DEBUG_PSEUDO_INVERSE
      cout << "Finished setting up the desired result vector " << endl; 
      cout << "Index of focus  =  " << li_index_of_focus  << endl;
      cout << "True position  of focus (not necessarily the same as what you input)  " << ld_actual_position_of_focus[0]  << " " << ld_actual_position_of_focus[1] << " " << ld_actual_position_of_focus[2] << " " << endl;     
      cout << "Requested position  of focus  " <<  cd_position_of_focus[0] << " " << cd_position_of_focus[1] << " " << cd_position_of_focus[2] << " " << endl;
      cout << "Requested value at focus " <<    cd_value_at_focus   << endl;
      cout << "Dist of this point to focus = " << ld_min_dist_to_focus<< endl;
#endif   
      
   /* One thing that must be looked at: we use doubles, not complexes. So we have to map the matrix problem Ap = r to the problem 
   (A_r   -A_i) (p_r) = (r_r)
   (A_i   A_r ) (p_i)   (r_i)
   Also, the matrix A must be stored in column-major order, etc. A_r has  li_number_of_points_in_space rows and ci_number_of_elements columns. The offset to the first entry in A_i is therefore  li_number_of_points_in_space, that to the first of -A_i is 2*li_number_of_points_in_space * ci_number_of_elements and to the second A_r is 2*li_number_of_points_in_space * ci_number_of_elements + li_number_of_points_in_space. Write these 3 as a set of offsets and then increment the offset by  2*li_number_of_points_in_space for each transducer element.
   */
   int li_offset_top_left = 0;
   int li_offset_bottom_left = li_number_of_points_in_space;
   int li_offset_top_right = 2*li_number_of_points_in_space * ci_number_of_elements ;
   int li_offset_bottom_right = 2*li_number_of_points_in_space * ci_number_of_elements + li_number_of_points_in_space;
   
   for (int li_ind = 0; li_ind < ci_number_of_elements; li_ind++)
   {
      li_counter = 0;
      for (int li_x = 0; li_x <(ci_x_resolution_optimization+1) ; li_x++)
      {
         for (int li_y = 0; li_y <(ci_y_resolution_optimization+1) ; li_y++)
         {
            for (int li_z = 0; li_z <(ci_z_resolution_optimization+1) ; li_z++)
            {  // First, get the value of the output that you would have if the input pressure were 1 for the current transducer element and 0 for all the others. 
               // These are set up in Setup_Pressures()               
               ld_real_pressure = ca_real_array_for_optimization[li_ind].Output_Value(li_x,li_y,li_z);
               ld_imag_pressure = ca_imag_array_for_optimization[li_ind].Output_Value(li_x,li_y,li_z);
               // Next, input this pressure into the appropriate place in the big array. 
               cpd_array_for_transfer_matrix[li_offset_top_left + li_counter]  =  ld_real_pressure;
               cpd_array_for_transfer_matrix[li_offset_bottom_left + li_counter]  = ld_imag_pressure;
               cpd_array_for_transfer_matrix[li_offset_top_right + li_counter]  = -ld_imag_pressure;
               cpd_array_for_transfer_matrix[li_offset_bottom_right + li_counter]  =  ld_real_pressure;
               ++li_counter; // increment the counter variable.
#ifdef DEBUG_PSEUDO_INVERSE
      //cout << "Pressures at point " << li_counter << ": " << ld_real_pressure << ", "  << ld_imag_pressure<< endl;
      //cin >> lc_tmp; 
#endif   

               //  
            }
         }
      }
      // Move to the next pair of columns in the big transfer matrix. 
      li_offset_top_left += 2*li_number_of_points_in_space;
      li_offset_bottom_left += 2*li_number_of_points_in_space;
      li_offset_top_right += 2*li_number_of_points_in_space;
      li_offset_bottom_right += 2*li_number_of_points_in_space;      
   }
#ifdef DEBUG_PSEUDO_INVERSE
      cout << "Finished setting up the transfer matrix " << endl; 
#endif   

   
   // Now run the pseudoinverse. 
   // I'm breaking with my usual notation convention to have the same as the lapack convention.
   // ASSUME THAT M > N. i.e. that we have more points in space of interest than we have pressures. 
   int M = 2*li_number_of_points_in_space;
   int N = 2*ci_number_of_elements;
   if (N> M) std::cout << " *       ***************************************  WARNING ************************************ overdetermined problem, I'm not sure what to do " << std::endl;
   int nrhs = 1;
   int lda = M;
   int ldb = M; 
   double *S = new double[N];
   double rcond = -1.0; // This means that we use machine precision for deciding when singular values can be treated as 0. 
   int  rank ;
   int lwork = 5*M; // Need a pretty big workspace. This should (I guess ) fit into memory and be big enough. Even with M = 10^6, this is only 40 Mb, so we're not running low on memory.
   double *work = new double[lwork];
   int info;  
   
#ifdef DEBUG_PSEUDO_INVERSE
      cout << "About to use LAPACK routine " << endl; 
      cout << "M = " << M << endl;
      cout << "N = " << N << endl;
      cout << "lwork = " << lwork << endl;
    /*  for (int li_ind = 0; li_ind < M; li_ind++)
      {
         cout << li_ind << "  " << cpd_desired_result_of_optimization[li_ind] << ",    ";
      }
      cout << endl;
    */
#endif   
	 dgelss_(&M,&N,&nrhs, cpd_array_for_transfer_matrix, &lda, cpd_desired_result_of_optimization , &ldb, S, &rcond, &rank,  work, &lwork, &info);
   
#ifdef DEBUG_PSEUDO_INVERSE
      cout << "Finished with lapack " << endl;
   /*    for (int li_ind = 0; li_ind < N; li_ind++)
       {
          cout << li_ind << "  " << S[li_ind] << endl;
       }
   */  
#endif   
   
   if(info !=0) std::cout << " Trouble with dgelss, info = " << info << std::endl;
   else
   {
      for (int li_ind =  0; li_ind < 2*ci_number_of_elements; li_ind++)
      {
         apd_pressures[li_ind] = cpd_desired_result_of_optimization[li_ind]; // This is where the result is stored. 
      }
   } 

   delete[] cpd_desired_result_of_optimization;
   delete[] cpd_optimized_pressure;
   delete[] cpd_array_for_transfer_matrix;
   delete[] S;
   delete[] work;
}


/*
 ADDITIONAL STUFF
 1) Ian and I agree that the psuedoinverse is quite good, but lacks finesse. So we're going to try something else. I'm not making enough changes to justify a version change. In principle, some of the routines I add could be used to clean up the pseudoinverse code a little. But then I would need a new version number, since you don't change working code without adding a number. So I'm not going to bother. 
 2) The new method will involve solving the psuedoinverse problem in a more direct way, i.e. by writing it as a problem containing a happy square matrix, and inverting this. The advantage of this is that I can add in a secondary constraint to enforce maximizing the pressure at the transducers. This is helpful because we've previously found that several of the transducers switch off to give a marginal tightening of the focus. This is cool, but a bit dumb. Secondary constraints can't be done using the pseudoinverse. 
 3) Additionally, I want to experiment with penalizing secondary maxima. The smoothest way I can see how to do this is to premultiply with an exponential that exaggerates the importance of deviations from the desired result far from the focus. That way, bits near the focus aren't going to be mistaken for secondary foci (and we don't need to do a ridiculously slow derivative test!)
 4) So, I'm introducing a method that generates the transfer matrix, another that generates the desired result (vector of zeros with a nonzero at the focus) and a vector of spatial positions (actually 3 vectors, x,y,z). All useful methods, honestly!  
 
*/

void Optimize_Transducer::Generate_Desired_Output_Vector_OLD(double *aa_desired_output_vector)
{ // THIS FUNCTION USES A MAPPING FROM 3-D TO 1-D THAT IS ALSO USED IN GENERATING THE TRANSFER MATRIX AND THE VECTOR OF POSITIONS. IT IS THEREFORE ESSENTIAL THAT YOU DON'T FIDDLE WITH IT UNNECESSARILY. ANY CHANGES TO THE ORDERING USED HERE MUST BE REFLECTED IN THE OTHER FUNCTIONS!!!!! 
    // Generates the desired output, which has a nonzero entry in the region near the focus and is zero elsewhere. Remember that this must have length  2*li_number_of_points_in_space, because it has real and imaginary parts separate. 
      // SIZE OF ARRAYS MUST BE 2*li_number_of_points_in_space
   int li_number_of_points_in_space = (ci_x_resolution_optimization+1)*(ci_y_resolution_optimization+1)*(ci_z_resolution_optimization+1) ; // Note the +1 is quite important. Derives from my way of describing a resolution, i.e. resolution is number of divisions, not number of regions.
    
   // Now we step through the grid where we want to optimize the solution. A lot of the following apparent garbage is just for the sake of figuring out where we are in space so that we can look for the grid point that corresponds to the focus. 
   double ld_dx = (cd_xrange_of_optimized_region[1]-cd_xrange_of_optimized_region[0])/double(ci_x_resolution_optimization);
   double ld_dy = (cd_yrange_of_optimized_region[1]-cd_yrange_of_optimized_region[0])/double(ci_y_resolution_optimization);
   double ld_dz = (cd_zrange_of_optimized_region[1]-cd_zrange_of_optimized_region[0])/double(ci_z_resolution_optimization);
   double ld_x, ld_y, ld_z; 
   
   double ld_dist_to_focus ; // just a temporary starting point. 
   double ld_min_dist_to_focus = 10000.0; // just a temporary starting point.
   double ld_actual_position_of_focus[3] ; // a temporary variable just intended for output to screen. Basically, if the user uses a grid on which the focus doesn't appear, then we have to choose a nearby point at which to position the peak. This variable contains that point.  
   
   int li_counter, li_index_of_focus;
   li_counter = 0;
   ld_x = cd_xrange_of_optimized_region[0];
   // NOTE REALLY IMPORTANT: The order of x,y,z in the following loop defines the map from 3-d space to 1-d array. So it MUST be the same as that used in generating the transfer matrix. Otherwise you get garbage. You don't want garbage. 
   for (int li_x = 0; li_x <(ci_x_resolution_optimization+1) ; li_x++)
   {
      ld_y = cd_yrange_of_optimized_region[0];
      for (int li_y = 0; li_y <(ci_y_resolution_optimization+1) ; li_y++)
      {
         ld_z = cd_zrange_of_optimized_region[0];
         for (int li_z = 0; li_z <(ci_z_resolution_optimization+1) ; li_z++)
         {
            aa_desired_output_vector[li_counter] = 0.0; 
            aa_desired_output_vector[li_counter + li_number_of_points_in_space] = 0.0; 
            ld_dist_to_focus = (ld_x-cd_position_of_focus[0])*(ld_x-cd_position_of_focus[0])+(ld_y-cd_position_of_focus[1])*(ld_y-cd_position_of_focus[1])+(ld_z-cd_position_of_focus[2])*(ld_z- cd_position_of_focus[2]); // OK, this is the squared distance, but it's good enough.
            if(ld_dist_to_focus < ld_min_dist_to_focus) 
            {
               ld_min_dist_to_focus = ld_dist_to_focus ; //update the minimum.
               li_index_of_focus = li_counter; // get the appropriate index.
               ld_actual_position_of_focus[0] = ld_x;
               ld_actual_position_of_focus[1] = ld_y;
               ld_actual_position_of_focus[2] = ld_z;
#ifdef DEBUG_PSEUDO_INVERSE
                   //     cout << "Current position  " <<  ld_x << " " << ld_y << " " << ld_z << " " << endl;
                    //    cout << "Dist of this point to focus = " << ld_min_dist_to_focus<< endl;
                    //    cout << "Current index " << li_counter << endl;
#endif   
            } 
            
            ++li_counter;
            ld_z +=ld_dz; 
         }
         ld_y += ld_dy;
      }
      ld_x +=ld_dx;
   }
   aa_desired_output_vector[li_index_of_focus] = cd_value_at_focus;
   // And that's that! All nicely set up.  
}
            

void Optimize_Transducer::Generate_Desired_Output_Vector(double *aa_desired_output_vector)
{   // THIS FUNCTION USES A MAPPING FROM 3-D TO 1-D THAT IS ALSO USED IN GENERATING THE TRANSFER MATRIX AND THE VECTOR OF POSITIONS. IT IS THEREFORE ESSENTIAL THAT YOU DON'T FIDDLE WITH IT UNNECESSARILY. ANY CHANGES TO THE ORDERING USED HERE MUST BE REFLECTED IN THE OTHER FUNCTIONS!!!!! 
    // Generates the desired output, which has a nonzero entry in the region near the focus and is zero elsewhere. Remember that this must have length  2*li_number_of_points_in_space, because it has real and imaginary parts separate. 
    // SIZE OF ARRAYS MUST BE 2*li_number_of_points_in_space
   int li_number_of_points_in_space = (ci_x_resolution_optimization+1)*(ci_y_resolution_optimization+1)*(ci_z_resolution_optimization+1) ; 
   // Note the +1 is quite important. Derives from my way of describing a resolution, i.e. resolution is number of divisions, not number of regions.
    
   // Now we step through the grid where we want to optimize the solution. A lot of the following apparent garbage is just for the sake of figuring out where we are in space so that we can look for the grid point that corresponds to the focus. 
   double ld_dx = (cd_xrange_of_optimized_region[1]-cd_xrange_of_optimized_region[0])/double(ci_x_resolution_optimization);
   double ld_dy = (cd_yrange_of_optimized_region[1]-cd_yrange_of_optimized_region[0])/double(ci_y_resolution_optimization);
   double ld_dz = (cd_zrange_of_optimized_region[1]-cd_zrange_of_optimized_region[0])/double(ci_z_resolution_optimization);
   double ld_x, ld_y, ld_z; 
   
   double ld_dist_to_focus ; // just a temporary starting point. 
   int li_counter, li_index_of_focus;
   double ld_dist_to_geometric_focus;
   li_counter = 0;
   ld_x = cd_xrange_of_optimized_region[0];
   // NOTE REALLY IMPORTANT: The order of x,y,z in the following loop defines the map from 3-d space to 1-d array. So it MUST be the same as that used in generating the transfer matrix. Otherwise you get garbage. You don't want garbage. 
   for (int li_x = 0; li_x <(ci_x_resolution_optimization+1) ; li_x++)
   {
      ld_y = cd_yrange_of_optimized_region[0];
      for (int li_y = 0; li_y <(ci_y_resolution_optimization+1) ; li_y++)
      {
         ld_z = cd_zrange_of_optimized_region[0];
         for (int li_z = 0; li_z <(ci_z_resolution_optimization+1) ; li_z++)
         {
            aa_desired_output_vector[li_counter] = 0.0; 
            aa_desired_output_vector[li_counter + li_number_of_points_in_space] = 0.0; 
            ld_dist_to_focus = (ld_x-cd_position_of_focus[0])*(ld_x-cd_position_of_focus[0])+(ld_y-cd_position_of_focus[1])*(ld_y-cd_position_of_focus[1])+(ld_z-cd_position_of_focus[2])*(ld_z- cd_position_of_focus[2]); // OK, this is the squared distance, but it's good enough.
            if(ld_dist_to_focus < cd_desired_width_of_focus) 
            {
               ld_dist_to_geometric_focus = ld_x*ld_x +ld_y*ld_y +ld_z*ld_z ; 
               aa_desired_output_vector[li_counter] = cd_value_at_focus*exp(-cd_scale_factor_for_wide_focus*ld_dist_to_focus/cd_desired_width_of_focus)*exp(cd_bias_factor_for_wide_focus*ld_dist_to_geometric_focus);
               std::cout << "Setting up output: " << ld_x << " " <<ld_y << " " <<ld_z << " " <<cd_value_at_focus*exp(-cd_scale_factor_for_wide_focus*ld_dist_to_focus/cd_desired_width_of_focus)*exp(cd_bias_factor_for_wide_focus*ld_dist_to_geometric_focus) << endl;                  
            }             
            ++li_counter;
            ld_z +=ld_dz; 
         }
         ld_y += ld_dy;
      }
      ld_x +=ld_dx;
   }     
}

void Optimize_Transducer::Generate_Spatial_Vector_For_Optimization(double *apd_spatial_vector_x, double *apd_spatial_vector_y, double *apd_spatial_vector_z, double *apd_position_of_focus)
{   // THIS FUNCTION USES A MAPPING FROM 3-D TO 1-D THAT IS ALSO USED IN GENERATING THE TRANSFER MATRIX AND THE VECTOR OF DESIRED VALUES. IT IS THEREFORE ESSENTIAL THAT YOU DON'T FIDDLE WITH IT UNNECESSARILY. ANY CHANGES TO THE ORDERING USED HERE MUST BE REFLECTED IN THE OTHER FUNCTIONS!!!!! 
    // Generates the spatial vectors, which is useful because then we can weight non-focal peaks more than weighting near-focal peaks. Remember that these must have length  2*li_number_of_points_in_space, because all other vectors have real and imaginary parts separate, so we have a position for both reals and imaginaries. Also outputs the position of the true focus, i.e. the point on the grid closest to the desired focus, which we shall use as the real desired focus. 
   // SIZE OF ARRAYS MUST BE 2*li_number_of_points_in_space
   int li_number_of_points_in_space = (ci_x_resolution_optimization+1)*(ci_y_resolution_optimization+1)*(ci_z_resolution_optimization+1) ; // Note the +1 is quite important. Derives from my way of describing a resolution, i.e. resolution is number of divisions, not number of regions.
    
   // Now we step through the grid where we want to optimize the solution. A lot of the following apparent garbage is just for the sake of figuring out where we are in space so that we can look for the grid point that corresponds to the focus. 
   double ld_dx = (cd_xrange_of_optimized_region[1]-cd_xrange_of_optimized_region[0])/double(ci_x_resolution_optimization);
   double ld_dy = (cd_yrange_of_optimized_region[1]-cd_yrange_of_optimized_region[0])/double(ci_y_resolution_optimization);
   double ld_dz = (cd_zrange_of_optimized_region[1]-cd_zrange_of_optimized_region[0])/double(ci_z_resolution_optimization);
   double ld_x, ld_y, ld_z; 
   
   double ld_dist_to_focus ; // just a temporary starting point. 
   double ld_min_dist_to_focus = 10000.0; // just a temporary starting point.
    
   int li_counter, li_index_of_focus;
   li_counter = 0;
   ld_x = cd_xrange_of_optimized_region[0];
   // NOTE REALLY IMPORTANT: The order of x,y,z in the following loop defines the map from 3-d space to 1-d array. So it MUST be the same as that used in generating the transfer matrix. Otherwise you get garbage. You don't want garbage. 
   for (int li_x = 0; li_x <(ci_x_resolution_optimization+1) ; li_x++)
   {
      ld_y = cd_yrange_of_optimized_region[0];
      for (int li_y = 0; li_y <(ci_y_resolution_optimization+1) ; li_y++)
      {
         ld_z = cd_zrange_of_optimized_region[0];
         for (int li_z = 0; li_z <(ci_z_resolution_optimization+1) ; li_z++)
         {
            apd_spatial_vector_x[li_counter] = ld_x; 
            apd_spatial_vector_x[li_counter + li_number_of_points_in_space] = ld_x; 
            apd_spatial_vector_y[li_counter] = ld_y; 
            apd_spatial_vector_y[li_counter + li_number_of_points_in_space] = ld_y; 
            apd_spatial_vector_z[li_counter] = ld_z; 
            apd_spatial_vector_z[li_counter + li_number_of_points_in_space] = ld_z; 
            ld_dist_to_focus = (ld_x-cd_position_of_focus[0])*(ld_x-cd_position_of_focus[0])+(ld_y-cd_position_of_focus[1])*(ld_y-cd_position_of_focus[1])+(ld_z-cd_position_of_focus[2])*(ld_z- cd_position_of_focus[2]); // OK, this is the squared distance, but it's good enough.
            if(ld_dist_to_focus < ld_min_dist_to_focus) 
            {
               ld_min_dist_to_focus = ld_dist_to_focus ; //update the minimum.
               li_index_of_focus = li_counter; // get the appropriate index.
               apd_position_of_focus[0] = ld_x;
               apd_position_of_focus[1] = ld_y;
               apd_position_of_focus[2] = ld_z;
            } 
            
            ++li_counter;
            ld_z +=ld_dz; 
         }
         ld_y += ld_dy;
      }
      ld_x +=ld_dx;
   }
   /*  I think it makes sense to force the "desired focus" to be equal to the actual focus. */
   cd_position_of_focus[0] = apd_position_of_focus[0];  
   cd_position_of_focus[1] = apd_position_of_focus[1];  
   cd_position_of_focus[2] = apd_position_of_focus[2];  

}





void Optimize_Transducer::Generate_Transfer_Matrix(double *apd_transfer_matrix)
{ //  // THIS FUNCTION USES A MAPPING FROM 3-D TO 1-D THAT IS ALSO USED IN GENERATING THE SPATIAL VECTOR AND THE VECTOR OF DESIRED VALUES. IT IS THEREFORE ESSENTIAL THAT YOU DON'T FIDDLE WITH IT UNNECESSARILY. ANY CHANGES TO THE ORDERING USED HERE MUST BE REFLECTED IN THE OTHER FUNCTIONS!!!!! 
    // Generates the transfer matrix IN COLUMN MAJOR ORDER. 
    /* We use doubles, not complexes. So we have to map the matrix problem Ap = r to the problem 
   (A_r   -A_i) (p_r) = (r_r)
   (A_i   A_r ) (p_i)   (r_i)
   Also, the matrix A must be stored in column-major order, etc. A_r has  li_number_of_points_in_space rows and ci_number_of_elements columns. The offset to the first entry in A_i is therefore  li_number_of_points_in_space, that to the first of -A_i is 2*li_number_of_points_in_space * ci_number_of_elements and to the second A_r is 2*li_number_of_points_in_space * ci_number_of_elements + li_number_of_points_in_space. Write these 3 as a set of offsets and then increment the offset by  2*li_number_of_points_in_space for each transducer element.
   SIZE OF ARRAY MUST BE:  4*ci_number_of_elements*li_number_of_points_in_space
   */
     
   int li_number_of_points_in_space = (ci_x_resolution_optimization+1)*(ci_y_resolution_optimization+1)*(ci_z_resolution_optimization+1) ; // Note the +1 is quite important. Derives from my way of describing a resolution, i.e. resolution is number of divisions, not number of regions.
   int li_counter; // usual useful counter index. 
   
   // These are the 4 offsets to help me keep to column major order.
   int li_offset_top_left = 0;
   int li_offset_bottom_left = li_number_of_points_in_space;
   int li_offset_top_right = 2*li_number_of_points_in_space * ci_number_of_elements ;
   int li_offset_bottom_right = 2*li_number_of_points_in_space * ci_number_of_elements + li_number_of_points_in_space;
   
   double ld_real_pressure, ld_imag_pressure; 
   
   for (int li_ind = 0; li_ind < ci_number_of_elements; li_ind++)
   {
      li_counter = 0;
      for (int li_x = 0; li_x <(ci_x_resolution_optimization+1) ; li_x++)
      {
         for (int li_y = 0; li_y <(ci_y_resolution_optimization+1) ; li_y++)
         {
            for (int li_z = 0; li_z <(ci_z_resolution_optimization+1) ; li_z++)
            { // First, get the value of the output that you would have if the input pressure were 1 for the current transducer element and 0 for all the others. 
               // These are set up in Setup_Pressures()               
               ld_real_pressure = ca_real_array_for_optimization[li_ind].Output_Value(li_x,li_y,li_z);
               ld_imag_pressure = ca_imag_array_for_optimization[li_ind].Output_Value(li_x,li_y,li_z);
              // Next, input this pressure into the appropriate place in the big array. 
               apd_transfer_matrix[li_offset_top_left + li_counter]  =  ld_real_pressure;
               apd_transfer_matrix[li_offset_bottom_left + li_counter]  = ld_imag_pressure;
               apd_transfer_matrix[li_offset_top_right + li_counter]  = -ld_imag_pressure;
               apd_transfer_matrix[li_offset_bottom_right + li_counter]  =  ld_real_pressure;
               ++li_counter; // increment the counter variable.
            }
         }
      }
      // Move to the next pair of columns in the big transfer matrix. 
      li_offset_top_left += 2*li_number_of_points_in_space;
      li_offset_bottom_left += 2*li_number_of_points_in_space;
      li_offset_top_right += 2*li_number_of_points_in_space;
      li_offset_bottom_right += 2*li_number_of_points_in_space;      
   }
}





void Optimize_Transducer::Generate_Pressures_Using_Fancy_Method1(double *apd_pressures)
{ // Uses an attempt at a sophisticated method to generate the optimal pressures.  Writes the corresponding pressures into the array, real first then imaginary. Array must be twice as long as the number of transducer elements.  
   /*
   This particular fancy method considers the overdetermined problem Ap=u, aggravates problems far from the focus by multiplying with an exponential (with parameters fixed in the Optimize_Transducer class), call this operation T,  then minimizes |TAp - Tu|^2 - lambda p*p. The latter ensures that the pressures are as high as possible. lambda needs to be chosen such that the optimization routine is solvable and yields a minimum (rather than a non-minimal extremum). When explicitly doing the minimization, we have to solve the matrix equation (B^T B - lambda I)p = B^T u, where B = TA, and ^T means transpose (everything is real, so this is the same as Hermitian conjugate). lambda must be smaller than the smallest eigenvalue of C = B^T B. Really, we want C-lambda I to be well-conditioned, so we need to set up an acceptable condition number. This takes some guessing!
   Anyway, first attempt will be to say that lambda is such that the condition number is exactly equal to whatever condition number I've stuck into the given optimization file. In future, I may change this, maybe have lambda a parameter that the user can choose. Who knows?!
    
   */
   

   // First, get the required matrices and vectors. 
   int li_number_of_points_in_space = (ci_x_resolution_optimization+1)*(ci_y_resolution_optimization+1)*(ci_z_resolution_optimization+1) ; // Note the +1 is quite important. Derives from my way of describing a resolution, i.e. resolution is number of divisions, not number of regions.
   
   double *cpd_desired_result_of_optimization; 
   double *cpd_optimized_pressure;
   double *cpd_transfer_matrix; 
   double *cpd_spatial_vector_x; 
   double *cpd_spatial_vector_y; 
   double *cpd_spatial_vector_z; 
   double *cpd_actual_position_of_focus; // The focus that was chosen might not be on a grid point. So we force it to a grid point. 
   
   
   try{
   	 	cpd_desired_result_of_optimization = new double[2*li_number_of_points_in_space]; 
   		cpd_optimized_pressure = new double[2*ci_number_of_elements];
   		cpd_transfer_matrix = new double[4*ci_number_of_elements*li_number_of_points_in_space]; 
   		cpd_spatial_vector_x = new double[2*li_number_of_points_in_space]; 
   		cpd_spatial_vector_y = new double[2*li_number_of_points_in_space]; 
   		cpd_spatial_vector_z = new double[2*li_number_of_points_in_space]; 
   		cpd_actual_position_of_focus = new double[3]; // The focus that was chosen might not be on a grid point. So we force it to a grid point. 
   }
   catch (bad_alloc&)
   {
	   cerr << "Error allocating memory in Generate_Pressures_Using_Fancy_Method1." << endl;
	   exit(-2);
	   
   }
   this->Generate_Transfer_Matrix(cpd_transfer_matrix);
   this->Generate_Spatial_Vector_For_Optimization(cpd_spatial_vector_x, cpd_spatial_vector_y, cpd_spatial_vector_z, cpd_actual_position_of_focus);
   this->Generate_Desired_Output_Vector(cpd_desired_result_of_optimization);
   
   std::cout << "Actual position of focus is " <<  cpd_actual_position_of_focus[0] << ", " <<  cpd_actual_position_of_focus[1] << ", "<<  cpd_actual_position_of_focus[2] << std::endl;
   // OK, those are the basic elements of the problem. But the problem is far more exciting than just this stuff. 
   // We must multiply the transfer matrix  by an exponential function. In fact, I believe it must be absolute exponential, i.e. exp(|x|), so that we symmetrically exaggerate any non-focal peaks. 
   // Remember that we are using column-major order in the transfer matrix. 
   double ld_premultiplier, mu1,dist; 
   mu1  = cpd_optimization_parameters[1]; // This is the inverse width of the exponential. I've labelled it mu1, because if you want to introduce even fancier things, they should be called mu2, mu3, ... 
   for (int li_row = 0; li_row < li_number_of_points_in_space; li_row++)
   { // Each row corresponds to a different point in space. Actually, that's not quite true, because real and imaginary parts are obviously at the same points in space.
      dist = fabs(cpd_spatial_vector_x[li_row]-cpd_actual_position_of_focus[0]) + fabs(cpd_spatial_vector_y[li_row]-cpd_actual_position_of_focus[1])+  fabs(cpd_spatial_vector_z[li_row]-cpd_actual_position_of_focus[2]); 
      ld_premultiplier = exp(mu1*dist);           
      for (int li_col = 0; li_col < 2*ci_number_of_elements; li_col++)
      {
         cpd_transfer_matrix[li_row + 2*li_number_of_points_in_space*li_col] *= ld_premultiplier;  // This is the top half of the matrix, i.e. the bit that generates the real part of the solution 
         cpd_transfer_matrix[li_row + li_number_of_points_in_space + 2*li_number_of_points_in_space*li_col] *= ld_premultiplier; // This is the bottom half of the matrix. 
      }
   }
   // Now, if B is the premultiplied transfer matrix, we need to solve the matrix problem (B^T B - lambda I) p = B^T u, where p is the array of pressures, u is the output and lambda is equal to (as a first attempt) lambda = (max_allowable_condition_number*smallest_eigenvalue_of_B^T*B - largest_eigenvalue_of_B^T*B)/(max_allowable_condition_number - 1). 
   
   // First steps now are to get C = B^T B and x = B^T u. 
   
   /* However, it would be good to establish some sort of notation right up front. I'll use Lapack's notation. We are always solving problems with 1 right-hand-side, and we are going to need to refer to the rows and columns of all sorts of matrices and vectors. So let's have easy names for them. Specifically, M is the number of rows in matrix B, which is 2*li_number_of_points_in_space. N is the number of columns of B, i.e. 2*ci_number_of_elements. NRHS is the number of right-hand-sides, i.e. 1. INFO is an info integer, used to check that things worked. If info returns nonzero, quit with an error message, don't try to repair the damage. 
   I'm going to assume that LAPACK is happy with me submitting N as LDA, etc. i.e. repeated use of the same memory address within one function call, for input-only arguments. We'll see if that flies, if not, add in the extra variables.
   I am going to do the following steps:
   1) Generate C and x. 
   2) Factorize C into QTQ^T, where Q is orthogonal and T is tridiagonal (and symmetric and positive, because C is). T will be represented by two vectors, D and E (D main diag, E first superdiag). Q is written into C, as a weird set of rotation vectors. There is also a vector TAU that contains the scale factors of the rotations. D and E need to be copied, because they are destroyed a little later. 
   3) Compute eigenvalues of C, which are the same as the eigenvalues of T. Check the eigenvalues are all positive. Check the condition number. Assign lambda to be related to the condition number of T such that (T-lambda I) still has reasonable condition number. (Note T-lambda I is worse conditioned than T, because T's eigenvalues are positive and lambda brings the smallest eigenvalue closer to zero.)   
   4) Note that we now want to solve Q(T-lambda I)Q^T p = x, i.e. (T-lambda I)Q^T p = Q^T x. Generate  Q^T x and T-lambda I. 
   5) Solve (T-lambda I) r = Q^T x.
   6) Compute p = Q r. And we're done! 
   */
   
   // ASSUME THAT M > N. i.e. that we have more points in space of interest than we have pressures. 
   int M = 2*li_number_of_points_in_space;
   int N = 2*ci_number_of_elements;
   if (N> M) std::cout << " *       ***************************************  WARNING ************************************ overdetermined problem, I'm not sure what to do " << std::endl;
   int NRHS = 1;
   int INFO;  
   
   // 1)
      // First, generate C = B^T B using BLAS 3, DSYRK:
   char UPLO='U'; // i.e. upper part of matrix will be stored. 
   char TRANS='T';
   double alpha = 1.0; 
   double beta = 0.0;
   double *C = new double[N*N]; // space for matrix C. 
   dsyrk_(&UPLO, &TRANS, &N, &M, &alpha, cpd_transfer_matrix, &M, &beta, C, &N);
      //
      // Next, generate x = B^T u, using BLAS 2, DGEMV: 
   alpha = 1.0; beta = 0.0; TRANS='T';    // Just to be sure that there's no cunning way of the machine changing my parameters:
   int INC = 1; // increment between consecutive elements. Not really an issue here.
   double *x = new double[N]; 
   dgemv_(&TRANS,&M, &N, &alpha, cpd_transfer_matrix, &M,cpd_desired_result_of_optimization, &INC, &beta,  x, &INC);
  // for (int ind = 0; ind < N; ind++) std::cout << ind << ", " << x[ind] << std::endl; 
   //std::cout << "************************************************" << std::endl;
    
   
   // 2) 
      // Factorize C into Q T Q^T using DSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
   UPLO = 'U'; // C is stored in upper part of matrix.
   double *D = new double[N]; // main diagonal of tridiag matrix.
   double *E = new double[N-1]; // first superdiagonal of tridiag  matrix.
   double *TAU = new double [N-1]; // The scalar coefficients of the rotations that make up Q. 
   int LWORK = -1; // This is the size of the work vector. I'm going to use LAPACK's builtin stuff to tell me what size this should be.
   double *WORK = new double[1]; // For now, just need 1 entry. 
   dsytrd_(&UPLO, &N, C, &N, D, E, TAU, WORK, &LWORK, &INFO); // do it with LWORK = -1 to get a work-space size query. WORK[0] is now the optimum size. 
   if (INFO) {      std::cout << "Failed to get ideal workspace" << INFO << std::endl;      return;   }
   LWORK = int(WORK[0]+0.5); // round off upward, just in case it's messing me around.
   delete WORK; // free up the memory WORK points to.
   WORK =  new double[LWORK]; // Grab some more memory for WORK.
   dsytrd_(&UPLO, &N, C, &N, D, E, TAU, WORK, &LWORK, &INFO); // This time, it's for real. D and E are our tridiag bits, C is overwritten with Q but in a weird symbolic fashion.
   if (INFO) {      std::cout << "Failed to reduce matrix C" << INFO << std::endl;      return;   }
   
   // 3)
      // First, backup D and E. Or rather, just copy them into scrap vectors:
   double *D2 = new double[N];
   double *E2 = new double[N-1];
   int N_minus_1 = N-1; // Stupid, yes, but recall that FORTRAN insists on pass-by-reference. So we have no choice. 
   dcopy_(&N,D,&INC,D2,&INC); // copy D into D2.
   dcopy_(&N_minus_1,E,&INC,E2,&INC); // copy E into E2.
   // Now get the eigenvalues of T with LAPACK:  DSTERF( N, D, E, INFO )
   dsterf_(&N, D2, E2, &INFO);
   if (INFO) {      std::cout << "Couldn't get eigenvalues. Sorry. " << INFO << std::endl;      return;   }
   
   // 4) 
      // First, get lambda
   double lambda;
   /*  I've realised that if I try to maximize the pressure squared at the transducer, I fall prey to one or another element getting a huge value and leaving the others switched off. Therefore I want to minimize the differences between different elements. However, this involves using a difference matrix, which doesn't commute with Q, and then everything needs to be redone. And I can't use superfast tridiagonal routines. So instead, let's MINIMIZE the pressure squared. Sound stupid? Yes. But if you minimize a positive quantity, you penalize the largest elements more than the smallest. It suddenly becomes preferable for the system to distribute things as evenly as possible. And that's what I want! Remember to normalize the result and you should be apples.
   Anyway, what this means is that lambda must be negative. From Numerical Recipes, a reasonable value is Tr(C^T C)/Tr(1^2) = sum of squares of eigenvalues of C divided by N. 
   
   DEPRECATED BIT, MAXIMIZES PRESSURE SQUARED.
   double eps0 = D2[0]; // smallest eigenvalue 
   double epsN = D2[N-1];  // largest eigenvalue. 
   if (eps0 <= 0.0)    {      std::cout << "I thought this was supposed to be positive!!"  << std::endl;      std::cout << "Smallest e-value is "  << eps0 << std::endl;
      std::cout << "Largest e-value is "  << epsN << std::endl;      return;   }
   if (epsN/eps0 >= cd_condition_number)
   { // Then the problem is badly-enough conditioned that sticking in an extra lambda would mess everything up.
      lambda = 0.0; 
   }
   else 
   {
      lambda = (cd_condition_number*eps0 - epsN)/(cd_condition_number - 1.0); 
   }
   //lambda=  -lambda;
   END OF DEPRECATED BIT, MAXIMIZES PRESSURE SQUARED.   
   */
   
  
	lambda = - ddot_(&N, D2, &INC, D2, &INC)/N; // get the sum of squares of eigenvalues. 
   
   lambda *= cpd_optimization_parameters[0]; // Introduce the ability to scale lambda without recompiling.
   
   
   std::cout << "Lambda is " << lambda << std::endl;
   //std::cout << cd_condition_number << "  "  << eps0  << "  " <<  epsN << std::endl;
   std::cout << "Smallest eigenvalue  is " << D2[0] << std::endl;
   //for (int ind = 0; ind < N; ind++) std::cout << ind << ", " << D2[ind] << ", " << D[ind] << std::endl;
      // Now get T-lambda I:
   for (int ind = 0; ind < N; ind++)
   {
      D[ind] -= lambda;
   } 
      // Now get Q^T x  using LAPACK:  DORMTR( SIDE, UPLO, TRANS, M, N, A, LDA, TAU, C, LDC, WORK, LWORK, INFO ) 
      // This function uses Q in the form returned by dsytrd, so we just pass in the matrix C and it should all work.
      // Overwrites the result of Q^T x into the vector x.  
   char SIDE = 'L';
   UPLO = 'U'; TRANS = 'T'; // i.e. Q is stored in upper diagonal, and we are calculating Q^T x.    
   std::cout << "Inc is " << INC << std::endl;
   dormtr_(&SIDE, &UPLO, &TRANS, &N, &NRHS, C, &N, TAU, x, &N, WORK, &LWORK, &INFO); // Note that I'm using the LWORK previously optimized. This routine actually uses a much smaller workspace, but there's no point in deleting the big vector and getting a little one; I'm not so low on memory. 
   if (INFO) {      std::cout << "Couldn't construct Q^T x.  Sorry. " << INFO << std::endl;      return;   }

   
   // 5) Solve (T-lambda I) r = Q^T x.
      // This uses  DPTSV( N, NRHS, D, E, B, LDB, INFO ). T-lambda I should be positive, symmetric, tridiagonal. If it's not, this is not going to work.
   dptsv_( &N, &NRHS, D, E, x, &N, &INFO ); 
   if (INFO) {  std::cout << " Couldn't invert the tridiagonal matrix.  Sorry. " << INFO << std::endl;      return;   }
   // result is copied into x. 

   // 6) Finally, get the pressure p = Q r, using DORMTR again. 
   SIDE = 'L'; UPLO = 'U'; TRANS = 'N'; // Multiply from the left, matrix is stored in upper half, no transpose this time. 
   dormtr_(&SIDE, &UPLO, &TRANS, &N, &NRHS, C, &N, TAU, x, &N, WORK, &LWORK, &INFO);  
   if (INFO) {  std::cout << " Well, you failed at the final hurdle. Couldn't invert the orthogonal matrix Q.  Sorry. " << INFO << std::endl;      return;   }   
    

   // Finally, copy the result across. 
   dcopy_(&N,x,&INC,apd_pressures,&INC);

   delete[] cpd_desired_result_of_optimization;
   delete[] cpd_optimized_pressure;
   delete[] cpd_transfer_matrix;
   delete[] cpd_spatial_vector_x;
   delete[] cpd_spatial_vector_y;
   delete[] cpd_spatial_vector_z;
   delete[] cpd_actual_position_of_focus;
   delete[] C;
   delete[] x;
   delete[] D;
   delete[] E;
   delete[] D2;
   delete[] E2;
   delete[] TAU;
   delete[] WORK;
   
}












void Optimize_Transducer::Generate_Pressures_Using_Fancy_Method2(double *apd_pressures)
{ // Uses an attempt at a sophisticated method to generate the optimal pressures.  Writes the corresponding pressures into the array, real first then imaginary. Array must be twice as long as the number of transducer elements.  
   /*
   This particular fancy method considers the overdetermined problem Ap=u, aggravates problems far from the focus by multiplying with an exponential (with parameters fixed in the Optimize_Transducer class), call this operation T,  then minimizes |TAp - Tu|^2 + lambda p^T*H*p, where H is the matrix with N-1 on the diagonal (N = number of elements) and -1 everywhere else; actually, this is not quite true, because of the whole real-imaginary divide, and H is as described except with zeros coupling real to imaginary parts. The latter ensures that the variation in pressure on the transducer is as low as possible. lambda needs to be chosen such that the optimization routine is solvable and yields a minimum (rather than a non-minimal extremum), and it can be any positive value. When explicitly doing the minimization, we have to solve the matrix equation (B^T B + lambda H)p = B^T u, where B = TA, and ^T means transpose (everything is real, so this is the same as Hermitian conjugate). Lambda is unfortunately kind-of a free parameter, as long as it is positive. So let the user choose it. A reasonable guess would be Tr(H^T H)/Tr(C^T C) where C = B^T B; this will tend to balance the desire to minimize the least squares distance to desired result and to minimize the pressure variation. 
    
   */
   

   // First, get the required matrices and vectors. 
   int li_number_of_points_in_space = (ci_x_resolution_optimization+1)*(ci_y_resolution_optimization+1)*(ci_z_resolution_optimization+1) ; // Note the +1 is quite important. Derives from my way of describing a resolution, i.e. resolution is number of divisions, not number of regions.
   double *cpd_desired_result_of_optimization = new double[2*li_number_of_points_in_space]; 
   double *cpd_optimized_pressure = new double[2*ci_number_of_elements];
   double *cpd_transfer_matrix = new double[4*ci_number_of_elements*li_number_of_points_in_space]; 
   double *cpd_spatial_vector_x = new double[2*li_number_of_points_in_space]; 
   double *cpd_spatial_vector_y = new double[2*li_number_of_points_in_space]; 
   double *cpd_spatial_vector_z = new double[2*li_number_of_points_in_space]; 
   double *cpd_actual_position_of_focus = new double[3]; // The focus that was chosen might not be on a grid point. So we force it to a grid point. 
   
   this->Generate_Transfer_Matrix(cpd_transfer_matrix);
   this->Generate_Spatial_Vector_For_Optimization(cpd_spatial_vector_x, cpd_spatial_vector_y, cpd_spatial_vector_z, cpd_actual_position_of_focus);
   this->Generate_Desired_Output_Vector(cpd_desired_result_of_optimization);
   
   std::cout << "Actual position of focus is " <<  cpd_actual_position_of_focus[0] << ", " <<  cpd_actual_position_of_focus[1] << ", "<<  cpd_actual_position_of_focus[2] << std::endl;
   // OK, those are the basic elements of the problem. But the problem is far more exciting than just this stuff. 
   // We must multiply the transfer matrix  by an exponential function. In fact, I believe it must be absolute exponential, i.e. exp(|x|), so that we symmetrically exaggerate any non-focal peaks. 
   // Remember that we are using column-major order in the transfer matrix. 
   double ld_premultiplier, mu1,dist; 
   mu1  = cpd_optimization_parameters[1]; // This is the inverse width of the exponential. I've labelled it mu1, because if you want to introduce even fancier things, they should be called mu2, mu3, ... 
   for (int li_row = 0; li_row < li_number_of_points_in_space; li_row++)
   { // Each row corresponds to a different point in space. Actually, that's not quite true, because real and imaginary parts are obviously at the same points in space.
      dist = fabs(cpd_spatial_vector_x[li_row]-cpd_actual_position_of_focus[0]) + fabs(cpd_spatial_vector_y[li_row]-cpd_actual_position_of_focus[1])+  fabs(cpd_spatial_vector_z[li_row]-cpd_actual_position_of_focus[2]); 
      ld_premultiplier = (mu1*dist*dist) + 1 ;
      // cout << "Spatial vectors: " << cpd_spatial_vector_x[li_row] << " " <<           cpd_spatial_vector_y[li_row] << " " <<cpd_spatial_vector_z[li_row] << " " <<dist << endl;
      for (int li_col = 0; li_col < 2*ci_number_of_elements; li_col++)
      {
         cpd_transfer_matrix[li_row + 2*li_number_of_points_in_space*li_col] *= ld_premultiplier;  // This is the top half of the matrix, i.e. the bit that generates the real part of the solution 
         cpd_transfer_matrix[li_row + li_number_of_points_in_space + 2*li_number_of_points_in_space*li_col] *= ld_premultiplier; // This is the bottom half of the matrix. 
         //cout << "premultiplied by " << ld_premultiplier << endl; 
      }
   }
   // Now, if B is the premultiplied transfer matrix, we need to solve the matrix problem (B^T B + lambda H) p = B^T u, where p is the array of pressures and u is the output. 
   // First steps now are to get C = B^T B and x = B^T u. 
   /* However, it would be good to establish some sort of notation right up front. I'll use Lapack's notation. We are always solving problems with 1 right-hand-side, and we are going to need to refer to the rows and columns of all sorts of matrices and vectors. So let's have easy names for them. Specifically, M is the number of rows in matrix B, which is 2*li_number_of_points_in_space. N is the number of columns of B, i.e. 2*ci_number_of_elements. NRHS is the number of right-hand-sides, i.e. 1. INFO is an info integer, used to check that things worked. If info returns nonzero, quit with an error message, don't try to repair the damage. 
   I'm going to assume that LAPACK is happy with me submitting N as LDA, etc. i.e. repeated use of the same memory address within one function call, for input-only arguments. We'll see if that flies, if not, add in the extra variables.
   I am going to do the following steps:
   1) Generate C and x. 
   2) Generate C+lambda H. 
   3) Solve (C+lambda H) p = x.
   */
   // ASSUME THAT M > N. i.e. that we have more points in space of interest than we have pressures. 
   int M = 2*li_number_of_points_in_space;
   int N = 2*ci_number_of_elements;
   if (N> M) std::cout << " *       ***************************************  WARNING ************************************ overdetermined problem, I'm not sure what to do " << std::endl;
   int NRHS = 1;
   int INFO;  
   
   // 1)
      // First, generate C = B^T B using BLAS 3, DSYRK:
   char UPLO='U'; // i.e. upper part of matrix will be stored. 
   char TRANS='T';
   double alpha = 1.0; 
   double beta = 0.0;
   double *C = new double[N*N]; // space for matrix C. 
   dsyrk_(&UPLO, &TRANS, &N, &M, &alpha, cpd_transfer_matrix, &M, &beta, C, &N);
   //
      // Next, generate x = B^T u, using BLAS 2, DGEMV: 
   alpha = 1.0; beta = 0.0; TRANS='T';    // Just to be sure that there's no cunning way of the machine changing my parameters:
   int INC = 1; // increment between consecutive elements. Not really an issue here.
   double *x = new double[N]; 
   dgemv_(&TRANS,&M, &N, &alpha, cpd_transfer_matrix, &M,cpd_desired_result_of_optimization, &INC, &beta,  x, &INC);
  // for (int ind = 0; ind < N; ind++) std::cout << ind << ", " << x[ind] << std::endl; 
   //std::cout << "************************************************" << std::endl;
    
   double lambda = 0.0;
   
   for (int ind1 = 0; ind1 < N; ind1++)
   {
      lambda += C[ind1*N+ind1]; // Get the trace of C. 
   }
   lambda = fabs(lambda/(2.0*ci_number_of_elements*(ci_number_of_elements-1))); // Now lambda is Tr(C)/Tr(H). Seems appropriate.
   
   
   lambda *= cpd_optimization_parameters[0]; // Scale lambda as you like.  
   std::cout << "Lambda is " << lambda << std::endl;
      

   // 2) 
      // Now get C + lambda H: Only the upper triangular part of C is relevant (because I chose UPLO = "U" above) and since H is symmetric, the same holds for the result. So only bother with above axis things. 
   // Recall H is:  number_of_elements-1's on the diagonal,  -1's everywhere else in the real-real or imaginary-imaginary blocks and zero elsewhere. 
   for (int ind1 = 0; ind1 < N; ind1++)
   {
      C[ind1*N+ind1] += lambda*(ci_number_of_elements-1); // ind1*N+ind1 is the main diagonal. 
   } 
   for (int ind1 = 0; ind1 < ci_number_of_elements; ind1++)
   {
      for (int ind2 = 0; ind2 < ind1; ind2++)
      {
         C[ind1*N+ind2] -= lambda; // This is the top left block of the matrix. 
         C[N*ci_number_of_elements + ci_number_of_elements + ind1*N+ind2] -= lambda; // This is the bottom right block of the matrix. 
      } 
   } 
   
   
   
   // 3) Solve (C+lambda H) p = x.
      // This uses DPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO ), because C+lambda H is symmetric and positive definite. If it's not positive definite, we're not going to find a minimum, but a useless extremum. So we can rely on this positivity. Could have a zero, I guess. Oh well, see how it goes.
   UPLO= 'U'; NRHS=1; 
   dposv_(&UPLO, &N, &NRHS, C, &N, x, &N, &INFO ); // Pretty easy, yes? 
   if (INFO) {  std::cout << " Couldn't invert the  matrix.  Sorry. " << INFO << std::endl;      return;   }
   // result is copied into x. 
  
   // Finally, copy the result across. 
   dcopy_(&N,x,&INC,apd_pressures,&INC);


    
   
   
   delete[] cpd_desired_result_of_optimization;
   delete[] cpd_optimized_pressure;
   delete[] cpd_transfer_matrix;
   delete[] cpd_spatial_vector_x;
   delete[] cpd_spatial_vector_y;
   delete[] cpd_spatial_vector_z;
   delete[] cpd_actual_position_of_focus;
   delete[] C;
   delete[] x;
   
}








void Optimize_Transducer::Generate_Pressures_Using_Fancy_Method3(double *apd_pressures)
{ // Uses an attempt at a sophisticated method to generate the optimal pressures.  Writes the corresponding pressures into the array, real first then imaginary. Array must be twice as long as the number of transducer elements.  
   /*
	This particular fancy method considers the overdetermined problem Ap=u. I want to have the fluxes on each element limited to < 10W/cm^2, while the flux at the focus should be 50 000 W/cm^2. My attempt at constraining the fluxes of the elements is to minimize e^(|p|^2 - 10). This does prefer small intensities, but strongly penalizes intensities above 10. Possibly I should use 5, not 10, or some other reduced value. But see how it goes. There will be a Lagrange parameter lambda giving the weighting of the flux-limiter relative to the pressure distribution minimizer, and there will be a parameter that is somewhere around 10, limiting the value of the flux. Choosing these parameters may be a pain. 
    
	Note that this is a nonlinear problem. Therefore iterate it and see what happens.
	NOTE: First iteration attempt was to put all nonlinear stuff on the right. i.e.
		B^T B p - B^T phi + lambda f = 0, 
		where f_i = p_i exp(|p_i|^2 - mu).
		Therefore, iterate:
		B^T B p^n = B^T phi - lambda f^{n-1}.
	This failed horribly - needs extremely small lambda to have a hope of converging (lambda = 1e-12 is much too big...) So it will be difficult to find a lambda giving stability and accuracy.
	 
	Next attempt: Solve 
		(B^T B + lambda exp(|p_i|^2 - mu) delta_ij) p = B^T phi = x. 
		Now the pressure is contained in the matrix, which should stabilize things, but it means we have to keep regenerating matrices. 
		We introduce a matrix C that contains B^T B and a vector x that contains B^T phi. Then we just solve it. 
   */
    
   // First, get the required matrices and vectors. 
	int li_number_of_points_in_space = (ci_x_resolution_optimization+1)*(ci_y_resolution_optimization+1)*(ci_z_resolution_optimization+1) ; // Note the +1 is quite important. Derives from my way of describing a resolution, i.e. resolution is number of divisions, not number of regions.
	double *cpd_desired_result_of_optimization = new double[2*li_number_of_points_in_space]; 
	double *lpd_optimized_pressure = new double[2*ci_number_of_elements];
	double *cpd_transfer_matrix = new double[4*ci_number_of_elements*li_number_of_points_in_space]; 
	double *cpd_spatial_vector_x = new double[2*li_number_of_points_in_space]; 
	double *cpd_spatial_vector_y = new double[2*li_number_of_points_in_space]; 
	double *cpd_spatial_vector_z = new double[2*li_number_of_points_in_space]; 
	double *cpd_actual_position_of_focus = new double[3]; // The focus that was chosen might not be on a grid point. So we force it to a grid point. 
   
	this->Generate_Transfer_Matrix(cpd_transfer_matrix);
	this->Generate_Spatial_Vector_For_Optimization(cpd_spatial_vector_x, cpd_spatial_vector_y, cpd_spatial_vector_z, cpd_actual_position_of_focus);
	this->Generate_Desired_Output_Vector(cpd_desired_result_of_optimization);
   
	std::cout << "Actual position of focus is " <<  cpd_actual_position_of_focus[0] << ", " <<  cpd_actual_position_of_focus[1] << ", "<<  cpd_actual_position_of_focus[2] << std::endl;
   
    
   // ASSUME THAT M > N. i.e. that we have more points in space of interest than we have pressures. 
	int M = 2*li_number_of_points_in_space;
	int N = 2*ci_number_of_elements;
	if (N> M) std::cout << " *       ***************************************  WARNING ************************************ overdetermined problem, I'm not sure what to do " << std::endl;
	const int NRHS = 1;
	int INFO;  
   
   // 1)
      // First, generate C = B^T B using BLAS 3, DSYRK:
	char UPLO='U'; // i.e. upper part of matrix will be stored. 
	char TRANS='T';
	double alpha = 1.0; 
	double beta = 0.0;
	double *C = new double[N*N]; // space for matrix C. 
	dsyrk_(&UPLO, &TRANS, &N, &M, &alpha, cpd_transfer_matrix, &M, &beta, C, &N);
	//
      // Next, generate x = B^T u, using BLAS 2, DGEMV: 
   alpha = 1.0; beta = 0.0; TRANS='T';    // Just to be sure that there's no cunning way of the machine changing my parameters:
   const int INC = 1; // increment between consecutive elements. Not really an issue here.
   double *x = new double[N]; 
   dgemv_(&TRANS,&M, &N, &alpha, cpd_transfer_matrix, &M,cpd_desired_result_of_optimization, &INC, &beta,  x, &INC);
  // for (int ind = 0; ind < N; ind++) std::cout << ind << ", " << x[ind] << std::endl; 
   //std::cout << "************************************************" << std::endl;
    
   //  C and x are backups; don't work directly with them. Instead of x, work with lpd_optimized_pressure, and instead of C lets introduce a matrix A.
   double *A = new double[N*N]; // space for matrix A.
   double lambda = cpd_optimization_parameters[0];
   double mu = cpd_optimization_parameters[1];
   double mixing = cpd_optimization_parameters[2];
   double one_minus_mixing = 1.0 - mixing; 
   std::cout << "lambda " << lambda << std::endl;
   std::cout << "mu " << mu << std::endl;
   std::cout << "mixing " << mixing << std::endl;
   double minus_lambda = -lambda; // Dumb? Have to pass by reference!
   double *lpd_pressure_old = new double[2*ci_number_of_elements]; // Keep track of previous result; use for iterations.
   double diff = 100.0; // Needs to start somewhere.
   const int N2 = N*N; // Useful for copying matrices, because we can just copy them as vectors. Because C is better than Fortran...
   double mod_pressure;
   double *WORK = new double[1]; // For now, just need 1 entry. 
   int *IPIV = new int[N];
   int LWORK = -1; // Gotta get correct workspace size. 
   dsysv_(&UPLO, &N, &INC, A, &N, IPIV, lpd_optimized_pressure, &N, WORK, &LWORK, &INFO); // Doesn't matter that A and lpd_optimized_pressure are empty. All this does is generate the ideal size for WORK.
   LWORK = int(WORK[0]+0.5); // round off upward, just in case it's messing me around.
   delete WORK; // free up the memory WORK points to.
   WORK =  new double[LWORK]; // Grab some more memory for WORK.
   memset(lpd_pressure_old,0,2*ci_number_of_elements*sizeof(double)); // Clear old pressure.
   char tmp; 
   double alph; // A desymmetrizing factor that makes |p|^2=mu desirable, but would rather undershoot than overshoot.
   // Now start iteratin'
   while (diff>1.e-5){ // diff measures the change between two successive iterations.
	   dcopy_(&N2,C,&INC,A,&INC); // copy C into A.
	   for (int p = 0; p < ci_number_of_elements; p++){
		   mod_pressure = lpd_pressure_old[p]*lpd_pressure_old[p] + lpd_pressure_old[p+ci_number_of_elements]*lpd_pressure_old[p+ci_number_of_elements];
		   // Recall that your matrix is twice the usual size because I've split real and imaginary parts. So there's a bit of fiddlin' to be done.
		   alph = (mod_pressure>mu) ? 1.0: -0.05; // If pressure is small, don't penalize too much.		    
		   A[p + p*2*ci_number_of_elements] += lambda*exp(alph*(mod_pressure - mu));
		   A[p + p*2*ci_number_of_elements+ ci_number_of_elements*(2*ci_number_of_elements+1)] += lambda*exp(alph*(mod_pressure - mu));
		   lpd_optimized_pressure[p] = x[p];		   
		   lpd_optimized_pressure[p+ ci_number_of_elements] = x[p+ ci_number_of_elements];		   
	   }
	   // Matrix A  and vector lpd_optimized_pressure are now ready for action.
	   UPLO = 'U';
	   dsysv_(&UPLO, &N, &INC, A, &N, IPIV, lpd_optimized_pressure, &N, WORK, &LWORK, &INFO); // solve the problem. 
	   if (INFO) std::cerr << "Trouble, can't solve the problem " << std::endl;
	   diff = 0.0;
	   for (int p = 0 ;p < 2*ci_number_of_elements; p++ ) {
		   diff += (lpd_optimized_pressure[p]-lpd_pressure_old[p])*(lpd_optimized_pressure[p]-lpd_pressure_old[p]);
		   /*if ((lpd_optimized_pressure[p]-lpd_pressure_old[p])*(lpd_optimized_pressure[p]-lpd_pressure_old[p])>1.e-8){
			   std::cout << p <<std::endl;
			   std::cout << lpd_optimized_pressure[p] <<std::endl;
			   std::cout << lpd_pressure_old[p] <<std::endl;
			   
		   }
		   */	   
	   }
	   // Don't do a full copy. Do a mixed copy . DEPRECATED dcopy_(&N,lpd_optimized_pressure,&INC,lpd_pressure_old,&INC); // Must update old one
	   
	   dscal_(&N, &one_minus_mixing, lpd_pressure_old, &INC); // Multiply old pressure by 1-mixing.
	   daxpy_(&N, &mixing, lpd_optimized_pressure, &INC, lpd_pressure_old, &INC); // Add in mixing*new pressure.  
	   std::cout << "Difference " << diff << std::endl;
	   //std::cin >> tmp;  
   }

   
   
   
   
   
   
   
   /* DEPRECATED
   
   
   // 2) 
      // Factorize C into Q T Q^T using DSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
   UPLO = 'U'; // C is stored in upper part of matrix.
   double *D = new double[N]; // main diagonal of tridiag matrix.
   double *E = new double[N-1]; // first superdiagonal of tridiag  matrix.
   double *TAU = new double [N-1]; // The scalar coefficients of the rotations that make up Q. 
    // This is the size of the work vector. I'm going to use LAPACK's builtin stuff to tell me what size this should be.
   
   dsytrd_(&UPLO, &N, C, &N, D, E, TAU, WORK, &LWORK, &INFO); // do it with LWORK = -1 to get a work-space size query. WORK[0] is now the optimum size. 
   if (INFO) {      std::cout << "Failed to get ideal workspace" << INFO << std::endl;      return;   }
   LWORK = int(WORK[0]+0.5); // round off upward, just in case it's messing me around.
   delete WORK; // free up the memory WORK points to.
   WORK =  new double[LWORK]; // Grab some more memory for WORK.
   dsytrd_(&UPLO, &N, C, &N, D, E, TAU, WORK, &LWORK, &INFO); // This time, it's for real. D and E are our tridiag bits, C is overwritten with Q but in a weird symbolic fashion.
   if (INFO) {      std::cout << "Failed to reduce matrix C" << INFO << std::endl;      return;   }
   
   // 3)
      // First, backup D and E. Or rather, just copy them into scrap vectors:
   double *D2 = new double[N];
   double *E2 = new double[N-1];
   int N_minus_1 = N-1; // Stupid, yes, but recall that FORTRAN insists on pass-by-reference. So we have no choice. 
   dcopy_(&N,D,&INC,D2,&INC); // copy D into D2.
   dcopy_(&N_minus_1,E,&INC,E2,&INC); // copy E into E2.
   
   // Get eigenvalues.
   dsterf_(&N, D, E, &INFO);
   std::cout << "Eigenvalue 0 is " << D[0] << std::endl;
   std::cout << "Eigenvalue 1 is " << D[1] << std::endl;
   std::cout << "Eigenvalue N-1 is " << D[N-1] << std::endl;
   // Restore vectors
   dcopy_(&N,D2,&INC,D,&INC); // copy D into D2.
   dcopy_(&N_minus_1,E2,&INC,E,&INC); // copy E into E2.
   
   
   
    
         // Now get Q^T x  using LAPACK:  DORMTR( SIDE, UPLO, TRANS, M, N, A, LDA, TAU, C, LDC, WORK, LWORK, INFO ) 
      // This function uses Q in the form returned by dsytrd, so we just pass in the matrix C and it should all work.
      // Overwrites the result of Q^T x into the vector x.  
   char SIDE = 'L';
   UPLO = 'U'; TRANS = 'T'; // i.e. Q is stored in upper diagonal, and we are calculating Q^T x.    
   //std::cout << "Inc is " << INC << std::endl;
   dormtr_(&SIDE, &UPLO, &TRANS, &N, &NRHS, C, &N, TAU, x, &N, WORK, &LWORK, &INFO); // Note that I'm using the LWORK previously optimized. This routine actually uses a much smaller workspace, but there's no point in deleting the big vector and getting a little one; I'm not so low on memory. 
   if (INFO) {      std::cout << "Couldn't construct Q^T x.  Sorry. " << INFO << std::endl;      return;   }
   */
   
   /*
   OK, at this point it is worth refreshing things a bit:
   C contains the matrix Q in a weird symbolic fashion.
   x contains the vector Q^T x. 
   D2 and D contain the main diagonal of T.
   E2 and E contain the superdiagonal of T.
   
   I shall use
   cpd_optimized_pressure to contain the new pressure values (basically it is a work variable)
   cpd_pressure_old  contains the pressure values from the previous iteration.
   */
   /* DEPRECATED 
  	   // Generate f^n
	   for (int p = 0; p < ci_number_of_elements; p++){
		   mod_pressure = lpd_pressure_old[p]*lpd_pressure_old[p] + lpd_pressure_old[p+ci_number_of_elements]*lpd_pressure_old[p+ci_number_of_elements]; 
		   lpd_optimized_pressure[p] = lpd_pressure_old[p]*exp(mod_pressure - mu);
		   lpd_optimized_pressure[p+ci_number_of_elements] = lpd_pressure_old[p+ci_number_of_elements]*exp(mod_pressure - mu);
	   }
	   //Generate Q^T f^n
	   SIDE = 'L'; UPLO = 'U'; TRANS = 'T'; // i.e. Q is stored in upper diagonal, and we are calculating Q^T x.    
       dormtr_(&SIDE, &UPLO, &TRANS, &N, &NRHS, C, &N, TAU, lpd_optimized_pressure, &N, WORK, &LWORK, &INFO);
	   
	   // Get - lambda Q^T f^n.
	   dscal_(&N, &minus_lambda,lpd_optimized_pressure, &INC);
	   // Get x - lambda Q^T f^n
	   daxpy_(&N, &one, x, &INC, lpd_optimized_pressure, &INC);
	   // Solve Tr = x-lambda Q^T f^n
	   dptsv_( &N, &NRHS, D, E, lpd_optimized_pressure , &N, &INFO );
	   dcopy_(&N,D2,&INC,D,&INC); // D has been overwritten. Replace.
	   dcopy_(&N_minus_1,E2,&INC,E,&INC); // E has been overwritten. Replace.
	   // Finally, the new pressure is computed as Qr:
	   SIDE = 'L'; UPLO = 'U'; TRANS = 'N'; // Multiply from the left, matrix is stored in upper half, no transpose this time. 
	   dormtr_(&SIDE, &UPLO, &TRANS, &N, &NRHS, C, &N, TAU, lpd_optimized_pressure , &N, WORK, &LWORK, &INFO);  
	   if (INFO) {  std::cout << " Well, you failed at the final hurdle. Couldn't invert the orthogonal matrix Q.  Sorry. " << INFO << std::endl;      return;   }
	   // Now we compute the difference between new and old and start again.
	   diff =0.0;
	   for (int p = 0 ;p < 2*ci_number_of_elements; p++ ) {
		   diff += (lpd_optimized_pressure[p]-lpd_pressure_old[p])*(lpd_optimized_pressure[p]-lpd_pressure_old[p]);
		   if ((lpd_optimized_pressure[p]-lpd_pressure_old[p])*(lpd_optimized_pressure[p]-lpd_pressure_old[p])>1.e-8){
			   std::cout << p <<std::endl;
			   std::cout << lpd_optimized_pressure[p] <<std::endl;
			   std::cout << lpd_pressure_old[p] <<std::endl;
			   
		   }
		   
		   
	   }
	   dcopy_(&N,lpd_optimized_pressure,&INC,lpd_pressure_old,&INC); // Must update old one
	   std::cout << "Difference " << diff << std::endl;
	  // std::cin >> tmp;
   }
   */
   // Finished. Copy the result across. 
   dcopy_(&N,lpd_optimized_pressure,&INC,apd_pressures,&INC);


    
   
   
   delete[] cpd_desired_result_of_optimization;
   delete[] lpd_optimized_pressure;
   delete[] lpd_pressure_old;
   delete[] cpd_transfer_matrix;
   delete[] cpd_spatial_vector_x;
   delete[] cpd_spatial_vector_y;
   delete[] cpd_spatial_vector_z;
   delete[] cpd_actual_position_of_focus;
   delete[] C;
   delete[] A;
   
   delete[] x;
   delete[] WORK;
   
}




void Optimize_Transducer::Generate_Pressures_Using_Piston_Time_Reversal(double *apd_pressures)
{ // 7 optimization parameters must be available, specifying the centre, radius and orientation of the piston being targeted. Using a piston means that we can consider only unidirectional waves. Which is a nice restriction.
	/* So the plan is: 
	1) set up a piston.
	2) Calculate the field at each computational element of the main transducer 
	3) Output these if you want them (if you pass an 8th parameter and it is nonzero)
	4) Average the pressures over each transducer element. 
	5) Reverse the phase and set these as the correct pressures. I think we have to reverse the phase. Makes sense to me, which is never such a good sign.
	 
	cpd_optimization_parameters[0]  = centre_X
	cpd_optimization_parameters[1]  = centre_Y
	cpd_optimization_parameters[2]  = centre_Z
	cpd_optimization_parameters[3]  = radius
	cpd_optimization_parameters[4]  = normal_X
	cpd_optimization_parameters[5]  = normal_Y
	cpd_optimization_parameters[6]  = normal_Z
	< optionally > cpd_optimization_parameters[7]  = output option - nonzero implies output values at each computational element on the main transducer.
	
	*/  
	if (ci_number_of_optimization_parameters < 7) {
		std::cerr << "Too few parameters for time reversal. You only input " << ci_number_of_optimization_parameters << std::endl;
		exit(-1);
	}
	// NOTE that we use the same number of radials and angles in the fictitious piston as in the nonfictitious ones. 
	Circular_Piston_Element time_reverser(cpd_optimization_parameters[3] ,  ci_number_of_radial_bits, ci_number_of_angles, 1.0, 0.0, cpd_optimization_parameters[0],  cpd_optimization_parameters[1], cpd_optimization_parameters[2], cpd_optimization_parameters[4], cpd_optimization_parameters[5],cpd_optimization_parameters[6]);  // Use nice constructor
	
	double total_real, total_imag,tmp_real, tmp_imag, posX,posY,posZ;
	int num_rad, num_angle;	
	// Now run over all transducer elements, and run over each of their computational elements and calculate the pressures on these.
	for (int p = 0; p <ci_number_of_elements; p++){
		// reset the pressures. 
		total_real=0.0;
		total_imag=0.0;
		cte_array_of_elements[p]->Output_Element_Numbers(num_rad, num_angle); // Get the actual number of radial and angular bits, not just expected number.		  
		for (int q = 0; q < num_rad*num_angle; q++){
			cte_array_of_elements[p]->Output_Position_Of_Computational_Element(posX,posY,posZ,q); // Get the position of the element.
			time_reverser.Calculate_Pressure(posX,posY,posZ, cd_omega_div_c, tmp_real,tmp_imag); // Get the pressure at the position due to the time-reverser. 
			total_real += tmp_real; 
			total_imag += tmp_imag;
			if ((ci_number_of_optimization_parameters > 7)&&(cpd_optimization_parameters[7])) {// Then we output results. Send them to stderr so that redirection is easier.
				std::cerr << posX << ", "<< posY << ", "<< posZ << ", "<< tmp_real << ", "<< tmp_imag << std::endl; 
			} 
		}
		apd_pressures[p] = total_real; 
		apd_pressures[p+ci_number_of_elements] = -total_imag; 
	} 
}


