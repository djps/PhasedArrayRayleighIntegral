/***************************************************************************
 *   Copyright (C) 2009 by Simon Woodford                                  *
 *   simon.woodford@icr.ac.uk                                              *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

// NOTE This program has been significantly modified since its original incarnation. Instead of searching for local maxima on a full set of Cartesian data, it has been adapted to search on a filtered set of Cartesian data. So I still assume a nice grid, not an arbitrary unstructured mesh!

/*
This program analyzes 3-d data. Basically, you feed it a file name and 9 parameters. 
The file contains (at least) 4 columns, consisting of the x,y,z coordinates and the value at these coords. I assume that the data is on a full rectangular grid, i.e. no points have been dropped from the file. 
The first 3 parameters are the 3 sizes in each direction, i.e. dx,dy,dz
Next 1 parameter tells how many data points are on each line. I assume comma-separation. 
Anyway, the next 4 parameters tell you which positions in the line are of interest. So if we have 1,2,3,6 then x is column 1, y is col 2, z col 3 and value of interest is col 6. This is the typical input for my acoustic files (as of April 2009). Note that I assume that the file is written such that the x coordinate changes first, then the y, then the z. If this is not the case, simply change the parameters around to make everything reasonable.  
The 9th parameter tells how many local maxima we look for. 
OK, once the program has loaded the file, it searches for the largest local maxima. It then goes in each direction from these maxima  up to the point where the value has dropped below half of this, and then to 1/4 (useful if you want the place where pressure drops to half and all you have is pressure squared). It outputs these data in a pretty format. .
There are a lot of data to output. We output the position and value of each maximum and the square root of the maximum (useful for same reason as the 1/4 data are useful), the FWHM in each direction and the associated coordinates for each of these (useful for checking whether we are at the edge of the domain, also good for skewness estimates),  and the FWQM in each direction and associated coords.     

Compile is easy, with no additonal dependencies:
  icc -o GetDampedFWHMAtNLocalMaxima GetDampedFWHMAtNLocalMaxima.cpp
  chmod +x GetDampedFWHMAtNLocalMaxima_a

*/
/* modified by RST 2011-06-16:  This does simple attenuation correction in Z from the origin. Additional command line parameters specify: (i) atten_constant, the frequency-independent attenuation in Nepers/distance unit and (ii) atten_power, the power-of-frequency-in-MHz attenuation parameter. The frequency is read (in kHz) from the value following _Freq_ in the input data file name.  The attenuation is applied in the Z-axis only, from the origin.  The intensity attenuation value is 2 * atten_constant * freq_MHz ^ atten_pow. The output intensity I_out = I_in * exp(-z * atten_intensity). */
/* modified by RST 2011-06-17:  This now reports just the first-found point on a mazimum plateau. */
/* modified by RST 2012-07-02:  The attenuation Freq parameter may be taken from the command line, or the attenuaton parameters may be excluded. */  

#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstring>
#include <vector>
#include <iomanip>
#include <map>
#include <cstdlib>
#include <cmath>

  
  
using namespace std;

// declarations --------------------------------------------
int Find_Local_Maxima(int number_of_local_maxima, int* xpos, int* ypos, int* zpos,  double* max,  const map<pair<int,pair<int,int> >, double >  &data);
void FindLinearDimensionsOfHalfMax(const map<pair<int,pair<int,int> >, double >  &data, int &xmin, int &xmax,int &ymin, int &ymax,int &zmin, int &zmax, int startx, int starty, int startz, double max);
void FindLinearDimensionsOfQuarterMax(const map<pair<int,pair<int,int> >, double >  &data, int &xmin, int &xmax,int &ymin, int &ymax,int &zmin, int &zmax,  int startx, int starty, int startz, double max);
int Find_Damped_Local_Maxima(int number_of_local_maxima, int* xpos, int* ypos, int* zpos,  double* max,  const map<pair<int,pair<int,int> >, double >  &data, double dz, double absorption);


// function definitions -----------------------------------

int main(int argc, char *argv[])
{
  // First things first, declare and define the parameters.
	double dx,dy,dz;
  int nperline;

   if (argc < 11)
   {
      cout << "Wrong number of parameters" << endl;
      return -1;
   }
 
   ifstream input;
   input.open(argv[1]);
   dx = atof(argv[2]);
   dy = atof(argv[3]);
   dz = atof(argv[4]);
   nperline = atoi(argv[5]);
   int n1,n2,n3,n4;
   n1 = atoi(argv[6]);
   n2 = atoi(argv[7]);
   n3  = atoi(argv[8]);
   n4 = atoi(argv[9]);
  // cerr << n1 << ", "<< n2 << ", "<< n3 << ", "<< n4 << ", "<< nperline << endl;
   int number_of_local_maxima = atoi(argv[10]);
   
   double atten_constant = 0.0;
   double atten_power = 1.0;
   double freq_MHz = 1.7;
   if( argc > 11 ) {
   	 atten_constant = atof(argv[11]);
     atten_power = atof(argv[12]);
     if( argc > 13 ) {
       freq_MHz = atof(argv[13] );
     } else {    	
       char *pcTemp = strstr(argv[1], "_Freq_");
       if( !pcTemp ){
   	     cout << "Could not find the '_Freq_' marker in the input filename." << endl;
   	     return -2;
   	   }
       int freq_kHz;
       sscanf(pcTemp+6, "%d" , &freq_kHz);
       double freq_MHz = double( freq_kHz ) / 1000.0;
     }
   } 
   
   double atten_intensity = 2.0 * atten_constant * pow(freq_MHz, atten_power);
   
   double *temp = new double [nperline];    
   char Comma;
   double x,y,z;
   string DiscardMe;
   map<pair<int,pair<int,int> >, double >  data; // The 3 indices are x,y and z in integer form. The double is the actual data. 
   pair<int,pair<int,int> > tempPair;
   
   long nLinesLoaded = 0;
   while(input.good()){ 
   // Load the data!
   		for (int indperline = 0; indperline < (nperline-1); indperline++){
           input >> temp[indperline] >> Comma;
		//   cerr << temp[indperline] << ", ";
        }
		input >> temp[nperline-1];
	//	cerr << temp[nperline-1] << endl;
		getline(input,DiscardMe); // get rid of rest of line.
		
		x = temp[n1-1];
		y = temp[n2-1];
		z = temp[n3-1];
//		cerr << x << ", "<< y << ", "<< z << endl;
		x/= dx;
		y/= dy; 
		z/= dz;
		//cerr << x << ", "<< y << ", "<< z << endl;
		tempPair.first = floor(x+0.5); // A good way to round off almost-integer doubles. 
		tempPair.second.first = floor(y+0.5);
		tempPair.second.second = floor(z+0.5);
		data.insert(pair<pair<int,pair<int,int> >, double>(tempPair, temp[n4-1]));
		//cerr << data.size()<<endl;
		/*for (map<pair<int,pair<int,int> >, double >::iterator it = data.begin(); it!=data.end(); it++){
			cerr << it->first.first << ", "<< it->first.second.first << ", "<< it->first.second.second << ", "<< it->second << endl;
   }*/
		//cerr << temp[0] << ", " << temp[1] << ", " << temp[2] << ", " << temp[3] << ", " << temp[4] << ", " << temp[5] << endl;
		//if( 0 == (data.size()%100000) ){
		//	 cout << "Loaded " << data.size() << endl;
		//} 
   }   
   input.close();
   
   if( data.size() == 0 ){
   	 cout << "Oops! No data loaded from " << argv[1];
   	 exit(2);
   }    
   //cout << "Data loaded... " << endl;
   // Right, the data should be loaded now. Hooray!
   int xminhalf,yminhalf,zminhalf;
   int xmaxhalf,ymaxhalf,zmaxhalf;
   int xminquarter,yminquarter,zminquarter;
   int xmaxquarter,ymaxquarter,zmaxquarter;
   double *max = new double [number_of_local_maxima]; 
   int *xpos = new int [number_of_local_maxima] ;
   int *ypos = new int [number_of_local_maxima] ;
   int *zpos = new int [number_of_local_maxima] ;
   
   //cout << "Entering local minima finder. " << endl;
   int foundNumberOfMaxima = Find_Local_Maxima(number_of_local_maxima, xpos, ypos,zpos, max, data);
   if (foundNumberOfMaxima != number_of_local_maxima) {
	   cerr<< "Not enough maxima in data. You searched for " << number_of_local_maxima << " and there were only " << foundNumberOfMaxima << endl;
	   cerr<< "Continuing with reduced number of maxima." << endl;
	   number_of_local_maxima = foundNumberOfMaxima;
   }
   //cout << "done with that. " << endl; 
   //cout << max[0] << ", " << xpos[0] << ", " << ypos[0] << ", " << zpos[0] << endl ;
    
   for (int pp =  0; pp < number_of_local_maxima ; pp++)
   {   
      FindLinearDimensionsOfHalfMax(data,xminhalf,xmaxhalf,yminhalf,ymaxhalf,zminhalf,zmaxhalf,xpos[pp],ypos[pp],zpos[pp],max[pp]);
      FindLinearDimensionsOfQuarterMax(data,xminquarter,xmaxquarter,yminquarter,ymaxquarter,zminquarter,zmaxquarter,xpos[pp],ypos[pp],zpos[pp],max[pp]);
      cout << max[pp] << ", " << sqrt(max[pp]) << ", " << dx*xpos[pp] << ", " << dy*ypos[pp] << ", " << dz*zpos[pp] << ", " ;  // position and maximum.       
      cout << dx*xmaxhalf-dx*xminhalf << ", " << dy*ymaxhalf-dy*yminhalf << ", "  << dz*zmaxhalf-dz*zminhalf <<  ", "; // FWHMs. 
      cout <<  dx*xminhalf << ", " << dx*xmaxhalf   ;   // skewness.
      cout <<  ", " <<  dy*yminhalf << ", " << dy*ymaxhalf ;   // skewness
      cout <<  ", " <<  dz*zminhalf << ", " << dz*zmaxhalf ;// skewness
      cout    << ", "  << dx*xmaxquarter-dx*xminquarter << ", "  << dy*ymaxquarter-dy*yminquarter << ", "  << dz*zmaxquarter-dz*zminquarter ; // FWQMs
      cout <<  ", " <<  dx*xminquarter << ", " << dx*xmaxquarter  ;  // skewness 
      cout <<  ", " <<  dy*yminquarter << ", " << dy*ymaxquarter;   // skewness
      cout <<  ", " <<  dz*zminquarter << ", " << dz*zmaxquarter ;   // skewness
      cout << "; " ; // make place for the next maximum.
   }
   foundNumberOfMaxima = Find_Damped_Local_Maxima(number_of_local_maxima, xpos, ypos, zpos, max, data, dz, atten_intensity); //
   if (foundNumberOfMaxima != number_of_local_maxima) {
	   cerr<< "Not enough maxima in data. You searched for " << number_of_local_maxima << " and there were only " << foundNumberOfMaxima << endl;
	   cerr<< "Continuing with reduced number of maxima." << endl;
	   number_of_local_maxima = foundNumberOfMaxima;
   }
   for (int pp =  0; pp < number_of_local_maxima ; pp++){
	   cout << max[pp] << ", " << sqrt(max[pp]) << ", " << dx*xpos[pp] << ", " << dy*ypos[pp] << ", " << dz*zpos[pp] << "; " ;  // position and maximum.       
   } 
    
   // In principle, it would be a good idea to put endls at the end of each of these. But often we shall use this program to analyze multiple data files, spit the results into a file and then load that file in a spreadsheet. In such a case, it makes much more sense to have all the data for each file on a single line. So my default will be to spit out all this stuff on one line.
   
   cout << endl;  // finish by flushing and ending the line.

   /*  The following makes a nice pretty set of data. But that's not what I want - I want to save this data to a file. So output all of it on a line as CSV's. 
   cout << "Max value in " << n4 <<"th column of array was  " << max << endl;
   cout << "If you're interested, its square root is " << sqrt(max) << endl;
   cout << "The position of this maximum was at " << XXX[xpos] << ", " <<YYY[ypos] << ", " <<ZZZ[zpos]  << endl;
   cout << "Might be a good idea to compare that to where you told the focus to go. " << endl;
   cout << "**************************   3 dB data    ********************************** " << endl;
   cout << "X positions: Left " << XXX[xminhalf] << "  Right " << XXX[xmaxhalf] << " Distance  "  << XXX[xmaxhalf]-XXX[xminhalf] << endl;   
   cout << "Y positions: Left " << YYY[yminhalf] << "  Right " << YYY[ymaxhalf] << " Distance  "  << YYY[ymaxhalf]-YYY[yminhalf] << endl;   
   cout << "Z positions: Left " << ZZZ[zminhalf] << "  Right " << ZZZ[zmaxhalf] << " Distance  "  << ZZZ[zmaxhalf]-ZZZ[zminhalf] << endl;   
   cout << "**************************   6 dB data    ********************************** " << endl;
   cout << "X positions: Left " << XXX[xminquarter] << "  Right " << XXX[xmaxquarter] << " Distance  "  << XXX[xmaxquarter]-XXX[xminquarter] << endl;   
   cout << "Y positions: Left " << YYY[yminquarter] << "  Right " << YYY[ymaxquarter] << " Distance  "  << YYY[ymaxquarter]-YYY[yminquarter] << endl;   
   cout << "Z positions: Left " << ZZZ[zminquarter] << "  Right " << ZZZ[zmaxquarter] << " Distance  "  << ZZZ[zmaxquarter]-ZZZ[zminquarter] << endl;   
   cout << endl;
   */
   /*
   cout << max << ", " << sqrt(max) << ", " << XXX[xpos] << ", " << YYY[ypos] << ", " << ZZZ[zpos] << ", " ;
   cout <<  XXX[xminhalf] << ", " << XXX[xmaxhalf] << ", "  << XXX[xmaxhalf]-XXX[xminhalf] ;   
   cout <<  ", " <<  YYY[yminhalf] << ", " << YYY[ymaxhalf] << ", "  << YYY[ymaxhalf]-YYY[yminhalf] ;   
   cout <<  ", " <<  ZZZ[zminhalf] << ", " << ZZZ[zmaxhalf] << ", "  << ZZZ[zmaxhalf]-ZZZ[zminhalf] ;   
   cout <<  ", " <<  XXX[xminquarter] << ", " << XXX[xmaxquarter] << ", "  << XXX[xmaxquarter]-XXX[xminquarter] ;   
   cout <<  ", " <<  YYY[yminquarter] << ", " << YYY[ymaxquarter] << ", "  << YYY[ymaxquarter]-YYY[yminquarter] ;   
   cout <<  ", " <<  ZZZ[zminquarter] << ", " << ZZZ[zmaxquarter] << ", "  << ZZZ[zmaxquarter]-ZZZ[zminquarter] ;   
   cout << endl;
   */
   
   return 0;
  
} // main()

bool compare(pair<pair<int,pair<int,int> >, double> a, pair<pair<int,pair<int,int> >, double> b){ // The sort function sorts elements such that the left-most element yields a "true" when compared to all others. Thus when < is used, the list ends up in ascending order. I want descending order, so I use >:
	return a.second > b.second; 	
}

int Find_Local_Maxima(int number_of_local_maxima, int* xpos, int* ypos, int* zpos,  double* max,  const map<pair<int,pair<int,int> >, double >  &data){ // Find the desired number of local maxima. Return the number of local maxima found, which allows you to confirm that everything went well.
	vector<pair<pair<int,pair<int,int> >, double> > maxima; // Yes, this is a ridiculous construction. But I need to use a vector if I want to use sort(). Anyway, all local maxima are pushed to this vector.
	double neighbours[26]; // The 26 neighbours of any point.
	int count;
	int isMax ;
	pair<int,pair<int,int> > position; // Useful to have one of these.
	for (map<pair<int,pair<int,int> >, double >::const_iterator it = data.begin(); it != data.end(); it++){ // Test each point to see whether it is a local maximum.
		count = 0;
		position = it->first; // start at the position of the test point. 
		// Now get all neighbouring values.
		for (int dx = -1; dx <=1; dx++){
			position.first += dx;
			for (int dy = -1; dy <=1; dy++){
				position.second.first += dy;				
				for (int dz = -1; dz <=1; dz++){
					position.second.second+=dz;
					if (dx || dy || dz) neighbours[count++] = (data.find(position) != data.end()) ? (data.find(position))->second : (it->second)-1.0 ;
					position.second.second -= dz;
					// Got all the neighbours.
					
				}
				position.second.first -= dy;			
			} 
			position.first -=dx;
		}
		isMax = 1;
		for (int p = 0 ; p < 26; p++) if (neighbours[p] > it->second) isMax = 0; // if any neighbours are larger, this point is not a max.
		if (isMax) {// We got a local maximum!
                        // look for adjoining maximum
                        bool bJoins = false;
                        for( vector<pair<pair<int,pair<int,int> >, double> >::const_iterator pMaxima = maxima.begin(); pMaxima != maxima.end(); pMaxima++ ){
                            if( (abs(pMaxima->first.first - position.first) < 2) && (abs(pMaxima->first.second.first - position.second.first) < 2) && (abs(pMaxima->first.second.second - position.second.second) < 2) ) {
                              bJoins = true;
                              break;
                            }
                        }
                        if( !bJoins ) {
			  maxima.push_back(pair<pair<int,pair<int,int> >, double> (position,it->second));
                            //if (it->second > 1000.0) cerr << "Found a local max: " << lt <<", "<< rt <<", "<< tp <<", "<< bt <<", "<< ft <<", "<< bk <<", " << it->second << endl;
                            //if (it->second > 1000.0) cerr << "Hello Found a local max: " << lt <<", "<< rt <<", "<< tp <<", "<< bt <<", "<< ft <<", "<< bk <<", " << it->second << endl;
                        }
		}
	}
	sort(maxima.begin(), maxima.end()-1, compare);
	/*
	cerr << maxima[0].second << ", " << maxima[0].first.first <<", " << maxima[0].first.second.first <<", " << maxima[0].first.second.second << endl;
	cerr << maxima[1].second << ", " << maxima[1].first.first <<", " << maxima[1].first.second.first <<", " << maxima[1].first.second.second << endl;
	cerr << maxima[2].second << ", " << maxima[2].first.first <<", " << maxima[2].first.second.first <<", " << maxima[2].first.second.second << endl;
	cerr << maxima[3].second << ", " << maxima[3].first.first <<", " << maxima[3].first.second.first <<", " << maxima[3].first.second.second << endl;
	*/
	int returnedMaxima = number_of_local_maxima;
	if (maxima.size()<number_of_local_maxima) returnedMaxima=maxima.size(); // if there are insufficient maxima, return fewer than requested.
	// Get the memory
	/*if (xpos) delete[] xpos;
	if (ypos) delete[] ypos;
	if (zpos) delete[] zpos;
	if (max) delete[] max;	
	xpos = new int [returnedMaxima];
	ypos = new int [returnedMaxima];
	zpos = new int [returnedMaxima];
	max = new double [returnedMaxima];
	*/
	for (int p = 0; p < returnedMaxima; p++){
		xpos[p] = maxima[p].first.first;
		//cerr << xpos[p]<<endl;
		ypos[p] = maxima[p].first.second.first;
		zpos[p] = maxima[p].first.second.second;
		max[p] = maxima[p].second;		
	}
	return returnedMaxima;
}

void FindLinearDimensionsOfHalfMax(const map<pair<int,pair<int,int> >, double > &data, int &xmin, int &xmax,int &ymin, int &ymax,int &zmin, int &zmax, int posmax_x,int posmax_y,int posmax_z,  double max)
{ // OK, this does most of the legwork.
   // Basically, the aim of this function is to start at a given "max" and move away from it until the value in the array dips below half of this maximum. Of course, the "max" need not be a true maximum, so this will be useful for info re local maxima.
   double cutoff;
   cutoff = 0.5*max;
   double tmp;
   // initialize things
   xmin = posmax_x;
   ymin = posmax_y;
   zmin = posmax_z;
   xmax = posmax_x;
   ymax = posmax_y;
   zmax = posmax_z;
   tmp = max;
   map<pair<int,pair<int,int> >, double >::const_iterator  it,it2;
   pair<int,pair<int,int> > position;
   position.first =  posmax_x;
   position.second.first =  posmax_y;
   position.second.second =  posmax_z;
   
   it = data.find(position); // point to the maximum in question.
   it2=it; // this one will be varied. 
   while ((it2!=data.end())&&((tmp=it2->second)>cutoff)){ // step in x, starting from the maximum position and working downwards. Note that I check that it2 exists before assigning tmp to its value. 
	   xmin--; // decrement xmin.      
	   position.first = xmin; // change position.
	   it2 = data.find(position); // Find the value at the new position, or not if it doesn't exist.      
   }
   position.first = posmax_x; // reset position.
   it2 = data.find(position);
   while ((it2!=data.end())&&((tmp=it2->second)>cutoff)){ // step in y, starting from the maximum position and working downwards. Note that I check that it2 exists before assigning tmp to its value. 
	   ymin--; // decrement ymin.      
	   position.second.first = ymin; // change position.
	   it2 = data.find(position); // Find the value at the new position, or not if it doesn't exist.      
   }
   position.second.first = posmax_y; // reset position.
   it2 = data.find(position);
   while ((it2!=data.end())&&((tmp=it2->second)>cutoff)){ // step in z, starting from the maximum position and working downwards. Note that I check that it2 exists before assigning tmp to its value. 
	   zmin--; // decrement zmin.      
	   position.second.second = zmin; // change position.
	   it2 = data.find(position); // Find the value at the new position, or not if it doesn't exist.      
   }
   position.second.second = posmax_z; // reset position.
   it2 = data.find(position);
   while ((it2!=data.end())&&((tmp=it2->second)>cutoff)){ // step in x, starting from the maximum position and working downwards. Note that I check that it2 exists before assigning tmp to its value. 
	   xmax++; // decrement xmin.      
	   position.first = xmax; // change position.
	   it2 = data.find(position); // Find the value at the new position, or not if it doesn't exist.      
   }
   position.first = posmax_x; // reset position.
   it2 = data.find(position);
   while ((it2!=data.end())&&((tmp=it2->second)>cutoff)){ // step in y, starting from the maximum position and working downwards. Note that I check that it2 exists before assigning tmp to its value. 
	   ymax++; // decrement ymin.      
	   position.second.first = ymax; // change position.
	   it2 = data.find(position); // Find the value at the new position, or not if it doesn't exist.      
   }
   position.second.first = posmax_y; // reset position.
   it2 = data.find(position);
   while ((it2!=data.end())&&((tmp=it2->second)>cutoff)){ // step in z, starting from the maximum position and working downwards. Note that I check that it2 exists before assigning tmp to its value. 
	   zmax++; // decrement zmin.      
	   position.second.second = zmax; // change position.
	   it2 = data.find(position); // Find the value at the new position, or not if it doesn't exist.      
   }
   position.second.second = posmax_z; // reset position.   
}



void FindLinearDimensionsOfQuarterMax(const map<pair<int,pair<int,int> >, double > &data, int &xmin, int &xmax,int &ymin, int &ymax,int &zmin, int &zmax, int posmax_x,int posmax_y,int posmax_z,  double max)
{ // OK, this does most of the legwork.
   // Basically, the aim of this function is to start at a given "max" and move away from it until the value in the array dips below quarter of this maximum. Of course, the "max" need not be a true maximum, so this will be useful for info re local maxima.
	double cutoff;
	cutoff = 0.25*max;
	double tmp;
   // initialize things
	xmin = posmax_x;
	ymin = posmax_y;
	zmin = posmax_z;
	xmax = posmax_x;
	ymax = posmax_y;
	zmax = posmax_z;
	tmp = max;
	map<pair<int,pair<int,int> >, double >::const_iterator  it,it2;
	pair<int,pair<int,int> > position;
	position.first =  posmax_x;
	position.second.first =  posmax_y;
	position.second.second =  posmax_z;
   
	it = data.find(position); // point to the maximum in question.
	it2=it; // this one will be varied. 
	while ((it2!=data.end())&&((tmp=it2->second)>cutoff)){ // step in x, starting from the maximum position and working downwards. Note that I check that it2 exists before assigning tmp to its value. 
		xmin--; // decrement xmin.      
		position.first = xmin; // change position.
		it2 = data.find(position); // Find the value at the new position, or not if it doesn't exist.      
	}
	position.first = posmax_x; // reset position.
	it2 = data.find(position);
	while ((it2!=data.end())&&((tmp=it2->second)>cutoff)){ // step in y, starting from the maximum position and working downwards. Note that I check that it2 exists before assigning tmp to its value. 
		ymin--; // decrement ymin.      
		position.second.first = ymin; // change position.
		it2 = data.find(position); // Find the value at the new position, or not if it doesn't exist.      
	}
	position.second.first = posmax_y; // reset position.
	it2 = data.find(position);
	while ((it2!=data.end())&&((tmp=it2->second)>cutoff)){ // step in z, starting from the maximum position and working downwards. Note that I check that it2 exists before assigning tmp to its value. 
		zmin--; // decrement zmin.      
		position.second.second = zmin; // change position.
		it2 = data.find(position); // Find the value at the new position, or not if it doesn't exist.      
	}
	position.second.second = posmax_z; // reset position.
	it2 = data.find(position);
	while ((it2!=data.end())&&((tmp=it2->second)>cutoff)){ // step in x, starting from the maximum position and working downwards. Note that I check that it2 exists before assigning tmp to its value. 
		xmax++; // decrement xmin.      
		position.first = xmax; // change position.
		it2 = data.find(position); // Find the value at the new position, or not if it doesn't exist.      
	}
	position.first = posmax_x; // reset position.
	it2 = data.find(position);
	while ((it2!=data.end())&&((tmp=it2->second)>cutoff)){ // step in y, starting from the maximum position and working downwards. Note that I check that it2 exists before assigning tmp to its value. 
		ymax++; // decrement ymin.      
		position.second.first = ymax; // change position.
		it2 = data.find(position); // Find the value at the new position, or not if it doesn't exist.      
	}
	position.second.first = posmax_y; // reset position.
	it2 = data.find(position);
	while ((it2!=data.end())&&((tmp=it2->second)>cutoff)){ // step in z, starting from the maximum position and working downwards. Note that I check that it2 exists before assigning tmp to its value. 
		zmax++; // decrement zmin.      
		position.second.second = zmax; // change position.
		it2 = data.find(position); // Find the value at the new position, or not if it doesn't exist.      
	}
	position.second.second = posmax_z; // reset position.   
}


int Find_Damped_Local_Maxima(int number_of_local_maxima, int* xpos, int* ypos, int* zpos,  double* max,  const map<pair<int,pair<int,int> >, double >  &am_data, double dz, double attenuation){ // Applies exponential damping to the field, then finds the desired number of local maxima. Return the number of local maxima found, which allows you to confirm that everything went well.
	vector<pair<pair<int,pair<int,int> >, double> > maxima; // Yes, this is a ridiculous construction. But I need to use a vector if I want to use sort(). Anyway, all local maxima are pushed to this vector.
	double neighbours[26]; // The 26 neighbours of any point.
	map<pair<int,pair<int,int> >, double >  data; // This will hold the exponentially damped elements.
	pair<int,pair<int,int> > pos;  
	for (map<pair<int,pair<int,int> >, double >::const_iterator it = am_data.begin(); it != am_data.end(); it++){ // copy data across.
		pos=it->first;
		data.insert(pair<pair<int,pair<int,int> >, double> (pos, (it->second)*exp(-dz*attenuation*pos.second.second)));//Note that the geometric focus at z=0 is considered the point of zero damping.
	}
	int count;
	int isMax ;
	pair<int,pair<int,int> > position; // Useful to have one of these.
	for (map<pair<int,pair<int,int> >, double >::iterator it = data.begin(); it != data.end(); it++){ // Test each point to see whether it is a local maximum.
		count = 0;
		position = it->first; // start at the position of the test point. 
		// Now get all neighbouring values.
		for (int dx = -1; dx <=1; dx++){
			position.first += dx;
			for (int dy = -1; dy <=1; dy++){
				position.second.first += dy;
				for (int dz = -1; dz <=1; dz++){
					position.second.second+=dz;
					if (dx || dy || dz) neighbours[count++] = (data.find(position) != data.end()) ? (data.find(position))->second : (it->second)-1.0 ;
					position.second.second -= dz;
					// Got all the neighbours.

				}
				position.second.first -= dy;
			}
			position.first -=dx;
		}
		isMax = 1;
		for (int p = 0 ; p < 26; p++) if (neighbours[p] > it->second) isMax = 0; // if any neighbours are larger, this point is not a max.
		if (isMax) {// We got a local maximum!
                       // look for adjoining maximum
                        bool bJoins = false;
                        for( vector<pair<pair<int,pair<int,int> >, double> >::const_iterator pMaxima = maxima.begin(); pMaxima != maxima.end(); pMaxima++ ){
                            if( (abs(pMaxima->first.first - position.first) < 2) && (abs(pMaxima->first.second.first - position.second.first) < 2) && (abs(pMaxima->first.second.second - position.second.second) < 2) ) {
                              bJoins = true;
                              break;
                            }
                        }
                        if( !bJoins ) {
			  maxima.push_back(pair<pair<int,pair<int,int> >, double> (position,it->second));
                            //if (it->second > 1000.0) cerr << "Found a local max: " << lt <<", "<< rt <<", "<< tp <<", "<< bt <<", "<< ft <<", "<< bk <<", " << it->second << endl;
                            //if (it->second > 1000.0) cerr << "Hello Found a local max: " << lt <<", "<< rt <<", "<< tp <<", "<< bt <<", "<< ft <<", "<< bk <<", " << it->second << endl;
                        }
		}
	}
	sort(maxima.begin(), maxima.end()-1, compare);
	/*
	cerr << maxima[0].second << ", " << maxima[0].first.first <<", " << maxima[0].first.second.first <<", " << maxima[0].first.second.second << endl;
	cerr << maxima[1].second << ", " << maxima[1].first.first <<", " << maxima[1].first.second.first <<", " << maxima[1].first.second.second << endl;
	cerr << maxima[2].second << ", " << maxima[2].first.first <<", " << maxima[2].first.second.first <<", " << maxima[2].first.second.second << endl;
	cerr << maxima[3].second << ", " << maxima[3].first.first <<", " << maxima[3].first.second.first <<", " << maxima[3].first.second.second << endl;
	*/
	int returnedMaxima = number_of_local_maxima;
	if (maxima.size()<number_of_local_maxima) returnedMaxima=maxima.size(); // if there are insufficient maxima, return fewer than requested.
	// Get the memory
	/*if (xpos) delete[] xpos;
	if (ypos) delete[] ypos;
	if (zpos) delete[] zpos;
	if (max) delete[] max;	
	xpos = new int [returnedMaxima];
	ypos = new int [returnedMaxima];
	zpos = new int [returnedMaxima];
	max = new double [returnedMaxima];
	*/
	for (int p = 0; p < returnedMaxima; p++){
		xpos[p] = maxima[p].first.first;
		//cerr << xpos[p]<<endl;
		ypos[p] = maxima[p].first.second.first;
		zpos[p] = maxima[p].first.second.second;
		max[p] = maxima[p].second;		
	}
	return returnedMaxima;
}




void FindEqualNeighbours( const map<pair<int,pair<int,int> >, double >  &data, map<pair<int,pair<int,int> >, double >& neighbours)
{
  /* given the sparse 3D array 'data' and the first datum in 'neighbours' (me), find the neighbouring values in 'data' that have the same value as me.  Add them to neighbours */


   pair< pair<int,pair<int,int> >, double>   me       = *neighbours.begin();
         pair<int,pair<int,int> >            position = me.first;;

   for (int dx = -1; dx <=1; dx++){
     position.first += dx;
     for (int dy = -1; dy <=1; dy++){
       position.second.first += dy;
       for (int dz = -1; dz <=1; dz++){
	 position.second.second+=dz;

	 if( (dx || dy || dz) && (data.find(position) != data.end()) && ((data.find(position))->second == me.second) ){

           // this is an equal-value neighbour.

             // add it to the map of neighbours (if it is not there already)
           pair< pair<int,pair<int,int> >, double>   neighbour( position, me.second);
             // very ugly construction. See Stroustroup 3rd edition pg. 488. Make a pair of an iterator to the map and a bool
           pair< map<pair<int,pair<int,int> >, double >::const_iterator, bool > p = neighbours.insert(neighbour);
           if( p.second ) {
             // this is a new neighbour.  recurse to try to find more...
             FindEqualNeighbours( data, neighbours);
           }

         }  // if( new equal-value-neighbour )

         position.second.second -= dz;
       } // for dz
       position.second.first -= dy;
     } // for dy
     position.first -=dx;
   } // for dx

}

