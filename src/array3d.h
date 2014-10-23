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
 
#ifndef ARRAY3D_H
#define ARRAY3D_H

//#define DEBUG_ARRAY

/*
@author Simon Woodford <simon.woodford@icr.ac.uk>
 This is an implementation of a 3-d array. The idea is that if we solve the heat equation in a 
 ADI FD method, we get tridiagonal matrices. However, if we force the 3-d array into a 1-d array, 
 then we no longer have true tridiagonal matrices because everything gets mixed up. SO: we need to 
 implement fancy tridiag multiplication schemes and solvers that know about this internal trouble. 
 Should be relatively easy, although a pain to parallelize. 
*/

#include <iostream>

//typedef double T; 
// Not sure how to use copy constructors with templates, so I'll just use a typedef for now.
template <class T>

class Array3D{

	public:

		Array3D();

		~Array3D();
		
		// Number of elements in each direction. 
		long int nx, ny, nz; 
		
		// stores the largest of nx, ny, nz.
		long int nmax; 
		
		// This is the huge array that holds everything.  
		// Note that the matrix will be stored as a long array in the order T_0,0,0, T_0,0,1, T_0,0,2 ... T_0,0,NZ-1, T_0,1,0...
		T* array; 
		
		// This stores the distance between the succesive x, y or z values. stride_z =1, stride_y = nz, stride_x = ny*nz. 
		// Of course, it's dumb to store them, but it's more dumb to calculate them every time you want to do a transposition...
		long int stride_x, stride_y, stride_z,number_of_elements; 
		
		// Useful to grab the memory for our work vectors only once, and hold onto it until the bitter end. 
		T *work; 
		T *work2; 

		// Default constructors, etc. are good idea.
		Array3D(long int nx_in, long int ny_in, long int nz_in);
		Array3D(long int nx_in, long int ny_in, long int nz_in, T *in);
		Array3D(Array3D<T> const &CopyFromMe);

		// Resizing arrays is sometimes needed:
		void Resize_Array(long int nx_in, long int ny_in, long int nz_in);

		// Useful addition member
		void operator += (Array3D<T> const &AddMeToYou);

		// Useful member function that outputs the value of a particular element:
		T Output_Value(int x, int y, int z) const { return array[x*stride_x+y*stride_y+z*stride_z]; };
		void Input_Value(int x, int y, int z, T value){ array[x*stride_x+y*stride_y+z*stride_z] = value; };
		void Output_Array_Of_Values(int x, int y, int z, T *output, int option) const;
		void Add_Value(int x, int y, int z, T value){ array[x*stride_x+y*stride_y+z*stride_z] += value; };
		void Subtract_Value(int x, int y, int z, T value){ array[x*stride_x+y*stride_y+z*stride_z] -= value; };

		// Some more useful functions for searching, etc.

		// find the position of the maximum value. Assumes that a less than operator has been defined on T. 
		void Find_Max_And_Position(int &posx, int &posy, int &posz, T &max_val) const; 

		// Finds the position and values of the largest "number_of_maxima" local maxima. Puts them into the given arrays in descending order. 
		void Find_Local_Maxima(int number_of_maxima, int *posx, int *posy, int *posz, T *values) const; 
		  
		int FindInsertionPostion(T insertme, T *values, int N) const;

		// Overload this for debugging. 
		void Input_Value(int x, int y, int z, T value, int option){ 
			if (option){
				std::cout << "Input  " << x*stride_x+y*stride_y+z*stride_z << ", " << number_of_elements << std::endl; 
				array[x*stride_x+y*stride_y+z*stride_z] = value;
			}
		};

		// Nontraditional matrix multiply in the sense that we can choose which index to multiply along. 
		// That is, we calculate M_ij T_klj or M_ij T_kjl or M_ij T_jkl, depending on the index chosen.  
		// Note that a is the vector below the main diagonal, b is the diagonal, c is the vector above the 
		// main diagonal. All three are assumed to have length N, so a[0] = 0 and c[N-1] = 0 by assumption. 
		// The result of this multiplication is written over array. Note that  Index is either 1 (x), 2 (y) or 3 (z) 
		// and N must correspond to the appropriate number of elements. Otherwise an error is returned. 
		// This allows easy testing - just make nx,ny,nz different. 
		void TriDiagMultiply(T *a, T *b, T *c, int Index, int N, int index1, int index2);

		// Now, solve a tridiagonal system. Again, a,b,c are the sub, main and superdiagonals, Index is 1, 2 or 3 
		// and N is the appropriate number of elements. Depending on index, we are either solving M_ij x_klj = d_kli, 
		// M_ij x_kjl = d_kil or M_ij x_jkl = d_ikl. At the call of this function, array contains the right-hand-side d. 
		// After the call, array contains the solution x. Some of the input vectors are overwritten. Sorry about this. 
		// Also note that we feed in two other indices. These are the values of the indices that do not take part in 
		// the matrix equation. 
		// THIS ROUTINE ONLY SOLVES THE MATRIX EQUATION FOR A SINGLE SET OF THESE OTHER INDICES!!! 
		// IT DOES NOT SOLVE IT FOR ALL INDICES AND IT WIPES OUT THE VALUES OF THE INPUT VECTORS. 
		// The reason for this is that the inhomogeneous heat equation does not break up nicely into independent 1-d problems. 
		// If it did, then each matrix inversion would be independent of the uninvolved indices. But it doesn't and the matrix 
		// inversion depends on these indices. 
		void TriDiagSolve(T *a, T *b, T *c, int Index, int N, int index1, int index2); 

		T CalculateDifference(Array3D<T> const &CompareMe)
		{ 
			// Straightforward routine to calculate the sum of the squared differences between the current 
			// array and the one that is passed in. Useful in testing things.
			T temp = 0;
			for (long int i = 0; i<number_of_elements; ++i) temp += (array[i] -CompareMe.array[i])*(array[i] -CompareMe.array[i]);
			return temp;
		}

		T SumOfElements(void)
		{ 
		// Sums up all elements in array. 
			T temp = 0;
			for (long int i = 0; i<number_of_elements; ++i) temp +=  array[i];
			return temp;
		}

		T SumOfElementsProduct(Array3D<T> const &MultiplyMe)
		{ 
			// Sums up all elements in array. 
			T temp = 0;
			for (long int i = 0; i<number_of_elements; ++i) temp +=  array[i]*MultiplyMe.array[i];
			return temp;
		}

};

#endif
