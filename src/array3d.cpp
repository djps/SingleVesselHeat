/**
 @file array3d.cpp
 */

#ifndef ARRAY3D_CPP
#define ARRAY3D_CPP

#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>

#include "array3d.h"


// Note that you cannot use templates in the same way as using usual classes.
// All of this needs to be #included in other files, not just the header file.
// This is because we need to recompile all of this when we instantiate the
// templates. So use #include "array3d.cpp" not array3d.h

template <class T> Array3D<T>::Array3D()
{ // default. Must assign array or it won't make sense to destroy it in the destructor.
  nx=1;
  ny=1;
  nz=1;
  number_of_elements=nx*ny*nz;
  stride_x = ny*nz;
  stride_y = nz;
  stride_z = 1;
  array = new T[number_of_elements];
  memset(array,0,number_of_elements*sizeof(T));

  if ((nx>ny)&&(nx>nz)) nmax = nx;
  else if (ny>nz) nmax = ny;
  else nmax = nz;

  work = new T[nmax];
  work2 = new T[nmax];

}

template <class T> Array3D<T>::~Array3D()
{
  delete[] array;
  delete[] work;
  delete[] work2;
}

template <class T> Array3D<T>::Array3D(long int nx_in, long int ny_in, long int nz_in)
{
  nx=nx_in;
  ny=ny_in;
  nz=nz_in;
  number_of_elements=nx*ny*nz;
  stride_x = ny*nz;
  stride_y = nz;
  stride_z = 1;
  array = new T [number_of_elements];
  //for (long int i = 0; i< number_of_elements; ++i) array[i] =0;
  memset(array,0,number_of_elements*sizeof(T)); // Set to 0 efficiently.

  if ((nx>ny)&&(nx>nz)) nmax = nx;
  else if (ny>nz) nmax = ny;
  else nmax = nz;

  work = new T[nmax];
  work2 = new T[nmax];
}

template <class T> Array3D<T>::Array3D(long int nx_in, long int ny_in, long int nz_in, T *in)
{
  nx=nx_in;
  ny=ny_in;
  nz=nz_in;
  number_of_elements=nx*ny*nz;
  stride_x = ny*nz;
  stride_y = nz;
  stride_z = 1;
  array = new T [number_of_elements];

  if ((nx>ny)&&(nx>nz)) nmax = nx;
  else if (ny>nz) nmax = ny;
  else nmax = nz;

  work = new T[nmax];
  work2 = new T[nmax];

  //for (int i = 0; i< number_of_elements; ++i) array[i] = in[i];
  memcpy((void*) array,(void*) in,number_of_elements*sizeof(T));
}

template <class T> Array3D<T>::Array3D(Array3D<T> const&  CopyFromMe)
{
  nx = CopyFromMe.nx;
  ny = CopyFromMe.ny;
  nz = CopyFromMe.nz;
  number_of_elements=CopyFromMe.number_of_elements ;
  stride_x = CopyFromMe.stride_x;
  stride_y = CopyFromMe.stride_y;
  stride_z = CopyFromMe.stride_z;
  array = new T [number_of_elements];

  //for (int i = 0; i<number_of_elements; ++i) array[i] = CopyFromMe.array[i];
  memcpy((void*) array,(void*) CopyFromMe.array , number_of_elements*sizeof(T));

  if ((nx>ny)&&(nx>nz)) nmax = nx;
  else if (ny>nz) nmax = ny;
  else nmax = nz;

  work = new T[nmax];
  work2 = new T[nmax];

}

template <class T> Array3D<T>&  Array3D<T>::operator=(Array3D<T> const&  CopyFromMe)
{ // assignment operator. Slightly different to copy.
	nx = CopyFromMe.nx;
	ny = CopyFromMe.ny;
	nz = CopyFromMe.nz;
	number_of_elements=CopyFromMe.number_of_elements;
	stride_x = CopyFromMe.stride_x;
	stride_y = CopyFromMe.stride_y;
	stride_z = CopyFromMe.stride_z;
	delete array;
	array = new T [number_of_elements];
  //for (int i = 0; i<number_of_elements; ++i) array[i] = CopyFromMe.array[i];
	memcpy((void*) array,(void*) CopyFromMe.array , number_of_elements*sizeof(T));
	if ((nx>ny)&&(nx>nz)) nmax = nx;
	else if (ny>nz) nmax = ny;
	else nmax = nz;

	work = new T[nmax];
	work2 = new T[nmax];
	return *this; // Technically, you should return a reference. That way, you can string them together.
}



template <class T> void Array3D<T>::operator+=(Array3D<T> const& AddMeToYou)
{// Adds the contents of AddMeToYou.array to array. Check conformity first.
   if ((nx==AddMeToYou.nx)&&(ny==AddMeToYou.ny)&&(nz==AddMeToYou.nz))
   {
     for (long int i = 0; i<number_of_elements; ++i) array[i] += AddMeToYou.array[i];
   }
   else {std::cerr << "Trouble adding the two vectors" << std::endl;}
}


template <class T> void Array3D<T>::AddScalarTimesArray(T scalar, Array3D<T> const& AddMeToYou){
	// Also very useful. Adds scalar *AddMeToYou.array to array. Check conformity first.
	if ((nx==AddMeToYou.nx)&&(ny==AddMeToYou.ny)&&(nz==AddMeToYou.nz))
	{
		for (long int i = 0; i<number_of_elements; ++i) array[i] += scalar*AddMeToYou.array[i];
	}
	else {std::cerr << "Trouble adding the two vectors" << std::endl;}
}


template <class T> void Array3D<T>::TriDiagMultiply(T *a, T* b, T* c, int Index, int N, int index1, int index2)
{// This does a multiplication of the tridiag matrix [a,b,c] (a below main diag, b main, c above) along the index Index.
// Index is either 1, 2, or 3, corresponding to x, y, z. N must be the appropriate number of points.
// Note that it only multiplies for one set of index1 and index2, which correspond to the indices that aren't involved
// in the matrix multiplication.

#ifdef DEBUG_ARRAY
   // for (int p = 0;p<N;++p) std::cout << a[p] << " " << b[p] << " " << c[p] << " " << std::endl;
  //  int temp;
  //  std::cin >> temp ;
#endif

  long int offset, offset2;

  if(Index==1)
  { // In this case, multiply along the x index, for the given  combination  of the y and z indices.
    if (N!= nx) std::cout << "Trouble. You've input the wrong length vector. " << std::endl;
    if (index1 >= ny) std::cout << "Trouble. You're trying to solve the equation for an invalid y index " << std::endl;
    if (index2 >= nz) std::cout << "Trouble. You're trying to solve the equation for an invalid z index " << std::endl;
    offset = index1*stride_y+index2*stride_z;
    offset2 = 0;
    work[0] = b[0]*array[offset+offset2]+c[0]*array[offset+offset2+stride_x]; // The top left part of the matrix multiplication.
    for (int x=1; x<(nx-1); ++x)
    {
	  offset2 += stride_x;
	  work[x] = a[x]*array[offset+offset2-stride_x] + b[x]*array[offset+offset2] + c[x]*array[offset+offset2 +stride_x]; // The main chunk of the  matrix multiplication.
    }
    work[nx-1] = a[nx-1]*array[offset+offset2]+b[nx-1]*array[offset+offset2+stride_x]; // The last part of the matrix multiplication.
	// Now work contains the new vector. Insert it back into the array:
	offset2 = 0;
	for (int x=0; x<nx; ++x)
	{
	  array[offset+offset2] = work[x]; // The main chunk of the  matrix multiplication.
	  offset2 += stride_x;
	 }

  }

  else if(Index==2)
  { // This time, go along the y index.
    if (N!= ny) std::cout << "Trouble. You've input the wrong length vector. " << std::endl;
    if (index1 >= nx) std::cout << "Trouble. You're trying to solve the equation for an invalid x index " << std::endl;
    if (index2 >= nz) std::cout << "Trouble. You're trying to solve the equation for an invalid z index " << std::endl;

        offset = index1*stride_x+index2*stride_z;
        offset2 = 0;
        work[0] = b[0]*array[offset+offset2]+c[0]*array[offset+offset2+stride_y]; // The top left part of the matrix multiplication.
        for (int y=1; y<(ny-1); ++y)
        {
          offset2 += stride_y;
          work[y] = a[y]*array[offset+offset2-stride_y] + b[y]*array[offset+offset2] + c[y]*array[offset+offset2 +stride_y]; // The main chunk of the  matrix multiplication.
        }
        work[ny-1] = a[ny-1]*array[offset+offset2]+b[ny-1]*array[offset+offset2+stride_y]; // The last part of the matrix multiplication.
        // Now work contains the new vector. Insert it back into the array:
        offset2 = 0;
        for (int y=0; y<ny; ++y)
        {
          array[offset+offset2] = work[y]; // The main chunk of the  matrix multiplication.
          offset2 += stride_y;
        }

  }

  else if(Index==3)
  { // Finally we go along the Z axis.
    if (N!= nz) std::cout << "Trouble. You've input the wrong length vector. " << std::endl;
    if (index1 >= nx) std::cout << "Trouble. You're trying to solve the equation for an invalid x index " << std::endl;
    if (index2 >= ny) std::cout << "Trouble. You're trying to solve the equation for an invalid y index " << std::endl;

    offset = index1*stride_x+index2*stride_y;
        offset2 = 0;
        work[0] = b[0]*array[offset+offset2]+c[0]*array[offset+offset2+stride_z]; // The top left part of the matrix multiplication.
        for (int z=1; z<(nz-1); ++z)
        {
          offset2 += stride_z;
          work[z] = a[z]*array[offset+offset2-stride_z] + b[z]*array[offset+offset2] + c[z]*array[offset+offset2 +stride_z]; // The main chunk of the  matrix multiplication.
        }
        work[nz-1] = a[nz-1]*array[offset+offset2]+b[nz-1]*array[offset+offset2+stride_z]; // The last part of the matrix multiplication.
        // Now work contains the new vector. Insert it back into the array:
        offset2 = 0;
        for (int z=0; z<nz; ++z)
        {
          array[offset+offset2] = work[z]; // The main chunk of the  matrix multiplication.
          offset2 += stride_z;
        }

  }
  else{std::cout<< "Index not supported, Multiply routine" <<  Index << std::endl;}
}

/* Tridiagonal solver when there's no problem about stepping through memory. Copied from Wikipedia.
Fills solution into x. Warning: will modify c and d! */
template <class T> void TridiagonalSolve(const T *a, const T *b, T *c, T *d, T *x, unsigned int n)
{
        int i;
#ifdef DEBUG_ARRAY
        if (b[0] == 0) std::cout << b[0] << "  If this is zero, there's trouble " ;
#endif
        /* Modify the coefficients. */
        c[0] /= b[0];                           /* Division by zero risk. */
        d[0] /= b[0];                           /* Division by zero would imply a singular matrix. */
        for(i = 1; i < n; i++){
                T id = (b[i] - c[i-1] * a[i]);     /* Division by zero risk. */
#ifdef DEBUG_ARRAY
               if (id == 0) std::cout << i << " " << id << "  If this is zero, there's trouble " ;
#endif
                c[i] /= id;                             /* Last value calculated is redundant. */
                d[i] = (d[i] - d[i-1] * a[i])/id;
        }

        /* Now back substitute. */
        x[n - 1] = d[n - 1];
        for(i = n - 2; i >= 0; i--)
                x[i] = d[i] - c[i] * x[i + 1];

#ifdef DEBUG_ARRAY
        /*std::cout << std::endl;
        int temp;
        std::cin >> temp ;
		*/
#endif

}


template <class T> void Array3D<T>::TriDiagSolve(T *a, T* b, T* c, int Index, int N, int index1, int index2)
{
// Basically, we copy the desired right-hand-side vector into a work vector, then call the traditional tridiagonal solver.
// Bit of a clunky method, but it's almost certainly faster than anything else I can think of. Note that we need two work
// vectors, since the tridiagonal routine doesn't copy solution over the right-hand-side vector.
  long int offset, offset2;
  if (Index == 1)
  {
	// In this case, solve along the x index, index1 refers to y, index2 to z.
    if (N != nx) std::cout << "Trouble. You've input the wrong length vector. " << std::endl;
    if (index1 >= ny) std::cout << "Trouble. You're trying to solve the equation for an invalid y index " << std::endl;
    if (index2 >= nz) std::cout << "Trouble. You're trying to solve the equation for an invalid z index " << std::endl;

    offset = index1*stride_y+index2*stride_z;
    offset2 = 0;
    for (int x=0; x<nx; ++x)
    {
		work[x] = array[offset+offset2];       // copy the 3-d array into the work vector.
		offset2 += stride_x;
    }
    // Now solve the system.
    TridiagonalSolve( a,  b,  c, work, work2, nx);
    // Solution is copied into work_2_x.

	// Now copy the solution into the array where it belongs.
    offset2 = 0;
    for (int x=0; x<nx; ++x)
    {
		array[offset+offset2] = work2[x]; // The main chunk of the  matrix multiplication.
		offset2 += stride_x;
    }
  }

  else if(Index==2)
  {
	// In this case, solve along the y index, index1 is x, index2 is z.
	  if (N!= ny) std::cout << "Trouble. You've input the wrong length vector. " << std::endl;
	  if (index1 >= nx) std::cout << "Trouble. You're trying to solve the equation for an invalid x index " << std::endl;
	  if (index2 >= nz) std::cout << "Trouble. You're trying to solve the equation for an invalid z index " << std::endl;

	  offset = index1*stride_x+index2*stride_z;
	  offset2 = 0;
	  for (int y=0; y<ny; ++y)
	  {
		  work[y] = array[offset+offset2];       // copy the 3-d array into the work vector.
		  offset2 += stride_y;
	  }
     // Now solve the system.
	  TridiagonalSolve(a,  b,  c, work, work2, ny);
     // Solution is copied into work_2.

 // Now copy the solution into the array where it belongs.
	  offset2 = 0;
	  for (int y=0; y<ny; ++y)
	  {
		  array[offset+offset2] = work2[y]; // The main chunk of the  matrix multiplication.
		  offset2 += stride_y;
	  }
  }

  else if(Index==3)
  { // In this case, solve along the z index. index1 is x, index2 is y.
	  if (N!= nz) std::cout << "Trouble. You've input the wrong length vector. " << std::endl;
	  if (index1 >= nx) std::cout << "Trouble. You're trying to solve the equation for an invalid x index " << std::endl;
	  if (index2 >= ny) std::cout << "Trouble. You're trying to solve the equation for an invalid y index " << std::endl;


	  offset = index1*stride_x+index2*stride_y;
	  offset2 = 0;
	  for (int z=0; z<nz; ++z)
	  {
		  work[z] = array[offset+offset2];  // copy the 3-d array into the work vector.
		  offset2 += stride_z;
	  }
     // Now solve the system.
	  TridiagonalSolve(a,  b, c, work, work2, nz);
     // Solution is copied into work_2.

 // Now copy the solution into the array where it belongs.
	  offset2 = 0;
	  for (int z=0; z<nz; ++z)
	  {
		  array[offset+offset2] = work2[z]; // The main chunk of the  matrix multiplication.
		  offset2 += stride_z;
	  }
  }
  else{std::cout<< "Index not supported, Solver routine " <<  Index << std::endl;}
}



template <class T> void Array3D<T>::Output_Array_Of_Values(int x, int y, int z, T *output, int option) const
{
// I think this may make things faster. It sends the whole stream of numbers along a cerain direction into the
// output array. direction is x if option ==1, y if opt =2, z if 3.
    if (option==1)
    {
        long int offset = y*stride_y + z*stride_z;
        int count = 0;
        for (int xx = 0; xx<nx; ++xx)
        {
            output[xx] = array[count + offset];
            count += stride_x;
        }
    }
    else if (option==2)
    {
        long int offset = x*stride_x + z*stride_z;
        long int count = 0;
        for (int yy = 0; yy<ny; ++yy)
        {
            output[yy] = array[count + offset];
            count += stride_y;
        }
    }
    else if (option==3)
    {
        long int offset = y*stride_y + x*stride_x;
        long int count = 0;
        for (int zz = 0; zz<nz; ++zz)
        {
            output[zz] = array[count + offset];
            count += stride_z;
        }
    }

}

template <class T> void Array3D<T>::Input_Array_Of_Values(int x, int y, int z, T *input, int option)
{
// I think this may make things faster. It sends the whole stream of numbers from the input array into the
// main 3D array. direction is x if option ==1, y if opt =2, z if 3.
	if (option==1)
	{
		long int offset = y*stride_y + z*stride_z;
		int count = 0;
		for (int xx = 0; xx<nx; ++xx)
		{
			array[count + offset] = input[xx];
			count += stride_x;
		}
	}
	else if (option==2)
	{
		long int offset = x*stride_x + z*stride_z;
		long int count = 0;
		for (int yy = 0; yy<ny; ++yy)
		{
			array[count + offset] = input[yy] ;
			count += stride_y;
		}
	}
	else if (option==3)
	{
		long int offset = y*stride_y + x*stride_x;
		long int count = 0;
		for (int zz = 0; zz<nz; ++zz)
		{
			array[count + offset] = input[zz] ;
			count += stride_z;
		}
	}

}


/**
 @brief Resizes the current array.
 @param [in] nx_in x-value
 @param [in] ny_in y-value
 @param [in] nz_in z-value
*/
template <class T> void Array3D<T>::resize(long int nx_in, long int ny_in, long int nz_in)
{
    nx = nx_in;
    ny = ny_in;
    nz = nz_in;
    stride_z = 1L;
    stride_y = nz;
    stride_x = ny*nz;
    number_of_elements = nx*ny*nz;
    nmax = (nx>ny) ? nx:ny;
    nmax = (nmax>nz)?nmax:nz; // Get maximum.
    delete[] work;
    delete[] work2;
    delete[] array;
    array = new T[number_of_elements];
    work = new T[nmax];
    work2 = new T[nmax];
}


#endif
