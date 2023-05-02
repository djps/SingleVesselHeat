//----------------------------------------------------------------------------------

/**
 * \file    array3d.cpp
 * \author  David Sinden
 * \date    2016-09-21
 * \brief   Methods of Array3D class
 */

//----------------------------------------------------------------------------------
#ifndef ARRAY3D_CPP
#define ARRAY3D_CPP

#include "array3d.h"

// Note that you cannot use templates in the same way as using usual classes. All of this needs to be #included in other files, not just the header file. This is because we need to recompile all of this when we instantiate the templates. So use #include "array3d.cpp" not array3d.h

template <class T>
Array3D<T>::Array3D()
{
   nx = 1;
   ny = 1;
   nz = 1;
   number_of_elements = nx * ny * nz;
   stride_x = ny * nz;
   stride_y = nz;
   stride_z = 1;
   array = new T[number_of_elements];
   for (long int i = 0; i < number_of_elements; ++i){
      array[i] = 0.0;
   }

   if ((nx > ny) && (nx > nz)) nmax = nx;
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


template <class T>
Array3D<T>::Array3D(long int nx_in, long int ny_in, long int nz_in)
{
   nx = nx_in;
   ny = ny_in;
   nz = nz_in;
   number_of_elements = nx * ny * nz;
   stride_x = ny * nz;
   stride_y = nz;
   stride_z = 1;
   array = new T[number_of_elements];
   for (long int i = 0; i < number_of_elements; ++i){
      array[i] = 0.0;
   }

   if ((nx > ny) && (nx > nz)) nmax = nx;
   else if (ny>nz) nmax = ny;
   else nmax = nz;

   work = new T[nmax];
   work2 = new T[nmax];
}


template <class T>
void Array3D<T>::Resize_Array(long int nx_in, long int ny_in, long int nz_in)
{
   nx = nx_in;
   ny = ny_in;
   nz = nz_in;
   number_of_elements = nx * ny * nz;
   stride_x = ny * nz;
   stride_y = nz;
   stride_z = 1;
   delete[] array;
   array = new T [number_of_elements];
   for (long int i = 0; i < number_of_elements; ++i){
      array[i] = 0;
   }

   if ((nx > ny) && (nx > nz)) nmax = nx;
   else if (ny>nz) nmax = ny;
   else nmax = nz;

   work = new T[nmax];
   work2 = new T[nmax];
}


template <class T>
Array3D<T>::Array3D(long int nx_in, long int ny_in, long int nz_in, T *in)
{
   nx = nx_in;
   ny = ny_in;
   nz = nz_in;
   number_of_elements = nx * ny * nz;
   stride_x = ny *nz;
   stride_y = nz;
   stride_z = 1;
   array = new T [number_of_elements];

   if ((nx > ny) && (nx > nz)) nmax = nx;
   else if (ny > nz) nmax = ny;
   else nmax = nz;

   work = new T[nmax];
   work2 = new T[nmax];

   for (int i = 0; i < number_of_elements; ++i) array[i] = in[i];
}


template <class T>
Array3D<T>::Array3D(Array3D<T> const& CopyFromMe)
{
   nx = CopyFromMe.nx;
   ny = CopyFromMe.ny;
   nz = CopyFromMe.nz;
   number_of_elements = nx * ny * nz;
   stride_x = ny * nz;
   stride_y = nz;
   stride_z = 1;
   array = new T [number_of_elements];

   for (int i = 0; i < number_of_elements; ++i) array[i] = CopyFromMe.array[i];

   if ((nx > ny) && (nx > nz)){nmax = nx;}
   else if (ny>nz) nmax = ny;
   else nmax = nz;

   work = new T[nmax];
   work2 = new T[nmax];
}


template <class T>
void Array3D<T>::operator+=(Array3D<T> const& AddMeToYou)
{
   // Adds the contents of AddMeToYou.array to array. Check conformity first.
   if ((nx == AddMeToYou.nx) && (ny == AddMeToYou.ny) && (nz == AddMeToYou.nz))
   {
      for (long int i = 0; i < number_of_elements; ++i) array[i] += AddMeToYou.array[i];
   }
   else {std::cout << "Trouble adding the two vectors" << std::endl;}
}


template <class T>
void Array3D<T>::Find_Max_And_Position(int &posx, int &posy, int &posz, T &max_val) const
{
   // Absolutely standard search for maximum. Could make it more efficient in terms of memory usage, but I don't feel like it.
   T tmpmax = array[0];
   // std::cout << "Current max " << tmpmax << std::endl;

   for (int indz = 0; indz < nz; ++indz)
   {
      // As I said, we assume that x is incremented first, so it is on the internal loop.
      for (int indy = 0; indy < ny; ++indy)
      {
         for (int indx = 0; indx < nx; ++indx)
         {
            if (array[indx * stride_x + indy * stride_y + indz * stride_z] > tmpmax)
            {
               tmpmax = array[indx * stride_x + indy * stride_y + indz * stride_z];
               // std::cout << "Current max " << tmpmax << std::endl;
               posx = int(indx);
               posy = int(indy);
               posz = int(indz);
            }
         }
      }
   }
   max_val = tmpmax;
}


template <class T> void Array3D<T>::Find_Local_Maxima(int number_of_maxima, int *posx , int *posy, int *posz, T *values) const
{
   T grid_of_points[3][3][3];
   int indx, indy, indz, insertion_position;
   for (indx = 0; indx < number_of_maxima; indx++)
   {
      values[indx] = array[0];  // initialize everything to the upper-left of the array.
      posx[indx] = 0;
      posy[indx] = 0;
      posz[indx] = 0;
   }

   for (indx = 1; indx < (nx-1); indx++)
   {
      //std::cout << indx << ", " << std::endl;
      for (indy = 1; indy < (ny-1); indy++)
      {
         // Initialize the cube to be centered at indx,indy,z = 0. Note that we center at z = 0 since the first step in the next routine is to roll the thing forward. Of course, it makes no sense to center it at z = 0, since the -1 plane isn't defined. So we only bother with the z = 0 and z=1 planes, which go into the 1 and 2 of the grid, while the 0 of the grid remains empty.   .
         for (int indx2 = 0; indx2 < 3; indx2++)
         {
            for (int indy2 = 0; indy2 < 3; indy2++)
            {
               for (int indz2 = 0; indz2 < 2; indz2++)
               {
                  grid_of_points[indx2][indy2][indz2+1] = array[(indx-1+indx2)*stride_x + (indy-1+indy2)*stride_y + indz2*stride_z];
               }
            }
         }
         for (indz = 1; indz < (nz-1); indz++)
         {  // At each step of the loop, we update the grid, check for local maximumness, then if appropriate,  insert the value into the array. I've written a secondary routine that finds the position at which insertion should occur, assuming the array of values is already ordered.
            // Move the data cube.
            for (int indx2=0;indx2<3;indx2++)
            {
               for (int indy2=0;indy2<3;indy2++)
               {
                  grid_of_points[indx2][indy2][0] = grid_of_points[indx2][indy2][1]; // roll the grid forward.
                  grid_of_points[indx2][indy2][1] = grid_of_points[indx2][indy2][2];
                  grid_of_points[indx2][indy2][2] = array[(indx-1+indx2)*stride_x+(indy-1+indy2)*stride_y+(indz+1)*stride_z] ; // get the next plane of values.
               }
            }

            if ((grid_of_points[1][1][1]>grid_of_points[0][0][0]) &&(grid_of_points[1][1][1]>grid_of_points[0][0][1]) &&(grid_of_points[1][1][1]>grid_of_points[0][0][2]) &&(grid_of_points[1][1][1]>grid_of_points[0][1][0]) &&(grid_of_points[1][1][1]>grid_of_points[0][1][1]) &&(grid_of_points[1][1][1]>grid_of_points[0][1][2]) &&(grid_of_points[1][1][1]>grid_of_points[0][2][0]) &&(grid_of_points[1][1][1]>grid_of_points[0][2][1]) &&(grid_of_points[1][1][1]>grid_of_points[0][2][2]) &&(grid_of_points[1][1][1]>grid_of_points[1][0][0]) &&(grid_of_points[1][1][1]>grid_of_points[1][0][1]) &&(grid_of_points[1][1][1]>grid_of_points[1][0][2]) &&(grid_of_points[1][1][1]>grid_of_points[1][1][0])                                                     &&(grid_of_points[1][1][1]>grid_of_points[1][1][2]) &&(grid_of_points[1][1][1]>grid_of_points[1][2][0]) &&(grid_of_points[1][1][1]>grid_of_points[1][2][1]) &&(grid_of_points[1][1][1]>grid_of_points[1][2][2]) &&(grid_of_points[1][1][1]>grid_of_points[2][0][0]) &&(grid_of_points[1][1][1]>grid_of_points[2][0][1]) &&(grid_of_points[1][1][1]>grid_of_points[2][0][2]) &&(grid_of_points[1][1][1]>grid_of_points[2][1][0]) &&(grid_of_points[1][1][1]>grid_of_points[2][1][1]) &&(grid_of_points[1][1][1]>grid_of_points[2][1][2]) &&(grid_of_points[1][1][1]>grid_of_points[2][2][0]) &&(grid_of_points[1][1][1]>grid_of_points[2][2][1]) && (grid_of_points[1][1][1]>grid_of_points[2][2][2])
                 && (grid_of_points[1][1][1] > values[number_of_maxima-1]))
            { // In this case, we have a local maximum that is larger than the smallest of the existing list of local maxima. So it must be inserted.
             //  std::cout << " Overwriting previous local max with value at point  " << indx << ", " << indy << ", " << indz << std::endl ;
              // std::cout << "Old value " << values[0]  << ", new value " << grid_of_points[1][1][1] << std::endl ;

               insertion_position = FindInsertionPostion(grid_of_points[1][1][1],values,number_of_maxima); // in principle, this can return insertion_position = number_of_maxima, i.e. out of range. But because we only enter this loop if the value must be inserted, this problem is avoided.
               for (int pp = (number_of_maxima-1); pp>insertion_position; pp--)
               {
                  values[pp] = values[pp-1];  // shift everything to the right by one place.
                  posx[pp] = posx[pp-1];
                  posy[pp] = posy[pp-1];
                  posz[pp] = posz[pp-1];
               }
               values[insertion_position] = grid_of_points[1][1][1];
               posx[insertion_position] = indx;
               posy[insertion_position] = indy;
               posz[insertion_position] = indz;
               // Debug use only
       //        for (int pp = 0; pp < number_of_maxima ; pp++) std::cout << values[pp] << ",  " ;
     //          std::cout << std::endl;
            }
         }
      }
   }


}



template <class T>
int Array3D<T>::FindInsertionPostion(T insertme, T *values, int N) const
{
   int position, min, max;

   // left-hand side of range.
   min = 0;

   // right-hand side of range.
   max = N-1;

   // failure -  the inserted value is too small. By putting this here, we forbid the loop to exit unsuccesfully.
   if (values[N-1] > insertme) position = N;
   else {
      while (min <= max)
      {
         // The following uses position as a dummy variable.
         // look at the middle of the range. Note integer division truncates remainders.
         position = (max + min) / 2;
         // Move the left of the range.
         if (values[position] > insertme) min = position+1;
         // move the right.
         else max = position;
         // test if we're done.
         if (insertme >= values[min])
         {
            // then we're done, x is inserted at min.
            // this is the position for insertion.
            position = min;
            // end the while loop.
            max = min-1;
         }
      }
   }
   //std::cout << ".........DONE" << std::endl;

   /* Uncomment this if you want to avoid insertion at equality:
   while(values[position]>=insertme) {position--;}
   */
   return position;
}



template <class T>
void Array3D<T>::TriDiagMultiply(T *a, T* b, T* c, int Index, int N, int index1, int index2)
{
   #ifdef DEBUG_ARRAY
   for (int p = 0; p < N; ++p) std::cout << a[p] << " " << b[p] << " " << c[p] << " " << std::endl;
   int temp;
   std::cin >> temp;
   #endif

   long int offset, offset2;
   if (Index==1)
   {
      // In this case, multiply along the x index, for the given  combination  of the y and z indices.
      if (N != nx) std::cout << "Trouble. You've input the wrong length vector. " << std::endl;
      if (index1 >= ny) std::cout << "Trouble. You're trying to solve the equation for an invalid y index " << std::endl;
      if (index2 >= nz) std::cout << "Trouble. You're trying to solve the equation for an invalid z index " << std::endl;
      offset = index1 * stride_y + index2 * stride_z;
      offset2 = 0;
      work[0] = b[0] * array[offset + offset2] + c[0]*array[offset + offset2 + stride_x]; // The top left part of the matrix multiplication.
      for (int x = 1; x < (nx-1); ++x)
      {
         offset2 += stride_x;
         work[x] = a[x] * array[offset+offset2-stride_x] + b[x] * array[offset + offset2] + c[x] * array[offset + offset2 + stride_x]; // The main chunk of the  matrix multiplication.
      }
      work[nx-1] = a[nx-1] * array[offset + offset2] + b[nx-1] * array[offset + offset2 + stride_x]; // The last part of the matrix multiplication.
      // Now work contains the new vector. Insert it back into the array:
      offset2 = 0;
      for (int x = 0; x < nx; ++x)
      {
         array[offset + offset2] = work[x]; // The main chunk of the  matrix multiplication.
         offset2 += stride_x;
      }
   }

  else if (Index==2)
  { // This time, go along the y index.
      if (N != ny) std::cout << "Trouble. You've input the wrong length vector. " << std::endl;
      if (index1 >= nx) std::cout << "Trouble. You're trying to solve the equation for an invalid x index " << std::endl;
      if (index2 >= nz) std::cout << "Trouble. You're trying to solve the equation for an invalid z index " << std::endl;

      offset = index1*stride_x + index2*stride_z;
      offset2 = 0;
      work[0] = b[0]*array[offset+offset2] + c[0]*array[offset+offset2+stride_y]; // The top left part of the matrix multiplication.
      for (int y = 1; y < (ny-1); ++y)
      {
         offset2 += stride_y;
         work[y] = a[y]*array[offset+offset2-stride_y] + b[y]*array[offset+offset2] + c[y]*array[offset+offset2 +stride_y]; // The main chunk of the  matrix multiplication.
      }
      work[ny-1] = a[ny-1]*array[offset+offset2]+b[ny-1]*array[offset+offset2+stride_y]; // The last part of the matrix multiplication.
      // Now work contains the new vector. Insert it back into the array:
      offset2 = 0;
      for (int y = 0; y < ny; ++y)
      {
         array[offset+offset2] = work[y]; // The main chunk of the  matrix multiplication.
         offset2 += stride_y;
      }
  }

   else if (Index==3)
   {
      // Finally we go along the Z axis.
      if (N != nz) std::cout << "Trouble. You've input the wrong length vector. " << std::endl;
      if (index1 >= nx) std::cout << "Trouble. You're trying to solve the equation for an invalid x index " << std::endl;
      if (index2 >= ny) std::cout << "Trouble. You're trying to solve the equation for an invalid y index " << std::endl;
      offset = index1*stride_x + index2*stride_y;
      offset2 = 0;
      work[0] = b[0]*array[offset+offset2] + c[0]*array[offset+offset2+stride_z]; // The top left part of the matrix multiplication.
      for (int z = 1; z < (nz-1); ++z)
      {
         offset2 += stride_z;
         work[z] = a[z]*array[offset+offset2-stride_z] + b[z]*array[offset+offset2] + c[z]*array[offset+offset2 +stride_z]; // The main chunk of the  matrix multiplication.
      }
      work[nz-1] = a[nz-1]*array[offset+offset2] + b[nz-1]*array[offset+offset2+stride_z]; // The last part of the matrix multiplication.
      // Now work contains the new vector. Insert it back into the array:
      offset2 = 0;
      for (int z = 0; z < nz; ++z)
      {
         array[offset + offset2] = work[z]; // The main chunk of the  matrix multiplication.
         offset2 += stride_z;
      }
   }

   else
   {
      std::cout<< "Index not supported, Multiply routine" <<  Index << std::endl;
   }
}


/**
 * \brief Tridiagonal solver when there's no problem about stepping through memory - copied from Wikipedia
 * \param[in] a Major diagonal
 * \param[in] b lower diagonal
 * \param[in, out] c upper diagonal
 * \param[in, out] d vector
 * \param[in, out] x vector
 * \param[in] n size of system
 * \warning  Will modify c and d
 * \note Not a member of Array3D class, but called by Array3D::TriDiagSolve
 */
template <class T>
void TridiagonalSolve(const T *a, const T *b, T *c, T *d, T *x, const unsigned int n)
{
   unsigned int i;

#ifdef DEBUG_ARRAY
   std::cout << b[0] << "  If this is zero, there's trouble ";
#endif

   // Modify the coefficients.
   c[0] /= b[0];                           /* Division by zero risk. */
   d[0] /= b[0];                           /* Division by zero would imply a singular matrix. */
   for (i = 1; i < n; i++){
      T id = (b[i] - c[i-1] * a[i]);     /* Division by zero risk. */

#ifdef DEBUG_ARRAY
      std::cout << i << " " << id << "  If this is zero, there's trouble ";
#endif

      c[i] /= id;                             /* Last value calculated is redundant. */
      d[i] = (d[i] - d[i-1] * a[i]) / id;
   }

   /* Now back substitute. */
   x[n - 1] = d[n - 1];
   for (i = n - 2; i >= 0; i--){
      x[i] = d[i] - c[i] * x[i + 1];
   }

#ifdef DEBUG_ARRAY
   std::cout << std::endl;
   int temp;
   std::cin >> temp;
#endif

}


template <class T>
void Array3D<T>::TriDiagSolve(T *a, T* b, T* c, int Index, int N, int index1, int index2)
{
   long int offset, offset2;
   if (Index==1)
   { // In this case, solve along the x index, index1 refers to y, index2 to z.
      if (N!= nx) std::cout << "Trouble. You've input the wrong length vector. " << std::endl;
      if (index1 >= ny) std::cout << "Trouble. You're trying to solve the equation for an invalid y index " << std::endl;
      if (index2 >= nz) std::cout << "Trouble. You're trying to solve the equation for an invalid z index " << std::endl;

      offset = index1 * stride_y + index2 * stride_z;
      offset2 = 0;
      for (int x = 0; x < nx; ++x)
      {
         work[x] = array[offset+offset2];       // copy the 3-d array into the work vector.
         offset2 += stride_x;
      }

      // Now solve the system.
      TridiagonalSolve(a,  b,  c, work, work2, nx);
      // Solution is copied into work2.

      // Now copy the solution into the array where it belongs.
      offset2 = 0;
      for (int x = 0; x < nx; ++x)
      {
         array[offset + offset2] = work2[x]; // The main chunk of the  matrix multiplication.
         offset2 += stride_x;
      }
   }

   else if (Index==2)
   { // In this case, solve along the y index, index1 is x, index2 is z.
      if (N != ny) std::cout << "Trouble. You've input the wrong length vector. " << std::endl;
      if (index1 >= nx) std::cout << "Trouble. You're trying to solve the equation for an invalid x index " << std::endl;
      if (index2 >= nz) std::cout << "Trouble. You're trying to solve the equation for an invalid z index " << std::endl;

      offset = index1 * stride_x + index2 * stride_z;
      offset2 = 0;
      for (int y = 0; y < ny; ++y)
      {
         work[y] = array[offset + offset2];       // copy the 3-d array into the work vector.
         offset2 += stride_y;
      }
      // Now solve the system.
      TridiagonalSolve(a,  b,  c, work, work2, ny);
      // Solution is copied into work_2.

      // Now copy the solution into the array where it belongs.
      offset2 = 0;
      for (int y = 0; y < ny; ++y)
      {
         array[offset + offset2] = work2[y]; // The main chunk of the  matrix multiplication.
         offset2 += stride_y;
      }
   }

   else if (Index==3)
   { // In this case, solve along the z index. index1 is x, index2 is y.
      if (N != nz) std::cout << "Trouble. You've input the wrong length vector. " << std::endl;
      if (index1 >= nx) std::cout << "Trouble. You're trying to solve the equation for an invalid x index " << std::endl;
      if (index2 >= ny) std::cout << "Trouble. You're trying to solve the equation for an invalid y index " << std::endl;

      offset = index1 * stride_x + index2 * stride_y;
      offset2 = 0;
      for (int z = 0; z < nz; ++z)
      {
         work[z] = array[offset+offset2];       // copy the 3-d array into the work vector.
         offset2 += stride_z;
      }
      // Now solve the system.
      TridiagonalSolve(a,  b, c, work, work2, nz);

      offset2 = 0;
      for (int z = 0; z < nz; ++z)
      {
         array[offset+offset2] = work2[z]; // The main chunk of the  matrix multiplication.
         offset2 += stride_z;
      }
   }

   else{
      std::cout<< "Index not supported, Solver routine " << Index << std::endl;
   }
}



template <class T>
void Array3D<T>::Output_Array_Of_Values(int x, int y, int z, T *output, int option) const
{
   if (option == 1)
   {
      long int offset = y * stride_y + z * stride_z;
      int count = 0;
      for (int xx = 0; xx < nx; ++xx)
      {
         output[xx] = array[count + offset];
         count += stride_x;
      }
   }
   else if (option == 2)
   {
      long int offset = x * stride_x + z * stride_z;
      long int count = 0;
      for (int yy = 0; yy < ny; ++yy)
      {
         output[yy] = array[count + offset];
         count += stride_y;
      }
   }
   else if (option == 3)
   {
      long int offset = y * stride_y + x * stride_x;
      long int count = 0;
      for (int zz = 0; zz < nz; ++zz)
      {
         output[zz] = array[count + offset];
         count += stride_z;
      }
   }
}


#endif
