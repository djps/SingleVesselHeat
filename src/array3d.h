//----------------------------------------------------------------------------------

/**
 * \file    array3d.h
 * \author  David Sinden
 * \date    2016-09-21
 * \brief   The class Array3D and some inline functions
*/

//----------------------------------------------------------------------------------
#ifndef ARRAY3D_H
#define ARRAY3D_H

//#define DEBUG_ARRAY
#include <iostream>

//typedef double T; // Not sure how to use copy constructors with templates, so I'll just use a typedef for now.
template <class T>


/**
  \class   Array3D
  \brief   An implementation of a 3-d array.
  \details The idea is that if we solve the heat equation in a ADI FD method, we get tridiagonal matrices. However, if we force the 3-d array into a 1-d array, then we no longer have true tridiagonal matrices because everything gets mixed up. SO: we need to implement fancy tridiag multiplication schemes and solvers that know about this internal trouble. Should be relatively easy, although a pain to parallelize.
  \note The array will be stored as a long array in the order \f$ T_{0,0,0}, T_{0,0,1}, T_{0,0,2} \ldots T_{0,0,n_{z-1}}, T_{0,1,0} \ldots \f$
*/
class Array3D{

  public:

    /**
     * \brief default constructor: must assign array or it won't make sense to destroy it in the destructor
     */
    Array3D();


    /**
     * \brief Destructor
     */
    ~Array3D();


    //! Number of elements in \f$x\f$ direction
    long int nx;


    //! Number of elements in \f$y\f$ direction
    long int ny;


    //! Number of elements in \f$z\f$ direction
    long int nz;


    //! stores the largest of `nx`, `ny`, `nz`
    long int nmax;


    //! This is the huge array that holds everything
    T* array;


    //! distance in memory between succesive \f$x\f$ values, stride_x = ny * nz
    long int stride_x;


    //! distance in memory between succesive \f$y\f$ values, stride_y = nz
    long int stride_y;


    //! distance in memory between succesive \f$z\f$ values, stride_z = 1
    long int stride_z;


    //! total number of elements
    long int number_of_elements;


    //! Useful to grab the memory for our work vectors only once, and hold onto it until the bitter end - but why is it public
    T *work;


    //! Useful to grab the memory for our work vectors only once, and hold onto it until the bitter end
    T *work2;


    //! Default constructor with sizes
    Array3D(long int nx_in, long int ny_in, long int nz_in);


    //! Default constructos with sizes and data
    Array3D(long int nx_in, long int ny_in, long int nz_in, T *in);


    //! Default constructor copying from another array3d
    Array3D(Array3D<T> const& CopyFromMe);


    /**
     * \brief Resizing arrays is sometimes needed
     */
    void Resize_Array(long int nx_in, long int ny_in, long int nz_in);


    /**
     * \brief Useful addition operator
     */
    void operator += (Array3D<T> const& AddMeToYou);


    /**
     * \brief Useful member function that outputs the value of a particular element
     */
    T Output_Value(int x, int y, int z) const {
      return array[x * stride_x + y*stride_y + z*stride_z];
    };


    /**
     * \brief Useful member function that sets the value of a particular element:
     */
    void Input_Value(int x, int y, int z, T value){
      array[x * stride_x + y * stride_y + z * stride_z] = value;
    };


    /**
     * \brief I think this may make things faster: it sends the whole stream of numbers along a cerain direction into the output array, where direction is x if `option` == 1, y if `option` is 2, z if òption` is 3
     */
    void Output_Array_Of_Values(int x, int y, int z, T *output, int option) const;


    /**
     * \brief Useful member function that adds to the value of a particular element by a given ammount, `value`
     */
    void Add_Value(int x, int y, int z, T value){
      array[x * stride_x + y * stride_y + z * stride_z] += value;
    };


    /**
     * \brief Useful member function that substracts the value of a particular element by a given ammount, `value`
     */
    void Subtract_Value(int x, int y, int z, T value){
      array[x * stride_x + y * stride_y + z * stride_z] -= value;
    };


    /**
     * \brief Find the position of the maximum value, assumes that a less than operator has been defined on type `T`
     */
    void Find_Max_And_Position(int &posx, int &posy, int &posz, T &max_val) const;


    /**
     * \brief  Finds the position and values of the largest number_of_maxima local maxima then puts them into the given arrays in descending order
     * \details I'm not sure of the most efficient way to hunt for local maxima. My feeling is that the best way is to have a moving cube of 27 data points and run this cube over the whole array. At each step, check whether the central point is larger than all the others. An alternative is to have a moving set of nearest neighbours (in 3-d, this would be a set of 6 points plus one central point). However, this has obvious difficulties since it is easy to think of a function that is diagonal to the Cartesian grid and would present false local maxima. So rather go with a full grid of 27 unless it proves to be absurdly slow.
     * I'm not sure how to check the borders efficiently. So I think it's easiest if I simply assume that the borders definitely do not contain local maxima. We can ensure this by editing the array when it is loaded, but that's not the problem of this routine.
     */
    void Find_Local_Maxima(int number_of_maxima, int *posx , int *posy, int *posz, T *values) const;


    /**
     * \brief Finds the insertion point
     * \details Straightforward search to find where "insertme" method fits into the array of values which is in descending order.
     * Perhaps not too straightforward, though. Basically, it is a _binary insert_ where the range of possible positions is halved each time. Insertion will always occur at the left of the range, so we must ensure that anything to the left (or rather, with lower index) than the left of the range is definitely larger than the value we are inserting. The right-hand side of the range is moved to ensure that the iterations end, but is otherwise not really involved.
     * Note that equality implies insertion, i.e. things are replaced by a new element of equal value. This is simply easier to code. It is easily changed by adding a quick check at the end, but that will slow things down.
     */
    int FindInsertionPostion(T insertme, T *values, int N) const;


    /**
     * \brief Overloaded function used for debugging
     */
    void Input_Value(int x, int y, int z, T value, int option)
    {
      if (option) std::cout << "Input  " << x * stride_x + y * stride_y + z * stride_z << ", " << number_of_elements << std::endl;
      array[x * stride_x + y * stride_y + z * stride_z] = value;
    };


    /**
     * \brief    Nontraditional matrix multiply in the sense that we can choose which index to multiply along
     * \details  That is, we calculate \f$ M_{ij} T_{klj}\f$ or \f$M_{ij} T_{kjl}\f$ or \f$M_{ij} T_{jkl}\f$, depending on the index chosen.
     * Note that `a` is the vector below the main diagonal, `b` is the diagonal, `c` is the vector above the main diagonal. All three are assumed to have length `N`, so `a[0] = 0` and `c[N-1] = 0` by assumption. The result of this multiplication writes over class member Array3D::array
     * Note that `Index` is either
     *  -# x,
     *  -# y or
     *  -# z
     * and `N` must correspond to the appropriate number of elements, otherwise an error is returned. This allows easy testing - just make `nx`, `ny`, `nz` different.
     */
    void TriDiagMultiply(T *a, T* b, T* c, int Index, int N, int index1, int index2);


    /*
    * \brief   Solves a tridiagonal system which calls TridiagonalSolve
    * \details Again, a, b, c are the sub, main and superdiagonals, Index is 1, 2 or 3 and N is the appropriate number of elements.
    * Depending on index, we are either solving \f$ M_{ij} x_{klj} = d_{kli}\f$, \f$ M_{ij} x_{kjl} = d_{kil}\f$ or \f$ M_{ij} x_{jkl} = d_{ikl}\f$. At the call of this function, array contains the right-hand-side \f$ d\f$ . After the call, the member array contains the solution \f$ x\f$ .
    * Some of the input vectors are overwritten. Sorry about this.
    * Basically, we copy the desired right-hand-side vector into a work vector, then call the traditional tridiagonal solver. Bit of a clunky method, but it's almost certainly faster than anything else I can think of. Note that we need two work vectors, since the tridiagonal routine doesn't copy solution over the right-hand-side vector.
    * Also note that we feed in two other indices. These are the values of the indices that do not take part in the matrix equation.
    * \note THIS ROUTINE ONLY SOLVES THE MATRIX EQUATION FOR A SINGLE SET OF THESE OTHER INDICES. IT DOES NOT SOLVE IT FOR ALL INDICES AND IT WIPES OUT THE VALUES OF THE INPUT VECTORS. The reason for this is that the inhomogeneous heat equation does not break up nicely into independent 1-d problems. If it did, then each matrix inversion would be independent of the uninvolved indices. But it doesn't and the matrix inversion depends on these indices.
    */
    void TriDiagSolve(T *a, T* b, T* c, int Index, int N, int index1, int index2);


    /**
     * \brief Sum of squares difference between two Array3D objects
     */
    T CalculateDifference(Array3D<T> const& CompareMe)
    {
      T temp = 0;
      for (long int i = 0; i < number_of_elements; ++i) temp += (array[i] - CompareMe.array[i]) * (array[i] - CompareMe.array[i]);
      return temp;
    }


    /**
     * \brief Sums up all elements in array
     */
    T SumOfElements(void)
    {
      T temp = 0;
      for (long int i = 0; i < number_of_elements; ++i) temp +=  array[i] ;
      return temp;
    }


    /**
     * \brief Sums product with another array
     */
    T SumOfElementsProduct(Array3D<T> const& MultiplyMe)
    {
      T temp = 0;
      for (long int i = 0; i < number_of_elements; ++i) temp +=  array[i] * MultiplyMe.array[i];
      return temp;
    }

};

#endif
