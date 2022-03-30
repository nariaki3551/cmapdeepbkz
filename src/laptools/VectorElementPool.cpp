/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* MIT License                                                                     */
/*                                                                                 */
/* Copyright (c) 2022 Nariaki Tateiwa <n-tateiwa@kyudai.jp>                        */
/*                                                                                 */
/* Permission is hereby granted, free of charge, to any person obtaining a copy    */
/* of this software and associated documentation files (the "Software"), to deal   */
/* in the Software without restriction, including without limitation the rights    */
/* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell       */
/* copies of the Software, and to permit persons to whom the Software is           */
/* furnished to do so, subject to the following conditions:                        */
/*                                                                                 */
/* The above copyright notice and this permission notice shall be included in all  */
/* copies or substantial portions of the Software.                                 */
/*                                                                                 */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR      */
/* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        */
/* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     */
/* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          */
/* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,   */
/* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE   */
/* SOFTWARE.                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    VectorElementPool.cpp
 * @brief   Base class for lattice vector container.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "VectorElementPool.h"

#include <algorithm>
#include <iterator>
#include <vector>
#include <float.h>
#include "VectorElement.h"
#include "eigen3/Eigen/src/Core/DenseCoeffsBase.h"

namespace LapTools
{



///
/// @brief test whether the vector is inseartable
/// @param[in] vectorElement lattice vector
/// @return false if it cannot be inserted else true
///
template<typename BasisFloat>
bool
VectorElementBasePool<BasisFloat>::simpledInsertableTest(
      VectorElementPoolType<BasisFloat> vectorElement
      )
{
   if( limitUsageOfPool < 0 || static_cast<int>(ascendingPool.size()) < limitUsageOfPool )
   {
      return true;
   }
   else if( vectorElement->squaredNorm() > (*std::prev(ascendingPool.end()))->squaredNorm() )
   {
      return false;
   }
   return true;
}


///
/// @brief insert vector
/// @param[in] vectorElement lattice vector
/// @return true if it is inserted else false
/// @note delete VectorElement pointer if it is not inserted into this pool
///
template<typename BasisFloat>
bool
VectorElementBasePool<BasisFloat>::insert(
      VectorElementPoolType<BasisFloat> vectorElement
      )
{
   assert( vectorElement->vector(0) >= 0 );
   bool inserted = true;
   auto findVec = ascendingPool.find(vectorElement);
   if( findVec != ascendingPool.end() )
   {
      // there exists a vector in the pool whose norm is equal to the new vector
      inserted = false;
   }
   else
   {
      // the new vector is not in the pool
      ascendingPool.insert(vectorElement);
      if( limitUsageOfPool > 0 && static_cast<int>(ascendingPool.size()) > limitUsageOfPool )
      {
         auto i = std::prev(ascendingPool.end());
         if( vectorElement == (*i) ) inserted = false;
         ascendingPool.erase(i);
      }
   }
   return inserted;
}


///
/// @brief reset size of vector element pool
/// @param[in] inLimit limit of number of vectors in pool
///
template<typename BasisFloat>
void VectorElementBasePool<BasisFloat>::resetLimit(
         int inLimit
         )
{
   if( static_cast<int>(size()) > inLimit )
   {
      auto iter = getIterBegin();
      for( int i = 0; i < inLimit; ++i )
      {
         ++iter;
      }
      auto endIter = getIterEnd();
      while( iter != endIter )
      {
         erase(iter);
      }
   }
   limitUsageOfPool = inLimit;
}


///
/// @brief get quantile sqnorm value
/// @param[in] GH estimated shortest non-zero vector's norm of lattice
///
template<typename BasisFloat>
VectorElementPoolData
VectorElementBasePool<BasisFloat>::getPoolData(
      double GH
      )
{
   VectorElementPoolData data;
   data.size = size();
   data.max_value = 0.0;
   data.min_value = DBL_MAX;
   data.nLessThanAFOnePointFive = 0;
   data.nLessThanAFOnePointFour = 0;
   data.nLessThanAFOnePointThree = 0;
   data.nLessThanAFOnePointTwo = 0;
   data.nLessThanAFOnePointOne = 0;
   data.nLessThanAFOnePointZeroSeven = 0;
   data.nLessThanAFOnePointZeroFive = 0;
   data.nLessThanAFOnePointZeroTwo = 0;
   data.nLessThanAFOnePointZeroOne = 0;

   std::vector<double> sqnorms(size(), 0);
   int i = -1;
   for( const auto &p : ascendingPool )
   {
      sqnorms[++i] = p->squaredNorm();
      if( sqnorms[i] > data.max_value){ data.max_value = sqnorms[i]; }
      if( sqnorms[i] < data.min_value){ data.min_value = sqnorms[i]; }
      if( sqnorms[i] < GH * GH * 1.5 * 1.5 )
         data.nLessThanAFOnePointFive++;
      if( sqnorms[i] < GH * GH * 1.4 * 1.4 )
         data.nLessThanAFOnePointFour++;
      if( sqnorms[i] < GH * GH * 1.3 * 1.3 )
         data.nLessThanAFOnePointThree++;
      if( sqnorms[i] < GH * GH * 1.2 * 1.2 )
         data.nLessThanAFOnePointTwo++;
      if( sqnorms[i] < GH * GH * 1.1 * 1.1 )
         data.nLessThanAFOnePointOne++;
      if( sqnorms[i] < GH * GH * 1.07 * 1.07 )
         data.nLessThanAFOnePointZeroSeven++;
      if( sqnorms[i] < GH * GH * 1.05 * 1.05 )
         data.nLessThanAFOnePointZeroFive++;
      if( sqnorms[i] < GH * GH * 1.02 * 1.02 )
         data.nLessThanAFOnePointZeroTwo++;
      if( sqnorms[i] < GH * GH * 1.01 * 1.01 )
         data.nLessThanAFOnePointZeroOne++;
   }
   auto Q1 = sqnorms.size() / 4;
   auto Q2 = sqnorms.size() / 2;
   auto Q3 = Q1 + Q2;
   std::nth_element(sqnorms.begin(),          sqnorms.begin() + Q1, sqnorms.end());
   std::nth_element(sqnorms.begin() + Q1 + 1, sqnorms.begin() + Q2, sqnorms.end());
   std::nth_element(sqnorms.begin() + Q2 + 1, sqnorms.begin() + Q3, sqnorms.end());
   data.Q1_value     = sqnorms[Q1];
   data.median_value = sqnorms[Q2];
   data.Q3_value     = sqnorms[Q3];
   return data;
}


///
/// @brief get quantile sqnorm value
/// @param[in] GH estimated shortest non-zero vector's norm of lattice
///
template<typename BasisFloat>
const std::string
VectorElementBasePool<BasisFloat>::toStatString(
      double GH,
      std::string delimiter
      )
{
   std::ostringstream s;
   auto poolData = getPoolData(GH);
   s << delimiter << poolData.size
     << delimiter << poolData.min_value
     << delimiter << poolData.Q1_value
     << delimiter << poolData.median_value
     << delimiter << poolData.Q3_value
     << delimiter << poolData.max_value
     << delimiter << poolData.nLessThanAFOnePointFive
     << delimiter << poolData.nLessThanAFOnePointFour
     << delimiter << poolData.nLessThanAFOnePointThree
     << delimiter << poolData.nLessThanAFOnePointTwo
     << delimiter << poolData.nLessThanAFOnePointOne
     << delimiter << poolData.nLessThanAFOnePointZeroSeven
     << delimiter << poolData.nLessThanAFOnePointZeroFive
     << delimiter << poolData.nLessThanAFOnePointZeroTwo
     << delimiter << poolData.nLessThanAFOnePointZeroOne;
   return s.str();
}



///
/// instantiation
///
template class VectorElementBasePool<int>;
template class VectorElementBasePool<long int>;


}  // namespace LapTools
