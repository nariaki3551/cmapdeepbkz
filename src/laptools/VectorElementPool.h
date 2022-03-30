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

/**@file    VectorElementPool.h
 * @brief   Base class for lattice vector container.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __LAPTOOLS_VECTOR_ELEMENT_POOL_H__
#define __LAPTOOLS_VECTOR_ELEMENT_POOL_H__

#include <ostream>
#include <new>
#include <assert.h>
#include <stddef.h>
#include <set>
#include <memory>
#include "VectorElement.h"


namespace LapTools
{


template<typename BasisFloat=int>
using VectorElementPoolType = std::shared_ptr<VectorElementBase<BasisFloat>>;


template<typename BasisFloat=int>
struct VectorElementBaseSortCriterion
{
   bool operator(
         )
      (
         const VectorElementPoolType<BasisFloat> &n1,
         const VectorElementPoolType<BasisFloat> &n2
         ) const
   {
      // 1. compare using only squared norm
      // return ( n1->squaredNorm < n2->squaredNorm );

      // 2. full compare
      if( n1->squaredNorm() < n2->squaredNorm() ) return true;
      if( n1->squaredNorm() > n2->squaredNorm() ) return false;
      auto a = n1->getVector();
      auto b = n2->getVector();
      int a_size = static_cast<int>(a.size());
      for( int i = 0; i < a_size; ++i )
      {
         if( a[i] < b[i] ) return true;
         if( a[i] > b[i] ) return false;
      }
      return false;
   }
};


typedef struct VectorElementPoolData_
{
    size_t size;                             ///< size of pool
    double min_value;                        ///< minimum value
    double Q1_value;                         ///< first quantile value
    double median_value;                     ///< median value
    double Q3_value;                         ///< third quantile value
    double max_value;                        ///< maximum value
    long int nLessThanAFOnePointFive;        ///< number of vectors whose norm < 1.5
    long int nLessThanAFOnePointFour;        ///< number of vectors whose norm < 1.4
    long int nLessThanAFOnePointThree;       ///< number of vectors whose norm < 1.3
    long int nLessThanAFOnePointTwo;         ///< number of vectors whose norm < 1.2
    long int nLessThanAFOnePointOne;         ///< number of vectors whose norm < 1.1
    long int nLessThanAFOnePointZeroSeven;   ///< number of vectors whose norm < 1.07
    long int nLessThanAFOnePointZeroFive;    ///< number of vectors whose norm < 1.05
    long int nLessThanAFOnePointZeroTwo;     ///< number of vectors whose norm < 1.02
    long int nLessThanAFOnePointZeroOne;     ///< number of vectors whose norm < 1.01
} VectorElementPoolData;


template<typename BasisFloat=int>
class VectorElementBasePool
{

using AscendingPool = std::set<
   VectorElementPoolType<BasisFloat>,
   VectorElementBaseSortCriterion<BasisFloat>
   >;

protected:

   int limitUsageOfPool;         ///< limit usage of this pool, if negative, this is ignore
   AscendingPool ascendingPool;  ///<


public:

   ///
   /// @brief default constructor
   ///
   VectorElementBasePool(){ limitUsageOfPool = -1; }

   ///
   /// @brief constructor
   /// @param[in] limit max element size of this pool
   ///
   VectorElementBasePool(int limit){ limitUsageOfPool = limit; }

   ///
   /// @brief deconstructor
   ///
   ~VectorElementBasePool()
   {
   }


   ///
   /// @brief test whether the vector is inseartable
   /// @param[in] vectorElement lattice vector
   /// @return false if it cannot be inserted else true
   ///
   bool simpledInsertableTest(
         VectorElementPoolType<BasisFloat> vectorElement
         );


   ///
   /// @brief insert vector
   /// @param[in] vectorElement lattice vector
   /// @return true if it is inserted else false
   /// @note delete VectorElement pointer if it is not inserted into this pool
   ///
   bool insert(
         VectorElementPoolType<BasisFloat> vectorElement
         );


   ///
   /// @return true if this pool is empty else false
   ///
   bool isEmpty(
         )
   {
      return ( ascendingPool.size() == 0 );
   }


   ///
   /// @return size of pool
   ///
   size_t size(
         )
   {
      return ascendingPool.size();
   }


   ///
   /// @brief extract VectorElement
   /// @warning this function should be called when thie pool is not empty
   ///
   VectorElementPoolType<BasisFloat> extract(
         )
   {
      assert( ascendingPool.size() > 0 );
      auto ve = *ascendingPool.begin();
      ascendingPool.erase(ascendingPool.begin());
      return ve;
   }


   ///
   /// @brief reset size of vector element pool
   /// @param[in] inLimit limit of number of vectors in pool
   ///
   void resetLimit(
         int inLimit
         );


   ///
   /// get begin of iterator
   ///
   typename AscendingPool::iterator getIterBegin(
         )
   {
      return ascendingPool.begin();
   }


   ///
   /// get end of iterator
   ///
   typename AscendingPool::iterator getIterEnd(
         )
   {
      return ascendingPool.end();
   }


   ///
   /// erase vectorElement from ascendingPool
   ///
   void erase(
         typename AscendingPool::iterator &p
         )
   {
      ascendingPool.erase(p++);
   }


   ///
   /// @brief get quantile sqnorm value
   /// @param[in] GH estimated shortest non-zero vector's norm of lattice
   ///
   VectorElementPoolData getPoolData(
         double GH
         );


   ///
   /// output statistics infomation of vector pool
   /// @param[in] GH Gaussian Heuristic of lattice
   /// @param[in] delimiter
   ///
   const std::string toStatString(
         double GH,
         std::string delimiter=","
         );


   ///
   /// output vectorElements in ascendingPool
   ///
   const std::string toString(
         )
   {
      std::ostringstream s;
      s << std::endl << std::endl;
      s << "===== Vector Pool (Size " << ascendingPool.size() << ") =====" << std::endl;
      int index = 0;
      for( const auto &p : ascendingPool )
      {
         s << "index " << index << std::endl;
         s << p->toSimpleString();
         s << std::endl;
         index++;
      }
      s << "===========================" << std::endl;
      return s.str();
   }

};


}  // namespace LapTools


#endif // __LAPTOOLS_VECTOR_ELEMENT_POOL_H__
