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

/**@file    VectorElement.h
 * @brief   Base class for lattice vector.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __LAPTOOLS_VECTOR_ELEMENT_H
#define __LAPTOOLS_VECTOR_ELEMENT_H

#include <iostream>

#include "Def.h"


namespace LapTools
{


template<typename BasisFloat=int>
class VectorElementBase
{


public:

   LatticeVector<BasisFloat> vector;  ///< vector (first coefficients is not negative)
   double __squaredNorm;              ///< squared norm of vector

   ///
   /// @brief constrctor
   ///
   VectorElementBase(){}
   VectorElementBase(
         LatticeVector<BasisFloat> &inVector
         )
   {
      init(inVector);
   }
   VectorElementBase(VectorElementBase<BasisFloat> *ve)
   {
      vector = LatticeVector<BasisFloat>(ve->vector);
      setSquaredNorm();
   }
   void init(LatticeVector<BasisFloat> &inVector)
   {
      vector = LatticeVector<BasisFloat>(inVector);
      if( vector.coeff(0) < 0 ) vector = -vector;
      setSquaredNorm();
   }


   ///
   /// @brief deconstructor
   ///
   virtual ~VectorElementBase(){}


   ///
   /// @return size of vector
   ///
   int size(
         )
   {
      return vector.size();
   }


   ///
   /// @return norm of vector
   ///
   double norm(
         )
   {
      return std::sqrt(__squaredNorm);
   }


   ///
   /// setters
   ///
   void setSquaredNorm(
         )
   {
      __squaredNorm = 0.0;
      for( int i = vector.size() - 1; i > -1; --i )
         __squaredNorm += vector.coeff(i) * vector.coeff(i);
   }
   void setSquaredNorm(double inSquaredNorm){ __squaredNorm = inSquaredNorm; }
   void setVector(LatticeVector<BasisFloat> &inVector){ vector = inVector; }


   ///
   /// @breif getter of vector
   ///
   LatticeVector<BasisFloat> &
   getVector(
         )
   {
      return vector;
   }

   ///
   /// @breif getter of squaredNorm
   ///
   double &
   squaredNorm(
         )
   {
      return __squaredNorm;
   }


   ///
   /// @brief Stringfy Vector Element
   ///
   virtual std::string toString(
         )
   {
      std::ostringstream s;
      s << "vector: "      << vector.transpose() << std::endl;
      s << "squaredNorm: " << __squaredNorm        << std::endl;
      s << "norm: "        << std::sqrt(__squaredNorm)  << std::endl;
      s << std::endl;
      return s.str();
   }

   virtual std::string toSimpleString(
         )
   {
      std::ostringstream s;
      s << "vector: "      << vector.head(5).transpose() << "..." << std::endl;
      s << "squaredNorm: " << __squaredNorm       << std::endl;
      s << "norm: "        << std::sqrt(__squaredNorm) << std::endl;
      s << std::endl;
      return s.str();
   }


}; // class VectorElementBase


}  // namespace LapTools


#endif  // __LAPTOOLS_VECTOR_ELEMENT_H
