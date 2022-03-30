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

/**@file    KleinSampler.cpp
 * @brief   Base class for Klein sampler.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "KleinSampler.h"

#include <cstdlib>
#include "Def.h"
#include "Lattice.h"


namespace LapTools
{


///
/// @brief sample lattice vectors
/// @param[in] simpleReduce if it is true, simplify vector using basis
///
template<typename BasisFloat, typename GSFloat>
LatticeVector<BasisFloat> &
KleinSampler<BasisFloat, GSFloat>::sample(
      bool simpleReduce
      )
{
   coeff.setZero();
   for( int i = L->m-1; i > -1; --i )
   {
      coeff(i) = sampleZ(coeff(i), sqSPrimes(i));
      for( int j = 0; j < i; ++j ){ coeff(j) -= coeff(i) * L->mu(i,j); }
   }

   v.setZero();
   for( int i = 0; i < L->m; ++i )
   {
      for( int j = 0; j < L->n; ++j )
         v.coeffRef(j) += static_cast<BasisFloat>(coeff.coeff(i) * L->basis.coeff(i, j));
   }

   if( simpleReduce ){ reduce(v); }
   return v;
}


///
/// @breif return uniform random int
///
template<typename BasisFloat, typename GSFloat>
long
KleinSampler<BasisFloat, GSFloat>::sampleZ(
      double c,
      double sSquare
      )
{
   s__ = std::sqrt(sSquare);
   minC__ = std::floor( c - s__ * t );
   maxC__ = std::ceil ( c + s__ * t );

   uniformX = std::uniform_int_distribution<long>{minC__, maxC__};
   while( true )
   {
      x__ = uniformX(mt);
      tmp__ = x__ - c;
      rho__ = std::exp( - M_PI * tmp__ * tmp__ / sSquare );
      if( uniform(mt) <= rho__ )
      {
         return x__;
      }
   }
}


///
/// @breif reduce vector using basis of Lattice
/// @param[in] vector lattice vector
///
template<typename BasisFloat, typename GSFloat>
void
KleinSampler<BasisFloat, GSFloat>::reduce(
      LatticeVector<BasisFloat> &vector
      )
{
   int i, q;
   BasisFloat dot;
   bool reduced = true;
   while( reduced )
   {
      reduced = false;
      for( int k = L->m-1; k > -1; --k )
      {
         dot = 0;
         for( i = L->n-1; i > -1; --i )
            dot += vector.coeff(i) * L->basis.coeff(k, i);

         if( std::abs(dot) >= sqNorms.coeff(k) * 0.5001 )
         {
            q = std::round(dot / sqNorms.coeff(k) );  // q = round( <v, bi> / |bi|^2 )
            for( i = L->n-1; i > -1; --i )
               vector.coeffRef(i) -= q * L->basis.coeff(k, i); // v - q*w
            reduced = true;
         }
      }
   }
}


///
/// instantiation
///
template class KleinSampler<int, double>;
template class KleinSampler<int, long double>;
template class KleinSampler<long int, double>;
template class KleinSampler<long int, long double>;


}  // namespace LapTools
