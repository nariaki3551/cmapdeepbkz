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

/**@file    KleinSampler.h
 * @brief   Base class for Klein sampler.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __LAPTOOLS_KLEIN_SAMPLER_H__
#define __LAPTOOLS_KLEIN_SAMPLER_H__

#include <cmath>
#include <memory>
#include <random>
#include "Def.h"
#include "Lattice.h"
#include "VectorElement.h"
#include "Config.h"


namespace LapTools
{


///
/// @class KleinSampler
///
/// @tparam BasisFloat floating point of basis matrix
/// @tparam GSFloat floating point of Gram-Schmidt matrix
///
template<typename BasisFloat=int, typename GSFloat=double>
class KleinSampler
{

using VectorElementType = std::shared_ptr<VectorElementBase<BasisFloat>>;
using LatticePtr = std::shared_ptr<Lattice<BasisFloat, GSFloat>>;


protected:

   LatticePtr  L;                ///< lattice
   Config   config;              ///< hyper parameters

   double t;                     ///< log(L->m)
   VectorMu<GSFloat> sqSPrimes;  ///< |s'_i|^2
   VectorMu<GSFloat> sqNorms;    ///< squared norms of b_i
   VectorMu<GSFloat> coeff;      ///< for sample
   LatticeVector<BasisFloat> v;  ///< generated vector

   std::mt19937 mt;              ///< randomize generartor

   double s__;                   ///< allocated memory for sampleZ
   long minC__;                  ///< allocated memory for sampleZ
   long maxC__;                  ///< allocated memory for sampleZ
   long x__;                     ///< allocated memory for sampleZ
   double rho__;                 ///< allocated memory for sampleZ
   double tmp__;                 ///< allocated memory for sampleZ
   std::uniform_int_distribution<long> uniformX;
   std::uniform_real_distribution<double> uniform;


public:

   ///
   /// @brief constructor
   ///
   KleinSampler(){};


   ///
   /// @brief constructor
   /// @param[in] inL         Lattice
   ///
   KleinSampler(
         LatticePtr inL
         )
   {
      L = inL;

      // mt.seed(0);
      std::random_device rnd;
      mt.seed(rnd());
      uniform = std::uniform_real_distribution<double>{0.0, 1.0};


      sqSPrimes.setZero(L->m);
      sqNorms  .setZero(L->m);
      coeff    .setZero(L->m);
      v        .setZero(L->n);

      t = std::log(L->m);
      GSFloat sSquare = L->B.maxCoeff() * t;
      for( int i = 0; i < L->m; ++i ){ sqSPrimes(i) = sSquare / L->B(i); }
      for( int i = 0; i < L->m; ++i ){ sqNorms(i) = L->basis.row(i).squaredNorm(); }
   }


   ///
   /// @brief sample lattice vectors
   /// @param[in] simpleReduce if it is true, simplify vector using basis
   ///
   LatticeVector<BasisFloat> &
   sample(
         bool simpleReduce=false
         );


   ///
   /// @breif return uniform random int
   ///
   long sampleZ(
         double c,
         double sSquare
         );


   ///
   /// @breif reduce vector using basis of Lattice
   /// @param[in] vector lattice vector
   ///
   void
   reduce(
         LatticeVector<BasisFloat> &vector
         );


}; // class KleinSampler


}  // namespace LapTools


#endif  // __LAPTOOLS_KLEIN_SAMPLER_H__
