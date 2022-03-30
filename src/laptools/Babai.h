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

/**@file    Babai.h
 * @brief   Base class for Babai.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __LAPTOOLS_BABAI_H__
#define __LAPTOOLS_BABAI_H__

#include <cmath>
#include <algorithm>

#include "Lattice.h"


namespace LapTools
{

///
/// @class Babai
///
/// @tparam BasisFloat floating point of basis matrix
/// @tparam GSFloat floating point of Gram-Schmidt matrix
/// @tparam EnumGSFloat floating point of Gram-Schmidt matrix for Enumeration
///
template<typename BasisFloat=int, typename GSFloat=double>
class Babai
{

using LatticePtr = std::shared_ptr<Lattice<BasisFloat, GSFloat>>;


private:

   LatticePtr  L;          ///< lattice
   Config      config;     ///< hyper parameters
   int         verbose;    ///< verbose <= 0: none, 1: light, 2: medium, 3: heavy

	VectorMu<GSFloat> __nu; ///< allocated memory


public:

   ///
   /// @brief Constructor
   ///
   Babai(){}


   ///
   /// @brief Constructor
   /// @param[in] inL         Lattice
   /// @param[in] inVerbose   <= 0 : not, 1 : light, 2: heavy
   ///
   Babai(
         LatticePtr inL,
         int inVerbose=0
         )
      :
         verbose(inVerbose)
   {
      L = inL;
      config = L->config;
      __nu.setZero(L->m);
   }


   ///
   /// @brief Babai Algorithm
   /// @param[in, out] v lattice vector
   /// @param[in] coeffv coefficients v_i where v = sum(v_i b_i for i)
   /// @note Transform lattice vector v to v - w, where w is the closest vector to v from sublattice {b_0,...b_k-1}_Z.
   ///
   void reduce(
         LatticeVector<BasisFloat> &v,
         LatticeVector<BasisFloat> &coeffv
         )
   {
      int i, j;
      // nu(i) = <v, b_i^*> / |b_i^*|^2 = sum( coeffv(j) * mu(j, i) forall i <= j < m )
      __nu.setZero();
	   for( i = 0; i < L->m; ++i )
      {
         __nu.coeffRef(i) = static_cast<GSFloat>(coeffv.coeff(i));
         for( j = i + 1; j < L->m; ++j )
         {
            __nu.coeffRef(i) += static_cast<GSFloat>(coeffv.coeff(j)) * L->mu.coeff(j, i);
         }
      }

	   BasisFloat c = 0;
      GSFloat cc = 0.0;
      for( i = L->m; i > -1; --i )
      {
         if( std::fabsl(__nu.coeff(i)) > config.eta )
         {
            c = std::roundl(__nu.coeff(i));
            cc = static_cast<GSFloat>(c);
            for( j = 0; j < L->m; ++j )
            {
               v.coeffRef(j) -= c * L->basis.coeff(i, j);
               __nu.coeffRef(j) -= cc * L->mu.coeff(i, j);
            }
         }
      }
   }


}; // class Babai


}  // namespace LapTools


// /***********************************************************************************
// Babai Algorithm
// Calculate lattice vector v in Lattice(b0,...bn)
// such that v is the closest vector to vector w in Z^n,
// where k = begin, l = end
//
// Input
//    VectorXi &w: target vector
//
// Output
//    VectorXi  v: lattice vector
// ************************************************************************************/
// inline VectorXi
// Lattice::Babai(
//       VectorXi &w
//       )
// {
//    int begin = 0;
//    int end = m-1;
//    return Babai(begin,end,w);
// }

// /***********************************************************************************
// Babai Algorithm
// Calculate lattice vector v = sum( cj bj forall k <= j <= l ) (cj in Z)
// such that v is the closest vector to vector w in Z^n,
// where k = begin, l = end
//
// Input
//    int   begin: index of basis
//    int     end: index of basis
//    VectorXi &w: target vector
//
// Output
//    VectorXi  v: projected lattice vector
// ************************************************************************************/
// inline VectorXi
// Lattice::Babai(
//       int begin,
//       int end,
//       VectorXi &w
//       )
// {
//    if( begin > 0 || end < m-1 )
//    {
//       int blockSize = end - begin + 1;
//       MatrixXi subBasis = basis.block(begin,0,blockSize,n);
//       Lattice subL = Lattice(subBasis);
//       return subL.Babai(w);
//    }
//
//    int i, j;
//    // generate GSO matrix
//    MatrixXd GSO(m, n);
//    for( i = 0; i < m; i++ )
//    {
//       for( j = 0; j < n; j++ ) { GSO(i,j) = basis(i,j); }
//       for( j = 0; j < i; j++ ) { GSO.row(i) -= mu(i,j)*GSO.row(j); }
//    }
//
//    VectorXi b = w;
//    int c;
//    double dot;
//    for( i = end; i >= begin; i-- )
//    {
//       dot = 0;
//       for( j = 0; j < n; j++ ) dot += b(j) * GSO(i, j);
//       c = round( dot / sqnorm(i) );
//       b -= c*basis.row(i);
//    }
//    return w - b;
// }


#endif  // __LAPTOOLS_BABAI_H__
