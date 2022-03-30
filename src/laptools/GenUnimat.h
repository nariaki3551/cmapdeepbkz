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

/**@file    GenUnimat.h
 * @brief   Functions for generating randomized unimoduler matrix.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __LAPTOOLS_GEN_UNIMAT_H__
#define __LAPTOOLS_GEN_UNIMAT_H__

#include <iostream>
#include <random>
#include <vector>

#include <eigen3/Eigen/Core>


using namespace Eigen;

namespace LapTools
{


///
/// @brief shuffle rows of matrix
/// @param[in] M matrix
/// @param[in] seed randomize seed
///
inline bool
shuffleRows(
      MatrixXi &M,
      int seed
      )
{
   std::mt19937 engine(seed);

   int n = M.rows();
   std::vector<int> v(n);
   for( int i = 0; i < n; i++ ) v[i] = i;

   std::shuffle(v.begin(), v.end(), engine);

   MatrixXi P = MatrixXi::Zero(n, n);
   for( int i = 0; i < n; i++ ) P(i,v[i]) = 1;
   M = P*M;
   return true;
}


///
/// @brief generate unimodular matrix using LU decomposition
/// @param[in] M matrix
/// @param[in] seed randomize seed
/// @param[in] scale scale-parameter of unimodular matrix
/// @note 1. generate random integer matrix R whose diagonal is 1
///       2. LU-decompose R = L + U, and it holds det(L) = 1 and det(U) = 1
///       3. shuffle L and U's rows, it dosenot change determine of L and U
///       4. generate M = L * U, whose determine det(M) = det(L)det(U) = 1
///
inline bool
genUnimatLU(
      MatrixXi &M,
      int seed,
      double scale=1.2
      )
{
   int n = M.rows();
   std::srand(seed);
   MatrixXi R = (scale*MatrixXf::Random(n,n)).cast<int>();
   for(int i = 0; i < n; ++i ) R(i, i) = 1;
   MatrixXi L = R.triangularView<Lower>();
   MatrixXi U = R.triangularView<Upper>();

   shuffleRows(L,seed);
   shuffleRows(U,seed);

   M = L*U;
   return true;
}


///
/// @brief generate unimodular matrix using LU decomposition
/// @param[in] M matrix
/// @param[in] seed randomize seed
/// @param[in] scale scale-parameter of unimodular matrix
/// @note 1. generate identity matrix M
///       2. shuffle M's rows
///       3. triangler transformation with coefficnets in -1, 0, 1
///
inline bool
genUnimatRandomInsertions(
      MatrixXi &M,
      int seed,
      double scale=3
      )
{
   int n = M.rows();
   M = MatrixXi::Identity(n, n);
   std::uniform_int_distribution<int> uniform{0, n-1};
   std::mt19937 mt(seed);  ///< randomize generartor

   // 1. permute rows
   int nIter = 4 * n;
   for( int i = 0; i < nIter; ++i )
   {
      int a = uniform(mt);
      int b = a;
      while( b == a ){ b = uniform(mt); }
      M.row(a).swap(M.row(b));
   }

   // 2. triangular transformation matrix with coefficnets in -1, 0, 1
   std::uniform_int_distribution<int> flag{0, 1};
   for( int a = 0; a < n - 2; ++a )
   {
      std::uniform_int_distribution<int> uniform_b{a+1, n-1};
      for( int i = 0; i < scale; ++i )
      {
         int b = uniform_b(mt);
         if( flag(mt) ){ M.row(a) += M.row(b); }
         else          { M.row(a) -= M.row(b); }
      }
   }
   return true;
}


///
/// @brief generate unimodular matrix using LU decomposition
/// @param[in] M matrix
/// @param[in] seed randomize seed
/// @param[in] scale scale-parameter of unimodular matrix
///
inline bool
genUnimat(
      MatrixXi &M,
      int seed,
      std::string randomizeType,
      double scale
      )
{
   if( randomizeType == "LU" )
   {
      return genUnimatLU(M, seed, scale);
   }
   else if( randomizeType == "Swap" )
   {
      return genUnimatRandomInsertions(M, seed, 0.0);
   }
   else if( randomizeType == "Fplll" )
   {
      return genUnimatRandomInsertions(M, seed, 3.0);
   }
   else
   {
      std::cout << "randomizeType " << randomizeType << " is invalid" << std::endl;
      exit(1);
   }
}

}  // namespace LapTools


#endif // __LAPTOOLS_GEN_UNIMAT_H__
