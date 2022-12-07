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

/**@file    DeepLll.cpp
 * @brief   Base class for DeepLLL.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "DeepLll.h"

#include <fstream>
#include <new>
#include <algorithm>
#include <cmath>
#include "Config.h"
#include "Def.h"
#include "Lattice.h"
#include "Log.h"
#include "Timer.h"



namespace LapTools
{


///
/// @brief DeepLll reduction algorithm
///        for Matrix[begin:end] and algorithm start at "start index";
/// @param[in] start index of basis
/// @param[in] begin index of basis
/// @param[in] end   index of basis
/// @return bool true: normal termination, false: abnormal termination
///
template<typename BasisFloat, typename GSFloat>
bool
DeepLll<BasisFloat, GSFloat>::deeplll(
      int start,
      int begin,
      int end
      )
{
   if( end < 0 ) end = L->m - 1;
   double startTime = Timer::getElapsedTime();
   int i, j, k = std::max(1, start);
   GSFloat c;
   nIter = 0;
   LatticeBasis<BasisFloat> *basis = &(L->basis);
   MatrixMu<GSFloat> *mu = &(L->mu);
   VectorB<GSFloat> *B = &(L->B);

   bool shouldAbort = false;
   while( k <= end )
   {
      communicate(shouldAbort);
      if( shouldAbort ){ break; }

      ++nIter;
      L->sizeReduce(k, config.eta);
      c = basis->row(k).squaredNorm();
      for( i = 0; i < begin; ++i )
      {
         c -= mu->coeff(k,i) * mu->coeff(k,i) * B->coeff(i);  // c = norm(pi_i(bk))
      }
      for( i = begin; i < k; ++i )
      {
         if( c > config.delta * B->coeff(i) )
         {
            c -= mu->coeff(k,i) * mu->coeff(k,i) * B->coeff(i);
         }
         else
         {
            for( j = k; j > i; --j )
            {
               basis->row(j).swap(basis->row(j-1));
            }
            if( !GSOUpdateDeepLll(i, k) ) return false;
            if( i == 0 )
            {
               runningTime = Timer::getElapsedTime() - startTime;
               updateBestObjectiveValue('*');
            }
            k = i;
            // k = max(i, begin);
            break;
         }
      }
      k++;
   }
   runningTime = Timer::getElapsedTime() - startTime;
   return true;
}


///
/// @brief update GSO matrix with deep-insertion
/// @param[in] i
/// @param[in] k
/// @details Update GSO information after executing DeepLll from row i to row k
///          {...,b_i,...,b_k,...} -> {...,b_i-1,b_k,b_i,...,b_k-1,b_k+1,...}
/// @return true if B(i) >= B'(i) else false
///
template<typename BasisFloat, typename GSFloat>
bool
DeepLll<BasisFloat, GSFloat>::GSOUpdateDeepLll(
      int i,
      int k
      )
{
   int j, l;
   MatrixMu<GSFloat> *mu = &(L->mu);
   VectorB<GSFloat> *B = &(L->B);
   auto tmp0 = B->coeff(i);   /// for check
   S__.setZero();

   P__.coeffRef(k) = B->coeff(k);
   D__.coeffRef(k) = B->coeff(k);
   for( j = k - 1; j > i - 1; --j )
   {
      P__.coeffRef(j) = mu->coeff(k,j) * B->coeff(j);
      D__.coeffRef(j) = D__.coeff(j+1) + mu->coeff(k,j) * P__.coeff(j);
   }
   for( j = k; j > i; --j )
   {
      T__ = mu->coeff(k,j-1) / D__.coeff(j);
      for( l = L->m-1; l > k; --l )
      {
         S__.coeffRef(l-i-2) += mu->coeff(l,j) * P__.coeff(j);
         mu->coeffRef(l,j) = mu->coeff(l,j-1) - T__ * S__.coeff(l-i-2);
      }
      for( l = k; l > j; --l )
      {
         S__.coeffRef(l-i-2) += mu->coeff(l-1,j) * P__.coeff(j);
         mu->coeffRef(l,j) = mu->coeff(l-1,j-1) - T__ * S__.coeff(l-i-2);
      }
   }
   T__ = 1.0 / D__.coeff(i);
   for( l = L->m - 1; l > k; --l )
   {
      mu->coeffRef(l,i) = T__ * (S__.coeff(l-i-2) + mu->coeff(l,i) * P__.coeff(i));
   }
   for( l = k; l > i + 1; --l )
   {
      mu->coeffRef(l,i) = T__ * (S__.coeff(l-i-2) + mu->coeff(l-1,i) * P__.coeff(i));
   }
   mu->coeffRef(i+1,i) = T__ * P__.coeff(i);
   for( l = k; l > i; --l )
   {
      mu->block(l,0,1,i).swap(mu->block(l-1,0,1,i));
   }
   for( j = k; j > i; --j )
   {
      B->coeffRef(j) = D__.coeff(j) * B->coeff(j-1) / D__.coeff(j-1);
   }
   B->coeffRef(i) = D__.coeff(i);

   if( B->coeff(i) >= tmp0 )
   {
      std::cout << "Rank " << rank << " Thread " << threadId
         << " Insertion error at index " << i << " in DeepLll::GSOUpdateDeepLll" << std::endl;
      L->setGSO(0, L->m-1);
      return false;
   }
   return true;
}


///
/// @brief update bestObjectiveValue
/// @param[in] sigh line-header character
///
template<typename BasisFloat, typename GSFloat>
bool
DeepLll<BasisFloat, GSFloat>::updateBestObjectiveValue(
      char sigh
      )
{
   if( L->B(0) < bestObjectiveValue - SVP_MINEPSILON )
   {
      bestObjectiveValue = L->B(0);
      outputLog(sigh);
      return true;
   }
   return false;
}


///
/// @brief createt log
/// @param[in] sigh line-header character
///
template<typename BasisFloat, typename GSFloat>
void
DeepLll<BasisFloat, GSFloat>::outputLog(
   char sigh
   )
{
   int vthreas = 1;
   if( this->getVerbose() >= vthreas && ( !config.Quiet || osCsvLog ) )
   {
      std::string taskName       = "DeepLll";
      int         size           = static_cast<int>(config.delta*100);
      long int    iter           = nIter;
      double      progress       = 0.0;
      double      shortestNorm   = L->shortestNorm();
      double      logCost        = std::log( L->enumCost() );
      std::string appendix       = "";
      if( !config.Quiet )
      {
         std::cout
            << Log::getLog(
                  *L,sigh,rank,threadId,taskName,runningTime,
                  size,iter,progress,logCost,shortestNorm,appendix
                  )
            << std::endl;
      }
      if( osCsvLog )
      {
         *osCsvLog
            << Log::getCsvLog(
                  *L,sigh,rank,threadId,taskName,runningTime,
                  size,iter,progress,logCost,shortestNorm
                  )
            << std::endl;
      }
   }
}


///
/// instantiation
///
template class DeepLll<int, double>;
template class DeepLll<int, long double>;
template class DeepLll<long int, double>;
template class DeepLll<long int, long double>;


}  // namespace LapTools
