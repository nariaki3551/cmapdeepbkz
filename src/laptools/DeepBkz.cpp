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

/**@file    DeepBkz.cpp
 * @brief   Base class for DeepBKZ.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "DeepBkz.h"

#include <algorithm>
#include <assert.h>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <unordered_map>
#include "Config.h"
#include "DeepLll.h"
#include "Def.h"
#include "Enumeration.h"
#include "Lattice.h"
#include "Log.h"
#include "Timer.h"
#include "Utils.h"


namespace LapTools
{


///
/// @brief DeepBkz reduction algorithm
///        for Matrix[begin:end] and algorithm start at "start index"
/// @param[in] inBlocksize block size;
/// @param[in] inTimeLimit
/// @return bool true: normal terminate, false: abnormal one
/// @remark If lowerBound is positive value, algorithm terminates wheh it founds the lattice vector whose norm is less than lowerBound
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
bool
DeepBkz<BasisFloat, GSFloat, EnumGSFloat>::deepbkz(
      int inBlocksize,
      double inTimeLimit
      )
{
   // initialize
   assert( inBlocksize > 0 );
   blocksize = inBlocksize;
   if( inTimeLimit >= 0 ) timeLimit = inTimeLimit;

   if( verbose > 1 && rank == 1 && threadId == 0 )
      std::cout << Log::getLogHeader() << std::endl;

   outputLog('S');

   // main-loop
   double nextLogTime = runningTime + this->config.LogIntervalDeepBkz;
   bool shouldAbort = false;
   while( !hasReduced )
   {
      if( timeLimit > 0 && runningTime > timeLimit ){ break; }

      communicate(shouldAbort);
      if( shouldAbort ) { break; }

      preprocess(shouldAbort);
      if( shouldAbort ) { break; }

      tour(shouldAbort);
      if( shouldAbort ) { break; }

      if( nextLogTime < runningTime )
      {
         outputLog(' ');
         nextLogTime = runningTime + this->config.LogIntervalDeepBkz;
      }
   }

   outputLog('E');
   return true;
}


///
/// @brief preprocess of DeepBkz main loop
/// @param[out] shouldAbort true if it should abort else false
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
bool
DeepBkz<BasisFloat, GSFloat, EnumGSFloat>::preprocess(
      bool &shouldAbort
      )
{
   double startTime = Timer::getElapsedTime();
   lllObj.setVerbose(verbose);
   lllObj.setOsCsvLog(osCsvLog);
   lllObj.deeplll();
   lllObj.setVerbose(0);
   lllObj.setOsCsvLog(nullptr);
   runningTime += Timer::getElapsedTime() - startTime;
   updateBestObjectiveValue('*');
   return true;
}


///
/// @brief run tour
/// @param[out] shouldAbort true if it should abort else false
/// @remark runningTime is updated in this function
/// @remark autoAborted is updaetd in this function
/// @remark If lowerBound is positive value, wheh it founds the lattice vector whose norm is less than lowerBound
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
bool
DeepBkz<BasisFloat, GSFloat, EnumGSFloat>::tour(
      bool &shouldAbort
      )
{
   int begin = 0, end = L->m - 1;
   double startTime = Timer::getElapsedTime();
   bool fullTour = true;
   shouldAbort = false;
   nTour++;

   // int iniK = std::round(11.46+0.0757*L->n);
   int iniK = 0;
   int z = std::max(iniK, begin);   // counter
   int k = std::max(iniK, begin);   // index of the beginning of blocksize
   communicateInTour(shouldAbort);
   for( ; k < end && !shouldAbort; k++ )
   {
      bool enumSuccess = step(k, z, shouldAbort);
      if( !fullTour && enumSuccess ){ break; }
      communicateInTour(shouldAbort);
      if( lowerBound > 0 && L->shortestNorm() < lowerBound ){ shouldAbort = true; }
   }
   if( z == end ){ hasReduced = true; }   // not updated

   // check auto abort
   if( this->config.AutoAbortThreashold > 0 )
   {
      double slopeGSA = L->slopeGSA();
      if( slopeGSA < prevSlopeGSA ){ nSlopeNoDecrease++; }
      else{ nSlopeNoDecrease = 0; }
      prevSlopeGSA = slopeGSA;
      if( nSlopeNoDecrease >= this->config.AutoAbortThreashold )
      {
         outputLog('A');
         autoAborted = true;
         shouldAbort = true;
      }
   }

   // check cycle updation
   hashCounter.push(L->hash());
   if( hashCounter.mostFrequentCount() > 0.1 * hashCounterSize )
   {
      outputLog('C');
      shouldAbort = true;
   }
   runningTime += Timer::getElapsedTime() - startTime;
   return true;
}


///
/// @brief one loop of algorithm
/// @param[in] k
/// @param[out] z
/// @param[out] shouldAbort true if it should abort else false
/// @return bool return true if enumeration process is successed else false
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
bool
DeepBkz<BasisFloat, GSFloat, EnumGSFloat>::step(
      int k,
      int &z,
      bool &shouldAbort
      )
{
   int begin = 0, end = L->m - 1;
   int l = std::min(k+blocksize-1, end);
   int h = std::min(l+1, end);

   // preprocess
   if( timeLimit > 0 && runningTime >= timeLimit ){ outputLog('T'); shouldAbort = true; }
   if( shouldAbort ){ return true; }

   // sub-SVP
   enumObj.init(k, l, -1, (timeLimit > 0 ? timeLimit-runningTime : -1));
   bool enumSuccess = enumObj.projectedEnum(v, coeffv);

   // insert
   if( enumSuccess ){ postProcess(k, l); }

   // reduction for next sub-SVP
   if( enumSuccess )
   {
      lllObj.deeplll(k, begin, h);
      z = 0;
   }
   else
   {
      lllObj.deeplll(h-1, begin, h);
      z++;
   }

   updateBestObjectiveValue('*');
   return enumSuccess;
}


///
/// @breif MLLL process after Enum found vector v during DeepBkz, and
///        v is inserted in the basis
/// @param[in] k  begin index of block
/// @param[in] l  end index of block
/// @return
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
bool
DeepBkz<BasisFloat, GSFloat, EnumGSFloat>::postProcess(
      int k,
      int l
      )
{
   // r is tha last index whose coefficient is not zero
   int r = 0;
   for( r = l; r >= k; --r )
   {
      if( coeffv(r) != 0 ){ break; }
   }

   // if( begin > 0 && !(abs(coeffv(r))==1) ){ return false; }
   if( std::abs(coeffv(r)) == 1 )
   {
      if( coeffv(r) == -1 ){ v = -v; coeffv = -coeffv; }
      for( int i = r; i >= k+1 ; i-- )
         L->basis.row(i).swap(L->basis.row(i-1));
      L->basis.row(k) = v;
      GSOupdateBkz(k, r);
   }
   else
   {
      L->insertMlll(v, 0, k, l); // insert v between k-1 and k and reduce {b_0 ... b_l}
   }
   return true;
}


///
/// @brief update mu and sqnorm for basis matrix
///        {b_0,...,b_k-1,v,b_k,...,b_l, b_l+1}
/// @param[in] k index of basis
/// @param[in] r tha last index whose coefficient is not zero
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
bool
DeepBkz<BasisFloat, GSFloat, EnumGSFloat>::GSOupdateBkz(
      int k,
      int r
      )
{
   int i, j, s, l;
   MatrixMu<GSFloat> *mu = &(L->mu);
   VectorB<GSFloat> *B = &(L->B);
   auto tmp0 = B->coeff(k);
   S__.setZero();
   int begin = 0, end = L->m-1;

   for( i = 0; i <= r; ++i )
   {
      nu__.coeffRef(i) = coeffv.coeff(i);
      for( s = i+1; s <= r; ++s )
      {
         nu__.coeffRef(i) += mu->coeff(s,i) * coeffv.coeff(s);
      }
   }
   D__.coeffRef(r) = B->coeff(r);
   for( l = r-1; l >= k; --l )
   {
      D__.coeffRef(l) = D__.coeff(l+1) + nu__.coeff(l) * nu__.coeff(l) * B->coeff(l);
   }
   for( j = r; j >= k+1; --j )
   {
      T__ = nu__.coeff(j-1) / D__.coeff(j);
      for( i = L->m-1; i >=r+1; --i )
      {
         S__.coeffRef(i) += nu__.coeff(j) * mu->coeff(i,j) * B->coeff(j);
         mu->coeffRef(i,j) = mu->coeff(i,j-1) - T__ * S__.coeff(i);
      }
      for( i = r; i >= j+1; --i )
      {
         S__.coeffRef(i) += nu__.coeff(j) * mu->coeff(i-1,j) * B->coeff(j);
         mu->coeffRef(i,j) = mu->coeff(i-1,j-1) - T__ * S__.coeff(i);
      }
   }
   T__ = 1.0 / D__.coeff(k);
   for( i = L->m-1; i >= r+1; --i )
   {
      mu->coeffRef(i,k) = T__ * ( S__.coeff(i) + mu->coeff(i,k) * nu__.coeff(k) * B->coeff(k) );
   }
   for( i = r; i >= k+2; --i )
   {
      mu->coeffRef(i,k) = T__ * ( S__.coeff(i) + mu->coeff(i-1,k) * nu__.coeff(k) * B->coeff(k) );
   }
   mu->coeffRef(k+1,k) = T__ * nu__.coeff(k) * B->coeff(k);
   for( j = 0; j <= k-1; ++j )
   {
      for( i = r; i >= k+1; --i )
      {
         mu->coeffRef(i,j) = mu->coeff(i-1,j);
      }
      mu->coeffRef(k,j) = nu__.coeff(j);
   }
   for( i = r; i >=k+1; --i )
   {
      B->coeffRef(i) = D__.coeff(i) * B->coeff(i-1) / D__.coeff(i-1);
   }
   B->coeffRef(k) = D__.coeff(k);

   if( B->coeff(k) >= tmp0 )
   {
      std::cout << "Insertion error at index " << k << " in DeepBKZ::GSOupdateBkz" << std::endl;
      L->setGSO(begin, end);
      return true;
   }

   return true;
}


///
/// @brief update bestObjectiveValue
/// @param[in] sigh line-header character
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
bool
DeepBkz<BasisFloat, GSFloat, EnumGSFloat>::updateBestObjectiveValue(
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
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
void
DeepBkz<BasisFloat, GSFloat, EnumGSFloat>::outputLog(
   char sigh,
   std::ofstream *os
   )
{
   int vthreas = 1;
   if( verbose >= vthreas && ( !this->config.Quiet || osCsvLog || os ) )
   {
      std::string taskName       = "DeepBkz";
      int         size           = blocksize;
      long int    iter           = nTour;
      double      progress       = 0.0;
      double      shortestNorm   = L->shortestNorm();
      double      logCost        = std::log(L->enumCost());
      double      slopeGSA       = L->slopeGSA();
      double      halfSlopeGSA   = L->slopeGSA(std::round(L->m*0.5));
      double      logCostGH      = std::log(L->enumCost(L->GH));
      int nNotChangeRows = 0;
      for( int i = 0; i < L->m; i++ )
      {
         if( L->basis.row(i) != prevBasis.row(i) ) break;
         nNotChangeRows++;
      }
      std::ostringstream s_appendix;
      s_appendix
         << "  --  "
         << "rho "            << std::fixed << std::setprecision(5) << slopeGSA
         << ", half_rho "     << std::fixed << std::setprecision(5) << halfSlopeGSA
         << ", not_change "   << std::right << std::setw(2)         << nNotChangeRows
         << ", enumCost_GH "  << logCostGH;
      std::string appendix = s_appendix.str();
      if( !this->config.Quiet )
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
            << "," << slopeGSA
            << "," << halfSlopeGSA
            << "," << nNotChangeRows
            << "," << logCostGH
            << "," << mergeTime
            << std::endl;  // appendix
      }
      if( os )
      {
         *os
            << Log::getCsvLog(
                  *L,sigh,rank,threadId,taskName,runningTime,
                  size,iter,progress,logCost,shortestNorm
                  )
            << "," << slopeGSA
            << "," << halfSlopeGSA
            << "," << nNotChangeRows
            << "," << logCostGH
            << "," << mergeTime
            << std::endl;  // appendix
      }
      prevBasis = L->basis;
   }
}


///
/// @brief get header of csv formed log
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
std::string
DeepBkz<BasisFloat, GSFloat, EnumGSFloat>::getCsvLogHeader(
   )
{
   std::ostringstream s;
   s  << Log::getCsvLogHeader()
      << ",slopeGSA"
      << ",halfSlopeGSA"
      << ",nNotChangeRows"
      << ",logCostGH"
      << ",mergeTime";
   return s.str();
}


///
/// instantiation
///
template class DeepBkz<int, double, double>;
template class DeepBkz<int, long double, double>;
template class DeepBkz<int, long double, long double>;
template class DeepBkz<long int, double, double>;
template class DeepBkz<long int, long double, double>;
template class DeepBkz<long int, long double, long double>;


}  // namespace LapTools
