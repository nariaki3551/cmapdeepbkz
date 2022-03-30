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

/**@file    Enumeration.cpp
 * @brief   Base class for Enumeration.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "Enumeration.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <new>
#include <utility>
#include "Config.h"
#include "Def.h"
#include "Lattice.h"
#include "Log.h"
#include "Timer.h"


namespace LapTools
{


///
/// @brief set parameters of enumeration
/// @param[in] inBegin  begin index of basis for reduction
/// @param[in] inEnd    end   index of basis for reduction
/// @param[in] inLowerBound
/// @param[in] inTimeLimit
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
void
Enumeration<BasisFloat, GSFloat, EnumGSFloat>::init(
      int inBegin,
      int inEnd,
      double inLowerBound,
      double inTimeLimit
      )
{
   begin = inBegin;
   end = inEnd;
   if( end < 0 ){ end = L->m - 1; }
   lowerBound = inLowerBound;
   timeLimit = inTimeLimit;

   bestObjectiveValue = L->B(begin);
   sqRadius = bestObjectiveValue;
   if( inLowerBound > 0 )
      sqRadius = lowerBound * lowerBound;
   nSearch = 0;
   if( verbose > 0 ) updateEnumCost();

   // Pruning
   for( int i = 0; i < end-begin+1; ++i ){ sqR[i] = sqRadius; }

   // setting parameters for main process of Enum
   mu = Lattice<BasisFloat, EnumGSFloat>(L->basis).mu;
   sigma .setZero();
   r     .setZero();
   rho   .setZero();
   coeffs.setZero();
   c     .setZero();
   w     .setZero();
   bestCoeffv.setZero();
   k = begin;
   lastNonZero = 0;

   coeffs[begin] = 1;
   for( int i = 0; i < end + 2; ++i ){ r[i] = i; }
}


///
/// @brief setter of lowerBound
/// @param[in] inLowerBound lower bound of norm
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
void
Enumeration<BasisFloat, GSFloat, EnumGSFloat>::setLowerBound(
      double inLowerBound
      )
{
   assert( inLowerBound > 0 );
   sqRadius = std::min(sqRadius, static_cast<EnumGSFloat>(inLowerBound * inLowerBound));
   for( int i = 0; i < end-begin+1; ++i ){ sqR[i] = sqRadius; }
   if( verbose > 0 ) updateEnumCost();
}


///
/// set initial point of search
/// @param[in] fixedCoefs coefficient array w such that coeffs.last(w.size()) = w
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
void
Enumeration<BasisFloat, GSFloat, EnumGSFloat>::setFixedCoeffs(
      LatticeVector<BasisFloat> &fixedCoeffs
      )
{
   assert( fixedCoeffs.size() > 0 );
   if( (fixedCoeffs.array() == 0).all() )
   {
      numFixedCoeffs = fixedCoeffs.size();
      return;
   }
   numFixedCoeffs = fixedCoeffs.size();
   coeffs.setZero();
   coeffs.tail(numFixedCoeffs) = fixedCoeffs;
   k = end - numFixedCoeffs + 1;

   // pre-calculation(r, rho, sigma, c)
   for( int j = begin; j <= end; ++j ){ r[j] = end; }
   for( int j = end; j >= k; --j )
   {
      for( int i = end; i >= k; i-- )
         sigma(i, j) = sigma(i + 1, j) + coeffs(i) * L->mu(i, j);
      c(j) = -sigma(j+1, j);
      rho(j) = rho(j+1) + (coeffs(j) - c(j)) * (coeffs(j) - c(j)) * L->B(j);
   }

   // pre-calculation(lastNonZero)
   lastNonZero = end;
   for(int i = end; i >= k; --i )
   {
      if( coeffs[i] == 0 ) lastNonZero--;
      else break;
   }
}


///
/// @brief Enum algorithm
/// @details Vector resv satisifies norm(\pi_k(resv))^2 <= sqR_k for all k in be from begin to end
///          resv = sum( coeffs[i]*v_i for all begin <= i <= end )
///          Enumeration the lattice vector under to find the shortest vector
///          in partial lattice matrix from row begin to row end
/// @param[out] resv          vector found by Enum
/// @parem[out] coeffResv     coefficient of vector found by Enum
/// @param[in] updateLimit    terminate when Enum update solution updateLimit times
/// @param[in] inTimeLimit    time limit
/// @return bool true if Enum found a shorter vector else false
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
bool
Enumeration<BasisFloat, GSFloat, EnumGSFloat>::projectedEnum(
      LatticeVector<BasisFloat> &resv,
      LatticeVector<BasisFloat> &coeffResv,
      int updateLimit,
      double inTimeLimit
      )
{
   if( inTimeLimit > 0 ) timeLimit = inTimeLimit;
   int nUpdate = 0;
   bool success = false, step = false;
   bool shouldAbort = false;
   if( verbose > 0 )
   {
      startTime = Timer::getElapsedTime();
      nextLogTime = config.LogIntervalEnum;
      if( nSearch == 0 )
      {
         std::cout << Log::getLogHeader() << std::endl;
         outputLog('S');
      }
   }

   while( true )
   {
      // preprocess
      if( timeLimit > 0 || verbose > 0 )
      {
         runningTime = Timer::getElapsedTime() - startTime;
         if( timeLimit > 0 && runningTime > timeLimit ){ outputLog('T'); break; }
         if( verbose > 0 && runningTime > nextLogTime ){ outputLog(' '); nextLogTime += config.LogIntervalEnum; }
      }
      communicate(shouldAbort);
      if( shouldAbort ){ break; }

      // search enumeration tree
      tmpScaler__ = ( coeffs.coeff(k) - c.coeff(k) );
      rho.coeffRef(k) = rho.coeff(k+1) + tmpScaler__ * tmpScaler__ * L->B.coeff(k);
      ++nSearch;
      if( rho.coeff(k) < sqR[k-begin] - SVP_MINEPSILON )
      {
         if( k == begin )
         {
            /// for sending a ( incumbent ) short vector
            success = true;
            ++nUpdate;
            bestCoeffv = coeffs;
            updateBestObjectiveValue('*', rho.coeff(begin));

            // check for termination
            if( rho.coeff(begin) < lowerBound * lowerBound ){ outputLog('L'); stepProjectedEnum(); break; }
            if( updateLimit > 0 && nUpdate >= updateLimit ) { outputLog('I'); stepProjectedEnum(); break; }
         }
         else
         {
            if( verbose > 3 )
            {
               std::cout << "[R] iter = " << nSearch << "; k = " << k << "; coeffs = " << coeffs.transpose() << std::endl;
            }
            --k;
            r.coeffRef(k) = std::max(r.coeff(k), r.coeff(k+1));
            for( int i = r.coeff(k); i >= k; --i )
               sigma.coeffRef(i,k) = sigma.coeff(i+1,k) + coeffs.coeff(i) * mu.coeff(i,k);
            c.coeffRef(k) = -sigma.coeff(k+1,k);
            coeffs.coeffRef(k) = std::round( c.coeff(k) );
            w.coeffRef(k) = 1;
         }
         postProcess();
      }
      else
      {
         if( verbose > 3 )
         {
            std::cout << "[P] iter = " << nSearch << "; k = " << k << " coeffs = " << coeffs.transpose() << std::endl;
         }
         step = stepProjectedEnum();
         if( !step ){ outputLog('E'); hasFinished = true; break; }
      }
   }

   if( success )
   {
      resv = getBestVector();
      coeffResv = bestCoeffv;
   }

   return success;
}


///
/// @brief search just one node
/// @return false if all nodes are searched else true
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
bool
Enumeration<BasisFloat, GSFloat, EnumGSFloat>::stepProjectedEnum(
   )
{
   ++k;
   if( k >= end - numFixedCoeffs + 1 ){ return false; }

   r.coeffRef(k) = k;
   if( k >= lastNonZero )
   {
      lastNonZero = k;
      ++coeffs.coeffRef(k);
   }
   else
   {
      if( coeffs.coeff(k) > c.coeff(k) )
         coeffs.coeffRef(k) -= w.coeff(k);
      else
         coeffs.coeffRef(k) += w.coeff(k);
      ++w.coeffRef(k);
   }
   return true;
}


///
/// @brief Calculate norm(pi_begin(v))^2
/// @param[in] inBegin begin index of blockSize matrix
/// @param[in] inEnd last index of blockSize matrix
/// @param[in] inCoeff coefficients of vector
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
double
Enumeration<BasisFloat, GSFloat, EnumGSFloat>::getProjectedSquaradNorm(
      int inBegin,
      int inEnd,
      LatticeVector<BasisFloat> &inCoeff
      )
{
   double u, projSquaredNorm = 0;
   for( int i = inBegin; i <= inEnd; ++i )
   {
      u = inCoeff.coeff(i);
      for( int j = i+1; j <= inEnd; ++j )
         u += mu.coeff(j, i) * inCoeff.coeff(j);
      projSquaredNorm += u * u * L->B.coeff(i);
   }
   return projSquaredNorm;
}


///
/// @brief update bestObjectiveValue
/// @param[in] sigh line-header character
/// @param[in] objectiveValue new objectiveValue
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
bool
Enumeration<BasisFloat, GSFloat, EnumGSFloat>::updateBestObjectiveValue(
      char sigh,
      EnumGSFloat objectiveValue
      )
{
   if( objectiveValue < bestObjectiveValue - SVP_MINEPSILON )
   {
      bestObjectiveValue = objectiveValue;
      outputLog(sigh);
      sqRadius = std::min(sqRadius, bestObjectiveValue);
      for( int i = 0; i < end-begin+1; ++i ){ sqR[i] = sqRadius; }
      if( verbose > 0 ) updateEnumCost();
      return true;
   }
   return false;
}


///
/// @brief getter of best vector
/// @return lattice vector
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
LatticeVector<BasisFloat> &
Enumeration<BasisFloat, GSFloat, EnumGSFloat>::getBestVector(
      )
{
   bestVector.setZero(L->n);
   for( int i = begin; i <= end; ++i )
   {
      for( int j = 0; j < L->n; ++j )
      {
         bestVector.coeffRef(j) += bestCoeffv(i) * L->basis.coeff(i, j);
      }
   }
   return bestVector;
}


///
/// @brief createt log
/// @param[in] sigh line-header character
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
void
Enumeration<BasisFloat, GSFloat, EnumGSFloat>::outputLog(
   char sigh
   )
{
   int vthreas = 1;
   if( verbose >= vthreas && ( !config.Quiet || osCsvLog ) )
   {
      std::string taskName       = "Enum";
      int         size           = end - begin + 1;
      long int    iter           = nSearch;
      double      progress       = ( static_cast<double>(nSearch) / enumCost ) * 100.0;
      double      shortestNorm   = std::sqrt(bestObjectiveValue);
      double      logCost        = std::log(enumCost);
      std::ostringstream s_appendix;
      std::string appendix = s_appendix.str();
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
                  *L,sigh,rank,threadId,taskName,
                  runningTime,size,iter,progress,logCost,shortestNorm
                  )
            << std::endl;  // appendix
      }
   }
}


///
/// @brief divide enumeration tree
/// @param[in] N number of sub enumeration tree
/// @note Depending on the input, it may only return a number of divisions smaller than N.
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
std::deque<LatticeVector<BasisFloat>>
Enumeration<BasisFloat, GSFloat, EnumGSFloat>::divide(
      int N
      )
{
   std::deque<LatticeVector<BasisFloat>> q;
   LatticeVector<BasisFloat> root(0);
   q.push_back(root);
   int Lower = 0, Upper = 0;
   while( static_cast<int>(q.size()) < N )
   {
      LatticeVector<BasisFloat> fixedCoeffs(q.front());
      q.pop_front();
      range(fixedCoeffs, Lower, Upper);
      if( (fixedCoeffs.array() == 0).all() )
      {
         Upper = std::max(-Lower, Upper); Lower = 0;
      }
      if( Upper < Lower )
      {
         return q;
      }
      for( int coeff = Lower; coeff <= Upper; ++coeff )
      {
         LatticeVector<BasisFloat> x(fixedCoeffs.size() + 1);
         x(0) = coeff;
         x.tail(fixedCoeffs.size()) = fixedCoeffs;
         q.push_back(x);
      }
   }
   return q;
}


///
/// @brief calculation lower and upper bound of coefficient
/// @param[in] fixedCoeffs w = fixedCoeffs; coefficient array such that d-tail of coeffs equals w, where d = w.size()
/// @param[out] Lower lower bound of d+1-tail coefficients
/// @param[out] Upper upper bound of d+1-tail coefficients
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
void
Enumeration<BasisFloat, GSFloat, EnumGSFloat>::range(
      LatticeVector<BasisFloat> &fixedCoeffs,
      int &Lower,
      int &Upper
     )
{
   int l = end - fixedCoeffs.size();
   MatrixMu<EnumGSFloat> __sigma(L->m+1, L->m);
   __sigma.setZero();
   for( int j = end; j > l; --j )
   {
      for( int i = end; i >= j; --i )
      {
         __sigma(i, j) = __sigma(i + 1, j) + L->mu(i, j) * fixedCoeffs(i - l - 1);
      }
   }
   for( int i = end; i > l; --i )
   {
      __sigma(i, l) = __sigma(i + 1, l) + L->mu(i, l) * fixedCoeffs(i - l - 1);  // j = l
   }
   double tmp = 0.0;
   for( int j = l+1; j <= end; ++j )
      tmp += __sigma(j, j) * __sigma(j, j) * L->B(j);
   double Y = (sqR[l] - tmp) / L->B(l);
   Lower = std::ceil(- std::sqrt(Y) - __sigma(l+1, l));
   Upper = std::floor(std::sqrt(Y) - __sigma(l+1, l));
}




///
/// instantiation
///
template class Enumeration<int, double, double>;
template class Enumeration<int, long double, double>;
template class Enumeration<int, long double, long double>;
template class Enumeration<long int, double, double>;
template class Enumeration<long int, long double, double>;
template class Enumeration<long int, long double, long double>;


}  // namespace LapTools

