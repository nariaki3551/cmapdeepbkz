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

/**@file    GaussSieve.cpp
 * @brief   Base class for GaussSieve.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "GaussSieve.h"

#include <stdlib.h>
#include <cassert>
#include <cmath>
#include <fstream>
#include <new>
#include <utility>
#include "Def.h"
#include "Config.h"
#include "KleinSampler.h"
#include "Lattice.h"
#include "Log.h"
#include "Timer.h"
#include "VectorElement.h"
#include "VectorElementPool.h"


namespace LapTools
{


///
/// @brief Gauss Sieve Algorithm
/// @param[in] c           max collision
/// @param[out] shortestV  vector
/// @param[in] timeLimit   time limit
///
template<typename BasisFloat, typename GSFloat>
bool
GaussSieve<BasisFloat, GSFloat>::gaussSieve(
      LatticeVector<BasisFloat> &shortestV,
      double inTimeLimit
      )
{
   assert( hasInitialized );
   if( inTimeLimit > 0 ) timeLimit = inTimeLimit;
   if( shortestV.size() > 0 )
   {
      bestVector = shortestV;
      bestObjectiveValue = bestVector.squaredNorm();
   }
   nLoop = 0;
   nCollision = 0;
   bool success = false;
   VectorElementPtr v = nullptr;
   bool simpleReduce = false;
   bool shouldAbort = false;

   startTime = Timer::getElapsedTime();
   runningTime = 0.0;
   double nextLogTime = startTime + config.LogIntervalGaussSieve;


   // append basis-vector into List and Stack
   LatticeVector<BasisFloat> basisV = L->basis.row(0);
   List.insert(VectorElementPtr(new VectorElementType(basisV)));
   for( int i = 1; i < L->m; ++i )
   {
      basisV = L->basis.row(i);
      Stack.insert(VectorElementPtr(new VectorElementType(basisV)));
   }

   if( verbose > 1 ){ std::cout << Log::getLogHeader() << std::endl; }

   while( maxCollision < 0 || nCollision < maxCollision )
   {
      ++nLoop;
      runningTime = Timer::getElapsedTime() - startTime;
      // abort check
      if( timeLimit > 0 && runningTime > timeLimit ){ outputLog('T'); break; }
      if( nextLogTime < runningTime )
      {
         outputLog(' ');
         nextLogTime = runningTime + config.LogIntervalGaussSieve;
      }

      // communicate with LC
      communicate(shouldAbort);
      if( shouldAbort ){ break; }

      // sampling
      if( !Stack.isEmpty() )
      {
         v = Stack.extract();
         if( verbose > 3 )
            std::cout << "sampling frmo Stack #S: " << Stack.size() << ", #L: " << List.size() << ", v norm: " << v->norm() << ", v: " << v->vector.transpose() << std::endl;
      }
      else
      {
         v = VectorElementPtr(new VectorElementType(kleinSampler.sample(simpleReduce)));
         if( std::isnan(v->squaredNorm()) ){ continue; }
         if( verbose > 3 )
            std::cout << "sampling frmo Klein, #L: " << List.size() << ", v sqnorm: " << v->squaredNorm() << ", v: " << v->vector.transpose() << std::endl;
      }

      // reduce v
      if( gaussReduce(v) )
      {
         if( List.insert(v) )
         {
            if( verbose > 3 )
               std::cout << "insert #L: " << List.size() << ", norm: " << v->norm() << ", v: " << v->vector.transpose() << std::endl;
            if( v->squaredNorm() < bestObjectiveValue )
            {
               shortestV = v->vector;
               success = true;
               updateBestObjectiveValue('*', v);
            }
         }
         else
         {
            // v is not inserted into L, and v is freed when extract from L
            nCollision++;
         }
      }
      else
      {
         // v is becomes zero vector
         nCollision++;
      }
   }

   outputLog('E');
   return success;
}


///
/// @brief Gauss Reduce Algorithm
/// @details 1. vector p is reduced by v in List such that |v| <= |p|
///          2. vectors v in List are reduced by p such that |v| > |p|
/// @param[in] p lattice vector which will be reduced
/// @return false if p becomes 0-vector after reduce process else true
///
template<typename BasisFloat, typename GSFloat>
bool
GaussSieve<BasisFloat, GSFloat>::gaussReduce(
      VectorElementPtr p
      )
{
   VectorElementPtr v = nullptr;
   iterListBegin = List.getIterBegin();
   iterListEnd = List.getIterEnd();

   // reduce p by vectors in List
   bool reduced = true;
   while( reduced )
   {
      reduced = false;
      for( iterList = iterListBegin; iterList != iterListEnd; ++iterList )
      {
         v = *iterList;
         assert( v->vector.size() == L->n );
         if( v->squaredNorm() > p->squaredNorm() ) break;  // |v| > |p|
         if( reduce(p, v) )
         {
            if( verbose > 3 )
               std::cout << "reduce p, #L: " << List.size() << ", p norm: " << p->norm() << ", p: " << p->vector.transpose() << std::endl;
            reduced = true;
            break;
         }
      }
   }

   // reduceing 0-vector
   if( p->squaredNorm() == 0 ) return false;
   postProcessOfReduce(p);

   // move v from L to S such that |v| > |p| and |v-p| <= |v|
   while( iterList != iterListEnd )
   {
      v = *iterList;
      if( reduce(v, p) )
      {
         if( verbose > 3 )
            std::cout << "reduce v and remove from L, #L: " << List.size() << ", v norm: " << v->norm() << ", v: " << v->vector.transpose() << std::endl;
         List.erase(iterList);
         Stack.insert(v);
      }
      else
      {
         iterList++;
      }
   }

   return true;
}


///
/// @brief reduce v by w such that |w| < |v|
/// @param[in] v lattice vector which will be reduced
/// @param[in] w lattice vector
/// @return true if v is reduced by w else false
///
template<typename BasisFloat, typename GSFloat>
bool
GaussSieve<BasisFloat, GSFloat>::reduce(
      VectorElementPtr v,
      VectorElementPtr w
      )
{
   // |v-w| <= |v| eq |v-w|^2 <= |v|^2 eq |w|^2 <= 2 * <v, w>
   // |v-q*w| <= |v| for any q  eq  round( <v, w> / |w|^2 )
   int dot = v->vector.dot(w->vector);
   if( std::abs(dot) < w->squaredNorm() * config.eta ) return false;
   // int q = std::round(dot / w->squaredNorm());
   v->vector -= static_cast<int>(std::round(dot / w->squaredNorm())) * w->vector;
   v->setSquaredNorm();
   if( v->vector.coeff(0) < 0 ) v->vector *= -1;
   return true;
}


///
/// @brief update bestObjectiveValue
/// @param[in] sigh line-header character
/// @param[in] shortest vector element
///
template<typename BasisFloat, typename GSFloat>
bool
GaussSieve<BasisFloat, GSFloat>::updateBestObjectiveValue(
      char sigh,
      VectorElementPtr v
      )
{
   if( v->squaredNorm() < bestObjectiveValue - SVP_MINEPSILON )
   {
      bestObjectiveValue = v->squaredNorm();
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
GaussSieve<BasisFloat, GSFloat>::outputLog(
      char sigh
      )
{
   int vthreas = 1;
   if( verbose >= vthreas && ( !config.Quiet || osCsvLog ) )
   {
      std::string taskName      = "Sieve";
      int         size          = L->m;
      long int    iter          = nLoop;
      double      progress      = ( static_cast<double>(nCollision) / maxCollision ) * 100;
      double      shortestNorm  = std::sqrt( bestObjectiveValue );
      double      logCost       = std::log( List.size() );
      std::ostringstream s_appendix;
      s_appendix
         << "  --  "
         << "|L| " << List.size()
         << ", |S| " << Stack.size()
         << ", minS " << ( Stack.size() > 0 ? (*Stack.getIterBegin())->norm() : -1 );
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
                  *L,sigh,rank,threadId,taskName,runningTime,
                  size,iter,progress,logCost,shortestNorm
                  )
            << "," << List.size()
            << "," << Stack.size()
            << std::endl;  // appendix
      }
   }
}


///
/// @brief get header of csv formed log
///
template<typename BasisFloat, typename GSFloat>
std::string
GaussSieve<BasisFloat, GSFloat>::getCsvLogHeader(
      )
{
   std::ostringstream s;
   s  << Log::getCsvLogHeader()
      << ",ListSize"
      << ",StackSize";
   return s.str();
}


///
/// instantiation
///
template class GaussSieve<int, double>;
template class GaussSieve<int, long double>;
template class GaussSieve<long int, double>;
template class GaussSieve<long int, long double>;


}  // namespace LapTools
