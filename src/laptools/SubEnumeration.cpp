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

/**@file    SubEnumeration.cpp
 * @brief   Base class for SubEnumeration.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "SubEnumeration.h"

#include <fstream>
#include <iomanip>
#include <cmath>

#include "Def.h"
#include "Config.h"
#include "Enumeration.h"
#include "Lattice.h"
#include "Log.h"
#include "Timer.h"


namespace LapTools
{


///
/// @brief sub-Enumeration
/// @param[out] resv          vector found by Enum
/// @parem[out] coeffResv     coefficient of vector found by Enum
/// @param[in] inTimeLimit    time limit
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
bool
SubEnumeration<BasisFloat, GSFloat, EnumGSFloat>::subEnum(
      LatticeVector<BasisFloat> &resv,
      LatticeVector<BasisFloat> &coeffResv,
      double inTimeLimit
      )
{
   if( inTimeLimit > 0 ) timeLimit = inTimeLimit;
   bool success = false;
   double nextLogTime = 0.0;
   double startTime = Timer::getElapsedTime();

   // process 1
   // Root Enum of {b_begin, ..., b_{begin+projDim}
   Enumeration<BasisFloat, GSFloat, EnumGSFloat> rootEnum{L, rank, threadId, 0};
   rootEnum.init(begin, begin+projDim);
   success = rootEnum.projectedEnum(resv, coeffResv);
   if( success ){ updateBestObjectiveValue('R', resv); }

   // process 2
   // child enum of {b_{begin+projDim+1}, ..., b_{end}}
   int childBegin = begin + projDim;
   LatticeVector<BasisFloat> childResv, childCoeffResv, rootResv, rootCoeffResv;
   auto childL = std::make_shared<Lattice<BasisFloat, GSFloat>>(L->basis);
   childL->setConfig(L->config);
   childL->GH = childL->projectedGH(childBegin, end);

   Enumeration<BasisFloat, GSFloat, EnumGSFloat> childEnum{childL, rank, threadId, 0};
   double childLowerbound = std::sqrt(childL->B(childBegin) * 1.01);
   childEnum.init(childBegin, end, childLowerbound);

   while( childEnum.projectedEnum(childResv, childCoeffResv, 1) )
   {
      if( verbose > 1 )
      {
         std::cout << "SubEnumeration::subEnum Rank " << rank << " Thread " << threadId << " Child-Enum " << nChildEnum
            << " proj(v) " << childEnum.getProjectedSquaradNorm(childBegin, end, childCoeffResv)
            << " < " << childEnum.getLattice()->B(childBegin)
            << " -- " << std::fixed << std::setprecision(3) << childEnum.getProgress() * 100 << " % "
            << std::endl;
      }
      nChildEnum++;
      // rootEnum.getLattice()->insert(childBegin, childResv)
      rootEnum.getLattice()->basis.row(childBegin) = childResv;
      rootEnum.getLattice()->setGSO(childBegin, childBegin);
      rootEnum.getLattice()->sizeReduce(childBegin, 0.51);
      rootEnum.init(begin, childBegin);
      bool rootSucccess = rootEnum.projectedEnum(rootResv, rootCoeffResv);
      if( rootSucccess ){ updateBestObjectiveValue('*', resv); }

      runningTime = Timer::getElapsedTime() - startTime;
      if( nextLogTime < runningTime )
      {
         outputLog(' ');
         nextLogTime = runningTime + this->config.LogIntervalSubEnum;
      }
      if( timeLimit > 0 && runningTime >= timeLimit )
      {
         outputLog('T');
         break;
      }
   }
   outputLog('E');
   return success;
}



///
/// @brief update bestObjectiveValue
/// @param[in] sigh line-header character
/// @param[in] v incumbent vector
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
void
SubEnumeration<BasisFloat, GSFloat, EnumGSFloat>::updateBestObjectiveValue(
      char sigh,
      LatticeVector<BasisFloat> &v
      )
{
   double squaredNorm = v.squaredNorm();
   if( squaredNorm < bestObjectiveValue )
   {
      bestObjectiveValue = squaredNorm;
      outputLog(sigh);
   }
}


///
/// @brief createt log
/// @param[in] sigh line-header character
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
void
SubEnumeration<BasisFloat, GSFloat, EnumGSFloat>::outputLog(
   char sigh
   )
{
   int vthreas = 1;
   if( verbose >= vthreas && ( !config.Quiet || osCsvLog ) )
   {
      std::string taskName       = "SubEnum";
      int         size           = end - begin + 1;
      long int    iter           = nChildEnum;
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
template class SubEnumeration<int, double, double>;
template class SubEnumeration<int, long double, double>;
template class SubEnumeration<int, long double, long double>;
template class SubEnumeration<long int, double, double>;
template class SubEnumeration<long int, long double, double>;
template class SubEnumeration<long int, long double, long double>;


}  // namespace LapTools
