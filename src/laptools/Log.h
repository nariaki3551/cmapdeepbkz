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

/**@file    Log.h
 * @brief   Utilities for logging.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __LAPTOOLS_LOG_H__
#define __LAPTOOLS_LOG_H__

#include <iostream>
#include <iomanip>

#include "Lattice.h"
#include "Timer.h"


namespace LapTools
{


namespace Log
{


///
/// @brief convert time to be readable
/// @param[in, out] time
/// @param[out] unit
///
inline
void
convertTimeReadable(
      double &time,
      std::string &unit
      )
{
   if( time > 60*60*24*365 )
   {
      time /= 60*60*24*365.0;
      unit = "y";
   }
   else if( time > 60*60*24*10 )
   {
      time /= 60*60*24.0;
      unit = "d";
   }
   else if( time > 60*60*10 )
   {
      time /= 60*60.0;
      unit = "h";
   }
   else if( time > 9999 )
   {
      time /= 60.0;
      unit = "m";
   }
   else
   {
      unit = "s";
   }
}



///
/// @brief header of displayed log
/// @return header string
///
inline
std::string
getLogHeader()
{
   std::ostringstream s;
   s << "                                                  left    log                       " << std::endl;
   s << "   time rank  th     task size     iter    prog   time   cost     best  a.fac  h.fac";
   return s.str();
}


///
/// @brief log for display
/// @param[in] L              lattice
/// @param[in] sigh           header of log-line
/// @param[in] rank           solver-id
/// @param[in] threadId       solver-id
/// @param[in] taskName       algorithm name
/// @param[in] elapsedTime    elapsed time of task
/// @param[in] size           parameter of task
/// @param[in] iter           iterator of task
/// @param[in] progress       progress of task
/// @param[in] cost           logarithm metric of task cost using enum cost
/// @param[in] shortestNorm   shortestNorm found
/// @param[in] appendix       appendix info of task
///
template<typename BasisFloat, typename GSFloat>
inline
std::string
getLog(
      Lattice<BasisFloat, GSFloat> &L,
      char sigh,
      int rank,
      int threadId,
      std::string taskName,
      double elapsedTime,
      int size,
      long int iter,
      double progress,
      double cost,
      double shortestNorm,
      std::string appendix
      )
{

   // factors
   double approxFactor = L.approxFactor(shortestNorm);
   double hermiteFactor = L.hermiteFactor(shortestNorm);

   // time
   double leftTime = -1;
   if( cost != -1 )
   {
      // estimated time for compleation of enum search
      leftTime = std::exp(cost - std::log(L.config.EnumTraverseNodesPerSecond));
   }
   std::string elapsedTimeUnit, leftTimeUnit;
   convertTimeReadable(elapsedTime, elapsedTimeUnit);
   convertTimeReadable(leftTime, leftTimeUnit);

   std::ostringstream s_left_time;
   if( leftTime >= 1000 )
   {
      s_left_time << std::fixed << std::scientific << std::setprecision(0) << leftTime << leftTimeUnit;
   }
   else if( leftTime != -1 )
   {
      if( std::isnan(leftTime) )
         s_left_time << "inf";
      else
         s_left_time << static_cast<int>(leftTime) << leftTimeUnit;
   }

   // iteration
   std::ostringstream s_iter;
   if( iter >= 1e6 )
      s_iter << std::fixed << std::scientific << std::setprecision(1) << static_cast<long double>(iter);
   else
      s_iter << iter;

   // progress
   std::ostringstream s_prog;
   if( progress >= 999 )
      s_prog << std::fixed << std::scientific << std::setprecision(0) << progress << "%";
   else if( progress >= 0 )
      s_prog << std::fixed << std::setprecision(2) << progress << "%";

   // cost
   std::ostringstream s_cost;
   if( cost > 1e6 )
      s_cost << std::fixed << std::scientific << std::setprecision(2) << cost;
   else if( cost != -1 )
      s_cost << std::fixed << std::setprecision(2) << cost;

   std::ostringstream s;
   s << sigh;
   s << std::setw(5)  << int(elapsedTime) << elapsedTimeUnit;
   s << std::setw(5)  << rank;
   s << std::setw(4)  << threadId;
   s << std::setw(9)  << taskName;
   s << std::setw(5)  << size;
   s << std::setw(9)  << s_iter.str();
   s << std::setw(8)  << s_prog.str();
   s << std::setw(7)  << s_left_time.str();
   s << std::setw(7)  << s_cost.str();
   s << std::setw(9)  << std::fixed << std::setprecision(2) << shortestNorm;
   s << std::setw(7)  << std::fixed << std::setprecision(3) << approxFactor;
   s << std::setw(7)  << std::fixed << std::setprecision(3) << hermiteFactor;
   // s << std::setw(8)  << std::fixed << std::setprecision(2) << orthogonalFactor;
   s << appendix;
   return s.str();
}


///
/// @brief header of csv form
/// @return header string
///
inline
std::string
getCsvLogHeader()
{
   std::ostringstream s;
   s  << "time"
      << ",sigh"
      << ",taskTime"
      << ",rank"
      << ",threadId"
      << ",taskName"
      << ",size"
      << ",iter"
      << ",progress"
      << ",leftTime"
      << ",logCost"
      << ",shortestNorm"
      << ",approxFactor"
      << ",hermiteFactor"
      << ",orthogonalFactor";
   return s.str();
}


///
/// @brief log with csv form
/// @param[in] L              lattice
/// @param[in] sigh           header of log-line
/// @param[in] rank           solver-id
/// @param[in] threadId       solver-id
/// @param[in] taskName       algorithm name
/// @param[in] elapsedTime    elapsed time of task
/// @param[in] size           parameter of task
/// @param[in] iter           iterator of task
/// @param[in] progress       progress of task
/// @param[in] cost           logarithm metric of task cost using enum cost
/// @param[in] shortestNorm   shortestNorm found
/// @param[in] appendix       appendix info of task
///
template<typename BasisFloat, typename GSFloat>
inline
std::string
getCsvLog(
      Lattice<BasisFloat, GSFloat> &L,
      char sigh,
      int rank,
      int threadId,
      std::string taskName,
      double elapsedTime,
      int size,
      long int iter,
      double progress,
      double cost,
      double shortestNorm
      )
{
   // global time
   double time = Timer::getElapsedTime();

   // factors
   double approxFactor = L.approxFactor(shortestNorm);
   double hermiteFactor = L.hermiteFactor(shortestNorm);
   double orthogonalFactor = L.logOrthogonalityDefect();

   // time
   long double leftTime = -1;
   if( cost != -1 )
   {
      leftTime = std::exp(cost - std::log(L.config.EnumTraverseNodesPerSecond));
   }

   // sigh
   if( sigh == ' ') sigh = '-';

   std::ostringstream csvLog;
   csvLog << time          << ",";
   csvLog << sigh          << ",";
   csvLog << elapsedTime   << ",";
   csvLog << rank          << ",";
   csvLog << threadId      << ",";
   csvLog << taskName      << ",";
   csvLog << size          << ",";
   csvLog << iter          << ",";
   csvLog << progress      << ",";
   csvLog << leftTime      << ",";
   csvLog << cost          << ",";
   csvLog << shortestNorm  << ",";
   csvLog << approxFactor  << ",";
   csvLog << hermiteFactor << ",";
   csvLog << orthogonalFactor;
   return csvLog.str();
}


}  // namespace Log


}  // namespace LapTools


#endif   // __LAPTOOLS_LOG_H__
