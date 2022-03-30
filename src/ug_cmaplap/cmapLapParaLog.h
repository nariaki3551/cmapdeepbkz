/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*          This file is part of the program and software framework          */
/*  CMAP-LAP --- Configurable Massively Parallel Solver for Lattice Problems */
/*                                                                           */
/*  Copyright Written by Nariaki Tateiwa <n-tateiwa@kyudai.jp>,              */
/*                       Yuji Shinano <shinano@zib.de>,                      */
/*            Copyright (C) 2021 by Zuse Institute Berlin,                   */
/*            licensed under LGPL version 3 or later.                        */
/*            Commercial licenses are available through <licenses@zib.de>    */
/*                                                                           */
/* This code is free software; you can redistribute it and/or                */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>.     */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    cmapLapParaLog.h
 * @brief   Log tools
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_LOG_H__
#define __CMAP_LAP_PARA_LOG_H__

#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "cmapLapParaDef.h"

namespace ParaCMapLAP
{


namespace Logging
{


///
/// @return string of readable time
///
inline std::string
getReadableTimeString(
      double elapsedTime
      )
{
   std::ostringstream s;
   s << std::fixed << std::setprecision(0);
   if( elapsedTime > 60*60*24*365*999.0 )
   {
      s << std::scientific << elapsedTime / (60*60*24*365) << "y";
   }
   else if( elapsedTime > 60*60*24*999.0 )
   {
      s << elapsedTime / (60*60*24*365) << "y";
   }
   else if( elapsedTime > 60*60*999 )
   {
      s << elapsedTime / (60*60*24) << "d";
   }
   else if( elapsedTime > 60*999 )
   {
      s << elapsedTime / (60*60) << "h";
   }
   else if( elapsedTime > 9999 )
   {
      s << elapsedTime / 60 << "m";
   }
   else if( elapsedTime > 0 )
   {
      s << elapsedTime << "s";
   }
   return s.str();
}


///
/// @return string of solverType
///
inline std::string
getSolverTypeString(
      SolverType solverType
      )
{
   switch( solverType )
   {
      case DeepBkz:
         return "DeepBkz";
      case Enum:
         return "Enum";
      case Sieve:
         return "Sieve";
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaLog::getSolverTypeString: Invalid solver type = ", static_cast<int>(solverType));
   }
}


///
/// @return header of string to assignment table data
///
inline std::string
getCMapLapParaAssignmentTableHeader(
      )
{
   std::ostringstream s;
   s << "                                                  # of Active  # of Active  # of Active" << std::endl;
   s << "         Shortest Approx. # of Idle  # of Active      DeepBkz         Enum        Sieve       # of  # of Basis # of Vec." << std::endl;
   s << "   Time      Norm  Factor   Solvers      Solvers      Solvers      Solvers      Solvers      Tasks     in Pool   in Pool" << std::endl;
   return s.str();
}


///
/// @return header of string to assignment table data
///
inline std::string
getCMapLapParaCsvAssignmentTableHeader(
      )
{
   return "head,time,shortestNorm,approxFactor,hermiteFactor,rootHermiteFactor"
      ",numberOfIdleSolvers,numberOfActiveSolvers"
      ",numberOfActiveBkzSolvers,numberOfActiveEnumSolvers,numberOfActiveSieveSolvers"
      ",numberOfBkzTask,numberOfEnumTask,numberOfSieveTask"
      ",sizeOfInstancePool,sizeOfShareDataPool"
      "\n";
}


///
/// stringfy checkpoint data
/// @return string to checkpoint data
///
inline std::string
getCsvAssignmentTableString(
      char header,
      double time,
      double shortestNorm,
      double approxFactor,
      double hermiteFactor,
      double rootHermiteFactor,
      int numberOfIdleSolvers,
      int numberOfActiveSolvers,
      int numberOfActiveBkzSolvers,
      int numberOfActiveEnumSolvers,
      int numberOfActiveSieveSolvers,
      int numberOfBkzTask,
      int numberOfEnumTask,
      int numberOfSieveTask,
      int sizeOfInstancePool,
      int sizeOfShareDataPool,
      std::string delimiter=","
      )
{
   std::ostringstream s;
   s << header;
   s << delimiter; s << time;
   s << delimiter; s << shortestNorm;
   s << delimiter; s << approxFactor;
   s << delimiter; s << hermiteFactor;
   s << delimiter; s << rootHermiteFactor;
   s << delimiter; s << numberOfIdleSolvers;
   s << delimiter; s << numberOfActiveSolvers;
   s << delimiter; s << numberOfActiveBkzSolvers;
   s << delimiter; s << numberOfActiveEnumSolvers;
   s << delimiter; s << numberOfActiveSieveSolvers;
   s << delimiter; s << numberOfBkzTask;
   s << delimiter; s << numberOfEnumTask;
   s << delimiter; s << numberOfSieveTask;
   s << delimiter; s << sizeOfInstancePool;
   s << delimiter; s << sizeOfShareDataPool;
   return s.str();
}


///
/// @return header of statisticsFinalRun for csv
///
inline std::string
getCMapLapParaCsvStatisticsFinalRunHeader(
      )
{
   return "time,rank,threadId,interrupted,elapsedTime,totalComputingTime"
      ",totalIdleTime,idleTimeToFirstParaTask,idleTimeBetweenParaTasks,idleTimeAfterLastParaTask"
      ",idleTimeToWaitNotificationId,idleTimeToWaitAckCompletion,idleTimeToWaitToken"
      ",idleTimeToWaitSolverState,idleTimeToWaitPackedVector,idleTimeToWaitSolution,idleTimeToWaitBasis"
      ",idleTimeToWaitIsend"
      ",nParaTasksReceived,nParaTasksSolved"
      ",nParaTasksDeepBkzReceived,nParaTasksEnumReceived,nParaTasksSieveReceived"
      ",runningTimeDeepBkz,runningTimeEnum,runningTimeSieve"
      ",nVectorsReceivedDeepBkz,nVectorsReceivedEnum,nVectorsReceivedSieve"
      ",nVectorsSentDeepBkz,nVectorsSentEnum,nVectorsSentSieve"
      ",nBasesSentDeepBkz,nSolverStateSent"
      "\n";
}


///
/// @return header of string to solver state table data
///
inline std::string
getCMapLapParaCsvSolverStateHeader(
      std::string delimiter=","
      )
{
   std::ostringstream s;
   s << "sigh";
   s << delimiter << "time";
   s << delimiter << "rank";
   s << delimiter << "thread";
   s << delimiter << "taskTime";
   s << delimiter << "task";
   s << delimiter << "size";
   s << delimiter << "iter";
   s << delimiter << "prog";
   s << delimiter << "leftTime";
   s << delimiter << "logCost";
   s << delimiter << "shortestNorm";
   s << delimiter << "approxFactor";
   s << delimiter << "hermiteFactor";
   s << delimiter << "logOrthogonalFactor";
   s << delimiter << "meanMessageQueueSize";
   s << delimiter << "maxMessageQueueSize";
   s << std::endl;
   return s.str();
}


///
/// stringfy solver state data as csv format
/// @return string to solver state table data as csv format
///
inline std::string
toStringCsvLogBase(
      std::string taskName,
      double elapsedTime,
      int size,
      long int iter,
      double progress,
      double leftTime,
      double cost,
      double shortestNorm,
      double approxFactor,
      double hermiteFactor,
      double rootHermiteFactor,
      double orthogonalFactor,
      int meanMessageQueueSize,
      int maxMessageQueueSize,
      std::string delimiter=","
      )
{
   std::ostringstream s;
   if(elapsedTime > 0) s << elapsedTime;
   else                s << '-';
   s << delimiter << taskName;
   s << delimiter; (size                  >= 0) ? s << size                   : s << '-';
   s << delimiter; (iter                  >= 0) ? s << iter                   : s << '-';
   s << delimiter; (progress              >= 0) ? s << progress               : s << '-';
   s << delimiter; (leftTime              >= 0) ? s << leftTime               : s << '-';
   s << delimiter; (cost                  >= 0) ? s << cost                   : s << '-';
   s << delimiter; (shortestNorm          >= 0) ? s << shortestNorm           : s << '-';
   s << delimiter; (approxFactor          >= 0) ? s << approxFactor           : s << '-';
   s << delimiter; (hermiteFactor         >= 0) ? s << hermiteFactor          : s << '-';
   s << delimiter; (orthogonalFactor      >= 0) ? s << orthogonalFactor       : s << '-';
   s << delimiter; (meanMessageQueueSize  >= 0) ? s << meanMessageQueueSize   : s << "-";
   s << delimiter; (maxMessageQueueSize   >= 0) ? s << maxMessageQueueSize    : s << "-";
   return s.str();
}


///
/// @return header of string to solver state table data
///
inline std::string
getCMapLapParaSolverStateHeader(
      )
{
   std::ostringstream s;
   s  << "                    task                                   left    log                       " << std::endl;
   s  << "  time    rank  th  time     task size     iter    prog    time   cost     best  a.fac  h.fac" << std::endl;
   return s.str();
}


///
/// stringfy solver state data
/// @return string to solver state table data
///
inline std::string
toStringLogBase(
      std::string taskName,
      double      elapsedTime,
      int         size,
      long int    iter,
      double      progress,
      double      leftTime,
      double      cost,
      double      shortestNorm,
      double      approxFactor,
      double      hermiteFactor,
      double      rootHermiteFactor,
      double      orthogonalFactor,
      int         meanMessageQueueSize,
      int         maxMessageQueueSize,
      std::string appendix="",
      std::string delimiter=""
      )
{
   if( delimiter != "" )
   {
      return toStringCsvLogBase(taskName,elapsedTime,size,iter,progress,leftTime,
            cost,shortestNorm,approxFactor,hermiteFactor,rootHermiteFactor,orthogonalFactor,
            meanMessageQueueSize,maxMessageQueueSize,
            delimiter);
   }

   // iteration
   std::ostringstream s_iter;
   if( iter > static_cast<long int>(1e7) )
   {
      s_iter << std::fixed << std::scientific << std::setprecision(2) << static_cast<double>(iter);
   }
   else
   {
      s_iter << iter;
   }

   // cost
   std::ostringstream s_cost;
   if( cost > 1e6 )
   {
      s_cost << std::fixed << std::scientific << std::setprecision(2) << cost;
   }
   else if( cost > 0)
   {
      s_cost << std::fixed << std::setprecision(1) << cost;
   }

   // shortest norm
   std::ostringstream s_shortestNorm;
   if( shortestNorm > 0 )
   {
      s_shortestNorm << std::fixed << std::setprecision(2) << shortestNorm;
   }

   // approx factor
   std::ostringstream s_approxFactor;
   if( approxFactor > 0 )
   {
      s_approxFactor << std::fixed << std::setprecision(3) << approxFactor;
   }

   // hermite factor
   std::ostringstream s_hermiteFactor;
   if( hermiteFactor > 0 )
   {
      s_hermiteFactor << std::fixed << std::setprecision(3) << hermiteFactor;
   }

   std::ostringstream s;
   s << std::setw(6)  << getReadableTimeString(elapsedTime);
   s << std::setw(9)  << taskName;
   s << std::setw(5)  << size;
   s << std::setw(9)  << s_iter.str();
   s << std::setw(7)  << std::fixed << std::setprecision(2) << progress << "%";
   s << std::setw(8)  << getReadableTimeString(leftTime);
   s << std::setw(7)  << s_cost.str();
   s << std::setw(9)  << s_shortestNorm.str();
   s << std::setw(7)  << s_approxFactor.str();
   s << std::setw(7)  << s_hermiteFactor.str();
   s << appendix;
   return s.str();
}


///
/// stringfy solver state data
/// @return string to solver state table data
///
inline std::string
getSolverStateString(
      char sigh,
      double elapsedTime,
      int source,
      int threadId,
      SolverType solverType,
      double shortestNorm,
      std::string delimiter=""
   )
{
   std::string taskName             = Logging::getSolverTypeString(solverType);
   double      elapsedTaskTime      = -1;
   int         size                 = -1;
   long int    iter                 = -1;
   double      progress             = -1;
   double      leftTime             = -1;
   double      cost                 = -1;
   double      approxFactor         = -1;
   double      hermiteFactor        = -1;
   double      rootHermiteFactor    = -1;
   double      orthogonalFactor     = -1;
   int         meanMessageQueueSize = -1;
   int         maxMessageQueueSize  = -1;

   std::ostringstream s;
   std::string appendix = "";
   s << sigh;
   if( delimiter == "" )
   {
      s << std::setw(5)  << getReadableTimeString(elapsedTime);
      s << std::setw(8)  << source;
      s << std::setw(4)  << threadId;
      s << toStringLogBase(taskName,elapsedTaskTime,size,iter,progress,leftTime,
            cost,shortestNorm,approxFactor,hermiteFactor,rootHermiteFactor,orthogonalFactor,
            meanMessageQueueSize,maxMessageQueueSize,
            appendix);
   }
   else
   {
      s << delimiter << elapsedTime;
      s << delimiter << source;
      s << delimiter << threadId;
      s << delimiter << toStringCsvLogBase(taskName,elapsedTaskTime,size,iter,progress,leftTime,
            cost,shortestNorm,approxFactor,hermiteFactor,rootHermiteFactor,orthogonalFactor,
            meanMessageQueueSize,maxMessageQueueSize,
            delimiter);
   }
   return s.str();
}


///
/// stringfy solver state data
/// @return string to solver state table data
///
inline std::string
getSolverStateString(
      char sigh,
      double elapsedTime,
      int source,
      int threadId,
      std::string stateString,
      std::string delimiter=""
   )
{
   std::ostringstream s;
   s << sigh;
   if( delimiter == "" )
   {
      s << std::setw(5) << getReadableTimeString(elapsedTime);
      s << std::setw(8) << source;
      s << std::setw(4) << threadId;
      s << stateString;
   }
   else
   {
      s << delimiter << elapsedTime;
      s << delimiter << source;
      s << delimiter << threadId;
      s << delimiter << stateString;
   }
   return s.str();
}


///
/// @return header of string to checkpoint data
///
inline std::string
getCMapLapParaCsvCheckpointHeader(
      std::string delimiter=","
      )
{
   std::ostringstream s;
   s << "sigh";
   s << delimiter << "elapsedTime";
   s << delimiter << "numberOfActiveTasks";
   s << delimiter << "numberOfInactiveEnum";
   s << delimiter << "numberOfInactiveSieve";
   s << delimiter << "numberOfBasis";
   s << delimiter << "numberOfVector";
   s << delimiter << "checkpointTime";
   s << delimiter << "totalChechpointTime";
   s << delimiter << "idleTimeToPrube";
   s << delimiter << "idleTimeToWaitIsend";
   s << std::endl;
   return s.str();
}


///
/// stringfy checkpoint data
/// @return string to checkpoint data
///
inline std::string
getCsvCheckpointStateString(
      double elapsedTime,
      double checkpointTime,
      double totalCheckpointTime,
      double idleTime,
      double idleTimeToWaitIsend,
      int nActive,
      int nInactiveEnum,
      int nInactiveSieve,
      int nBasis,
      int nVector,
      std::string delimiter,
      char header='W'
      )
{
   std::ostringstream s;
   s << header;
   s << delimiter; (elapsedTime         >= 0) ? s << elapsedTime         : s << '-';
   s << delimiter; (nActive             >= 0) ? s << nActive             : s << '-';
   s << delimiter; (nInactiveEnum       >= 0) ? s << nInactiveEnum       : s << '-';
   s << delimiter; (nInactiveSieve      >= 0) ? s << nInactiveSieve      : s << '-';
   s << delimiter; (nBasis              >= 0) ? s << nBasis              : s << '-';
   s << delimiter; (nVector             >= 0) ? s << nVector             : s << '-';
   s << delimiter; (checkpointTime      >= 0) ? s << checkpointTime      : s << '-';
   s << delimiter; (totalCheckpointTime >= 0) ? s << totalCheckpointTime : s << "-";
   s << delimiter; (idleTime            >= 0) ? s << idleTime            : s << "-";
   s << delimiter; (idleTimeToWaitIsend >= 0) ? s << idleTimeToWaitIsend : s << "-";
   return s.str();
}


///
/// stringfy checkpoint data
/// @return string to checkpoint data
///
inline std::string
getCheckpointStateString(
      double elapsedTime,
      double checkpointTime,
      int nActive,
      int nInactiveEnum,
      int nInactiveSieve,
      int nBasis,
      int nVector,
      char header='W'
      )
{
   std::ostringstream s;
   s << std::endl;
   s << header;
   s << std::setw(5) << getReadableTimeString(elapsedTime);
   s << "   save ";
   s << (nActive+nInactiveEnum+nInactiveSieve) << " Tasks( ";
   s << nActive << " active, ";
   s << nInactiveEnum  << " inactive Enum, ";
   s << nInactiveSieve << " inactive Sieve), ";
   s << nBasis << " basis, ";
   s << nVector << " vectors,";
   s << " time " << std::fixed << std::setprecision(2) << checkpointTime << "s";
   s << std::endl << std::endl;
   return s.str();
}


///
/// @return header of string of osLogNotificationInterval
///
inline std::string
getNotificationIntervalLogHeader(
      )
{
   return "time,idleRatio,lowerIdleRatio,upperIdleRatio"
      ",notificationInterval,leftNotificationInterval";
}

///
/// @return stringfy of osLogNotificationInterval
///
inline std::string
getNotificationIntervalLog(
      double elapsedTime,
      double idleRatio,
      double lowerIdleRatio,
      double upperIdleRatio,
      double notificationInterval,
      double leftNotificationInterval,
      std::string delimiter=","
      )
{
   std::ostringstream s;
   s << elapsedTime;
   s << delimiter << idleRatio;
   s << delimiter << lowerIdleRatio;
   s << delimiter << upperIdleRatio;
   s << delimiter << notificationInterval;
   s << delimiter << leftNotificationInterval;
   return s.str();
}


} // namespace Logging


} // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_LOG_H__
