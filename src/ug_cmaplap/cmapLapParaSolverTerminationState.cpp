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

/**@file    cmapLapParaSolverTerminationState.cpp
 * @brief   This class contains solver termination state which is transferred form Solver to LC.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "cmapLapParaSolverTerminationState.h"
#include <ostream>
#include "ug/paraDef.h"
#ifdef UG_WITH_ZLIB
#include "ug/gzstream.h"
#endif
namespace UG { class ParaComm; }
namespace UG { class ParaInitiator; }


namespace ParaCMapLAP
{

///
/// stringfy CMapLapParaSolverTerminationState
///
std::string
CMapLapParaSolverTerminationState::toString(
      UG::ParaInitiator *initiator
      )
{
   std::ostringstream os;
   switch( interrupted )
   {
   case 0:
   {
      os << "######### Solver Rank = " << rank << " Thread = " << threadId << " is terminated. #########" << std::endl;
      break;
   }
   case 1:
   {
      os << "######### Solver Rank = " << rank << " Thread = " << threadId << "is interrupted. #########" << std::endl;
      break;
   }
   case 2:
   {
      os << "######### Solver Rank = " << rank << " Thread = " << threadId << "is at checkpoint. #########" << std::endl;
      break;
   }
   case 3:
   {
      os << "######### Solver Rank = " << rank << " Thread = " << threadId << "is at the end of racing process. #########" << std::endl;
      break;
   }
   default:
   {
      THROW_LOGICAL_ERROR1("invalid interrupted flag in CMapLapParaSolverTerminationState!");
   }
   }

   double totalIdleTime = idleTimeToFirstParaTask
                        + idleTimeBetweenParaTasks
                        + idleTimeAfterLastParaTask
                        + idleTimeToWaitNotificationId
                        + idleTimeToWaitAckCompletion
                        + idleTimeToWaitToken
                        + idleTimeToWaitSolverState
                        + idleTimeToWaitPackedVector
                        + idleTimeToWaitSolution
                        + idleTimeToWaitBasis
                        + idleTimeToWaitIsend;
   os << "#=== Elapsed time to terminate this Solver = " << runningTime << std::endl;
   os << "#=== Total computing time = " << runningTime - totalIdleTime << std::endl;
   os << "#=== Total idle time = " << totalIdleTime << std::endl;
   os << "#=== ( Idle time to start first Task = " << idleTimeToFirstParaTask
   << ", Idle time between Tasks = " << idleTimeBetweenParaTasks
   << ", Idle Time after last Task = " << idleTimeAfterLastParaTask
   << " )" << std::endl;
   os << "#=== ( Idle time to wait notification Id messages = " << idleTimeToWaitNotificationId << " )" << std::endl;
   os << "#=== ( Idle time to wait acknowledgment of completion = " << idleTimeToWaitAckCompletion << " )" << std::endl;
   os << "#=== ( Idle time to wait token = " << idleTimeToWaitToken << " )" << std::endl;
   os << "#=== ( Idle time to wait solver states = " << idleTimeToWaitSolverState << " )" << std::endl;
   os << "#=== ( Idle time to wait packed vectors = " << idleTimeToWaitPackedVector << " )" << std::endl;
   os << "#=== ( Idle time to wait solutions = " << idleTimeToWaitSolution << " )" << std::endl;
   os << "#=== ( Idle time to wait basis = " << idleTimeToWaitBasis << " )" << std::endl;
   os << "#=== ( Idle time to wait isend = " << idleTimeToWaitIsend << " )" << std::endl;
   os << "#=== The number of Tasks received in this solver = " << nParaTasksReceived << std::endl;
   os << "#=== The number of Tasks solved in this solver = " << nParaTasksSolved << std::endl;
   return os.str();
}

///
/// stringfy CMapLapParaSolverTerminationState
///
std::string
CMapLapParaSolverTerminationState::toCsvString(
      UG::ParaInitiator *initiator
      )
{
   std::ostringstream os;
   switch( interrupted )
   {
   case 0:
   {
      os << rank << "," << threadId << ",terminated";
      break;
   }
   case 1:
   {
      os << rank << "," << threadId << ",interrupted";
      break;
   }
   case 2:
   {
      os << rank << "," << threadId << ",checkpoint";
      break;
   }
   case 3:
   {
      os << rank << "," << threadId << ",end of racing process";
      break;
   }
   default:
   {
      THROW_LOGICAL_ERROR1("invalid interrupted flag in CMapLapParaSolverTerminationState!");
   }
   }

   double totalIdleTime = idleTimeToFirstParaTask
                        + idleTimeBetweenParaTasks
                        + idleTimeAfterLastParaTask
                        + idleTimeToWaitNotificationId
                        + idleTimeToWaitAckCompletion
                        + idleTimeToWaitToken
                        + idleTimeToWaitSolverState
                        + idleTimeToWaitPackedVector
                        + idleTimeToWaitSolution
                        + idleTimeToWaitBasis
                        + idleTimeToWaitIsend;
   os << "," << runningTime;                  // Elapsed time to terminate this Solver
   os << "," << runningTime - totalIdleTime;  // Total compting time
   os << "," << totalIdleTime;                // Total idle time
   os << "," << idleTimeToFirstParaTask;
   os << "," << idleTimeBetweenParaTasks;
   os << "," << idleTimeAfterLastParaTask;
   os << "," << idleTimeToWaitNotificationId;
   os << "," << idleTimeToWaitAckCompletion;
   os << "," << idleTimeToWaitToken;
   os << ',' << idleTimeToWaitSolverState;
   os << "," << idleTimeToWaitPackedVector;
   os << "," << idleTimeToWaitSolution;
   os << "," << idleTimeToWaitBasis;
   os << "," << idleTimeToWaitIsend;
   os << "," << nParaTasksReceived;
   os << "," << nParaTasksSolved;
   os << "," << nParaTasksDeepBkzReceived;
   os << "," << nParaTasksEnumReceived;
   os << "," << nParaTasksSieveReceived;
   os << "," << runningTimeDeepBkz;
   os << "," << runningTimeEnum;
   os << "," << runningTimeSieve;
   os << "," << nVectorsReceivedDeepBkz;
   os << "," << nVectorsReceivedEnum;
   os << "," << nVectorsReceivedSieve;
   os << "," << nVectorsSentDeepBkz;
   os << "," << nVectorsSentEnum;
   os << "," << nVectorsSentSieve;
   os << "," << nBasesSentDeepBkz;
   os << "," << nSolverStateSent;
   os << std::endl;
   return os.str();
}

#ifdef UG_WITH_ZLIB
void
CMapLapParaSolverTerminationState::write(
      gzstream::ogzstream &out
      )
{
   out.write(reinterpret_cast<char *>(&interrupted),                  sizeof(int));
   out.write(reinterpret_cast<char *>(&rank),                         sizeof(int));
   out.write(reinterpret_cast<char *>(&nParaTasksReceived),           sizeof(int));
   out.write(reinterpret_cast<char *>(&nParaTasksSolved),             sizeof(int));
   out.write(reinterpret_cast<char *>(&runningTime),                  sizeof(double));
   out.write(reinterpret_cast<char *>(&idleTimeToFirstParaTask),      sizeof(double));
   out.write(reinterpret_cast<char *>(&idleTimeBetweenParaTasks),     sizeof(double));
   out.write(reinterpret_cast<char *>(&idleTimeAfterLastParaTask),    sizeof(double));
   out.write(reinterpret_cast<char *>(&idleTimeToWaitNotificationId), sizeof(double));
   out.write(reinterpret_cast<char *>(&idleTimeToWaitToken),          sizeof(double));

   out.write(reinterpret_cast<char *>(&idleTimeToWaitSolverState),  sizeof(double));
   out.write(reinterpret_cast<char *>(&idleTimeToWaitPackedVector), sizeof(double));
   out.write(reinterpret_cast<char *>(&idleTimeToWaitSolution),     sizeof(double));
   out.write(reinterpret_cast<char *>(&idleTimeToWaitBasis),        sizeof(double));
   out.write(reinterpret_cast<char *>(&idleTimeToWaitIsend),        sizeof(double));
   out.write(reinterpret_cast<char *>(&nParaTasksDeepBkzReceived),  sizeof(int));
   out.write(reinterpret_cast<char *>(&nParaTasksEnumReceived),     sizeof(int));
   out.write(reinterpret_cast<char *>(&nParaTasksSieveReceived),    sizeof(int));
   out.write(reinterpret_cast<char *>(&runningTimeDeepBkz),         sizeof(double));
   out.write(reinterpret_cast<char *>(&runningTimeEnum),            sizeof(double));
   out.write(reinterpret_cast<char *>(&runningTimeSieve),           sizeof(double));
   out.write(reinterpret_cast<char *>(&nVectorsReceivedDeepBkz),    sizeof(int));
   out.write(reinterpret_cast<char *>(&nVectorsReceivedEnum),       sizeof(int));
   out.write(reinterpret_cast<char *>(&nVectorsReceivedSieve),      sizeof(int));
   out.write(reinterpret_cast<char *>(&nVectorsSentDeepBkz),        sizeof(int));
   out.write(reinterpret_cast<char *>(&nVectorsSentEnum),           sizeof(int));
   out.write(reinterpret_cast<char *>(&nVectorsSentSieve),          sizeof(int));
   out.write(reinterpret_cast<char *>(&nBasesSentDeepBkz),          sizeof(int));
   out.write(reinterpret_cast<char *>(&nSolverStateSent),           sizeof(int));
   // detTime and dualBound are not saved
}

bool
CMapLapParaSolverTerminationState::read(
      UG::ParaComm *comm,
      gzstream::igzstream &in
      )
{
   in.read(reinterpret_cast<char *>(&interrupted), sizeof(int));
   if( in.eof() ) return false;
   in.read(reinterpret_cast<char *>(&rank),                         sizeof(int));
   in.read(reinterpret_cast<char *>(&nParaTasksReceived),           sizeof(int));
   in.read(reinterpret_cast<char *>(&nParaTasksSolved),             sizeof(int));
   in.read(reinterpret_cast<char *>(&runningTime),                  sizeof(double));
   in.read(reinterpret_cast<char *>(&idleTimeToFirstParaTask),      sizeof(double));
   in.read(reinterpret_cast<char *>(&idleTimeBetweenParaTasks),     sizeof(double));
   in.read(reinterpret_cast<char *>(&idleTimeAfterLastParaTask),    sizeof(double));
   in.read(reinterpret_cast<char *>(&idleTimeToWaitNotificationId), sizeof(double));
   in.read(reinterpret_cast<char *>(&idleTimeToWaitToken),          sizeof(double));

   in.read(reinterpret_cast<char *>(&idleTimeToWaitSolverState),  sizeof(double));
   in.read(reinterpret_cast<char *>(&idleTimeToWaitPackedVector), sizeof(double));
   in.read(reinterpret_cast<char *>(&idleTimeToWaitSolution),     sizeof(double));
   in.read(reinterpret_cast<char *>(&idleTimeToWaitBasis),        sizeof(double));
   in.read(reinterpret_cast<char *>(&idleTimeToWaitIsend),        sizeof(double));
   in.read(reinterpret_cast<char *>(&nParaTasksDeepBkzReceived),  sizeof(int));
   in.read(reinterpret_cast<char *>(&nParaTasksEnumReceived),     sizeof(int));
   in.read(reinterpret_cast<char *>(&nParaTasksSieveReceived),    sizeof(int));
   in.read(reinterpret_cast<char *>(&runningTimeDeepBkz),         sizeof(double));
   in.read(reinterpret_cast<char *>(&runningTimeEnum),            sizeof(double));
   in.read(reinterpret_cast<char *>(&runningTimeSieve),           sizeof(double));
   in.read(reinterpret_cast<char *>(&nVectorsReceivedDeepBkz),    sizeof(int));
   in.read(reinterpret_cast<char *>(&nVectorsReceivedEnum),       sizeof(int));
   in.read(reinterpret_cast<char *>(&nVectorsReceivedSieve),      sizeof(int));
   in.read(reinterpret_cast<char *>(&nVectorsSentDeepBkz),        sizeof(int));
   in.read(reinterpret_cast<char *>(&nVectorsSentEnum),           sizeof(int));
   in.read(reinterpret_cast<char *>(&nVectorsSentSieve),          sizeof(int));
   in.read(reinterpret_cast<char *>(&nBasesSentDeepBkz),          sizeof(int));
   in.read(reinterpret_cast<char *>(&nSolverStateSent),           sizeof(int));

   // detTime and dualBound are not saved
   return true;
}

} // namespace ParaCMapLAP

#endif
