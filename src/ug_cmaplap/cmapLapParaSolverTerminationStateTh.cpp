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

/**@file    cmapLapParaSolverTerminationStateTh.cpp
 * @brief   CMapLapParaSolverTerminationState extension for threads communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "cmapLapParaSolverTerminationStateTh.h"
#include "cmapLapParaComm.h"
#ifdef _COMM_MPI_WORLD
#include "cmapLapParaSolverLocalComm.h"
#endif
#include <memory>
#include "ug/paraComm.h"


namespace ParaCMapLAP
{

CMapLapParaSolverTerminationStateTh *
CMapLapParaSolverTerminationStateTh::createDatatype(
      )
{
   return new CMapLapParaSolverTerminationStateTh(
         threadId,
         interrupted,
         rank,
         nParaTasksReceived,
         nParaTasksSolved,
         runningTime,
         idleTimeToFirstParaTask,
         idleTimeBetweenParaTasks,
         idleTimeAfterLastParaTask,
         idleTimeToWaitNotificationId,
         idleTimeToWaitAckCompletion,
         idleTimeToWaitToken,
         idleTimeToWaitSolverState,
         idleTimeToWaitPackedVector,
         idleTimeToWaitSolution,
         idleTimeToWaitBasis,
         idleTimeToWaitIsend,
         detTime,
         nParaTasksDeepBkzReceived,
         nParaTasksEnumReceived,
         nParaTasksSieveReceived,
         runningTimeDeepBkz,
         runningTimeEnum,
         runningTimeSieve,
         nVectorsReceivedDeepBkz,
         nVectorsReceivedEnum,
         nVectorsReceivedSieve,
         nVectorsSentDeepBkz,
         nVectorsSentEnum,
         nVectorsSentSieve,
         nBasesSentDeepBkz,
         nSolverStateSent
         );

}

void
CMapLapParaSolverTerminationStateTh::send(
      UG::ParaComm *comm,
      int destination,
      int tag
      )
{

#ifdef _COMM_CPP11
   CMapLapParaCommTh *commTh = dynamic_cast<CMapLapParaCommTh*>(comm);
   PARA_COMM_CALL(
      commTh->uTypeSend((void *)createDatatype(), UG::ParaSolverTerminationStateType, destination, tag)
   );
#else
   CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
   PARA_COMM_CALL(
         localComm->uTypeSend((void *)createDatatype(), ParaSolverTerminationStateType, destination, tag)
   );
#endif
}

void
CMapLapParaSolverTerminationStateTh::receive(
      UG::ParaComm *comm,
      int source,
      int tag
      )
{

   std::unique_ptr<CMapLapParaSolverTerminationStateTh> received;

#ifdef _COMM_CPP11
   CMapLapParaCommTh *commTh = dynamic_cast<CMapLapParaCommTh*>(comm);
   PARA_COMM_CALL(
      commTh->uTypeReceive((void **)&received, UG::ParaSolverTerminationStateType, source, tag)
   );

#else
   CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
   PARA_COMM_CALL(
      localComm->uTypeReceive((void **)&received, ParaSolverTerminationStateType, source, tag)
   );

#endif

   threadId                      = received->threadId;
   interrupted                   = received->interrupted;
   rank                          = received->rank;
   nParaTasksReceived            = received->nParaTasksReceived;
   nParaTasksSolved              = received->nParaTasksSolved;
   runningTime                   = received->runningTime;
   idleTimeToFirstParaTask       = received->idleTimeToFirstParaTask;
   idleTimeBetweenParaTasks      = received->idleTimeBetweenParaTasks;
   idleTimeAfterLastParaTask     = received->idleTimeAfterLastParaTask;
   idleTimeToWaitNotificationId  = received->idleTimeToWaitNotificationId;
   idleTimeToWaitAckCompletion   = received->idleTimeToWaitAckCompletion;
   idleTimeToWaitToken           = received->idleTimeToWaitToken;
   idleTimeToWaitSolverState     = received->idleTimeToWaitSolverState;
   idleTimeToWaitPackedVector    = received->idleTimeToWaitPackedVector;
   idleTimeToWaitSolution        = received->idleTimeToWaitSolution;
   idleTimeToWaitBasis           = received->idleTimeToWaitBasis;
   idleTimeToWaitIsend           = received->idleTimeToWaitIsend;
   detTime                       = received->detTime;
   nParaTasksDeepBkzReceived     = received->nParaTasksDeepBkzReceived;
   nParaTasksEnumReceived        = received->nParaTasksEnumReceived;
   nParaTasksSieveReceived       = received->nParaTasksSieveReceived;
   runningTimeDeepBkz            = received->runningTimeDeepBkz;
   runningTimeEnum               = received->runningTimeEnum;
   runningTimeSieve              = received->runningTimeSieve;
   nVectorsReceivedDeepBkz       = received->nVectorsReceivedDeepBkz;
   nVectorsReceivedEnum          = received->nVectorsReceivedEnum;
   nVectorsReceivedSieve         = received->nVectorsReceivedSieve;
   nVectorsSentDeepBkz           = received->nVectorsSentDeepBkz;
   nVectorsSentEnum              = received->nVectorsSentEnum;
   nVectorsSentSieve             = received->nVectorsSentSieve;
   nBasesSentDeepBkz             = received->nBasesSentDeepBkz;
   nSolverStateSent              = received->nSolverStateSent;
}

} // namespace ParaCMapLAP
