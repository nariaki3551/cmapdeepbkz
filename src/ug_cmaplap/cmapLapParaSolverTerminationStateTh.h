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

/**@file    cmapLapParaSolverTerminationStateTh.h
 * @brief   CMapLapParaSolverTerminationState extension for threads communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_SOLVER_TERMINATION_STATE_TH_H__
#define __CMAP_LAP_PARA_SOLVER_TERMINATION_STATE_TH_H__


#include "cmapLapParaSolverTerminationState.h"
namespace UG { class ParaComm; }


namespace ParaCMapLAP
{

///
/// class CMapLapParaSolverTerminationStateMpi
/// (Solver termination state in a ParaSolver communicated by thread)
///
class CMapLapParaSolverTerminationStateTh : public CMapLapParaSolverTerminationState
{

   ///
   /// create CMapLapParaSolverTerminationStateTh datatype
   /// @return MPI_Datatype for CMapLapParaSolverTerminationStateTh
   ///
   CMapLapParaSolverTerminationStateTh *createDatatype(
         );

public:

   ///
   /// default constructor
   ///
   CMapLapParaSolverTerminationStateTh(
         )
   {
   }

   ///
   /// constructor
   ///
   CMapLapParaSolverTerminationStateTh(
         int          inThreadId,                             ///< thread id
         int          inInterrupted,                          ///< indicate that this solver is interrupted or not.
                                                              ///< 0: not interrupted,
                                                              ///< 1: interrupted
                                                              ///< 2: checkpoint,
                                                              ///< 3: racing-ramp up
         int          inRank,                                 ///< rank of this solver
         int          inNParaNodesReceived,                   ///< number of ParaNodes received in this ParaSolver
         int          inNParaNodesSolved,                     ///< number of ParaNodes solved ( received ) in this ParaSolver
         double       inRunningTime,                          ///< this solver running time
         double       inIdleTimeToFirstParaNode,              ///< idle time to start solving the first ParaNode
         double       inIdleTimeBetweenParaNodes,             ///< idle time between ParaNodes processing
         double       inIdleTimeAfterLastParaNode,            ///< idle time after the last ParaNode was solved
         double       inIdleTimeToWaitNotificationId,         ///< idle time to wait notification Id messages
         double       inIdleTimeToWaitAckCompletion,          ///< idle time to wait ack completion message
         double       inIdleTimeToWaitToken,                  ///< idle time to wait token
         double       inIdleTimeToWaitSolverState,            ///< idle time to wait solver state
         double       inIdleTimeToWaitPackedVector,           ///< idle time to wait packed vector
         double       inIdleTimeToWaitSolution,               ///< idle time to wait solution
         double       inIdleTimeToWaitBasis,                  ///< idle time to wait basis
         double       inIdleTimeToWaitIsend,                  ///< idle time to wait Isend
         double       inDetTime,                              ///< deterministic time, -1: should be non-deterministic
         int          inNParaTasksDeepBkzReceived,            ///< number of DeepBkz ParaTasks received in this solver
         int          inNParaTasksEnumReceived,               ///< number of Enum ParaTasks received in this solver
         int          inNParaTasksSieveReceived,              ///< number of Sieve ParaTasks received in this solver
         double       inRunningTimeDeepBkz,                   ///< this solver running time of DeepBkz
         double       inRunningTimeEnum,                      ///< this solver running time of Enum
         double       inRunningTimeSieve,                     ///< this solver running time of Sieve
         int          inNVectorsReceivedDeepBkz,              ///< number of vectors received in DeepBkz
         int          inNVectorsReceivedEnum,                 ///< number of vectors received in Enum
         int          inNVectorsReceivedSieve,                ///< number of vectors received in Sieve
         int          inNVectorsSentDeepBkz,                  ///< number of vectors sent in DeepBkz
         int          inNVectorsSentEnum,                     ///< number of vectors sent in Enum
         int          inNVectorsSentSieve,                    ///< number of vectors sent in Sieve
         int          inNBasesSentDeepBkz,                    ///< number of vectors sent in DeepBkz
         int          inNSolverStateSent                      ///< number of solver states sent
         )
         : CMapLapParaSolverTerminationState(
               inThreadId,
               inInterrupted,
               inRank,
               inNParaNodesReceived,
               inNParaNodesSolved,
               inRunningTime,
               inIdleTimeToFirstParaNode,
               inIdleTimeBetweenParaNodes,
               inIdleTimeAfterLastParaNode,
               inIdleTimeToWaitNotificationId,
               inIdleTimeToWaitAckCompletion,
               inIdleTimeToWaitToken,
               inIdleTimeToWaitSolverState,
               inIdleTimeToWaitPackedVector,
               inIdleTimeToWaitSolution,
               inIdleTimeToWaitBasis,
               inIdleTimeToWaitIsend,
               inDetTime,
               inNParaTasksDeepBkzReceived,
               inNParaTasksEnumReceived,
               inNParaTasksSieveReceived,
               inRunningTimeDeepBkz,
               inRunningTimeEnum,
               inRunningTimeSieve,
               inNVectorsReceivedDeepBkz,
               inNVectorsReceivedEnum,
               inNVectorsReceivedSieve,
               inNVectorsSentDeepBkz,
               inNVectorsSentEnum,
               inNVectorsSentSieve,
               inNBasesSentDeepBkz,
               inNSolverStateSent
           )
   {
   }

   ///
   /// send this object
   /// @return always 0 (for future extensions)
   ///
   virtual void send(
         UG::ParaComm *comm,           ///< communicator used
         int destination,              ///< destination rank
         int tag                       ///< UG::TagTerminated
         );

   ///
   /// receive this object
   /// @return always 0 (for future extensions)
   ///
   virtual void receive(
         UG::ParaComm *comm,          ///< communicator used
         int source,                  ///< source rank
         int tag                      ///< UG::TagTerminated
         );

};

} // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_SOLVER_TERMINATION_STATE_TH_H__
