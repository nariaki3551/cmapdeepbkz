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


/**@file    cmapLapParaSolverTerminationState.h
 * @brief   ParaSolverTerminationState extension for CMAP-LAP
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_SOLVER_TERMINATION_STATE_H__
#define __CMAP_LAP_PARA_SOLVER_TERMINATION_STATE_H__

#include <memory>
#include "ug/paraSolverTerminationState.h"
#ifdef UG_WITH_ZLIB
#include "ug/gzstream.h"
#endif
namespace UG { class ParaComm; }
namespace UG { class ParaInitiator; }

namespace ParaCMapLAP
{

///
/// class CMapLapParaSolverTerminationState
/// (Solver termination state in a ParaSolver)
///
class CMapLapParaSolverTerminationState : public UG::ParaSolverTerminationState, public std::enable_shared_from_this<CMapLapParaSolverTerminationState>
{
protected:
   ///-------------------------------------
   /// Counters related to this ParaSolver
   ///-------------------------------------
   int    threadId;                                    ///< thread id
   double idleTimeToWaitSolverState;                   ///< idle time to wait solver state
   double idleTimeToWaitPackedVector;                  ///< idle time to wait packed vector
   double idleTimeToWaitSolution;                      ///< idle time to wait solution
   double idleTimeToWaitBasis;                         ///< idle time to wait basis
   double idleTimeToWaitIsend;                         ///< idle time to wait Isend
   int    nParaTasksDeepBkzReceived;                   ///< number of DeepBkz ParaTasks received in this solver
   int    nParaTasksEnumReceived;                      ///< number of Enum ParaTasks received in this solver
   int    nParaTasksSieveReceived;                     ///< number of Sieve ParaTasks received in this solver
   double runningTimeDeepBkz;                          ///< this solver running time of DeepBkz
   double runningTimeEnum;                             ///< this solver running time of Enum
   double runningTimeSieve;                            ///< this solver running time of Sieve
   int    nVectorsReceivedDeepBkz;                     ///< number of vectors received in DeepBkz
   int    nVectorsReceivedEnum;                        ///< number of vectors received in Enum
   int    nVectorsReceivedSieve;                       ///< number of vectors received in Sieve
   int    nVectorsSentDeepBkz;                         ///< number of vectors sent in DeepBkz
   int    nVectorsSentEnum;                            ///< number of vectors sent in Enum
   int    nVectorsSentSieve;                           ///< number of vectors sent in Sieve
   int    nBasesSentDeepBkz;                           ///< number of vectors sent in DeepBkz
   int    nSolverStateSent;                            ///< number of solver states sent


   ///-----------------------------
   ///  times for root node process
   ///-----------------------------

public:

   ///
   /// default constructor
   ///
   CMapLapParaSolverTerminationState(
         )
         : ParaSolverTerminationState(),
           threadId(-1),
           idleTimeToWaitSolverState(0.0),
           idleTimeToWaitPackedVector(0.0),
           idleTimeToWaitSolution(0.0),
           idleTimeToWaitBasis(0.0),
           idleTimeToWaitIsend(0.0),
           nParaTasksDeepBkzReceived(0),
           nParaTasksEnumReceived(0),
           nParaTasksSieveReceived(0),
           runningTimeDeepBkz(0.0),
           runningTimeEnum(0.0),
           runningTimeSieve(0.0),
           nVectorsReceivedDeepBkz(0),
           nVectorsReceivedEnum(0),
           nVectorsReceivedSieve(0),
           nVectorsSentDeepBkz(0),
           nVectorsSentEnum(0),
           nVectorsSentSieve(0),
           nBasesSentDeepBkz(0),
           nSolverStateSent(0)
   {
   }

   ///
   /// constructor
   ///
   CMapLapParaSolverTerminationState(
         int          inThreadId,                      ///< thread id
         int          inInterrupted,                   ///< indicate that this solver is interrupted or not.
                                                       ///< 0: not interrupted,
                                                       ///< 1: interrupted
                                                       ///< 2: checkpoint,
                                                       ///< 3: racing-ramp up
         int          inRank,                          ///< rank of this solver
         int          inNParaTasksReceived,            ///< number of ParaTasks received in this ParaSolver
         int          inNParaTasksSolved,              ///< number of ParaTasks solved ( received ) in this ParaSolver
         double       inRunningTime,                   ///< this solver running time
         double       inIdleTimeToFirstParaTask,       ///< idle time to start solving the first ParaTask
         double       inIdleTimeBetweenParaTasks,      ///< idle time between ParaTasks processing
         double       inIdleTimeAfterLastParaTask,     ///< idle time after the last ParaTask was solved
         double       inIdleTimeToWaitNotificationId,  ///< idle time to wait notification Id messages
         double       inIdleTimeToWaitAckCompletion,   ///< idle time to wait ack completion message
         double       inIdleTimeToWaitToken,           ///< idle time to wait token
         double       inIdleTimeToWaitSolverState,     ///< idle time to solver state
         double       inIdleTimeToWaitPackedVector,    ///< idle time to wait packed vector
         double       inIdleTimeToWaitSolution,        ///< idle time to wait solution
         double       inIdleTimeToWaitBasis,           ///< idle time to wait basis
         double       inIdleTimeToWaitIsend,           ///< idle time to wait Isend
         double       inDetTime,                       ///< deterministic time, -1: should be non-deterministic
         int          inNParaTasksDeepBkzReceived,     ///< number of DeepBkz ParaTasks received in this solver
         int          inNParaTasksEnumReceived,        ///< number of Enum ParaTasks received in this solver
         int          inNParaTasksSieveReceived,       ///< number of Sieve ParaTasks received in this solver
         double       inRunningTimeDeepBkz,            ///< this solver running time of DeepBkz
         double       inRunningTimeEnum,               ///< this solver running time of Enum
         double       inRunningTimeSieve,              ///< this solver running time of Sieve
         int          inNVectorsReceivedDeepBkz,       ///< number of vectors received in DeepBkz
         int          inNVectorsReceivedEnum,          ///< number of vectors received in Enum
         int          inNVectorsReceivedSieve,         ///< number of vectors received in Sieve
         int          inNVectorsSentDeepBkz,           ///< number of vectors sent in DeepBkz
         int          inNVectorsSentEnum,              ///< number of vectors sent in Enum
         int          inNVectorsSentSieve,             ///< number of vectors sent in Sieve
         int          inNBasesSentDeepBkz,             ///< number of vectors sent in DeepBkz
         int          inNSolverStateSent               ///< number of solver states sent
         )
         : ParaSolverTerminationState(inInterrupted, inRank, inNParaTasksReceived, inNParaTasksSolved,
                                      inRunningTime, inIdleTimeToFirstParaTask, inIdleTimeBetweenParaTasks, inIdleTimeAfterLastParaTask,
                                      inIdleTimeToWaitNotificationId, inIdleTimeToWaitAckCompletion, inIdleTimeToWaitToken, inDetTime),
         threadId(inThreadId)
   {
      idleTimeToWaitSolverState  = inIdleTimeToWaitSolverState;
      idleTimeToWaitPackedVector = inIdleTimeToWaitPackedVector;
      idleTimeToWaitSolution     = inIdleTimeToWaitSolution;
      idleTimeToWaitBasis        = inIdleTimeToWaitBasis;
      idleTimeToWaitIsend        = inIdleTimeToWaitIsend;
      nParaTasksDeepBkzReceived  = inNParaTasksDeepBkzReceived;
      nParaTasksEnumReceived     = inNParaTasksEnumReceived;
      nParaTasksSieveReceived    = inNParaTasksSieveReceived;
      runningTimeDeepBkz         = inRunningTimeDeepBkz;
      runningTimeEnum            = inRunningTimeEnum;
      runningTimeSieve           = inRunningTimeSieve;
      nVectorsReceivedDeepBkz    = inNVectorsReceivedDeepBkz;
      nVectorsReceivedEnum       = inNVectorsReceivedEnum;
      nVectorsReceivedSieve      = inNVectorsReceivedSieve;
      nVectorsSentDeepBkz        = inNVectorsSentDeepBkz;
      nVectorsSentEnum           = inNVectorsSentEnum;
      nVectorsSentSieve          = inNVectorsSentSieve;
      nBasesSentDeepBkz          = inNBasesSentDeepBkz;
      nSolverStateSent           = inNSolverStateSent;
   }

   ///
   /// destructor
   ///
   virtual ~CMapLapParaSolverTerminationState(
         )
   {
   }

   ///
   /// get thread ID
   /// @return thread ID
   ///
   virtual int getThreadId(
         )
   {
      return threadId;
   }

   ///
   /// get runningTime of DeepBkz
   /// @return runningTimeDeepBkz
   ///
   virtual double getRunningTimeDeepBkz(
         )
   {
      return runningTimeDeepBkz;
   }

   ///
   /// get runningTime of Enum
   /// @return runningTimeEnum
   ///
   virtual double getRunningTimeEnum(
         )
   {
      return runningTimeEnum;
   }

   ///
   /// get runningTime of Sieve
   /// @return runningTimeSieve
   ///
   virtual double getRunningTimeSieve(
         )
   {
      return runningTimeSieve;
   }

   /// stringfy CMapLapParaSolverTerminationState object
   /// @return string to show inside of CMapLapParaSolverTerminationState object
   ///
   virtual std::string toString(
         UG::ParaInitiator *initiator=nullptr   ///< pointer to ParaInitiator object
         );

   /// stringfy CMapLapParaSolverTerminationState object as csv format
   /// @return string to show inside of CMapLapParaSolverTerminationState object as csv format
   ///
   virtual std::string toCsvString(
         UG::ParaInitiator *initiator=nullptr   ///< pointer to ParaInitiator object
         );

#ifdef UG_WITH_ZLIB

   ///
   /// write CMapLapParaSolverTerminationState to checkpoint file
   ///
   virtual void write(
         gzstream::ogzstream &out      ///< gzstream to output
         );

   ///
   /// read CMapLapParaSolverTerminationState from checkpoint file
   ///
   virtual bool read(
         UG::ParaComm *comm,           ///< communicator used
         gzstream::igzstream &in       ///< gzstream to input
         );

#endif

};

} // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_SOLVER_TERMINATION_STATE_H__
