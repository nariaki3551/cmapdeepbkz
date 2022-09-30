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

/**@file    cmapLapParaLoadCoordinator.h
 * @brief   Load Coordinator.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_LOADCOORDINATOR_H__
#define __CMAP_LAP_PARA_LOADCOORDINATOR_H__

#include <iostream>
#include <deque>
#include <vector>
#include <memory>
#include "ug/paraComm.h"
#include "ug/paraTagDef.h"
#include "ug/paraParamSet.h"
#include "ug/paraSolverPool.h"
#include "ug/paraLoadCoordinator.h"
#include "cmapLapParaBasis.h"
#include "cmapLapParaLattice.h"
#include "cmapLapParaInitiator.h"
#include "cmapLapParaSolverPool.h"
#include "cmapLapParaInstancePool.h"
#include "cmapLapParaShareDataPool.h"
#include "cmapLapParaLoadCoordinatorTerminationState.h"
namespace UG { class ParaDeterministicTimer; }
namespace UG { class ParaInitiator; }
namespace UG { class PcmapDeepBkzParaLoadCoordinator; }
namespace UG { class ParaTimer; }
namespace ParaCMapLAP { class BasisElement; }
namespace ParaCMapLAP { class BasisElementQueue; }
namespace ParaCMapLAP { class CMapLapParaSolution; }
namespace ParaCMapLAP { class CMapLapParaSolverLocalComm; }
namespace ParaCMapLAP { class CMapLapParaSolverPool; }
namespace ParaCMapLAP { class CMapLapParaTask; }
namespace ParaCMapLAP { class CMapLapParaTaskPoolInAscendingOrder; }
namespace ParaCMapLAP { class CMapLapParaTaskQueue; }
namespace ParaCMapLAP { class ShareDataPool; }
namespace ParaCMapLAP { class VectorElementQueue; }


namespace ParaCMapLAP
{

struct LocalSolverThreadData_t {
   int                           argc;
   char                          **argv;
   CMapLapParaSolverLocalComm    *comm;
   UG::ParaParamSet              *paraParams;
   UG::ParaInstance              *instance;
   UG::ParaDeterministicTimer    *paraDetTimer;
   UG::ParaTimer                 *paraTimer;
   int                           rank;          ///< rank is always 1 in local solver
   int                           threadId;
};

typedef struct LocalSolverThreadData_t LocalSolverThreadData;


#ifdef UG_WITH_ZLIB

class CheckpointElement
{
public:
   std::shared_ptr<VectorElementQueue>    vectorElementQueue;
   std::shared_ptr<BasisElementQueue>     basisElementQueue;
   std::shared_ptr<CMapLapParaTaskQueue>  paraActiveTaskQueue;
   std::shared_ptr<CMapLapParaTaskQueue>  paraDeepBkzTaskQueue;
   std::shared_ptr<CMapLapParaTaskQueue>  paraEnumTaskQueue;
   std::shared_ptr<CMapLapParaTaskQueue>  paraSieveTaskQueue;
   std::shared_ptr<CMapLapParaSolution>   globalBestCMapLapSolution;
   std::shared_ptr<BasisElement>          incumbentBasis;
   char                                *probName;
   char                                *lastCheckpointTimeStr;
   std::string                         lcStatString;
   double                              writeCheckpointTime;

   CheckpointElement(
         )
      :
         vectorElementQueue(nullptr),
         basisElementQueue(nullptr),
         paraActiveTaskQueue(nullptr),
         paraDeepBkzTaskQueue(nullptr),
         paraEnumTaskQueue(nullptr),
         paraSieveTaskQueue(nullptr),
         globalBestCMapLapSolution(nullptr),
         incumbentBasis(nullptr),
         probName(nullptr),
         writeCheckpointTime(0)
   {
      lastCheckpointTimeStr = new char[256];
      probName = new char[256];
   }
   virtual ~CheckpointElement()
   {
      if( probName )                   delete [] probName;
      if( lastCheckpointTimeStr )      delete [] lastCheckpointTimeStr;
   }

   virtual void setLastCheckpointTimeStr(char *timeStr)
   {
      strcpy(lastCheckpointTimeStr, timeStr);
   }
   virtual void setProbName(const char *inProbName)
   {
      strcpy(probName, inProbName);
   }
};
#endif


class CMapLapParaLoadCoordinator : public UG::ParaLoadCoordinator
{

   typedef int(CMapLapParaLoadCoordinator::*CMapLapMessageHandlerFunctionPointer)(int, int);

   CMapLapMessageHandlerFunctionPointer *localMessageHandler;   ///< table for message handlers
   static CMapLapParaSolverLocalComm    *lcLocalComm;           ///< communicator for local solver
   static CMapLapParaSolverLocalComm    *lcCheckpointComm;      ///< communicator for checkpoint writer

   CMapLapParaInitiator *cmapLapParaInitiator;     ///< dynamic cast of ParaInitiator
   CMapLapParaSolverPool *cmapLapParaSolverPool;   ///< dynamic cast of ParaSolverPool

   ///
   /// Vector and Basis Pool
   ///
   std::shared_ptr<ShareDataPool> shareDataPool;
   std::shared_ptr<InstancePool> instancePool;

   ///
   /// Task Pool
   ///
   std::shared_ptr<CMapLapParaTaskPoolInAscendingOrder> paraDeepBkzTaskPool;
   std::shared_ptr<CMapLapParaTaskPoolInAscendingOrder> paraEnumTaskPool;
   std::shared_ptr<CMapLapParaTaskPoolInAscendingOrder> paraSieveTaskPool;

   ///
   /// Pools for local solvers
   ///
   std::shared_ptr<CMapLapParaSolverPool> paraLocalSolverPool;             ///< local solver pool
   std::shared_ptr<CMapLapParaTaskPoolInAscendingOrder> paraLocalTaskPool; ///< task pool for local solvers

   ///
   /// Objects to output statistics and progress state
   ///
   CMapLapParaLoadCoordinatorTerminationState lcts;      ///< LoadCoordinatorTerminationState: counters and times
   std::ofstream      ofsAssingmentTable;                ///< ofstream for assignment table
   std::ostream       *osAssignmentTable;                ///< ostream for assignment table to switch output location
   std::ofstream      ofsCsvAssingmentTable;             ///< ofstream for assignment table
   std::ostream       *osCsvAssignmentTable;             ///< ostream for assignment table to switch output location
   std::ofstream      ofsCsvLogSolvingStatus;            ///< ofstream for solving status
   std::ostream       *osCsvLogSolvingStatus;            ///< ostram for solving status to switch output location
   std::ofstream      ofsCsvStatisticsFinalRun;          ///< ofstream for statistics of the final run
   std::ostream       *osCsvStatisticsFinalRun;          ///< ostream for statistics of the final run
   std::ofstream      ofsCsvLogCheckpoint;               ///< ofstream for checkpoint log
   std::ostream       *osCsvLogCheckpoint;               ///< ostream for checkpoint log
   std::ofstream      ofsCsvLctsStat;                    ///< ofstream for checkpoint log
   std::ostream       *osCsvLctsStat;                    ///< ostream for checkpoint log
   std::ofstream      ofsLogShareDataPool;               ///< ofstream for vector pool log
   std::ostream       *osLogShareDataPool;               ///< ostream for vector pool log
   std::ofstream      ofsLogShareDataPoolStat;           ///< ofstream for vector pool log
   std::ostream       *osLogShareDataPoolStat;           ///< ostream for vector pool log
   std::ofstream      ofsLogSimilarity;                  ///< ofstream for shared lattice log
   std::ostream       *osLogSimilarity;                  ///< ostream for shared lattice log
   std::ofstream      ofsLogNotificationInterval;        ///< ofstream for notification interval log
   std::ostream       *osLogNotificationInterval;        ///< ostream for notification interval log


#ifdef UG_WITH_ZLIB
   ///
   /// Queue for Checkpoint Threading
   ///
   std::unique_ptr<CheckpointElement> currentCheckpointElement;
   std::unique_ptr<CheckpointElement> nextCheckpointElement;
#endif


   int                nSolvers;                           ///< number of all available solver threads
   int                nThreadsPerRank;                    ///< number of threads in a solver
   int                currentSeed;                        ///< maximum seed that has already been sent to solver
   int                nStatus;                            ///< number of status after output a log header
   double             previousAssingmentTableOutputTime;  ///< to keep assignment table output time

   bool               assigningParaTask;                  ///< true if lc is assigning paraTask else false
   std::vector<bool>  hasSentBasisRequest;                ///< manager of SendBasisRequest
   bool               lowerBoundIsReached;                ///< true if LC has the solution that is better than globalLowerBoundOfSquaredNorm
   double             globalLowerBoundOfSquaredNorm;      ///< global lattice vector's lowerbound of squared norm


   ///
   /// run a local solver thread
   ///
   static void * runSolverThread(
         void *threadData
         );

   ///////////////////////
   ///
   /// Message handlers
   ///
   ///////////////////////

   ///
   /// function to process TagNode message
   /// @return always 0 (for extension)
   ///
   virtual int processTagTask(
         int source,                                      ///< source solver rank
         int tag                                          ///< UG::TagTask
         )
   {
      return 0;
   }

   ///
   /// function to process UG::TagSolution message
   /// @return always 0 (for extension)
   ///
   virtual int processTagSolution(
         int source,                                      ///< source solver rank
         int tag                                          ///< UG::TagSolution
         );

   ///
   /// function to process UG::TagSolution message from local solver
   /// @return always 0 (for extension)
   ///
   virtual int processLocalTagSolution(
         int source,                                      ///< source solver rank
         int tag                                          ///< UG::TagSolution
         );

   ///
   /// function to process TagBasis message
   /// @return always 0 (for extension)
   ///
   virtual int processTagBasis(
         int source,                                      ///< source solver rank
         int tag                                          ///< TagBasis
         );

   ///
   /// function to process TagBasis message from local solver
   /// @return always 0 (for extension)
   ///
   virtual int processLocalTagBasis(
         int source,                                      ///< source solver rank
         int tag                                          ///< TagBasis
         );

   ///
   /// function to process TagVectorRequest message
   /// @return always 0 (for extension)
   ///
   virtual int processTagVectorRequest(
         int source,                                      ///< source solver rank
         int tag                                          ///< TagVectorRequest
         );

   ///
   /// function to process TagVectorRequest message from local solver
   /// @return always 0 (for extension)
   ///
   virtual int processLocalTagVectorRequest(
         int source,                                      ///< source solver rank
         int tag                                          ///< TagVectorRequest
         );

   ///
   /// function to process UG::TagSolverState message
   /// @return always 0 (for extension)
   ///
   virtual int processTagSolverState(
         int source,                                      ///< source solver rank
         int tag                                          ///< UG::TagSolverState
         );

   ///
   /// function to process TagCompletionOfCalculation message
   /// @return always 0 (for extension)
   ///
   virtual int processTagCompletionOfCalculation(
         int source,                                      ///< source solver rank
         int tag                                          ///< UG::TagCompletionOfCalculation
         );

   ///
   /// function to process TagCompletionOfCalculation message from local solver
   /// @return always 0 (for extension)
   ///
   virtual int processLocalTagCompletionOfCalculation(
         int source,                                      ///< source solver rank
         int tag                                          ///< UG::TagCompletionOfCalculation
         );

   ///
   /// function to process UG::TagTerminated message
   /// @return always 0 (for extension)
   ///
   virtual int processTagTerminated(
         int source,                                      ///< source solver rank
         int tag                                          ///< UG::TagTerminated
         );

   ///
   /// function to process TagCMapLapPackedVector message
   /// @return always 0 (for extension)
   ///
   virtual int processTagCMapLapPackedVector(
         int source,                                      ///< source solver rank
         int tag                                          ///< TagCMapLapPackedVector
         );

   ///
   /// function to process TagCMapLapPackedVector message from local solver
   /// @return always 0 (for extension)
   ///
   virtual int processLocalTagCMapLapPackedVector(
         int source,                                      ///< source solver rank
         int tag                                          ///< TagCMapLapPackedVector
         );

   ///
   /// function to process UG::TagTerminated message from local solver
   /// @return always 0 (for extension)
   ///
   virtual int processLocalTagTerminated(
         int source,                                      ///< source solver rank
         int tag                                          ///< UG::TagTerminated
         );

   ///
   /// function to process UG::TagHardTimeLimit message from local solver
   /// @return always 0 (for extension)
   ///
   virtual int processLocalTagHardTimeLimit(
         int source,                                      ///< source solver rank
         int tag                                          ///< UG::TagHardTimeLimit
         );

   ///////////////////////
   ///
   /// message handlers specialized for racing ramp-up
   ///
   ///////////////////////

   ///
   /// function to process UG::TagSolverState message in racing ramp-up stage
   /// @return always 0 (for extension)
   ///
   virtual int processRacingRampUpTagSolverState(
         int source,                                      ///< source solver rank
         int tag                                          ///< UG::TagSolverState
         )
   {
      return 0;
   }

   ///
   /// function to process TagCompletionOfCalculation message in racing ramp-up stage
   /// @return always 0 (for extension)
   ///
   virtual int processRacingRampUpTagCompletionOfCalculation(
         int source,                                      ///< source solver rank
         int tag                                          ///< UG::TagCompletionOfCalculation
         )
   {
      return 0;
   }


   ///
   /// for remove compile error of -Woverloaded-virtual
   ///
   using UG::ParaLoadCoordinator::run;

   ///
   /// run function to start main process
   ///
   virtual void run(
         );


   ///
   /// @return true if it should create and assign task
   ///
   virtual bool sendParaTasksToIdleSolvers(
         )
   {
      return 1;
   };

   ///
   /// @return true if it should create and assign task
   ///
   virtual bool shouldCreateAndAssignTask(
         SolverType &solverType,
         bool &hasSetInitialSolvers,
         bool &hasSetInitialDeepBkzSolvers
         );

   ///
   /// check the terminate condition
   /// @return true if it is satisfy the terminate condition else false
   ///
   virtual bool shouldBreakRun(
         );

   ///
   /// interrupt the required number of active solvers and return paraTask of solverType
   /// @return CMapLapParaTask pointer
   ///
   virtual std::shared_ptr<CMapLapParaTask> getParaTask(
         SolverType solverType,  ///< solver task should be created
         int begin,
         int end
         );

   ///
   /// create CMapLapParaTask DeepBkz and insert into DeepBkzTaskPool
   ///
   virtual int createParaTaskDeepBkz(
         int begin,
         int end
         );

   ///
   /// create paraTaskEnum using basis in InstancePool
   ///
   virtual int createParaTaskEnums(
         int begin,
         int end
         );

   // ///
   // /// create paraTaskEnum (SubEnum) using basis in InstancePool
   // ///
   // virtual int createParaTaskSubEnums(
   //       int begin,
   //       int end
   //       );

   ///
   /// create paraTaskSieve using basis in InstancePool
   ///
   virtual int createParaTaskSieve(
         int begin,
         int end
         );

   ///
   /// create a paraTask for local solvers in LoadCoordinator when a vector is inserted in vectorElement
   /// @return CMapLapParaTask*
   ///
   virtual CMapLapParaTask* createParaTaskForLocalSolver(
         );

   ///
   /// assign cmapLapParaTask
   /// @return true if success to activate cmapLapParaTask else false
   ///
   virtual bool assignParaTask(
         std::shared_ptr<CMapLapParaTask> cmapLapParaTask
         );

   ///
   /// assign cmapLapParaLocalTask
   /// @return true if success to assign cmapLapParaLocalTask else false
   ///
   virtual bool assignParaLocalTask(
         std::shared_ptr<CMapLapParaTask> cmapLapParaLocalTask
         );

   ///
   ///
   ///
   virtual void waitForAnyMessageFromAnyWhere(
         );

   ///
   ///
   ///
   virtual void waitForAnyMessageFromLocal(
         );


   ///
   /// send interrupt request to all solvers
   ///
   virtual void sendInterruptRequest(
        )
   {
      int exitSolverRequest = 0;    // do nothing
      if( paraSolverPool->getNumActiveSolvers() > 0 )
      {
         if( interruptIsRequested ) return;
         for( int i = 1; i < paraComm->getSize(); i++ )
         {
            PARA_COMM_CALL(
                  paraComm->send( &exitSolverRequest, 1, UG::ParaINT, i, UG::TagInterruptRequest )
            );
         }
      }
      interruptIsRequested = true;
   }

   ///
   /// update and notice NotificationInterval to solvers
   ///
   virtual void updateNotificationInterval(
         double termTime,                 ///< time of this term
         double &notificationInterval,    ///< current solver's notification interval
         double &leftNotificationInterval ///< previous solver's notification interval
         );


   ///
   /// output similarity of basis of each solver
   ///
   virtual void outputSimilarityOfBasisHeader(
         );

   ///
   /// output similarity of basis of each solver
   ///
   virtual void outputSimilarityOfBasis(
         std::deque<std::shared_ptr<LatticeBasis<int>>> &basisDeque
         );


   ///
   /// return seed and increment this
   ///
   virtual int seedGenerator(
         )
   {
      currentSeed++;
      return currentSeed;
   }


#ifdef UG_WITH_ZLIB

   ///
   /// function to update checkpoint files
   ///
   virtual void updateCheckpointFiles(
         );

   ///
   /// copy objects for checkpoint
   /// @return true if Load Coardinator chould send new task to Checkpoint Writer else false
   ///
   virtual bool copyCheckpointObjects(
         );

   ///
   /// run function to start checkpoint writer process
   ///
   static void * runCheckpointThread(
         void *threadData
         );

#endif

public:

   CMapLapParaLoadCoordinator(
      UG::ParaComm *inComm,
      UG::ParaParamSet *inParaParamSet,
      UG::ParaInitiator *inParaInitiator,
      bool *inRacingSolversExist,
      UG::ParaTimer *inParaTimer,
      UG::ParaDeterministicTimer *detTimer,
      int inNThreadsInSolver
      );

   virtual ~CMapLapParaLoadCoordinator(
         );

   ///
   /// increment nStatus
   ///
   virtual void incrementNStatus(
         )
   {
      if( nStatus == 0  )
      {
         *osLogSolvingStatus << std::endl << Logging::getCMapLapParaSolverStateHeader();
      }
      nStatus++;
      if( nStatus >= 20 ){ nStatus = 0; }
   }

   ///
   /// getter of lcLocalComm
   ///
   virtual CMapLapParaSolverLocalComm * getLcLocalComm(
         )
   {
      return lcLocalComm;
   }

   ///
   /// getter of cmapLapParaSolverPool
   ///
   virtual CMapLapParaSolverPool * getCmapLapParaSolverPool(
         )
   {
      return cmapLapParaSolverPool;
   }


   ///
   /// getter of paraDeepBkzTaskPool
   ///
   virtual std::shared_ptr<CMapLapParaTaskPoolInAscendingOrder> getParaDeepBkzTaskPool(
         )
   {
      return paraDeepBkzTaskPool;
   }

   ///
   /// getter of paraEnumTaskPool
   ///
   virtual std::shared_ptr<CMapLapParaTaskPoolInAscendingOrder> getParaEnumTaskPool(
         )
   {
      return paraEnumTaskPool;
   }

   ///
   /// getter of paraSieveTaskPool
   ///
   virtual std::shared_ptr<CMapLapParaTaskPoolInAscendingOrder> getParaSieveTaskPool(
         )
   {
      return paraSieveTaskPool;
   }

   ///
   /// getter of cmapLapParaParams
   ///
   virtual CMapLapParaParamSet * getCmapLapParaParams(
         )
   {
      return dynamic_cast<CMapLapParaParamSet *>(paraParams);
   }

   ///
   /// getter of cmapLapParaInitiator
   ///
   virtual CMapLapParaInitiator * getCmapLapParaInitiator(
         )
   {
      return cmapLapParaInitiator;
   }

   ///
   /// getter of instancePool
   ///
   virtual std::shared_ptr<InstancePool> getInstancePool(
         )
   {
      return instancePool;
   }

   ///
   /// getter of lcts
   ///
   virtual CMapLapParaLoadCoordinatorTerminationState & getLcts(
         )
   {
      return lcts;
   }

   ///
   /// getter of nSolvers
   ///
   virtual int getNSolvers(
         )
   {
      return nSolvers;
   }

   ///
   /// getter of osLogSolvingStatus
   ///
   virtual std::ostream * getOsLogSolvingStatus(
         )
   {
      return osLogSolvingStatus;
   }

   ///
   /// getter of osCsvLogSolvingStatus
   ///
   virtual std::ostream * getOsCsvLogSolvingStatus(
         )
   {
      return osCsvLogSolvingStatus;
   }

   ///
   /// output assignment table
   ///
   virtual void outputAssignmentTable(
         char shortest     ///< a charactor to show shortest or not ('*': indicate shortest)
         );

   ///
   /// output assignment table
   ///
   virtual void outputCsvAssignmentTable(
         char shortest     ///< a charactor to show shortest or not ('*': indicate shortest)
         );

   ///
   ///
   ///
   virtual bool isAssigningParaTask(
         )
   {
      return assigningParaTask;
   }

#ifdef UG_WITH_ZLIB

   ///
   /// warm start (restart)
   ///
   virtual void warmStart(
         );

#endif

};

} // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_LOADCOORDINATOR_H__
