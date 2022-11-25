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

/**@file    cmapLapParaLoadCoordinator.cpp
 * @brief   Load coordinator.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <string>
#include <time.h>
#include <memory>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <signal.h>
#include "cmapLapParaTask.h"
#include "cmapLapParaBasis.h"
#include "cmapLapParaTagDef.h"
#include "cmapLapParaCommTh.h"
#include "cmapLapParaSolver.h"
#include "cmapLapParaLattice.h"
#include "cmapLapParaInstance.h"
#include "cmapLapParaSolution.h"
#include "cmapLapParaParamSet.h"
#include "cmapLapParaTaskPool.h"
#include "cmapLapParaShareData.h"
#include "cmapLapParaInitiator.h"
#include "cmapLapParaSolverPool.h"
#include "cmapLapParaSimilarity.h"
#include "cmapLapParaInstanceTh.h"
#include "cmapLapParaSolverState.h"
#include "cmapLapParaPackedVector.h"
#include "cmapLapParaInstancePool.h"
#include "cmapLapParaShareDataPool.h"
#include "cmapLapParaSolverLocalComm.h"
#include "cmapLapParaLoadCoordinator.h"
#include "cmapLapParaCalculationState.h"
#include "cmapLapParaCheckpointWriter.h"
#include "cmapLapParaSolverTerminationState.h"
#include "cmapLapParaLoadCoordinatorTerminationState.h"


namespace ParaCMapLAP
{

CMapLapParaSolverLocalComm *CMapLapParaLoadCoordinator::lcLocalComm = nullptr;
CMapLapParaSolverLocalComm *CMapLapParaLoadCoordinator::lcCheckpointComm = nullptr;

///
/// Constructor
///
CMapLapParaLoadCoordinator::CMapLapParaLoadCoordinator(
      UG::ParaComm *inComm,
      UG::ParaParamSet *inParaParamSet,
      UG::ParaInitiator *inParaInitiator,
      bool *inRacingSolversExist,
      UG::ParaTimer *inParaTimer,
      UG::ParaDeterministicTimer *detTimer,
      int inNThreadsInSolver
      )
   : UG::ParaLoadCoordinator(
         N_CMAP_LAP_TAGS,
         inComm,
         inParaParamSet,
         inParaInitiator,
         inRacingSolversExist,
         inParaTimer,
         detTimer
         ),
   localMessageHandler(nullptr),
   currentCheckpointElement(nullptr),
   nextCheckpointElement(nullptr),
   nThreadsPerRank(inNThreadsInSolver),
   nStatus(0),
   previousAssingmentTableOutputTime(0.0),
   assigningParaTask(false),
   lowerBoundIsReached(false),
   globalLowerBoundOfSquaredNorm(-1.0)
{

   CMapLapParaParamSet *cmapLapParaParams = dynamic_cast<CMapLapParaParamSet *>(paraParams);
   cmapLapParaInitiator = dynamic_cast<CMapLapParaInitiator *>(paraInitiator);

   nSolvers = (paraComm->getSize()-1) * nThreadsPerRank;
   hasSentBasisRequest = std::vector<bool>(nSolvers, false);

   currentSeed = cmapLapParaInitiator->getRandomizeSeed();
   lcts.setProcessTimesOfMessageHandler(N_CMAP_LAP_TAGS);

   ///
   /// Setting of globalLowerBoundOfSquaredNorm
   ///
   double lowerBoundOfNorm = paraParams->getRealParamValue(LowerBoundOfNorm);
   if( lowerBoundOfNorm > 0.0 )
   {
      globalLowerBoundOfSquaredNorm = lowerBoundOfNorm * lowerBoundOfNorm;
      std::cout
         << "*** "
         << "CMapLap sets globalLowerBoundOfSquaredNorm to " << globalLowerBoundOfSquaredNorm
         << " ***"
         << std::endl;
   }
   double lowerBoundOfApproxFactor = paraParams->getRealParamValue(LowerBoundOfApproxFactor);
   if( lowerBoundOfApproxFactor > 0.0 )
   {
      lowerBoundOfNorm = getCmapLapParaInitiator()->getLattice().GH * lowerBoundOfApproxFactor;
      if( globalLowerBoundOfSquaredNorm < 0.0 )
      {
         globalLowerBoundOfSquaredNorm = lowerBoundOfNorm * lowerBoundOfNorm;
      }
      else
      {
         globalLowerBoundOfSquaredNorm = std::min(
               globalLowerBoundOfSquaredNorm,
               lowerBoundOfNorm * lowerBoundOfNorm
               );
      }
      std::cout
         << "*** "
         << "CMapLap sets globalLowerBoundOfSquaredNorm to " << globalLowerBoundOfSquaredNorm
         << " from lowerBoundOfApproxFactor " << lowerBoundOfApproxFactor
         << " ***"
         << std::endl;
   }

   ///
   /// Setting of Message Handler
   ///
   CMapLapMessageHandlerFunctionPointer *cmapLapMessageHandler = reinterpret_cast<CMapLapMessageHandlerFunctionPointer *>(messageHandler);
   cmapLapMessageHandler[UG::TagSolution] = &ParaCMapLAP::CMapLapParaLoadCoordinator::processTagSolution;
   cmapLapMessageHandler[TagBasis] = &ParaCMapLAP::CMapLapParaLoadCoordinator::processTagBasis;
   cmapLapMessageHandler[TagVectorRequest] = &ParaCMapLAP::CMapLapParaLoadCoordinator::processTagVectorRequest;
   cmapLapMessageHandler[UG::TagSolverState] = &ParaCMapLAP::CMapLapParaLoadCoordinator::processTagSolverState;
   cmapLapMessageHandler[UG::TagCompletionOfCalculation] = &ParaCMapLAP::CMapLapParaLoadCoordinator::processTagCompletionOfCalculation;
   cmapLapMessageHandler[UG::TagTerminated] = &ParaCMapLAP::CMapLapParaLoadCoordinator::processTagTerminated;
   cmapLapMessageHandler[UG::TagHardTimeLimit] = &ParaCMapLAP::CMapLapParaLoadCoordinator::processTagHardTimeLimit;
   cmapLapMessageHandler[TagCMapLapPackedVector] = &ParaCMapLAP::CMapLapParaLoadCoordinator::processTagCMapLapPackedVector;

   ///
   /// Setting of Pool and some parameters
   ///
   paraSolverPool = new CMapLapParaSolverPool(1, inComm, inParaParamSet, inParaTimer, nThreadsPerRank );
   cmapLapParaSolverPool = dynamic_cast<CMapLapParaSolverPool *>(paraSolverPool);
   paraDeepBkzTaskPool = std::make_shared<CMapLapParaTaskPoolInAscendingOrder>();
   paraEnumTaskPool = std::make_shared<CMapLapParaTaskPoolInAscendingOrder>();
   paraSieveTaskPool = std::make_shared<CMapLapParaTaskPoolInAscendingOrder>();
   shareDataPool = std::make_shared<ShareDataPool>(paraParams->getIntParamValue(ShareDataPoolSize));
   instancePool = std::make_shared<InstancePool>(paraParams->getIntParamValue(InstancePoolSize));
   auto instanceBasis = cmapLapParaInitiator->getBasis();
   instancePool->insert(DBL_MAX, instanceBasis);
   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm);
   LatticeVector<int> v = instanceBasis.row(0);
   std::shared_ptr<CMapLapParaSolution> initialSol{
      cmapLapParaComm->createCMapLapParaSolution(0, v, v.squaredNorm())};
   cmapLapParaInitiator->tryToSetIncumbentSolution(initialSol);


   if( cmapLapParaParams->getBoolParamValue(CheckpointThreading) )
   {
      int argc = 3;
      char *argv[4];
      char argv0 [] = "LC";
      char argv1 [] = "-ntpr";
      char argv2 [] = "2";
      argv[0] = argv0;
      argv[1] = argv1;
      argv[2] = argv2;
      argv[3] = NULL;

      lcCheckpointComm = new CMapLapParaSolverLocalComm();
      lcCheckpointComm->init(argc, argv);

      LocalSolverThreadData *solverThreadData = new LocalSolverThreadData();
      solverThreadData->paraParams = paraParams;
      solverThreadData->paraDetTimer = paraDetTimer;
      solverThreadData->paraTimer = paraTimer;
      solverThreadData->rank = 1;
      solverThreadData->threadId = 1;    // to make consistency in CMapLapParaSolver message handlers
                                         // for preventing to receive messages from LC
                                         // becase solvers whose threadID=0 are received messages from LC
      std::thread t(runCheckpointThread, solverThreadData);
      t.detach();  // to make indipendently checkpoint writer

      lcCheckpointComm->mainThreadInit(paraParams);
   }


   ///
   /// Setting of LocalSolver
   ///
   if( cmapLapParaParams->getIntParamValue(NumOfLocalSolversInLC) > 0 )
   {
      localMessageHandler = new CMapLapMessageHandlerFunctionPointer[N_CMAP_LAP_TAGS];
      for( int i = 0; i < N_CMAP_LAP_TAGS; i++ )
      {
         localMessageHandler[i] = 0;
      }
      // install message handlers
      localMessageHandler[UG::TagSolution] = &ParaCMapLAP::CMapLapParaLoadCoordinator::processLocalTagSolution;
      localMessageHandler[TagBasis] = &ParaCMapLAP::CMapLapParaLoadCoordinator::processLocalTagBasis;
      localMessageHandler[TagVectorRequest] = &ParaCMapLAP::CMapLapParaLoadCoordinator::processLocalTagVectorRequest;
      localMessageHandler[UG::TagCompletionOfCalculation] = &ParaCMapLAP::CMapLapParaLoadCoordinator::processLocalTagCompletionOfCalculation;
      localMessageHandler[UG::TagTerminated] = &ParaCMapLAP::CMapLapParaLoadCoordinator::processLocalTagTerminated;
      localMessageHandler[UG::TagHardTimeLimit] = &ParaCMapLAP::CMapLapParaLoadCoordinator::processLocalTagHardTimeLimit;
      localMessageHandler[TagCMapLapPackedVector] = &ParaCMapLAP::CMapLapParaLoadCoordinator::processLocalTagCMapLapPackedVector;

      char charNumOfLocalSolverInLC[256];
      sprintf(charNumOfLocalSolverInLC, "%d", cmapLapParaParams->getIntParamValue(NumOfLocalSolversInLC) + 1);
      int argc = 3;
      char *argv[4];
      char argv0 [] = "LC";
      char argv1 [] = "-ntpr";
      argv[0] = argv0;
      argv[1] = argv1;
      argv[2] = charNumOfLocalSolverInLC;
      argv[3] = NULL;
      lcLocalComm = new CMapLapParaSolverLocalComm();
      lcLocalComm->init(argc, argv);

      // Note that "rank" is thread Id in local solvers,
      // for each "rank" thread Id is always 1 in paraLocalSolverPool.
      // To ask something to paraLocalSolverPool, thread id have to set zero
      paraLocalSolverPool = std::make_shared<CMapLapParaSolverPool>(1, lcLocalComm, paraParams, inParaTimer, 1);
      paraLocalTaskPool = std::make_shared<CMapLapParaTaskPoolInAscendingOrder>();

      for( int j = 1; j <= cmapLapParaParams->getIntParamValue(NumOfLocalSolversInLC); j++ )
      {
         LocalSolverThreadData *solverThreadData = new LocalSolverThreadData();
         solverThreadData->argc = argc;
         solverThreadData->argv = argv;
         solverThreadData->paraParams = paraParams;
         solverThreadData->instance = paraInitiator->getParaInstance();
         solverThreadData->paraDetTimer = paraDetTimer;
         solverThreadData->paraTimer = paraTimer;
         solverThreadData->rank = j;
         solverThreadData->threadId = 1;    // to make consistency in CMapLapParaSolver message handlers
                                            // for preventing to receive messages from LC
                                            // becase solvers whose threadID=0 are received messages from LC
         std::thread t(runSolverThread, solverThreadData);
         t.detach();  // to make indipendently each solvers
      }

      lcLocalComm->mainThreadInit(paraParams);
   }


   ///
   /// Setting of osAssignmentTable
   ///
   if( paraParams->getBoolParamValue(UG::Quiet) )
   {
      osAssignmentTable = &std::cout;
   }
   else
   {
      std::ostringstream s;
      s << paraParams->getStringParamValue(UG::LogSolvingStatusFilePath)
         << paraInitiator->getParaInstance()->getProbName() << "_LC" << paraComm->getRank() << "_T.assign";
      ofsAssingmentTable.open(s.str().c_str(), std::ios::app );
      if( !ofsAssingmentTable )
      {
         std::cout << "Tabular assignment table file cannot open : file name = " << s.str() << std::endl;
         exit(1);
      }
      osAssignmentTable = &ofsAssingmentTable;
   }
   *osAssignmentTable << Logging::getCMapLapParaAssignmentTableHeader();


   ///
   /// Setting of osCsvAssignmentTable
   ///
   {
      std::ostringstream s;
      s << paraParams->getStringParamValue(UG::LogSolvingStatusFilePath)
         << paraInitiator->getParaInstance()->getProbName() << "_LC" << paraComm->getRank() << "_T.assign.csv";
      ofsCsvAssingmentTable.open(s.str().c_str(), std::ios::app );
      if( !ofsCsvAssingmentTable )
      {
         std::cout << "Tabular assignment table file cannot open : file name = " << s.str() << std::endl;
         exit(1);
      }
      osCsvAssignmentTable = &ofsCsvAssingmentTable;
   }
   *osCsvAssignmentTable << Logging::getCMapLapParaCsvAssignmentTableHeader();


   ///
   /// Setting of osCsvLogSolvingStatus
   ///
   if( !paraParams->getBoolParamValue(UG::LogSolvingStatus) )
   {
      osCsvLogSolvingStatus = &std::cout;
   }
   else
   {
      std::ostringstream s;
      s << paraParams->getStringParamValue(UG::LogSolvingStatusFilePath)
         << paraInitiator->getParaInstance()->getProbName() << "_LC" << paraComm->getRank() << ".status.csv";
      ofsCsvLogSolvingStatus.open(s.str().c_str(), std::ios::app );
      if( !ofsCsvLogSolvingStatus )
      {
         std::cout << "Solving status log csv file cannot open : file name = " << s.str() << std::endl;
         exit(1);
      }
      osCsvLogSolvingStatus = &ofsCsvLogSolvingStatus;
   }
   *osCsvLogSolvingStatus << Logging::getCMapLapParaCsvSolverStateHeader();

   ///
   /// Setting of osCsvStatisticsFinalRun
   ///
   if( paraParams->getBoolParamValue(UG::Quiet) )
   {
      osCsvStatisticsFinalRun = &std::cout;
   }
   else
   {
      std::ostringstream s;
      s << paraParams->getStringParamValue(UG::LogSolvingStatusFilePath)
        << paraInitiator->getParaInstance()->getProbName()
        << "_statistics_final_LC" << paraComm->getRank() << ".csv";
      ofsCsvStatisticsFinalRun.open(s.str().c_str(), std::ios::app );
      if( !ofsCsvStatisticsFinalRun )
      {
         std::cout << "Tabular assignment csv file cannot open : file name = " << s.str() << std::endl;
         exit(1);
      }
      osCsvStatisticsFinalRun = &ofsCsvStatisticsFinalRun;
   }
   *osCsvStatisticsFinalRun << Logging::getCMapLapParaCsvStatisticsFinalRunHeader();

   ///
   /// Setting of osCsvLogCheckpoint
   ///
   if( !paraParams->getBoolParamValue(UG::LogSolvingStatus) )
   {
      osCsvLogCheckpoint = &std::cout;
   }
   else
   {
      std::ostringstream s;
      s << paraParams->getStringParamValue(UG::LogSolvingStatusFilePath)
         << paraInitiator->getParaInstance()->getProbName() << "_LC" << paraComm->getRank() << ".checkpoint.csv";
      ofsCsvLogCheckpoint.open(s.str().c_str(), std::ios::app );
      if( !ofsCsvLogCheckpoint )
      {
         std::cout << "Tabular checkpoint log csv file cannot open : file name = " << s.str() << std::endl;
         exit(1);
      }
      osCsvLogCheckpoint = &ofsCsvLogCheckpoint;
   }
   *osCsvLogCheckpoint << Logging::getCMapLapParaCsvCheckpointHeader();

   ///
   /// Setting of osCsvLctsStat
   ///
   {
      std::ostringstream s;
      s << paraParams->getStringParamValue(UG::LogSolvingStatusFilePath)
         << paraInitiator->getParaInstance()->getProbName() << "_LC" << paraComm->getRank() << ".processTimes.csv";
      ofsCsvLctsStat.open(s.str().c_str(), std::ios::app );
      if( !ofsCsvLctsStat )
      {
         std::cout << "Tabular checkpoint log csv file cannot open : file name = " << s.str() << std::endl;
         exit(1);
      }
      osCsvLctsStat = &ofsCsvLctsStat;
   }
   *osCsvLctsStat
      << lcts.toStringProcessTimesOfMessageHandlerHeader(paraComm, ",")
      << std::endl;

   ///
   /// Setting of osLogShareDataPool
   ///
   if( !paraParams->getBoolParamValue(LogShareDataPoolAll) )
   {
      osLogShareDataPool = &std::cout;
   }
   else
   {
      std::ostringstream s;
      s << paraParams->getStringParamValue(UG::LogSolvingStatusFilePath)
         << paraInitiator->getParaInstance()->getProbName() << "_LC" << paraComm->getRank() << ".sharedatapool.all.txt";
      ofsLogShareDataPool.open(s.str().c_str(), std::ios::app );
      if( !ofsLogShareDataPool )
      {
         std::cout << "Tabular vector pool log file cannot open : file name = " << s.str() << std::endl;
         exit(1);
      }
      osLogShareDataPool = &ofsLogShareDataPool;
   }

   ///
   /// Setting of osLogShareDataPoolStat
   ///
   if( !paraParams->getBoolParamValue(LogShareDataPoolStat) )
   {
      osLogShareDataPoolStat = &std::cout;
   }
   else
   {
      std::ostringstream s;
      s << paraParams->getStringParamValue(UG::LogSolvingStatusFilePath)
         << paraInitiator->getParaInstance()->getProbName() << "_LC" << paraComm->getRank() << ".sharedatapool.stat.csv";
      ofsLogShareDataPoolStat.open(s.str().c_str(), std::ios::app );
      if( !ofsLogShareDataPoolStat )
      {
         std::cout << "Tabular vector pool log file cannot open : file name = " << s.str() << std::endl;
         exit(1);
      }
      osLogShareDataPoolStat = &ofsLogShareDataPoolStat;
      *osLogShareDataPoolStat << shareDataPool->toStatStringHeader() << std::endl;
   }

   ///
   /// Setting of osLogSimilarity
   ///
   {
      std::ostringstream s;
      s << paraParams->getStringParamValue(UG::LogSolvingStatusFilePath)
         << paraInitiator->getParaInstance()->getProbName() << "_LC" << paraComm->getRank() << ".similarity.csv";
      ofsLogSimilarity.open(s.str().c_str(), std::ios::app );
      if( !ofsLogSimilarity )
      {
         std::cout << "Tabular vector pool log file cannot open : file name = " << s.str() << std::endl;
         exit(1);
      }
      osLogSimilarity = &ofsLogSimilarity;
   }


   ///
   /// Setting of osLogNotificationInterval
   ///
   {
      std::ostringstream s;
      s << paraParams->getStringParamValue(UG::LogSolvingStatusFilePath)
         << paraInitiator->getParaInstance()->getProbName() << "_LC" << paraComm->getRank() << ".notificationInterval.csv";
      ofsLogNotificationInterval.open(s.str().c_str(), std::ios::app );
      if( !ofsLogNotificationInterval )
      {
         std::cout << "Tabular vector pool log file cannot open : file name = " << s.str() << std::endl;
         exit(1);
      }
      osLogNotificationInterval = &ofsLogNotificationInterval;
   }
   *osLogNotificationInterval
      << Logging::getNotificationIntervalLogHeader() << ","
      << lcts.toStringLocalProcessTimesOfMessageHandlerHeader(paraComm, ",")
      << std::endl;


}


void *
CMapLapParaLoadCoordinator::runSolverThread(
      void *threadData
      )
{
   std::unique_ptr<LocalSolverThreadData> solverThreadData{static_cast<LocalSolverThreadData *>(threadData)};

   assert( solverThreadData->rank >= 1 );
   UG::ParaParamSet *inParaParams =  solverThreadData->paraParams;
   // Note :: In the local solvers, (rank, threaID) = (threaID, 1)
   lcLocalComm->slaveThreadInit(solverThreadData->rank, inParaParams);

   int argc = solverThreadData->argc;
   char **argv = solverThreadData->argv;

   sigset_t sigsBlock;
   sigemptyset(&sigsBlock);
   sigaddset(&sigsBlock, SIGINT);
   pthread_sigmask(SIG_BLOCK, &sigsBlock, NULL);

   std::unique_ptr<UG::ParaSolver> paraSolver{new CMapLapParaSolver(
         argc, argv, 0,
         lcLocalComm,
         inParaParams,
         new CMapLapParaInstanceTh( *(dynamic_cast<CMapLapParaInstance *>(solverThreadData->instance)) ),
         solverThreadData->paraDetTimer,
         solverThreadData->paraTimer->getElapsedTime()
         )};

   if( inParaParams->getBoolParamValue(UG::StatisticsToStdout) )
   {
      lcLocalComm->lockApp();
      std::cout << "After Rank " << lcLocalComm->getRank() << "thread Id" << lcLocalComm->getStartTime()
            << " Solver initialized 1: " << solverThreadData->paraTimer->getElapsedTime() << std::endl;
      lcLocalComm->unlockApp();
   }

   paraSolver->run();

   return 0;
}



void *
CMapLapParaLoadCoordinator::runCheckpointThread(
      void *threadData
      )
{
   std::unique_ptr<LocalSolverThreadData> solverThreadData{static_cast<LocalSolverThreadData *>(threadData)};

   assert( solverThreadData->rank >= 1 );
   UG::ParaParamSet *inParaParams =  solverThreadData->paraParams;
   // Note :: In the checkpoint writer, (rank, threaID) = (threaID, 1)

   lcCheckpointComm->slaveThreadInit(solverThreadData->rank, inParaParams);

   sigset_t sigsBlock;
   sigemptyset(&sigsBlock);
   sigaddset(&sigsBlock, SIGINT);
   pthread_sigmask(SIG_BLOCK, &sigsBlock, NULL);

   std::unique_ptr<CMapLapParaCheckpointWriter> paraCheckpointWriter{new CMapLapParaCheckpointWriter(
         lcCheckpointComm,
         inParaParams
         )};

   paraCheckpointWriter->run();

   return 0;
}



int
CMapLapParaLoadCoordinator::processTagSolution(
      int source,                                      ///< source solver rank
      int tag                                          ///< UG::TagSolution
      )
{
   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm);

   std::shared_ptr<CMapLapParaSolution> sol{dynamic_cast<CMapLapParaSolution *>(paraComm->createParaSolution())};
   sol->receive(paraComm, source);

   SolverType solverType = cmapLapParaSolverPool->getSolverType(source, sol->getThreadId());
   if( cmapLapParaInitiator->tryToSetIncumbentSolution(sol) )
   {
      // save incumbent solution
      char solutionFileName[256];
      sprintf(solutionFileName,"%s%s_after_checkpointing_solution.gz",
            paraParams->getStringParamValue(UG::CheckpointFilePath),
            paraInitiator->getParaInstance()->getProbName());
      paraInitiator->writeCheckpointSolution(std::string(solutionFileName));

      if( logSolvingStatusFlag )
      {
         incrementNStatus();

         int threadId = sol->getThreadId();
         double shortestNorm = std::sqrt(sol->getObjectiveFunctionValue());
         *osLogSolvingStatus
            << Logging::getSolverStateString('*', paraTimer->getElapsedTime(), source, threadId, solverType, shortestNorm)
            << std::endl;
         *osCsvLogSolvingStatus
            << Logging::getSolverStateString('*', paraTimer->getElapsedTime(), source, threadId, solverType, shortestNorm, ",")
            << std::endl;
      }
      outputAssignmentTable('*');
      outputCsvAssignmentTable('*');
   }
   else
   {
      std::shared_ptr<CMapLapParaSolution> incumbent(
         dynamic_cast<CMapLapParaSolution *>(
            cmapLapParaInitiator->getGlobalBestIncumbentSolution()->clone(cmapLapParaComm)
      ));

      int threadId = sol->getThreadId();
      double shortestNorm = std::sqrt(incumbent->getObjectiveFunctionValue());
      if( logSolvingStatusFlag )
      {
         *osLogSolvingStatus
            << Logging::getSolverStateString('U', paraTimer->getElapsedTime(), source, threadId, solverType, shortestNorm)
            << std::endl;
         *osCsvLogSolvingStatus
            << Logging::getSolverStateString('U', paraTimer->getElapsedTime(), source, threadId, solverType, shortestNorm, ",")
            << std::endl;
      }
      incumbent->setThreadId(threadId);
      incumbent->send(cmapLapParaComm, source);
   }

   // insert ShareDataPool
   unsigned int solverThreadId = cmapLapParaSolverPool->getSolverThreadId(source, sol->getThreadId());
   std::shared_ptr<VectorElement> vectorElement(new VectorElement(sol->getVector(), solverThreadId, solverType));

   shareDataPool->insert(vectorElement);
   if( paraParams->getIntParamValue(NumOfLocalSolversInLC) > 0 )
   {
      std::shared_ptr<CMapLapParaTask> taskForLocalSolver(createParaTaskForLocalSolver());
      if( taskForLocalSolver )
      {
         paraLocalTaskPool->insert(taskForLocalSolver);
      }
   }

   return 0;
}


int
CMapLapParaLoadCoordinator::processLocalTagSolution(
      int source,                                      ///< source solver rank
      int tag                                          ///< UG::TagSolution
      )
{
   std::shared_ptr<CMapLapParaSolution> sol{dynamic_cast<CMapLapParaSolution *>(lcLocalComm->createParaSolution())};
   sol->receive(lcLocalComm, source);

   if( cmapLapParaInitiator->tryToSetIncumbentSolution(sol) )
   {
      // save incumbent solution
      char solutionFileName[256];
      sprintf(solutionFileName,"%s%s_after_checkpointing_solution.gz",
            paraParams->getStringParamValue(UG::CheckpointFilePath),
            paraInitiator->getParaInstance()->getProbName());
      paraInitiator->writeCheckpointSolution(std::string(solutionFileName));

      if( logSolvingStatusFlag )
      {
         incrementNStatus();

         // output log like solver state log
         assert( cmapLapParaSolverPool->getSolverType(source,0) == DeepBkz );
         int inSource = 0;
         int inThreadId = sol->getThreadId();
         double shortestNorm = std::sqrt(sol->getObjectiveFunctionValue());
         *osLogSolvingStatus
            << Logging::getSolverStateString('*', paraTimer->getElapsedTime(), inSource, inThreadId, DeepBkz, shortestNorm)
            << std::endl;
         *osCsvLogSolvingStatus
            << Logging::getSolverStateString('*', paraTimer->getElapsedTime(), inSource, inThreadId, DeepBkz, shortestNorm, ",")
            << std::endl;
      }
      outputAssignmentTable('*');
      outputCsvAssignmentTable('*');
   }

   // insert ShareDataPool
   unsigned int solverThreadId = cmapLapParaSolverPool->getSolverThreadId(source, sol->getThreadId());
   std::shared_ptr<VectorElement> vectorElement(new VectorElement(sol->getVector(), solverThreadId));

   shareDataPool->insert(vectorElement);

   if( paraParams->getIntParamValue(NumOfLocalSolversInLC) > 0 )
   {
      std::shared_ptr<CMapLapParaTask> taskForLocalSolver(createParaTaskForLocalSolver());
      if( taskForLocalSolver )
      {
         paraLocalTaskPool->insert(taskForLocalSolver);
      }
   }

   return 0;
}


int
CMapLapParaLoadCoordinator::processTagBasis(
      int source,                                      ///< source solver rank
      int tag                                          ///< TagBasis
      )
{
   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm);

   std::unique_ptr<CMapLapParaBasis> paraBasis{cmapLapParaComm->createCMapLapParaBasis()};
   paraBasis->receive(paraComm, source);
   instancePool->insert(paraBasis->getEnumCost(), paraBasis->getBasis());
   unsigned int solverThreadId = cmapLapParaSolverPool->getSolverThreadId(source, paraBasis->getThreadId());
   hasSentBasisRequest[solverThreadId] = false;
   return 0;
}


int
CMapLapParaLoadCoordinator::processLocalTagBasis(
      int source,                                      ///< source solver rank
      int tag                                          ///< TagBasis
      )
{
   std::unique_ptr<CMapLapParaBasis> paraBasis{lcLocalComm->createCMapLapParaBasis()};
   paraBasis->receive(lcLocalComm, source);
   instancePool->insert(paraBasis->getEnumCost(), paraBasis->getBasis());
   unsigned int solverThreadId = cmapLapParaSolverPool->getSolverThreadId(source, paraBasis->getThreadId());
   hasSentBasisRequest[solverThreadId] = false;
   return 0;
}


int
CMapLapParaLoadCoordinator::processTagVectorRequest(
      int source,                                      ///< source solver rank
      int tag                                          ///< TagVectorRequest
      )
{
   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm);

   std::vector<int> receiveVectorsData(2);
   PARA_COMM_CALL(
         paraComm->receive(&receiveVectorsData[0], 2, UG::ParaINT, source, TagVectorRequest)
         );
   int nSendVectors = receiveVectorsData[0];
   int threadId = receiveVectorsData[1];

   if( nSendVectors > 0 )
   {
      unsigned int solverThreadId = cmapLapParaSolverPool->getSolverThreadId(source, threadId);
      if( !shareDataPool->empty() )
      {
         std::shared_ptr<CMapLapParaPackedVector> packedVector = nullptr;
         std::vector<std::shared_ptr<VectorElement>> vectorElements(nSendVectors);
         int nCopiedElement = 0;
         if( lcLocalComm )
            nCopiedElement = shareDataPool->getCopiedElement(nSendVectors, solverThreadId, vectorElements, lcLocalComm);
         else if( lcCheckpointComm )
            nCopiedElement = shareDataPool->getCopiedElement(nSendVectors, solverThreadId, vectorElements, lcCheckpointComm);
         else
            nCopiedElement = shareDataPool->getCopiedElement(nSendVectors, solverThreadId, vectorElements, nullptr);
         if( nCopiedElement > 0 )
            packedVector = std::shared_ptr<CMapLapParaPackedVector>(
               cmapLapParaComm->createCMapLapParaPackedVector(threadId, vectorElements)
            );
         else
            packedVector = std::shared_ptr<CMapLapParaPackedVector>(
               cmapLapParaComm->createCMapLapParaPackedVector(threadId, 0)
            );
         packedVector->send(paraComm, source);
      }
   }

   return 0;
}


int
CMapLapParaLoadCoordinator::processLocalTagVectorRequest(
      int source,                                      ///< source solver rank
      int tag                                          ///< TagVectorRequest
      )
{
   std::vector<int> receiveVectorsData(2);
   PARA_COMM_CALL(
         lcLocalComm->receive(&receiveVectorsData[0], 2, UG::ParaINT, source, TagVectorRequest)
         );
   int nSendVectors = receiveVectorsData[0];
   int threadId = receiveVectorsData[1];

   if( nSendVectors > 0 )
   {
      unsigned int solverThreadId = cmapLapParaSolverPool->getSolverThreadId(source, threadId);
      if( !shareDataPool->empty() )
      {
         std::shared_ptr<CMapLapParaPackedVector> packedVector = nullptr;
         std::vector<std::shared_ptr<VectorElement>> vectorElements(nSendVectors);
         int nCopiedElement = 0;
         if( lcLocalComm )
            nCopiedElement = shareDataPool->getCopiedElement(nSendVectors, solverThreadId, vectorElements, lcLocalComm);
         else if( lcCheckpointComm )
            nCopiedElement = shareDataPool->getCopiedElement(nSendVectors, solverThreadId, vectorElements, lcCheckpointComm);
         else
            nCopiedElement = shareDataPool->getCopiedElement(nSendVectors, solverThreadId, vectorElements, nullptr);
         if( nCopiedElement > 0 )
            packedVector = std::shared_ptr<CMapLapParaPackedVector>(
               lcLocalComm->createCMapLapParaPackedVector(threadId, vectorElements)
            );
         else
            packedVector = std::shared_ptr<CMapLapParaPackedVector>(
               lcLocalComm->createCMapLapParaPackedVector(threadId, 0)
            );
         packedVector->send(lcLocalComm, source);
      }
   }

   return 0;
}


int
CMapLapParaLoadCoordinator::processTagSolverState(
      int source,                                      ///< source solver rank
      int tag                                          ///< UG::TagSolverState
      )
{
   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm);

   std::unique_ptr<CMapLapParaSolverState> solverState{dynamic_cast<CMapLapParaSolverState *>(cmapLapParaComm->createParaSolverState())};
   solverState->receive(cmapLapParaComm, source, tag);

   if( paraDetTimer
         && paraDetTimer->getElapsedTime() < solverState->getDeterministicTime() )
   {
      paraDetTimer->update( solverState->getDeterministicTime() - paraDetTimer->getElapsedTime() );
   }

   SolverType solverType = solverState->getSolverType();

   // for deepbkz
   int *basisArray;
   LatticeBasis<int> basis;

   // get basis
   if( solverType == DeepBkz )
   {
      std::shared_ptr<LatticeBasis<int>> prevBasis = cmapLapParaSolverPool->getBasis(
            source,
            solverState->getThreadId()
            );
      int n = cmapLapParaInitiator->getDimension();
      basisArray = solverState->getBasis();
      basis.setZero(n, n);
      int index = 0;
      for(int row = 0; row < n; ++row)
      {
         for(int col = 0; col < n; ++col)
         {
            basis(row, col) = basisArray[index++];
         }
      }
      int nNotChangeRows = 0;
      for( int i = 0; i < n; i++ )
      {
         if( basis.row(i) != prevBasis->row(i) ) break;
         nNotChangeRows++;
      }
      solverState->setNNotChangeRows(nNotChangeRows);
   }


   // update solver state
   switch( solverType )
   {
      case DeepBkz:
      {
         cmapLapParaSolverPool->updateSolverStatus(
               source,
               solverState->getThreadId(),
               solverState->getCurrentBlockSize(),
               solverState->getBasis(),
               solverState->getEnumCost(),
               solverState->getSlopeGSA(),
               solverState->getTopHalfSlopeGSA(),
               solverState->getOrthogonalFactor()
               );
         instancePool->insert(solverState->getEnumCost(), basis);
         break;
      }
      case Enum:
      {
         cmapLapParaSolverPool->updateSolverStatus(
               source,
               solverState->getThreadId(),
               solverState->getSearchNodeIndex(),
               solverState->getCoeffs(),
               solverState->getNumSearchedNodes()
               );
         break;
      }
      case Sieve:
      {
         break;
      }
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaLoadCoordinator::activateSolver: Invalid solver type = ",static_cast<int>(solverType));
   }

   if( cmapLapParaInitiator->getGlobalBestIncumbentSolution()
         && EPSLT(cmapLapParaInitiator->getGlobalBestIncumbentSolution()->getObjectiveFunctionValue(), solverState->getObjectiveFunctionValue(), DEFAULT_NUM_EPSILON)
            )
   {
      std::shared_ptr<CMapLapParaSolution> solution(
         dynamic_cast<CMapLapParaSolution *>(
            cmapLapParaInitiator->getGlobalBestIncumbentSolution()->clone(cmapLapParaComm)
            ));
      solution->setThreadId(solverState->getThreadId());
      solution->send(cmapLapParaComm, source);
   }


   if( logSolvingStatusFlag )
   {
      incrementNStatus();

      int threadId = solverState->getThreadId();
      *osLogSolvingStatus
         << Logging::getSolverStateString(' ', paraTimer->getElapsedTime(), source, threadId, solverState->toStringLog())
         << std::endl;
      *osCsvLogSolvingStatus
         << Logging::getSolverStateString(' ', paraTimer->getElapsedTime(), source, threadId, solverState->toStringLog(","), ",")
         << std::endl;
   }

   if( !paraParams->getBoolParamValue(NoWaitNotificationId) )
   {
      unsigned int notificationId[2];   /// notificationId[0]: notification ID
                                        /// notificationId[1]: thread ID
      notificationId[0] = solverState->getNotificaionId();
      notificationId[1] = solverState->getThreadId();
      PARA_COMM_CALL(
            cmapLapParaComm->send( notificationId, 2, UG::ParaUNSIGNED, source, UG::TagNotificationId)
            );
#ifdef _COMM_MPI_WORLD
      cmapLapParaComm->waitAllIsends();
#endif
   }
   return 0;

}


int
CMapLapParaLoadCoordinator::processTagCompletionOfCalculation(
      int source,
      int tag
      )
{
   std::unique_ptr<CMapLapParaCalculationState> calcState{dynamic_cast<CMapLapParaCalculationState *>(paraComm->createParaCalculationState())};
   calcState->receive(paraComm, source, tag);

   /// erase this terminated solver thread id from set of sent thread ids
   /// of all element in vector element pool
   unsigned int solverThreadId = cmapLapParaSolverPool->getSolverThreadId(source, calcState->getThreadId());
   shareDataPool->eraseSolverThreadId(solverThreadId, lcLocalComm);
   hasSentBasisRequest[solverThreadId] = false;

   if( logSolvingStatusFlag )
   {
      int threadId = calcState->getThreadId();
      *osLogSolvingStatus
         << Logging::getSolverStateString('<', paraTimer->getElapsedTime(), source, threadId, calcState->toStopString())
         << std::endl;
      *osCsvLogSolvingStatus
         << Logging::getSolverStateString('<', paraTimer->getElapsedTime(), source, threadId, calcState->toStopString(","), ",")
         << std::endl;
   }

   if( osStatisticsFinalRun )
   {
      *osStatisticsFinalRun << paraTimer->getElapsedTime()
            << " [S." << source << ", T." << calcState->getThreadId() << "] "
            << calcState->toStopString() << std::endl;
   }

   switch( calcState->getSolverType() )
   {
   case DeepBkz:
      cmapLapParaSolverPool->inactivateSolver(source, calcState->getThreadId(), paraDeepBkzTaskPool);
      break;
   case Enum:
      cmapLapParaSolverPool->inactivateSolver(source, calcState->getThreadId(), paraEnumTaskPool);
      break;
   case Sieve:
      cmapLapParaSolverPool->inactivateSolver(source, calcState->getThreadId(), paraSieveTaskPool);
      break;
   default:
      THROW_LOGICAL_ERROR2(
            "CMapLapParaLoadCoordinator::processTagCompletionOfCalculation: Invalid solver type = ",
            static_cast<int>(calcState->getSolverType()));
   }

   return 0;
}


///
/// processTagTerminated
///
int
CMapLapParaLoadCoordinator::processTagTerminated(
      int source,
      int tag
      )
{
   std::unique_ptr<CMapLapParaSolverTerminationState> termState{
      dynamic_cast<CMapLapParaSolverTerminationState *>(paraComm->createParaSolverTerminationState())};
   termState->receive(paraComm, source, tag);

   /// erase this terminated solver thread id from set of sent thread ids
   /// of all element in vector element pool
   unsigned int solverThreadId = cmapLapParaSolverPool->getSolverThreadId(source, termState->getThreadId());
   shareDataPool->eraseSolverThreadId(solverThreadId, lcLocalComm);

   if( paraDetTimer )
   {
      if( paraDetTimer->getElapsedTime() < termState->getDeterministicTime() )
      {
         paraDetTimer->update( termState->getDeterministicTime() - paraDetTimer->getElapsedTime() );
      }
      PARA_COMM_CALL(
            paraComm->send( NULL, 0, UG::ParaBYTE, source, UG::TagAckCompletion )
      );
   }

   if( osStatisticsFinalRun )
   {
      *osStatisticsFinalRun << termState->toString(paraInitiator);
      osStatisticsFinalRun->flush();
   }
   if( paraParams->getBoolParamValue(UG::StatisticsToStdout) )
   {
      std::cout << termState->toString(paraInitiator) << std::endl;
   }

   if( (!racingTermination) && termState->getInterruptedMode() == 1 )
   {
      computationIsInterrupted = true;
   }

   nTerminated++;
   return 0;
}

int
CMapLapParaLoadCoordinator::processLocalTagCompletionOfCalculation(
      int source,
      int tag
      )
{
   CMapLapParaCalculationState *calcState = dynamic_cast<CMapLapParaCalculationState *>(lcLocalComm->createParaCalculationState());
   calcState->receive(lcLocalComm, source, tag);

   if( logSolvingStatusFlag )
   {
      int inSource = 0;
      int inThreadId = calcState->getThreadId();
      *osLogSolvingStatus
         << Logging::getSolverStateString('^', paraTimer->getElapsedTime(), inSource, inThreadId, calcState->toStopString())
         << std::endl;
      *osCsvLogSolvingStatus
         << Logging::getSolverStateString('^', paraTimer->getElapsedTime(), inSource, inThreadId, calcState->toStopString(","), ",")
         << std::endl;
   }

   if( osStatisticsFinalRun )
   {
      *osStatisticsFinalRun << paraTimer->getElapsedTime()
            << " [S." << 0 << ", T." << calcState->getThreadId() << "] "
            << calcState->toStopString() << std::endl;
   }

   // paraLocalSolverPool->inactivateSolver(source, calcState->getThreadId(), paraLocalTaskPool);
   paraLocalSolverPool->inactivateSolver(source, 0, paraLocalTaskPool);

   return 0;
}

int
CMapLapParaLoadCoordinator::processTagCMapLapPackedVector(
      int source,
      int tag
      )
{
   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm);
   std::unique_ptr<CMapLapParaPackedVector> packedVector{cmapLapParaComm->createCMapLapParaPackedVector()};
   packedVector->receive(paraComm, source);
   SolverType solverType = cmapLapParaSolverPool->getSolverType(source, packedVector->getThreadId());
   while( true )
   {
      auto vectorElement = packedVector->extract();
      if( !vectorElement ) break;
      vectorElement->setSolverType(solverType);
      shareDataPool->insert(vectorElement);
   }
   return 0;
}

int
CMapLapParaLoadCoordinator::processLocalTagCMapLapPackedVector(
      int source,
      int tag
      )
{
   std::unique_ptr<CMapLapParaPackedVector> packedVector{lcLocalComm->createCMapLapParaPackedVector()};
   packedVector->receive(lcLocalComm, source);
   SolverType solverType = cmapLapParaSolverPool->getSolverType(source, packedVector->getThreadId());
   while( true )
   {
      auto vectorElement = packedVector->extract();
      if( !vectorElement ) break;
      vectorElement->setSolverType(solverType);
      shareDataPool->insert(vectorElement);
   }
   return 0;
}

int
CMapLapParaLoadCoordinator::processLocalTagTerminated(
      int source,
      int tag
      )
{
   std::unique_ptr<UG::ParaSolverTerminationState> paraSolverTerminationState{lcLocalComm->createParaSolverTerminationState()};
   paraSolverTerminationState->receive(lcLocalComm, source, tag);

   if( paraDetTimer )
   {
      if( paraDetTimer->getElapsedTime() < paraSolverTerminationState->getDeterministicTime() )
      {
         paraDetTimer->update( paraSolverTerminationState->getDeterministicTime() - paraDetTimer->getElapsedTime() );
      }
      PARA_COMM_CALL(
            lcLocalComm->send( NULL, 0, UG::ParaBYTE, source, UG::TagAckCompletion )
      );
   }

   if( osStatisticsFinalRun )
   {
      *osStatisticsFinalRun << paraSolverTerminationState->toString(paraInitiator);
      osStatisticsFinalRun->flush();
   }
   if( paraParams->getBoolParamValue(UG::StatisticsToStdout) )
   {
      std::cout << paraSolverTerminationState->toString(paraInitiator) << std::endl;
   }

   if( (!racingTermination) && paraSolverTerminationState->getInterruptedMode() == 1 )
   {
      computationIsInterrupted = true;
   }

   nTerminated++;

   return 0;
}

int
CMapLapParaLoadCoordinator::processLocalTagHardTimeLimit(
      int source,
      int tag
      )
{
   PARA_COMM_CALL(
         lcLocalComm->receive( NULL, 0, UG::ParaBYTE, source, UG::TagHardTimeLimit)
         );
   hardTimeLimitIsReached = true;
   return 0;
}


void
CMapLapParaLoadCoordinator::outputAssignmentTable(
      char shortest
      )
{
   if( shortest != ' ' ||
         paraTimer->getElapsedTime() - previousAssingmentTableOutputTime
         > paraParams->getRealParamValue(IntervalTimeOfAssignmentTableOutput) )
   {
      *osAssignmentTable << std::setw(1) << shortest;
      *osAssignmentTable << std::setw(6) << Logging::getReadableTimeString(paraTimer->getElapsedTime());
      if( cmapLapParaInitiator->getGlobalBestIncumbentSolution() )
      {
         *osAssignmentTable << std::setw(10) << std::right << std::setprecision(3) << std::fixed
               << std::sqrt(cmapLapParaInitiator->getGlobalBestIncumbentSolution()->getObjectiveFunctionValue());
      }
      else
      {
         *osAssignmentTable << std::setw(9) << std::right << std::setprecision(0) << std::fixed << 0;
      }
      if( cmapLapParaInitiator->getGlobalBestIncumbentSolution() )
      {
         *osAssignmentTable << std::setw(8) << std::right << std::setprecision(4) << std::fixed
               << std::sqrt(cmapLapParaInitiator->getGlobalBestIncumbentSolution()->getObjectiveFunctionValue())
                / cmapLapParaInitiator->getGH();
      }
      else
      {
         *osAssignmentTable << std::setw(8) << std::right << std::setprecision(0) << std::fixed << 0;
      }
      *osAssignmentTable << std::setw(10) << std::right << std::setprecision(0) << std::fixed << cmapLapParaSolverPool->getNumInactiveSolvers();
      *osAssignmentTable << std::setw(13) << std::right << std::setprecision(0) << std::fixed << cmapLapParaSolverPool->getNumActiveSolvers();
      *osAssignmentTable << std::setw(9) << std::right << std::setprecision(0) << std::fixed << cmapLapParaSolverPool->getNumActiveDeepBkzSolvers()
                         << "(+" << ( (paraLocalSolverPool) ? paraLocalSolverPool->getNumActiveDeepBkzSolvers() : 0 ) << ")";
      *osAssignmentTable << std::setw(9) << std::right << std::setprecision(0) << std::fixed << cmapLapParaSolverPool->getNumActiveEnumSolvers()
                         << "(+" << ( (paraLocalSolverPool) ? paraLocalSolverPool->getNumActiveEnumSolvers() : 0 ) << ")";
      *osAssignmentTable << std::setw(10) << std::right << std::setprecision(0) << std::fixed << cmapLapParaSolverPool->getNumActiveSieveSolvers()
                         << "[" << cmapLapParaSolverPool->getNumActiveSieveSolvers()*nThreadsPerRank << "]";
      int nTaskSetw = 7;
      if( paraLocalTaskPool )
      {
         if( paraLocalTaskPool->getNumOfTasks() >= 10 )   nTaskSetw -= 1;
         else if( paraLocalTaskPool->getNumOfTasks() >= 100 )  nTaskSetw -= 2;
         else if( paraLocalTaskPool->getNumOfTasks() >= 1000 ) nTaskSetw -= 3;
      }
      // number of tatsks ( DeepBkz + Enum ( Local Task ) )
      *osAssignmentTable << std::setw(nTaskSetw)  << std::right << std::setprecision(0) << std::fixed
         << paraDeepBkzTaskPool->getNumOfTasks() + paraEnumTaskPool->getNumOfTasks() + paraSieveTaskPool->getNumOfTasks()
         << "(+" << ( (paraLocalTaskPool) ? paraLocalTaskPool->getNumOfTasks() : 0 ) << ")";
      *osAssignmentTable << std::setw(12) << std::right << std::setprecision(0) << std::fixed << instancePool->size();
      *osAssignmentTable << std::setw(10) << std::right << std::setprecision(0) << std::fixed << shareDataPool->size();
      *osAssignmentTable << std::endl;
      previousAssingmentTableOutputTime = paraTimer->getElapsedTime();
   }
}


void
CMapLapParaLoadCoordinator::outputCsvAssignmentTable(
      char shortest
      )
{
   double shortestNorm = -1.0;
   if( cmapLapParaInitiator->getGlobalBestIncumbentSolution() )
   {
      shortestNorm = std::sqrt(cmapLapParaInitiator->getGlobalBestIncumbentSolution()->getObjectiveFunctionValue());
   }

   *osCsvAssignmentTable
      << Logging::getCsvAssignmentTableString(
         shortest,
         paraTimer->getElapsedTime(),
         shortestNorm,
         cmapLapParaInitiator->getLattice().approxFactor(shortestNorm),
         cmapLapParaInitiator->getLattice().hermiteFactor(shortestNorm),
         cmapLapParaInitiator->getLattice().rootHermiteFactor(shortestNorm),
         cmapLapParaSolverPool->getNumInactiveSolvers(),
         cmapLapParaSolverPool->getNumActiveSolvers(),
         cmapLapParaSolverPool->getNumActiveDeepBkzSolvers(),
         cmapLapParaSolverPool->getNumActiveEnumSolvers(),
         cmapLapParaSolverPool->getNumActiveSieveSolvers(),
         paraDeepBkzTaskPool->getNumOfTasks(),
         paraEnumTaskPool->getNumOfTasks(),
         paraSieveTaskPool->getNumOfTasks(),
         instancePool->size(),
         shareDataPool->size()
         )
      << std::endl;
}



///
/// @return true if it should create and assign task
///
bool
CMapLapParaLoadCoordinator::shouldCreateAndAssignTask(
      SolverType &solverType,      ///< solver-task should be created
      bool &hasSetInitialSolvers,  ///< if it is false, initial solver setting is not completed
      bool &hasSetInitialDeepBkzSolvers
      )
{
   if( nSolvers == 1 )
   {
      hasSetInitialDeepBkzSolvers = true;
   }
   if( hardTimeLimitIsReached || lowerBoundIsReached )
   {
      return false;
   }
   if( !hasSetInitialDeepBkzSolvers )
   {
      // note: before assign Sieve task, must assign DeepBkz or Enum task
      if( static_cast<int>(cmapLapParaSolverPool->getNumActiveDeepBkzSolvers()) < nSolvers )
      {
         solverType = DeepBkz;
         return true;
      }
      hasSetInitialDeepBkzSolvers = true;
      if( paraParams->getIntParamValue(NumOfInitialEnumSolvers) > 0 )
      {
         cmapLapParaSolverPool->interruptDeepBkzSolvers(
               paraParams->getIntParamValue(NumOfInitialEnumSolvers));
      }
   }

   int numOfInitialEnumSolvers    = paraParams->getIntParamValue(NumOfInitialEnumSolvers);
   int numOfInitialSieveSolvers   = paraParams->getIntParamValue(NumOfInitialSieveSolvers);
   if( numOfInitialEnumSolvers < 0 ){ numOfInitialEnumSolvers = nSolvers; };
   if( numOfInitialSieveSolvers < 0 ){ numOfInitialSieveSolvers = nSolvers; };
   int numOfActiveEnumSolvers     = cmapLapParaSolverPool->getNumActiveEnumSolvers();
   int numOfActiveSieveSolvers    = cmapLapParaSolverPool->getNumActiveSieveSolvers();

   if( numOfInitialEnumSolvers > numOfActiveEnumSolvers )
   {
      solverType = Enum;
      return true;
   }
   else if( numOfInitialSieveSolvers > numOfActiveSieveSolvers )
   {
      solverType = Sieve;
      return true;
   }
   else if( cmapLapParaSolverPool->getNumInactiveSolvers() > 0 )
   {
      solverType = DeepBkz;
      return true;
   }
   return false;
}


///
/// check the terminate condition
/// @return true if it is satisfy the terminate condition else false
///
bool
CMapLapParaLoadCoordinator::shouldBreakRun(
      )
{
   if( !isAssigningParaTask() && (!hardTimeLimitIsReached && !lowerBoundIsReached) )
      return false;
   else if( paraSolverPool->getNumActiveSolvers() > 0 )
      return false;
   else
      return true;
}


///
/// interrupt the required number of active solvers and return paraTask of solverType
/// @return CMapLapParaTask pointer
///
std::shared_ptr<CMapLapParaTask>
CMapLapParaLoadCoordinator::getParaTask(
      SolverType solverType,  ///< solver task should be created
      int begin,
      int end
      )
{
   switch( solverType )
   {
   case DeepBkz:
   {
      if( cmapLapParaSolverPool->getNumInactiveSolvers() == 0 )
      {
         if( !cmapLapParaSolverPool->isThereInterruptingDeepBkzSolver() )
         {
            int nInterruptSolvers = 1;
            cmapLapParaSolverPool->interruptDeepBkzSolvers(nInterruptSolvers);
         }
      }
      if( paraDeepBkzTaskPool->empty() )
      {
         createParaTaskDeepBkz(begin, end);
      }
      return paraDeepBkzTaskPool->extractTask();
   }
   case Enum:
   {
      if( cmapLapParaSolverPool->getNumInactiveSolvers() == 0 &&
            !cmapLapParaSolverPool->isThereInterruptingDeepBkzSolver() )
      {
         int nInterruptSolvers = 1;
         cmapLapParaSolverPool->interruptDeepBkzSolvers(nInterruptSolvers);
      }
      if( instancePool->size() == 0 )
      {
         cmapLapParaSolverPool->sendMessageForSolvers(lcLocalComm, TagBasisRequest, hasSentBasisRequest, DeepBkz);
         sleep(std::floor(paraParams->getRealParamValue(IReceiveInterval)));
         return nullptr;
      }
      if( paraEnumTaskPool->empty() )
      {
         createParaTaskEnums(begin, end);
      }
      return paraEnumTaskPool->extractTask();
   }
   case Sieve:
   {
      if( cmapLapParaSolverPool->getNumInactiveSolvers() == 0 &&
            !cmapLapParaSolverPool->isThereInterruptingDeepBkzSolver() )
      {
         cmapLapParaSolverPool->interruptDeepBkzAndEnumSolversToRunSieveSolver();
      }
      if( instancePool->size() == 0 )
      {
         cmapLapParaSolverPool->sendMessageForSolvers(lcLocalComm, TagBasisRequest, hasSentBasisRequest, DeepBkz);
         sleep(std::floor(paraParams->getRealParamValue(IReceiveInterval)));
         return nullptr;
      }
      if( paraSieveTaskPool->empty() )
      {
         createParaTaskSieve(begin, end);
      }
      return paraSieveTaskPool->extractTask();
   }
   default:
      THROW_LOGICAL_ERROR2("CMapLapParaLoadCoordinator::getParaTask: Invalid solver type = ", static_cast<int>(solverType));
   }
}

///
/// @return CMapLapParaTask DeepBkz pointer
///
int
CMapLapParaLoadCoordinator::createParaTaskDeepBkz(
      int begin,
      int end
      )
{
   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm );

   std::shared_ptr<LatticeBasis<int>> basis;
   if( !instancePool->empty() )
   {
      basis = instancePool->extractMinEnumCostBasis();
   }
   else
   {
      basis = std::shared_ptr<LatticeBasis<int>>(
            new LatticeBasis<int>(*(instancePool->getIncumbentBasis()->getBasis()))
            );
   }
   int u = paraParams->getIntParamValue(RandomizeRows);  // number of randomize rows
   int seed = seedGenerator();
   int blockSize = paraParams->getIntParamValue(DeepBkzStartBlockSize);
   double estimatedValue = -DBL_MAX;

   std::shared_ptr<CMapLapParaTask> cmapLapParaTask(
         dynamic_cast<CMapLapParaTask *>(
            cmapLapParaComm->createParaTask(
               UG::TaskId(), UG::TaskId(), estimatedValue, -1,
               begin, end, blockSize, u, seed, basis)));
   paraDeepBkzTaskPool->insert(cmapLapParaTask);
   return 1;
}


///
/// create paraTaskEnum using basis in InstancePool
///
int
CMapLapParaLoadCoordinator::createParaTaskEnums(
      int begin,
      int end
      )
{
   assert( instancePool->size() > 0 );

   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm );

   int nTaskEnums = 0;
   std::shared_ptr<LatticeBasis<int>> basis(instancePool->extractMinEnumCostBasis());

   // crate paraTaskEnum
   int start = begin;
   int last = end + 1;
   double prob = paraParams->getRealParamValue(EnumPruningParameter);
   auto coeffs = std::make_shared<LatticeVector<int>>(basis->rows());
   coeffs->setZero();

   Lattice<int, double> L(basis);
   double estimatedValue = L.enumCost(L.GH);

   std::shared_ptr<CMapLapParaTask> cmapLapParaTask(
         dynamic_cast<CMapLapParaTask *>(
            cmapLapParaComm->createParaTask(
               UG::TaskId(), UG::TaskId(), estimatedValue, -1,
               begin, end, start, last, prob, coeffs, basis)));

   // insert paraTask into paraTaskPool
   paraEnumTaskPool->insert(cmapLapParaTask);
   nTaskEnums++;
   return nTaskEnums;
}



///
/// create paraTaskSieve using basis in InstancePool
///
int
CMapLapParaLoadCoordinator::createParaTaskSieve(
      int begin,
      int end
      )
{
   assert( instancePool->size() > 0 );

   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm );

   std::shared_ptr<LatticeBasis<int>> basis(instancePool->extractMinEnumCostBasis());

   // create Sieve task
   double estimatedValue = -DBL_MAX;        // highest priority
   std::shared_ptr<CMapLapParaTask> cmapLapParaTask(
         dynamic_cast<CMapLapParaTask *>(
            cmapLapParaComm->createParaTask(
               UG::TaskId(), UG::TaskId(), estimatedValue, -1,
               begin, end, basis)));

   // insert paraTask into paraTaskPool
   paraSieveTaskPool->insert(cmapLapParaTask);
   return 1;
}

///
/// create a paraTask for local solvers in LoadCoordinator when a vector is inserted in vectorElement
/// @return CMapLapParaTask*
///
CMapLapParaTask*
CMapLapParaLoadCoordinator::createParaTaskForLocalSolver(
      )
{
   std::shared_ptr<LatticeBasis<int>> subBasis = shareDataPool->createSubBasis();
   if( subBasis->rows() < paraParams->getIntParamValue(LocalTaskNVectorLowerBound) ) return 0;

   int begin = 0;
   int end = subBasis->rows()-1;
   int u = 1;  // do not randomize
   int seed = seedGenerator();
   int blockSize = paraParams->getIntParamValue(BlockSizeForLocalSolver);
   double estimatedValue = -DBL_MAX;

   CMapLapParaTask* cmapLapParaTask = new CMapLapParaTaskLC(
         UG::TaskId(), UG::TaskId(), estimatedValue, -1,
         begin, end, blockSize, u, seed, subBasis);

   return cmapLapParaTask;
}

///
/// assign cmapLapParaTask
/// @return true if success to activate cmapLapParaTask else false
///
bool
CMapLapParaLoadCoordinator::assignParaTask(
      std::shared_ptr<CMapLapParaTask> cmapLapParaTask
      )
{
   assert( cmapLapParaTask );
   int threadId = -1;
   int destination = cmapLapParaSolverPool->activateSolver(cmapLapParaTask->getSolverType(), cmapLapParaTask, threadId );
   if( destination < 0 )
   {
      // cannot activate
      switch( cmapLapParaTask->getSolverType() )
      {
      case DeepBkz:
         paraDeepBkzTaskPool->insert(cmapLapParaTask);
         break;
      case Enum:
         paraEnumTaskPool->insert(cmapLapParaTask);
         break;
      case Sieve:
         // assert( assigningSieveSolver == true );
         paraSieveTaskPool->insert(cmapLapParaTask);
         break;
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaLoadCoordinator::run: Invalid solver type = ", static_cast<int>(cmapLapParaTask->getSolverType()));
      }
      return false;
   }
   else
   {
      lcts.nSent++;
      //  writeTransferLog(destination);
      if( logSolvingStatusFlag )
      {
         *osLogSolvingStatus
            << Logging::getSolverStateString('>', paraTimer->getElapsedTime(), destination, threadId, cmapLapParaTask->toStartString())
            << std::endl;
         *osCsvLogSolvingStatus
            << Logging::getSolverStateString('>', paraTimer->getElapsedTime(), destination, threadId, cmapLapParaTask->toStartString(","), ",")
            << std::endl;
      }
      assigningParaTask = false;
      return true;
   }
}

///
/// assign cmapLapParaLocalTask
/// @return true if success to assign cmapLapParaLocalTask else false
///
bool
CMapLapParaLoadCoordinator::assignParaLocalTask(
      std::shared_ptr<CMapLapParaTask> cmapLapParaLocalTask
      )
{
   int lcThreadId = -1;
   int lcDestination = paraLocalSolverPool->activateSolver(cmapLapParaLocalTask->getSolverType(), cmapLapParaLocalTask, lcThreadId );
   if( lcDestination < 0 )
   {
      paraLocalTaskPool->insert(cmapLapParaLocalTask);
      return false;
   }
   else
   {
      // lcts.nSent++;
      int source = 0;
      if( logSolvingStatusFlag )
      {
         *osLogSolvingStatus
            << Logging::getSolverStateString('v', paraTimer->getElapsedTime(), source, lcDestination, cmapLapParaLocalTask->toStartString())
            << std::endl;
         *osLogSolvingStatus
            << Logging::getSolverStateString('v', paraTimer->getElapsedTime(), source, lcDestination, cmapLapParaLocalTask->toStartString(","), ",")
            << std::endl;
      }
      return true;
   }
}

///
/// waiting for any message from anywhere
///
void
CMapLapParaLoadCoordinator::waitForAnyMessageFromAnyWhere(
      )
{
   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm );

   int source = -1;
   int tag = -1;
   double startTime = paraTimer->getElapsedTime();
   cmapLapParaComm->probe(&source, &tag);
   lcts.idleTime += ( paraTimer->getElapsedTime() - startTime );

   int status;
   if( messageHandler[tag] )
   {
      startTime = paraTimer->getElapsedTime();
      status = (this->*messageHandler[tag])(source, tag);
      if( status )
      {
         std::ostringstream s;
         s << "[ERROR RETURN form Message Hander]:" <<  __FILE__ <<  "] func = "
            << __func__ << ", line = " << __LINE__ << " - "
            << "process tag = " << tag << std::endl;
         abort();
      }
      lcts.addProcessTime(tag, paraTimer->getElapsedTime() - startTime);
   }
   else
   {
      THROW_LOGICAL_ERROR3( "No message hander for ", tag, " is not registered" );
   }
#ifdef _COMM_MPI_WORLD
   startTime = paraTimer->getElapsedTime();
   if( !hardTimeLimitIsReached && !lowerBoundIsReached )
   {
      cmapLapParaComm->testAllIsends();
   }
   else
   {
      cmapLapParaComm->waitAllIsends();
   }
   lcts.idleTimeToWaitIsend += ( paraTimer->getElapsedTime() - startTime );
#endif
}

///
/// waiting for any message from local
///
void
CMapLapParaLoadCoordinator::waitForAnyMessageFromLocal(
      )
{
   int lcSource = -1;
   int lcTag = -1;
   if( lcLocalComm->iProbe(&lcSource, &lcTag) )
   {
      int lcStatus;
      if( localMessageHandler[lcTag] )
      {
         lcStatus = (this->*localMessageHandler[lcTag])(lcSource, lcTag);
         if( lcStatus )
         {
            std::ostringstream s;
            s << "[ERROR RETURN form local Message Hander]:" <<  __FILE__ <<  "] func = "
               << __func__ << ", line = " << __LINE__ << " - "
               << "process tag = " << lcTag << std::endl;
            abort();
         }
      }
      else
      {
         THROW_LOGICAL_ERROR3( "No local message hander for ", lcTag, " is not registered" );
      }
   }
}


///
/// run function to start main process
///
void
CMapLapParaLoadCoordinator::run(
      )
{
   /// Always full SVP Algorithm
   int begin = 0;
   int end   = cmapLapParaInitiator->getDimension() - 1;

   if( paraParams->getBoolParamValue(OutputSimilarityOfBasis) )
      outputSimilarityOfBasisHeader();

   double previousLogShareDataPoolTime = 0.0;
   double previousOutputSimilarityOfBasisTime = 0.0;
   double previousLoadCoordinatorStatisticsStdoutTime = 0.0;
   double nextOutputLctsStatTime = paraTimer->getElapsedTime();

   // for adjust solver's notification interval
   double previousTermUpdateNotificationIntervalTime = paraTimer->getElapsedTime();
   double notificationInterval = paraParams->getRealParamValue(UG::NotificationInterval);
   double leftNotificationInterval = notificationInterval;

   bool hasSetInitialSolvers = false;
   bool hasSetInitialDeepBkzSolvers = false;
   assigningParaTask = false;

   SolverType solverType = Undefined;
   outputAssignmentTable(' ');
   outputCsvAssignmentTable(' ');


   for(;;)
   {
      // Assign Task to idle solver
      while( shouldCreateAndAssignTask(solverType, hasSetInitialSolvers, hasSetInitialDeepBkzSolvers) )
      {
         std::shared_ptr<CMapLapParaTask> cmapLapParaTask = getParaTask(solverType, begin, end);
         if( !cmapLapParaTask ) break;
         assigningParaTask = assignParaTask(cmapLapParaTask);
         if( !assigningParaTask ) break;
      }

      // Wait message from solvers
      waitForAnyMessageFromAnyWhere();

      // Check to see if the time limit has been reached
      if( !hardTimeLimitIsReached &&
            !lowerBoundIsReached &&
            paraParams->getRealParamValue(UG::TimeLimit) > 0 &&
            paraTimer->getElapsedTime() > paraParams->getRealParamValue(UG::TimeLimit) )
      {
         hardTimeLimitIsReached = true;
         interruptIsRequested = true;
         cmapLapParaSolverPool->interruptAllSolvers(lcLocalComm, TagTimeLimitRequest);
         if( paraLocalSolverPool ) paraLocalSolverPool->interruptAllSolvers(lcLocalComm, TagTimeLimitRequest);
#ifdef UG_WITH_ZLIB
         updateCheckpointFiles();
#endif
      }

      // Check to see if the lower bound has been reached
      if( !hardTimeLimitIsReached &&
            !lowerBoundIsReached &&
            globalLowerBoundOfSquaredNorm > 0 &&
            getCmapLapParaInitiator()->getGlobalBestIncumbentSolution()->getObjectiveFunctionValue()
               < globalLowerBoundOfSquaredNorm )
      {
         lowerBoundIsReached = true;
         interruptIsRequested = true;
         cmapLapParaSolverPool->interruptAllSolvers(lcLocalComm, TagTimeLimitRequest);
         if( paraLocalSolverPool ) paraLocalSolverPool->interruptAllSolvers(lcLocalComm, TagTimeLimitRequest);
#ifdef UG_WITH_ZLIB
         updateCheckpointFiles();
#endif
      }

      // Assign Task to idle local solver
      if( paraParams->getIntParamValue(NumOfLocalSolversInLC) > 0 )
      {
         waitForAnyMessageFromLocal();
         if( !hardTimeLimitIsReached &&
               !lowerBoundIsReached &&
               !paraLocalTaskPool->empty() &&
               paraLocalSolverPool->getNumInactiveSolvers() > 0 )
         {
            std::shared_ptr<CMapLapParaTask> cmapLapParaLocalTask(paraLocalTaskPool->extractTask());
            assignParaLocalTask(cmapLapParaLocalTask);
         }
      }

      // Check it should break from this main process loop
      if( shouldBreakRun() )
      {
         break;
      }


#ifdef UG_WITH_ZLIB
      // Create checkpoint files
      if( paraParams->getBoolParamValue(UG::Checkpoint) &&
            ( paraTimer->getElapsedTime() - previousCheckpointTime )
            > paraParams->getRealParamValue(UG::CheckpointInterval) )
      {
         if( !interruptIsRequested )
         {
            updateCheckpointFiles();
         }
      }
#endif

      // Logging : Assignment table
      outputAssignmentTable(' ');

      // Logging : Statistics of Load Coordinator
      if( paraTimer->getElapsedTime() > nextOutputLctsStatTime )
      {
         *osCsvLctsStat
            << lcts.toStringProcessTimesOfMessageHandler(paraComm, paraTimer->getElapsedTime(), ",")
            << std::endl;
         nextOutputLctsStatTime += paraParams->getRealParamValue(UG::CheckpointInterval) * 0.5;
      }

      // Logging : Statistics of Share-data pool
      if( (paraParams->getBoolParamValue(LogShareDataPoolAll) ||
               paraParams->getBoolParamValue(LogShareDataPoolStat)) &&
            ( paraTimer->getElapsedTime() - previousLogShareDataPoolTime )
            > paraParams->getRealParamValue(IntervalTimeOfLogShareDataPool)
            )
      {
         if( paraParams->getBoolParamValue(LogShareDataPoolAll) )
         {
            *osLogShareDataPool << paraTimer->getElapsedTime() << " s; ";
            shareDataPool->writeVectorNormsToOstream(osLogShareDataPool);
         }
         if( paraParams->getBoolParamValue(LogShareDataPoolStat) && shareDataPool->size() > 10 )
         {
            *osLogShareDataPoolStat << "ShareDataPoolStat," << paraTimer->getElapsedTime()
               << "," << shareDataPool->toStatString(cmapLapParaInitiator->getGH())
               << std::endl;
         }
         previousLogShareDataPoolTime = paraTimer->getElapsedTime();
      }

      // Logging : Similarity of solver's basis
      if( paraParams->getBoolParamValue(OutputSimilarityOfBasis) &&
            ( paraTimer->getElapsedTime() - previousOutputSimilarityOfBasisTime )
            > paraParams->getRealParamValue(IntervalTimeOfOutputSimilarityOfBasis)
            )
      {
         std::deque<std::shared_ptr<LatticeBasis<int>>> basisDeque = cmapLapParaSolverPool->getBasisOfSolvers();
         if( basisDeque.size() > 1 )
         {
            outputSimilarityOfBasis(basisDeque);
            previousOutputSimilarityOfBasisTime = paraTimer->getElapsedTime();
         }
      }

      // Adjust NotificationInterval for Load Coordinator's load
      if( (paraTimer->getElapsedTime() - previousTermUpdateNotificationIntervalTime)
            > paraParams->getRealParamValue(LCTermTimeUpdateNotificationInterval) )
      {
         updateNotificationInterval(
               paraTimer->getElapsedTime() - previousTermUpdateNotificationIntervalTime,
               notificationInterval,
               leftNotificationInterval
               );
         lcts.resetLocalProcessTimesOfMessageHandler();
         previousTermUpdateNotificationIntervalTime = paraTimer->getElapsedTime();
      }
      else if( paraTimer->getElapsedTime() - previousLoadCoordinatorStatisticsStdoutTime > 600 )
      {
         std::cout
            << "***"
            << " LoadCoordinator Statistics"
            << " Time " << previousLoadCoordinatorStatisticsStdoutTime << " - " << paraTimer->getElapsedTime()
            << " ***"
            << std::endl;
         double termTime = paraTimer->getElapsedTime() - previousLoadCoordinatorStatisticsStdoutTime;
         double idleRatio = 1.0 - lcts.totalLocalProcessTimesOfMessageHandler() / termTime;
         std::cout
            << "***    "
            << " LoadCoordinator Idle Ratio " << idleRatio
            << " ***"
            << std::endl;
         std::cout
            << "***    "
            << " LoadCoordinator handled " << lcts.localProcessCallsOfMessageHandler[UG::TagSolverState]
            << " solver statuses"
            << " ***"
            << std::endl;
         std::cout
            << "***    "
            << " LoadCoordinator handled " << lcts.localProcessCallsOfMessageHandler[UG::TagSolution]
            << " solutions"
            << " ***"
            << std::endl;
         lcts.resetLocalProcessTimesOfMessageHandler();
         previousLoadCoordinatorStatisticsStdoutTime = paraTimer->getElapsedTime();
      }
   }

   // send terminate request to CheckpointWriter
   if( paraParams->getBoolParamValue(CheckpointThreading) )
   {
      lcCheckpointComm->send( NULL, 0, UG::ParaBYTE, 1, UG::TagTerminateRequest ); // destination is CheckpointWriter( 1 )
   }

   if( paraLocalSolverPool ) paraLocalSolverPool->terminateAllSolvers(lcLocalComm);
   outputAssignmentTable('E');
   outputCsvAssignmentTable('E');

   // set final solver status
   if( hardTimeLimitIsReached )
   {
      cmapLapParaInitiator->setFinalSolverStatus(ParaCMapLAP::HardTimeLimitIsReached);
   }
   if( lowerBoundIsReached )
   {
      cmapLapParaInitiator->setFinalSolverStatus(ParaCMapLAP::LowerBoundIsReached);
   }

   // Logging : Statistics of Load Coordinator
   *osCsvLogCheckpoint
      << Logging::getCsvCheckpointStateString(
            paraTimer->getElapsedTime(),
            -1,lcts.totalCopyCheckpointTime+lcts.totalWriteCheckpointTime,
            lcts.idleTime,lcts.idleTimeToWaitIsend,
            -1,-1,-1,-1,-1,",",'T')
      << std::endl;
   *osCsvLctsStat
      << lcts.toStringProcessTimesOfMessageHandler(paraComm, paraTimer->getElapsedTime(), ",")
      << std::endl;

   // Logging : Statistics of Share-data pool
   if( paraParams->getBoolParamValue(LogShareDataPoolAll) )
   {
      *osLogShareDataPool << paraTimer->getElapsedTime() << " s; ";
      shareDataPool->writeVectorNormsToOstream(osLogShareDataPool);
   }
   if( paraParams->getBoolParamValue(LogShareDataPoolStat) && shareDataPool->size() > 10 )
   {
      *osLogShareDataPoolStat << "ShareDataPoolStat," << paraTimer->getElapsedTime()
         << "," << shareDataPool->toStatString(cmapLapParaInitiator->getGH())
         << std::endl;
   }

}


///
/// update and notice NotificationInterval to solvers
///
void
CMapLapParaLoadCoordinator::updateNotificationInterval(
      double termTime,                 ///< time of this term
      double &notificationInterval,    ///< current solver's notification interval
      double &leftNotificationInterval ///< previous solver's notification interval
      )
{
   double idleRatio = 1.0 - lcts.totalLocalProcessTimesOfMessageHandler() / termTime;
   if( !paraParams->getBoolParamValue(AutoAdjustmentNotificationInterval) )
   {
      return;
   }
   std::cout
      << "***"
      << " LoadCoordinator Statistics "
      << " Term at " << paraTimer->getElapsedTime()
      << " ***"
      << std::endl;
   std::cout
      << "***    "
      << " LoadCoordinator Idle Ratio " << idleRatio
      << " ***"
      << std::endl;
   std::cout
      << "***    "
      << " LoadCoordinator handled " << lcts.localProcessCallsOfMessageHandler[UG::TagSolverState]
      << " solver statuses"
      << " ***"
      << std::endl;
   std::cout
      << "***    "
      << " LoadCoordinator handled " << lcts.localProcessCallsOfMessageHandler[UG::TagSolution]
      << " solutions"
      << " ***"
      << std::endl;
   if( idleRatio < paraParams->getRealParamValue(LCLowerIdleRatio) )
   {
      double target = paraParams->getRealParamValue(LCUpperIdleRatio) * 0.99;
      std::cout
         << "***"
         << " Term at "     << paraTimer->getElapsedTime()
         << " idle ratio " << idleRatio << " <"
         << " lower "      << paraParams->getRealParamValue(LCLowerIdleRatio)
         << " : Change notification interval "
         << notificationInterval << " -> " << notificationInterval * std::min(target / idleRatio, 2.0)
         << " ***"
         << std::endl;
      leftNotificationInterval = notificationInterval;
      notificationInterval *= std::min(target / idleRatio, 2.0);
      for(int i = 1; i < paraComm->getSize(); i++ )
      {
         paraComm->send( &notificationInterval, 1, UG::ParaDOUBLE, i, TagUpdateNotificationInterval );
      }
   }
   else if( idleRatio > paraParams->getRealParamValue(LCUpperIdleRatio) )
   {
      double diff = (notificationInterval - leftNotificationInterval) * 0.5;
      if( diff > 1 )
      {
         std::cout
            << "*** "
            << "Term at "     << paraTimer->getElapsedTime()
            << " idle ratio " << idleRatio << " >"
            << " upper "      << paraParams->getRealParamValue(LCUpperIdleRatio)
            << " : Change notification interval "
            << notificationInterval << " -> " << notificationInterval - diff
            << " ***"
            << std::endl;
         notificationInterval -= diff;
         for(int i = 1; i < paraComm->getSize(); i++ )
         {
            paraComm->send( &notificationInterval, 1, UG::ParaDOUBLE, i, TagUpdateNotificationInterval );
         }
      }
      double target = paraParams->getRealParamValue(LCLowerIdleRatio);
      leftNotificationInterval = std::min(
            leftNotificationInterval * ( target / idleRatio ),
            leftNotificationInterval - diff * 0.5
            );
   }
   *osLogNotificationInterval
      << Logging::getNotificationIntervalLog(
            paraTimer->getElapsedTime(),
            idleRatio,
            paraParams->getRealParamValue(LCLowerIdleRatio),
            paraParams->getRealParamValue(LCUpperIdleRatio),
            notificationInterval,
            leftNotificationInterval
            ) << ","
      << lcts.toStringLocalProcessTimesOfMessageHandler(paraComm, ",")
      << std::endl;
}


///
/// deconstructor
///
CMapLapParaLoadCoordinator::~CMapLapParaLoadCoordinator(
      )
{
   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm );

   // for termination of threadId 1 -- nThreadsPerRank solvers
   for(int i = 1; i < paraComm->getSize(); i++ )
   {
      paraComm->send( NULL, 0, UG::ParaBYTE, i, UG::TagTerminateRequest );
   }

   double totalDeepBkzTime = 0;
   double totalEnumTime    = 0;
   double totalSieveTime   = 0;
   nTerminated = 0;
   size_t nSolversWithoutThreadZero =
      (static_cast<int>(cmapLapParaComm->getSize()) - 1) * static_cast<int>(nThreadsPerRank);
   while ( nTerminated < nSolversWithoutThreadZero )
   {
      int source = -1;
      int tag = -1;
      cmapLapParaComm->probe(&source, &tag);
      std::unique_ptr<CMapLapParaSolverTerminationState> termState{
         dynamic_cast<CMapLapParaSolverTerminationState *>(cmapLapParaComm->createParaSolverTerminationState())};
      termState->receive(cmapLapParaComm, source, tag);
      // for statistic output
      totalDeepBkzTime += termState->getRunningTimeDeepBkz();
      totalEnumTime    += termState->getRunningTimeEnum();
      totalSieveTime   += termState->getRunningTimeSieve();
      if( osStatisticsFinalRun )
      {
         *osStatisticsFinalRun << paraTimer->getElapsedTime()
               << " [S." << source << ", T." << termState->getThreadId() << "] "
               << termState->toString(paraInitiator) << std::endl;
         *osCsvStatisticsFinalRun << paraTimer->getElapsedTime()
               << "," << termState->toCsvString(paraInitiator);
      }
      nTerminated++;
   }

   // for termination of threadId 0 solvers
   for(int i = 1; i < paraComm->getSize(); i++ )
   {
      paraComm->send( NULL, 0, UG::ParaBYTE, i, UG::TagTerminateRequest );
   }

   cmapLapParaInitiator->outputFinalSolverStatistics(
         0, paraTimer->getElapsedTime(),
         totalDeepBkzTime, totalEnumTime, totalSieveTime
         );

   if( paraParams->getBoolParamValue(CheckpointThreading) )
   {
      // for termination of CheckpointWriter
      int checkpointSource = -1;
      int checkpointTag = -1;
      while( lcCheckpointComm->probe(&checkpointSource, &checkpointTag) )
      {
         lcCheckpointComm->receive( NULL, 0, UG::ParaBYTE, 1, checkpointTag );
         if( checkpointTag == UG::TagTerminated ) break;
      }
      assert( checkpointTag == UG::TagTerminated );
      std::cout
         << "*** "
         << "Time "
         << paraTimer->getElapsedTime()
         << " CMapLapParaLoadCoordinator::deconstructor"
         << " receive " << checkpointTag << " from CheckpointWriter "
         << " ***"
         << std::endl;
   }

   if( localMessageHandler ) delete[] localMessageHandler;
   if( lcLocalComm ) delete lcLocalComm;
   if( lcCheckpointComm ) delete lcCheckpointComm;

}


#ifdef UG_WITH_ZLIB

///
/// copy objects for checkpoint
/// @return true if Load Coardinator should send new task tag
///
bool
CMapLapParaLoadCoordinator::copyCheckpointObjects(
      )
{
   bool shouldSendNewTaskToChekcpointWriter = false;

   if( !currentCheckpointElement )
   {
      shouldSendNewTaskToChekcpointWriter = true;
   }
   else
   {
      int lcSource = -1;
      int lcTag = -1;
      if( lcCheckpointComm->iProbe(&lcSource, &lcTag) )
      {
         if( lcTag == UG::TagCompletionOfCalculation )
         {
            // CheckpointWriter finishes to write the previous checkpoint files
            lcCheckpointComm->receive( NULL, 0, UG::ParaBYTE, 1, lcTag );
            strcpy(lastCheckpointTimeStr, currentCheckpointElement->lastCheckpointTimeStr);
            shouldSendNewTaskToChekcpointWriter = true;
            if( logSolvingStatusFlag )
            {
               *osLogSolvingStatus
                  << paraTimer->getElapsedTime()
                  << " CMapLapParaLoadCoordinator::copyCheckpointObjects"
                  << " CheckpointWriter has finised to write; "
                  << " lastCheckpointTimeStr " << lastCheckpointTimeStr
                  << std::endl;

               *osLogSolvingStatus
                  << Logging::getCheckpointStateString(
                        paraTimer->getElapsedTime(),
                        currentCheckpointElement->writeCheckpointTime,
                        -1,-1,-1,-1,-1,'W')
                  << std::endl;
               *osCsvLogCheckpoint
                  << Logging::getCsvCheckpointStateString(
                        paraTimer->getElapsedTime(),
                        currentCheckpointElement->writeCheckpointTime,
                        -1,-1,-1,-1,-1,-1,-1,-1,",",'W')
                  << std::endl;
            }

            currentCheckpointElement = std::move(nextCheckpointElement);
            if( currentCheckpointElement )
            {
               currentCheckpointElement->setLastCheckpointTimeStr(lastCheckpointTimeStr);
            }
         }
         else
         {
            THROW_LOGICAL_ERROR2( lcTag, " is invalid tag for checkpoint writer" );
         }
      }
      else
      {
         // CheckpointWriter does not finish to write previous checkpoint files
         if( logSolvingStatusFlag )
         {
            *osLogSolvingStatus
               << paraTimer->getElapsedTime()
               << " CMapLapParaLoadCoordinator::copyCheckpointObjects"
               << " CheckpointWriter does not finish to write the previous checkpoint files"
               << std::endl;
         }
      }
   }

   if( !paraParams->getBoolParamValue(CheckpointReserving) &&
         currentCheckpointElement )
   {
      return shouldSendNewTaskToChekcpointWriter;
   }

   //
   // create checkpointElement
   //
   double startTime = paraTimer->getElapsedTime();
   std::unique_ptr<CheckpointElement> checkpointElement{new CheckpointElement()};
   checkpointElement->vectorElementQueue = std::shared_ptr<VectorElementQueue>(
         shareDataPool->getVectorElementQueue(
            paraComm, paraParams->getIntParamValue(WriteSizeShareDataPool)
            )
         );
   checkpointElement->basisElementQueue = std::shared_ptr<BasisElementQueue>(
         instancePool->getBasisElementQueue(paraComm, -1)
         );
   checkpointElement->paraActiveTaskQueue = std::shared_ptr<CMapLapParaTaskQueue>(
         cmapLapParaSolverPool->getParaActiveTaskQueue(paraComm)
         );
   checkpointElement->paraDeepBkzTaskQueue = paraDeepBkzTaskPool->getParaTaskQueue(paraComm);
   checkpointElement->paraEnumTaskQueue = paraEnumTaskPool->getParaTaskQueue(paraComm);
   checkpointElement->paraSieveTaskQueue = paraSieveTaskPool->getParaTaskQueue(paraComm);
   checkpointElement->globalBestCMapLapSolution = std::shared_ptr<CMapLapParaSolution>(
         dynamic_cast<CMapLapParaSolution *>(
            cmapLapParaInitiator->getGlobalBestIncumbentSolution()->clone(paraComm)
            )
         );
   checkpointElement->incumbentBasis = std::shared_ptr<BasisElement>(
         instancePool->getIncumbentBasis()->clone()
         );
   lcts.isCheckpointState = true;
   checkpointElement->lcStatString              = lcts.toString();
   checkpointElement->setLastCheckpointTimeStr(lastCheckpointTimeStr);
   checkpointElement->setProbName(cmapLapParaInitiator->getParaInstance()->getProbName());
   double checkpointTime = paraTimer->getElapsedTime() - startTime;


   //
   // write log
   //
   lcts.nCopyCheckpointing++;
   lcts.totalCopyCheckpointTime += checkpointTime;
   if( logSolvingStatusFlag )
   {
      *osLogSolvingStatus
         << Logging::getCheckpointStateString(
               paraTimer->getElapsedTime(),checkpointTime,
               checkpointElement->paraActiveTaskQueue->size(),       // nActive
               checkpointElement->paraDeepBkzTaskQueue->size()       // nInactive
                  + checkpointElement->paraEnumTaskQueue->size(),
               checkpointElement->paraSieveTaskQueue->size(),        // nInactive(Sieve)
               checkpointElement->basisElementQueue->size(),         // nBasis
               checkpointElement->vectorElementQueue->size(),        // nVector
               'C')
         << std::endl;
      *osCsvLogCheckpoint
         << Logging::getCsvCheckpointStateString(
               paraTimer->getElapsedTime(),checkpointTime,
               lcts.totalCopyCheckpointTime+lcts.totalWriteCheckpointTime,
               lcts.idleTime,lcts.idleTimeToWaitIsend,
               checkpointElement->paraActiveTaskQueue->size(),       // nActive
               checkpointElement->paraDeepBkzTaskQueue->size()       // nInactive
                  + checkpointElement->paraEnumTaskQueue->size(),
               checkpointElement->paraSieveTaskQueue->size(),        // nInactive(Sieve)
               checkpointElement->basisElementQueue->size(),         // nBasis
               checkpointElement->vectorElementQueue->size(),        // nVector
               ",", 'C')
         << std::endl;
   }


   assert( !nextCheckpointElement || currentCheckpointElement );
   if( !currentCheckpointElement )
   {
      if( logSolvingStatusFlag )
      {
         *osLogSolvingStatus
            << paraTimer->getElapsedTime()
            << " CMapLapParaLoadCoordinator::copyCheckpointObjects"
            << " set currentCheckpointElement"
            << std::endl;
      }
      currentCheckpointElement = std::move(checkpointElement);
   }
   else
   {
      if( logSolvingStatusFlag )
      {
         if( nextCheckpointElement )
         {
            *osLogSolvingStatus
               << paraTimer->getElapsedTime()
               << " CMapLapParaLoadCoordinator::copyCheckpointObjects"
               << " remove old nextCheckpointElement"
               << std::endl;
         }
         *osLogSolvingStatus
            << paraTimer->getElapsedTime()
            << " CMapLapParaLoadCoordinator::copyCheckpointObjects"
            << " set nextCheckpointElement"
            << std::endl;
      }
      nextCheckpointElement = std::move(checkpointElement);
   }
   return shouldSendNewTaskToChekcpointWriter;
}



void
CMapLapParaLoadCoordinator::updateCheckpointFiles(
      )
{
   double startCheckpointTime = paraTimer->getElapsedTime();

   if( paraParams->getBoolParamValue(CheckpointThreading) )
   {
      bool shouldSendNewTaskToChekcpointWriter = copyCheckpointObjects();
      if( shouldSendNewTaskToChekcpointWriter )
      {
         assert( currentCheckpointElement );
         std::cout
            << "*** "
            << "Time "
            << paraTimer->getElapsedTime()
            << " CMapLapParaLoadCoordinator::updateCheckpointFiles"
            << " send currentCheckpointElement "
            << " ***"
            << std::endl;
         lcCheckpointComm->uTypeSend(
               (void *)currentCheckpointElement.get(),
               ParaCheckpointElementType,
               1, // CheckpointWriter Thread
               UG::TagTask);
      }
      previousCheckpointTime = paraTimer->getElapsedTime();
      return;
   }


   time_t timer;
   char timeStr[30];


   // get checkpoint
   time(&timer);
   // make checkpoint time string
#ifdef _MSC_VER
   int bufsize = 256;
   ctime_s(timeStr, bufsize, &timer);
#else
   ctime_r(&timer, timeStr);
#endif
   for( int i = 0; timeStr[i] != '\0' && i < 26; i++ )
   {
      if( timeStr[i] == ' ') timeStr[i] = '_';
      if( timeStr[i] == '\n' ) timeStr[i] = '\0';
   }
   char *newCheckpointTimeStr = &timeStr[4];    // remove a day of the week
   if( strcmp(newCheckpointTimeStr,lastCheckpointTimeStr) == 0 )
   {
      int l = strlen(newCheckpointTimeStr);
      newCheckpointTimeStr[l] = 'a';
      newCheckpointTimeStr[l+1] = '\0';
   }

   ///
   /// save tasks information
   ///
   char taskFileName[256];
   sprintf(taskFileName,"%s%s_%s_tasks.gz",
         paraParams->getStringParamValue(UG::CheckpointFilePath),
         paraInitiator->getParaInstance()->getProbName(),newCheckpointTimeStr);
   gzstream::ogzstream checkpointTaskStream;
   checkpointTaskStream.open(taskFileName, std::ios::out | std::ios::binary);
   if( !checkpointTaskStream )
   {
      std::cout << "Checkpoint file for tasks cannot open. file name = " << taskFileName << std::endl;
      exit(1);
   }
   int nActive = cmapLapParaSolverPool->writeActiveTasksToCheckpointFile(checkpointTaskStream);
   int nInactiveDeepBkz = paraDeepBkzTaskPool-> writeTasksToCheckPointFile(checkpointTaskStream);
   int nInactiveEnum = paraEnumTaskPool-> writeTasksToCheckPointFile(checkpointTaskStream);
   int nInactiveSieve = paraSieveTaskPool-> writeTasksToCheckPointFile(checkpointTaskStream);


   ///
   /// save incumbent solution
   ///
   char solutionFileName[256];
   sprintf(solutionFileName,"%s%s_%s_solution.gz",
         paraParams->getStringParamValue(UG::CheckpointFilePath),
         paraInitiator->getParaInstance()->getProbName(),newCheckpointTimeStr);
   gzstream::ogzstream checkpointSolutionStream;
   checkpointSolutionStream.open(solutionFileName, std::ios::out | std::ios::binary);
   cmapLapParaInitiator->writeCheckpointSolution(checkpointSolutionStream);
   checkpointSolutionStream.close();


   ///
   /// save InstancePool
   ///
   char instancePoolFileName[256];
   sprintf(instancePoolFileName,"%s%s_%s_instancePool.gz",
         paraParams->getStringParamValue(UG::CheckpointFilePath),
         paraInitiator->getParaInstance()->getProbName(),newCheckpointTimeStr);
   gzstream::ogzstream checkpointInstancePoolStream;
   checkpointInstancePoolStream.open(instancePoolFileName, std::ios::out | std::ios::binary);
   if( !checkpointInstancePoolStream )
   {
      std::cout << "Checkpoint file for instancePool cannot open. file name = " << instancePoolFileName << std::endl;
      exit(1);
   }
   int nBasis = instancePool->writeBasisElementsToCheckpointFile(checkpointInstancePoolStream, -1);


   ///
   /// save Incumbent Basis
   ///
   char basisTextFileName[256];
   sprintf(basisTextFileName,"%s%s_%s_basis.txt",
         paraParams->getStringParamValue(UG::CheckpointFilePath),
         paraInitiator->getParaInstance()->getProbName(),newCheckpointTimeStr);
   instancePool->writeIncumbentBasis(basisTextFileName);


   ///
   /// save ShareDataPool
   ///
   char shareDataPoolFileName[256];
   sprintf(shareDataPoolFileName,"%s%s_%s_shareDataPool.gz",
         paraParams->getStringParamValue(UG::CheckpointFilePath),
         paraInitiator->getParaInstance()->getProbName(),newCheckpointTimeStr);
   gzstream::ogzstream checkpointShareDataPoolStream;
   checkpointShareDataPoolStream.open(shareDataPoolFileName, std::ios::out | std::ios::binary);
   if( !checkpointShareDataPoolStream )
   {
      std::cout << "Checkpoint file for ShareDataPool cannot open. file name = " << shareDataPoolFileName << std::endl;
      exit(1);
   }
   int nVector = shareDataPool->writeVectorElementsToCheckpointFile(
         checkpointShareDataPoolStream,
         lcCheckpointComm,
         paraParams->getIntParamValue(WriteSizeShareDataPool));


   ///
   /// save LoadCoordinator statistics
   ///
   char loadCoordinatorStatisticsFileName[256];
   sprintf(loadCoordinatorStatisticsFileName,"%s%s_%s_loadCoordinatorStatistics.txt",
         paraParams->getStringParamValue(UG::CheckpointFilePath),
         paraInitiator->getParaInstance()->getProbName(),newCheckpointTimeStr);
   std::ofstream lcStaticFile(loadCoordinatorStatisticsFileName);
   lcts.isCheckpointState = true;
   lcStaticFile << lcts.toString();


   ///
   /// remove old checkfiles after second checkpointing
   ///
   if( lastCheckpointTimeStr[0] != ' ' )  //  first checkpoint, should not remove
   {
      sprintf(taskFileName,"%s%s_%s_tasks.gz",
            paraParams->getStringParamValue(UG::CheckpointFilePath),
            paraInitiator->getParaInstance()->getProbName(),lastCheckpointTimeStr);
      sprintf(solutionFileName,"%s%s_%s_solution.gz",
            paraParams->getStringParamValue(UG::CheckpointFilePath),
            paraInitiator->getParaInstance()->getProbName(),lastCheckpointTimeStr);
      sprintf(basisTextFileName,"%s%s_%s_basis.txt",
            paraParams->getStringParamValue(UG::CheckpointFilePath),
            paraInitiator->getParaInstance()->getProbName(),lastCheckpointTimeStr);
      sprintf(instancePoolFileName,"%s%s_%s_instancePool.gz",
            paraParams->getStringParamValue(UG::CheckpointFilePath),
            paraInitiator->getParaInstance()->getProbName(),lastCheckpointTimeStr);
      sprintf(shareDataPoolFileName,"%s%s_%s_shareDataPool.gz",
            paraParams->getStringParamValue(UG::CheckpointFilePath),
            paraInitiator->getParaInstance()->getProbName(),lastCheckpointTimeStr);
      sprintf(loadCoordinatorStatisticsFileName,"%s%s_%s_loadCoordinatorStatistics.txt",
            paraParams->getStringParamValue(UG::CheckpointFilePath),
            paraInitiator->getParaInstance()->getProbName(),lastCheckpointTimeStr);
      if( remove(taskFileName) )
      {
         std::cout << "checkpoint task file cannot be removed: errno = " << strerror(errno) << std::endl;
         exit(1);
      }
      if ( remove(solutionFileName) )
      {
         std::cout << "checkpoint solution file cannot be removed: errno = " << strerror(errno) << std::endl;
         exit(1);
      }
      if ( remove(basisTextFileName) )
      {
         std::cout << "checkpoint Basis text file cannot be removed: errno = " << strerror(errno) << std::endl;
         exit(1);
      }
      if ( remove(instancePoolFileName) )
      {
         std::cout << "checkpoint InstancePool file cannot be removed: errno = " << strerror(errno) << std::endl;
         exit(1);
      }
      if ( remove(shareDataPoolFileName) )
      {
         std::cout << "checkpoint ShareDataPool file cannot be removed: errno = " << strerror(errno) << std::endl;
         exit(1);
      }
      if ( remove(loadCoordinatorStatisticsFileName) )
      {
         std::cout << "checkpoint LoadCoordinatorStatistics file cannot be removed: errno = " << strerror(errno) << std::endl;
         exit(1);
      }
   }


   ///
   /// remove after_checkpoint_solution if it exists
   ///
   char afterCheckpointingSolutionFileName[256];
   sprintf(afterCheckpointingSolutionFileName,"%s%s_after_checkpointing_solution.gz",
         paraParams->getStringParamValue(UG::CheckpointFilePath),
         paraInitiator->getParaInstance()->getProbName() );
   gzstream::igzstream afterCheckpointingSolutionStream;
   afterCheckpointingSolutionStream.open(afterCheckpointingSolutionFileName, std::ios::in | std::ios::binary);
   if( afterCheckpointingSolutionStream  )
   {
      // afater checkpointing solution file exists
      afterCheckpointingSolutionStream.close();
      if ( remove(afterCheckpointingSolutionFileName) )
      {
         std::cout << "after checkpointing solution file cannot be removed: errno = " << strerror(errno) << std::endl;
         exit(1);
      }
   }


   double checkpointTime = paraTimer->getElapsedTime() - startCheckpointTime;
   lcts.totalWriteCheckpointTime += checkpointTime;
   lcts.nWriteCheckpointing++;

   if( logSolvingStatusFlag )
   {
      *osLogSolvingStatus
         << Logging::getCheckpointStateString(
               paraTimer->getElapsedTime(),checkpointTime,
               nActive,nInactiveDeepBkz+nInactiveEnum,nInactiveSieve,
               nBasis,nVector);
      *osCsvLogCheckpoint
         << Logging::getCsvCheckpointStateString(
               paraTimer->getElapsedTime(),checkpointTime,
               lcts.totalCopyCheckpointTime+lcts.totalWriteCheckpointTime,
               lcts.idleTime,lcts.idleTimeToWaitIsend,
               nActive,nInactiveDeepBkz+nInactiveEnum,nInactiveSieve,
               nBasis,nVector, ",");
   }

   // update last checkpoint time string
   strcpy(lastCheckpointTimeStr,newCheckpointTimeStr);
   previousCheckpointTime = paraTimer->getElapsedTime();
}


void
CMapLapParaLoadCoordinator::warmStart(
      )
{
   restarted = true;

   if( logSolvingStatusFlag )
   {
      *osLogSolvingStatus << "WarmStart from the checkpoint file " << std::string(cmapLapParaInitiator->getPrefixWarm()) << std::endl;
   }

   ///
   /// set global incumbet solution from _solution.gz and after_checkpointing_solution.gz
   ///
   char afterCheckpointingSolutionFileName[256];
   int prefixWarmParentPathIndex = std::string(cmapLapParaInitiator->getPrefixWarm()).rfind("/");
   std::string prefixWarmParentPath = std::string(cmapLapParaInitiator->getPrefixWarm()).substr(0, prefixWarmParentPathIndex + 1);
   sprintf(afterCheckpointingSolutionFileName, "%s%s_after_checkpointing_solution.gz",
         prefixWarmParentPath.c_str(),
         paraInitiator->getParaInstance()->getProbName() );
   gzstream::igzstream afterCheckpointingSolutionStream;
   afterCheckpointingSolutionStream.open(afterCheckpointingSolutionFileName, std::ios::in | std::ios::binary);
   if( afterCheckpointingSolutionStream.good() )
   {
      double shortestNorm = paraInitiator->readSolutionFromCheckpointFile(afterCheckpointingSolutionFileName);
      if( !paraParams->getBoolParamValue(UG::Quiet) )
      {
         std::ostringstream strs;
         strs << shortestNorm;
         cmapLapParaInitiator->writeSolution("[Warm started from "+std::string(cmapLapParaInitiator->getPrefixWarm())+" : the solution "+strs.str()+" from the checkpoint file]");
      }
      std::cout
         << "*** "
         << "CMapLapParaLoadCoordinator::warmStart load solution " << shortestNorm << " from after checkpoint solution " << afterCheckpointingSolutionFileName
         << " ***"
         << std::endl;
   }
   else
   {
      std::cout << "CMapLapParaLoadCoordinator::warmStart donot find " << afterCheckpointingSolutionFileName << std::endl;

      ///
      /// load incumbent solution
      ///
      char solutionFileName[256];
      sprintf(solutionFileName,"%s_solution.gz", cmapLapParaInitiator->getPrefixWarm());
      gzstream::igzstream checkpointSolutionStream;
      checkpointSolutionStream.open(solutionFileName, std::ios::in | std::ios::binary);
      if( !checkpointSolutionStream )
      {
         std::cout
            << "*** "
            << "Checkpoint file for solution cannot open. file name = " << solutionFileName
            << " ***"
            << std::endl;
         exit(1);
      }
      double shortestNorm = paraInitiator->readSolutionFromCheckpointFile(solutionFileName);
      if( !paraParams->getBoolParamValue(UG::Quiet) )
      {
         std::ostringstream strs;
         strs << shortestNorm;
         cmapLapParaInitiator->writeSolution("[Warm started from "+std::string(cmapLapParaInitiator->getPrefixWarm())+" : the solution "+strs.str()+" from the checkpoint file]");
      }
      std::cout
         << "*** "
         << "CMapLapParaLoadCoordinator::warmStart load solution " << shortestNorm << " from checkpoint solution " << solutionFileName
         << " ***"
         << std::endl;
      checkpointSolutionStream.close();
   }
   afterCheckpointingSolutionStream.close();

   if( !paraParams->getBoolParamValue(WarmStartOnlyPool) )
   {
      ///
      /// read tasks information
      ///
      char tasksFileName[256];
      sprintf(tasksFileName,"%s_tasks.gz", cmapLapParaInitiator->getPrefixWarm());
      gzstream::igzstream  checkpointTasksStream;  ///< gzstream for checkpoint tasks file
      checkpointTasksStream.open(tasksFileName, std::ios::in | std::ios::binary);
      if( !checkpointTasksStream.good() )
      {
         std::cerr << "ERROR: Opening file `" << tasksFileName << "' failed.\n";
         exit(1);
      }
      std::unique_ptr<CMapLapParaTaskPoolInAscendingOrder> tempTaskPool{new CMapLapParaTaskPoolInAscendingOrder()};
      tempTaskPool->readTasksFromCheckPointFile(paraComm, checkpointTasksStream);
      checkpointTasksStream.close();

      std::shared_ptr<CMapLapParaTask> paraTask(nullptr);

      if( logSolvingStatusFlag )
      {
         *osLogSolvingStatus
            << "warmStart: randomize basis; size = " << paraParams->getIntParamValue(RandomizeRows) << std::endl;
      }
      while( (paraTask=tempTaskPool->extractTask()) )
      {
         switch( paraTask->getSolverType() )
         {
         case DeepBkz:
            paraTask->setU(0); // do not randamize
            // paraTask->setU(paraParams->getIntParamValue(RandomizeRows));
            paraDeepBkzTaskPool->insert(paraTask);  // active DeepBkz tasks in the previous run
            break;
         case Enum:
            paraEnumTaskPool->insert(paraTask);
            break;
         case Sieve:
            paraSieveTaskPool->insert(paraTask);
            break;
         default:
            THROW_LOGICAL_ERROR2("CMapLapParaLoadCoordinator::warmStart: Invalid solver type = ", static_cast<int>(paraTask->getSolverType()));
         }
      }
      std::cout
         << "*** "
         << "CMapLapParaLoadCoordinator::warmStart: Invalid solver type warning is no problem during warm start."
         << " ***"
         << std::endl;
      if( logSolvingStatusFlag )
      {
         *osLogSolvingStatus << "warmStart read "
            << paraDeepBkzTaskPool->getNumOfTasks() << " DeepBkzTasks"
            << paraEnumTaskPool->getNumOfTasks()    << ", EnumTasks"
            << paraSieveTaskPool->getNumOfTasks()   << ", SieveTasks" << std::endl;
      }
   }


   ///
   /// load InstancePool
   ///
   char instancePoolFileName[256];
   sprintf(instancePoolFileName,"%s_instancePool.gz", cmapLapParaInitiator->getPrefixWarm());
   gzstream::igzstream checkpointInstancePoolStream;
   checkpointInstancePoolStream.open(instancePoolFileName, std::ios::in | std::ios::binary);
   if( !checkpointInstancePoolStream )
   {
      std::cout << "Checkpoint file for instancePool cannot open. file name = " << instancePoolFileName << std::endl;
      exit(1);
   }
   instancePool->readBasisElementsFromCheckpointFile(paraComm, checkpointInstancePoolStream);
   checkpointInstancePoolStream.close();


   bool clearSolverThreadIds = paraParams->getBoolParamValue(WarmStartOnlyPool);

   ///
   /// load ShareDataPool
   ///
   char shareDataPoolFileName[256];
   sprintf(shareDataPoolFileName,"%s_shareDataPool.gz", cmapLapParaInitiator->getPrefixWarm());
   gzstream::igzstream checkpointShareDataPoolStream;
   checkpointShareDataPoolStream.open(shareDataPoolFileName, std::ios::in | std::ios::binary);
   if( !checkpointShareDataPoolStream )
   {
      std::cout << "Checkpoint file for ShareDataPool cannot open. file name = " << shareDataPoolFileName << std::endl;
      exit(1);
   }
   shareDataPool->readVectorElementsFromCheckpointFile(paraComm, checkpointShareDataPoolStream, clearSolverThreadIds);
   checkpointShareDataPoolStream.close();


   /// If sieve solvers with highest priority exist, sieve solvers should be activeted.


   // run main process
   run();
}


void
CMapLapParaLoadCoordinator::outputSimilarityOfBasisHeader(
      )
{
   int interval = 10;
   int nCalcIndexes = std::ceil(cmapLapParaInitiator->getDimension() / static_cast<double>(interval));
   std::vector<int> grassmannIndexes(nCalcIndexes);
   for( int i = 0; i < nCalcIndexes; ++i )
      grassmannIndexes[i] = interval * i;
   BasisSimilarity::setGrassmannIndexes(grassmannIndexes);
   *osLogSimilarity
      << ParaCMapLAP::BasisSimilarity::outputSimilarityOfBasisHeader(
            cmapLapParaInitiator->getDimension()
            )
      << std::endl;
}


void
CMapLapParaLoadCoordinator::outputSimilarityOfBasis(
      std::deque<std::shared_ptr<LatticeBasis<int>>> &basisDeque
      )
{
   *osLogSimilarity
      << ParaCMapLAP::BasisSimilarity::outputSimilarityOfBasis(
            basisDeque,
            paraTimer->getElapsedTime(),
            cmapLapParaInitiator->getDimension(),
            paraParams->getIntParamValue(NumOfSamplingForSimilarityOfBasis)
            )
      << std::endl;
}

#endif // End of UG_WITH_ZLIB

} // namespace ParaCMapLAP
