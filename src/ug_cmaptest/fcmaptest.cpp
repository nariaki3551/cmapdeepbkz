/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*          This file is part of the program and software framework          */
/*                 CMAP-TEST --- Test configure for CMAP-LAP                 */
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

/**@file    fcmaptest.cpp
 * @brief   FiberCMapTest MAIN.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <assert.h>
#include <pthread.h>
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <sys/signal.h>
#include <unistd.h>
#include <iostream>
#include <string>
#include <thread>
#include <memory>
#include "ug/paraCommCPP11.h"
#include "ug/paraInstance.h"
#include "ug/paraTimer.h"
#include "cmapLapParaParamSet.h"
#include "cmapLapParaCommTh.h"
#include "cmapLapParaDef.h"
#include "cmapLapParaDeterministicTimer.h"
#include "cmapLapParaInitiator.h"
#include "cmapLapParaInstance.h"
#include "cmapLapParaLoadCoordinator.h"
#include "cmapLapParaSolverLocalComm.h"
#include "cmapTestParaSolver.h"

using namespace ParaCMapLAP;

static std::shared_ptr<CMapLapParaCommTh> comm(nullptr);
static std::shared_ptr<CMapLapParaParamSet> paraParams(nullptr);
static std::unique_ptr<CMapLapParaLoadCoordinator> paraLc(nullptr);
static CMapLapParaSolverLocalComm **localComm = nullptr;
static int nSolvers = 0;

struct SolverThreadData_t {
   int argc;
   char **argv;
   UG::ParaTimer *paraTimer;
   int rank;
   int threadId;
};

typedef struct SolverThreadData_t SolverThreadData;


///
/// interrupt handler for CTRL-C interrupts
///
static
void interruptHandler(
   int signum  ///< interrupt signal number
   )
{
   std::cout << "in interrupt hander" << std::endl;
   paraLc->interrupt();
   delete paraLc.release();
   _exit(1);
}


///
/// @brief output usage
///
void
outputCommandLineMessages(
      char **argv
      )
{
   std::cout << std::endl;
   std::cout << "syntax: " << argv[0]
             << " fcmaptest_param_file problem_file_name "
             << "[-sth <number>] [-w <prefix_warm>] [-fsol <solution_file>]"
             << "[-isol <initial_solution_file>]"
             << "[-o <number> <setting> ... ] [-seed <number>]"
             << std::endl;
   std::cout << "  -sth <number>         : number of threads " << std::endl;
   std::cout << "  -w <prefix_warm>      : warm start file prefix ( prefix_warm_nodes.gz and prefix_warm_solution.txt are read )" << std::endl;
   std::cout << "  -fsol <solution_file> : specify output solution file" << std::endl;
   std::cout << "  -isol <solution_file> : specify initial solution file" << std::endl;
   std::cout << "  -o <number> <setting_file> ... : number and path of setting fileis for overwrite parameter " << std::endl;
   std::cout << "  -seed <number>        : randomized seed" << std::endl;
   std::cout << "File names need to be specified by full path form." << std::endl;
}


///
/// @brief save parameter sets as .prm file
///
void
outputParaParamSet(
      CMapLapParaInitiator &paraInitiator
      )
{
   if( !paraParams->getBoolParamValue(UG::Quiet) )
   {
      std::ofstream ofsParamsOutputFile;
      std::ostringstream s;
      if( paraInitiator.getPrefixWarm() )
      {
         s << paraInitiator.getPrefixWarm();
      }
      else
      {
         s << paraParams->getStringParamValue(UG::LogSolvingStatusFilePath)
         << paraInitiator.getParaInstance()->getProbName();
      }
      s << ".prm";
      ofsParamsOutputFile.open(s.str().c_str());
      if( !ofsParamsOutputFile ){
         std::cout << "Cannot open ParaParams output file: file name = " << s.str() << std::endl;
         exit(1);
      }
      paraParams->write(&ofsParamsOutputFile);
      ofsParamsOutputFile.close();
   }
}


///
/// @brief create and run sovler threads
///
void *
runSolverThread(
      void *threadData
      )
{
   SolverThreadData *solverThreadData = static_cast<SolverThreadData *>(threadData);

   assert( solverThreadData->rank < comm->getSize() );

   if( solverThreadData->threadId == 0)
   {
      dynamic_cast<UG::ParaCommCPP11 *>(comm.get())->solverInit(solverThreadData->rank, paraParams.get());
      localComm[ solverThreadData->rank - 1 ]->mainThreadInit(paraParams.get());
   }
   else
   {
      comm->solverInit(solverThreadData->rank, paraParams.get());
      localComm[ solverThreadData->rank - 1 ]->slaveThreadInit(solverThreadData->threadId, paraParams.get());
   }

   int argc = solverThreadData->argc;
   char **argv = solverThreadData->argv;
   UG::ParaTimer *paraTimer = solverThreadData->paraTimer;
   std::unique_ptr<CMapLapParaDeterministicTimer> detTimer(nullptr);

   sigset_t sigsBlock;
   sigemptyset(&sigsBlock);
   sigaddset(&sigsBlock, SIGINT);
   pthread_sigmask(SIG_BLOCK, &sigsBlock, NULL);

   if( paraParams->getBoolParamValue(UG::Deterministic) )
   {
      detTimer.reset(new CMapLapParaDeterministicTimer());
   }

   UG::ParaInstance *paraInstance = comm->createParaInstance();
   if( solverThreadData->threadId == 0)
   {
      paraInstance->bcast(comm.get(), 0, ON_MEMORY_TRANSFER);
   }
   dynamic_cast<CMapLapParaInstance *>(paraInstance)->bcast(localComm[ solverThreadData->rank - 1 ], 0, ON_MEMORY_TRANSFER);

   CMapTestParaSolver paraSolver(argc, argv, comm.get(), localComm[ solverThreadData->rank - 1 ], paraParams.get(), paraInstance, detTimer.get(), paraTimer->getElapsedTime());

   if( paraParams->getBoolParamValue(UG::StatisticsToStdout) )
   {
      comm->lockApp();
      std::cout << "After Rank " << comm->getRank() << " Solver initialized 1: " << paraTimer->getElapsedTime() << std::endl;
      comm->unlockApp();
   }

   paraSolver.run();

   delete solverThreadData;

   return 0;
}


/**************************************************************************************
 *                                                                                    *
 * Command line see outputCommandLineMessages()                                       *
 *                                                                                    *
 **************************************************************************************/
int
main (
      int  argc,
      char **argv
     )
{
   bool racingSolversExist = false;

   if( argc < 3 )
   {
      outputCommandLineMessages(argv);
      return 1;
   }
   for( int i = 1; i < argc; i++ )
   {
      if( strcmp(argv[i], "-h") == 0 )
      {
         outputCommandLineMessages(argv);
         return 1;
      }
   }

   int nThreadsPerRank = 1;
   /// Defulat comSize is the number of cores in a compute node
// #ifdef _MSC_VER
//    SYSTEM_INFO sysinfo;
//    GetSystemInfo(&sysinfo);
//    nThreadsPerRank = sysinfo.dwNumberOfProcessors; //includes logical cpu
// #else
//    nThreadsPerRank = sysconf(_SC_NPROCESSORS_CONF);
// #endif

   int nOverwriteFilenames = 0;
   char **overwriteFilenames = 0;

   for( int i = 1; i < argc; i++ )
   {
      if( strcmp(argv[i], "-ntpr") == 0 )
      {
         i++;
         if( i < argc )
         {
            nThreadsPerRank = atoi(const_cast<const char*>(argv[i]));   // if -sth 0, then it is considered as use the number of cores system has
         }
         else
         {
            std::cerr << "missing the number of solver threads after parameter '-ntpr'" << std::endl;
            exit(1);
         }
      }
      if( strcmp(argv[i], "-o") == 0 )
      {
         i++;
         if( i < argc )
         {
            nOverwriteFilenames = atoi(argv[i]);
         }
         else
         {
            std::cerr << "missing the filename after parameter '-o'" << std::endl;
            exit(1);
         }
         overwriteFilenames = new char* [nOverwriteFilenames];
         for( int j = 0; j < nOverwriteFilenames; j++ )
         {
            i++;
            if ( i < argc )
            {
               overwriteFilenames[j] = argv[i];
            }
            else
            {
               std::cerr << "missing the filename after parameter '-o'" << std::endl;
               exit(1);
            }
         }
      }
   }

   ///
   /// catch SIGINT
   ///
#ifdef NO_SIGACTION
   (void)signal(SIGINT, interruptHandler);
#else
   struct sigaction newaction;

   // initialize new signal action
   newaction.sa_handler = interruptHandler;
   newaction.sa_flags = SA_RESTART | SA_NODEFER | SA_RESETHAND;
   (void)sigemptyset(&newaction.sa_mask);

   // set new signal action
   (void)sigaction(SIGINT, &newaction, NULL);
#endif

   comm.reset(new CMapLapParaCommTh());
   comm->init(argc,argv);

   std::unique_ptr<UG::ParaTimer> paraTimer(comm->createParaTimer());
   paraTimer->init(comm.get());

   nSolvers = comm->getSize() - 1;
   paraParams = std::shared_ptr<CMapLapParaParamSet>(new CMapLapParaParamSet());
   paraParams->read(comm.get(), argv[1]);
   for( int i = 0; i < nOverwriteFilenames; i++ )
      paraParams->read(comm.get(), overwriteFilenames[i]);
   comm->lcInit(paraParams.get());

   localComm = new CMapLapParaSolverLocalComm*[nSolvers];
   for( int i = 0; i < nSolvers; i++ )
   {
      localComm[i] = new CMapLapParaSolverLocalComm();
      localComm[i]->init(argc, argv);
   }

   CMapLapParaInitiator paraInitiator(comm.get(), paraTimer.get());
   paraInitiator.init(paraParams.get(), argc, argv);

   if( paraParams->getBoolParamValue(UG::StatisticsToStdout) )
   {
      std::cout << "After Initiator initialized: " << paraTimer->getElapsedTime() << std::endl;
   }

   UG::ParaInstance *paraInstance = paraInitiator.getParaInstance();
   if( paraParams->getIntParamValue(UG::OutputParaParams) > 0 )
   {
      outputParaParamSet(paraInitiator);
   }

   std::thread **threads = new std::thread*[nSolvers * nThreadsPerRank];
   int k = 0;
   for( int i = 0; i < nSolvers; i++ )
   {
      for( int j = 0; j < nThreadsPerRank; j++ )
      {
         SolverThreadData *solverThreadData = new SolverThreadData;
         solverThreadData->argc = argc;
         solverThreadData->argv = argv;
         solverThreadData->paraTimer = paraTimer.get();
         solverThreadData->rank = (i+1);
         solverThreadData->threadId = j;
         threads[k++] = new std::thread(runSolverThread, solverThreadData);
      }
   }

   std::unique_ptr<CMapLapParaDeterministicTimer> detTimer(nullptr);
   if( paraParams->getBoolParamValue(UG::Deterministic) )
   {
      detTimer.reset(new CMapLapParaDeterministicTimer());
   }

   paraInstance->bcast(comm.get(), 0, ON_MEMORY_TRANSFER);
   paraInitiator.sendSolverInitializationMessage();  // This messages should be received in constructor of the target Solver

   paraLc.reset(new CMapLapParaLoadCoordinator(comm.get(), paraParams.get(), &paraInitiator, &racingSolversExist, paraTimer.get(), detTimer.get(), nThreadsPerRank));

   if( paraInitiator.isWarmStarted() )
   {
#ifdef UG_WITH_ZLIB
      paraLc->warmStart();
#endif
   }
   else
   {
      paraLc->parallelDispatch();
   }

   delete paraLc.release();
   comm.reset();
   paraParams.reset();

   for( int i = 0; i < nSolvers * nThreadsPerRank; i++ )
   {
      threads[i]->join();
   }
   for( int i = 0; i < nSolvers * nThreadsPerRank; i++ )
   {
      delete threads[i];
   }
   delete [] threads;
   for( int i = 0; i < nSolvers; i++ )
   {
      delete localComm[i];
   }
   delete [] localComm;
   delete [] overwriteFilenames;
   return 0;

} /* END main */
