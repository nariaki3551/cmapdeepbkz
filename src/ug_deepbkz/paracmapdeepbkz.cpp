/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*          This file is part of the program and software framework          */
/*  CMAP-DeepBKZ --- Configurable Massively Parallel Solver for DeepBKZ      */
/*                                                                           */
/*  Copyright Written by Nariaki Tateiwa <n-tateiwa@kyudai.jp>,              */
/*                       Yuji Shinano <shinano@zib.de>,                      */
/*            Copyright (C) 2021 by Zuse Institute Berlin,                   */
/*            licensed under LGPL version 3 or later.                        */
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

/**@file    paracmapdeepbkz.cpp
 * @brief   CMAP-DeepBKZ MAIN.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <unistd.h>
#include <thread>
#include <signal.h>
#include "cmapLapParaParamSet.h"
#include "cmapLapParaSolverLocalComm.h"
#include "cmapLapParaCommMpi.h"
#include "cmapLapParaInstance.h"
#include "cmapLapParaDeterministicTimer.h"
#include "cmapLapParaInitiator.h"
#include "cmapDeepBkzParaSolver.h"
#include "cmapDeepBkzParaLoadCoordinator.h"

using namespace ParaCMapLAP;

static std::shared_ptr<CMapLapParaParamSet> paraParams(nullptr);
static std::unique_ptr<CMapDeepBkzParaLoadCoordinator> paraLc(nullptr);
static std::unique_ptr<CMapLapParaSolverLocalComm> localComm(nullptr);

struct SolverThreadData_t {
   int argc;
   char **argv;
   CMapLapParaParamSet *paraParams;
   CMapLapParaCommMpi *comm;
   UG::ParaTimer *paraTimer;
   int rank;
   int threadId;
};

typedef struct SolverThreadData_t SolverThreadData;

void
outputCommandLineMessages(
      char **argv
      )
{
   std::cout << std::endl;
   std::cout << "syntax: " << argv[0] << " #MPI_processes(#solvers + 1) paracmapdeepbkz_param_file problem_file_name "
             << "[-sl <settings>] [-s <settings>] [-sr <root_settings>] [-w <prefix_warm>] [-fsol <solution_file>]"
             << "[-isol <initial_solution_file>]"
             << "[-seed <number>]"
             << std::endl;
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

   CMapLapParaCommMpi *comm = solverThreadData->comm;

   paraParams.reset(solverThreadData->paraParams);

   assert( solverThreadData->rank < comm->getSize() );
   localComm->slaveThreadInit(solverThreadData->threadId, paraParams.get());

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

   dynamic_cast<CMapLapParaInstance *>(paraInstance)->bcast(localComm.get(), 0, ON_MEMORY_TRANSFER);

   CMapDeepBkzParaSolver paraSolver(argc, argv, comm, localComm.get(), paraParams.get(), paraInstance, detTimer.get(), paraTimer->getElapsedTime());

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

   UG::ParaSysTimer sysTimer;
   sysTimer.start();

#ifdef UG_WITH_UGS
   char *configFileName = 0;
   for( int i = 1; i < argc; ++i )
   {
      if ( strcmp(argv[i], "-ugsc") == 0 )
      {
         i++;
         if( i < argc )
         {
            configFileName = argv[i];
            break;
         }
         else
         {
            std::cerr << "missing file name after parameter '-ugsc" << std::endl;
            exit(1);
         }
      }
   }

   CMapLapParaCommMpi *comm = 0;
   UGS::UgsParaCommMpi *commUgs = 0;   // commUgs != 0 means ParaCMapLAP runs under UGS
   if( configFileName )
   {
      commUgs = new UGS::UgsParaCommMpi();
      commUgs->init(argc,argv);
      comm = new CMapLapParaCommMpi(commUgs->getSolverMPIComm());
      comm->setUgsComm(commUgs);
   }
   else
   {
      comm = new CMapLapParaCommMpi();
   }
#else
   CMapLapParaCommMpi *comm = new CMapLapParaCommMpi();
#endif

   comm->init(argc,argv);

   std::shared_ptr<UG::ParaTimer> paraTimer(new UG::ParaTimerMpi());
   paraTimer->init(comm);

   paraParams = std::shared_ptr<CMapLapParaParamSet>(new CMapLapParaParamSet());
   paraParams->read(comm, argv[1]);
   for( int i = 0; i < nOverwriteFilenames; i++ )
   {
      paraParams->read(comm, overwriteFilenames[i]);
   }

   if( comm->getRank() == 0 )
   {
      if( argc < 3 )
      {
         outputCommandLineMessages(argv);
         return 1;
      }
      paraParams->bcast(comm, 0);
      comm->lcInit(paraParams.get());
      CMapLapParaInitiator paraInitiator(comm, paraTimer.get());
      if( paraInitiator.init(paraParams.get(), argc, argv) )
      {
      }
      std::cout << "** Initiatior was initilized after " << paraTimer->getElapsedTime() << " sec." << std::endl;
      UG::ParaInstance *paraInstance = paraInitiator.getParaInstance();
      if( paraParams->getIntParamValue(UG::OutputParaParams) > 0 )
      {
         outputParaParamSet(paraInitiator);
      }
      paraInstance->bcast(comm, 0, 0);  // 0: root is LC, 0: On memory transfer
      paraInitiator.sendSolverInitializationMessage();  // This messages should be received in constructor of the target Solver
      std::cout << "** Instance data were sent to all solvers after " << paraTimer->getElapsedTime() << " sec." << std::endl;
      std::unique_ptr<CMapLapParaDeterministicTimer> detTimer(nullptr);
      if( paraParams->getBoolParamValue(UG::Deterministic) )
      {
         detTimer.reset(new CMapLapParaDeterministicTimer());
      }

      paraLc.reset(new CMapDeepBkzParaLoadCoordinator(comm, paraParams.get(), &paraInitiator, &racingSolversExist, paraTimer.get(), detTimer.get(), nThreadsPerRank));

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
   }
   else
   {
      if( argc < 3 )
      {
         return 1;
      }
      paraParams->bcast(comm, 0);
      comm->solverInit(paraParams.get());
      UG::ParaInstance *paraInstance = comm->createParaInstance();
      paraInstance->bcast(comm, 0, ON_MEMORY_TRANSFER);   // 0: on  memory
      std::unique_ptr<CMapLapParaDeterministicTimer> detTimer(nullptr);
      if( paraParams->getBoolParamValue(UG::Deterministic) )
      {
         detTimer.reset(new CMapLapParaDeterministicTimer());
      }

      localComm.reset(new CMapLapParaSolverLocalComm());
      localComm->init(argc, argv);

      for( int j = 1; j < nThreadsPerRank; j++ )
      {
         SolverThreadData *solverThreadData = new SolverThreadData;
         solverThreadData->argc = argc;
         solverThreadData->argv = argv;
         solverThreadData->paraParams = paraParams.get();
         solverThreadData->comm = comm;
         solverThreadData->paraTimer = paraTimer.get();
         solverThreadData->rank = comm->getRank();
         solverThreadData->threadId = j;
         std::thread t(runSolverThread, solverThreadData);
         t.detach();
      }

      localComm->mainThreadInit(paraParams.get());
      paraInstance->bcast(localComm.get(), 0, ON_MEMORY_TRANSFER);
      CMapDeepBkzParaSolver paraSolver(argc, argv, comm, localComm.get(), paraParams.get(), paraInstance, detTimer.get(), paraTimer->getElapsedTime());
      paraSolver.run();

      //if( paraInstance ) delete paraInstance;  /** deleted in paraSolver destructor */
   }

   sysTimer.stop();
   std::cout << "[ Rank: " << comm->getRank() << " ], UTime = " << sysTimer.getUTime()
         << ", STime = " << sysTimer.getSTime() << ", RTime = " << sysTimer.getRTime() << std::endl;

   if( racingSolversExist ) comm->abort();

   delete comm;

#ifdef UG_WITH_UGS
   if( commUgs ) delete commUgs;
#endif

   delete [] overwriteFilenames;
   return 0;

} /* END main */

