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

/**@file    cmapLapParaSolver.cpp
 * @brief   ParaSolver extension for CMAP-LAP.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <cfloat>
#include <cstring>
#include <cstdlib>
#include <string>
#include "ug/paraComm.h"
#include "ug/paraTask.h"
#include "ug/paraInstance.h"
#include "ug/paraSolver.h"
#include "ug/paraSolution.h"
#include "ug/paraParamSet.h"
#include "cmapLapParaTagDef.h"
#include "cmapLapParaDef.h"
#include "cmapLapParaSolution.h"
#include "cmapLapParaLattice.h"
#include "cmapLapParaInstance.h"
#include "cmapLapParaSolver.h"
#include "cmapLapParaBasis.h"
#include "cmapLapParaTask.h"
#include "cmapLapParaComm.h"
#include "cmapLapParaPackedVector.h"
#include "cmapLapParaParamSet.h"
#include "cmapLapParaSolverState.h"
#include "cmapLapParaCalculationState.h"


namespace ParaCMapLAP
{


///
/// run this Solver
/// (measure idleTimeToFirstParaTask, idleTimeBetweenParaTasks, idleTimeAfterLastParaTask)
///
void
CMapLapParaSolver::run(
      )
{
   for(;;)
   {
      /***************************************************
       *  Wait a new ParaNode from ParaLoadCoordinator   *
       *  If Termination message is received, then break *
       ***************************************************/
      if( terminationMode != UG::TimeLimitTerminationMode )
      {
         if( !currentTask )
         {
            if( receiveNewTaskAndReactivate() == false )
            {
               break;
            }
         }
      }
      else
      {
         iReceiveMessages();
         if( !stayAliveAfterInterrupt )
         {
            idleTimeAfterLastParaTask = paraTimer->getElapsedTime() - previousStopTime;
            break;
         }
         continue;
      }

      /** set start time and ilde times */
      paraTaskStartTime = paraTimer->getElapsedTime();
      if( previousStopTime < 0.0 )
      {
         idleTimeToFirstParaTask = paraTaskStartTime - (idleTimeToWaitToken - previousIdleTimeToWaitToken);
      }
      else
      {
         idleTimeBetweenParaTasks += ( paraTaskStartTime - previousStopTime - ( idleTimeToWaitToken - previousIdleTimeToWaitToken ) );
      }


      /******************
       * start solving  *
       ******************/
      assert(!newTask);
      solve();
      /*****************************************************
      * notify completion of a calculation of a ParaNode  *
      *****************************************************/
      previousStopTime = paraTimer->getElapsedTime();
      double compTime = previousStopTime - paraTaskStartTime;
      previousIdleTimeToWaitToken = idleTimeToWaitToken;

      /****************************************************************************
      * send completion of calculation and update counters and accumulation time  *
      ****************************************************************************/
      if( paraParams->getBoolParamValue(UG::Deterministic) )
      {
         do
         {
            iReceiveMessages();
         } while ( !waitToken(getRank()) );
      }

      iReceiveMessages();     /** Before sending completion state, receiving message should be checked.
                               *   When subproblem terminated with no branch, solver lost a timing for receiving new node */
      sendCompletionOfCalculation(compTime);

      if( paraParams->getBoolParamValue(UG::Deterministic) )
      {
         waitAckCompletion();
         // if( hasToken() ) passToken();
         paraDetTimer->update(1.0);
         previousCommTime = paraDetTimer->getElapsedTime();
#ifdef _DEBUG_DET
         std::cout << previousCommTime << " run2 R." << getRank() << ": token passed" << std::endl;
#endif
         passToken(getRank());
      }

      /**************************
      * update current ParaTask *
      ***************************/
      assert( currentTask );
      delete currentTask;
      if( newTask )
      {
         currentTask = newTask;
         newTask = 0;
      }
      else
      {
         currentTask = 0;
      }
      if( terminationMode && !stayAliveAfterInterrupt )
      {
         idleTimeAfterLastParaTask = paraTimer->getElapsedTime() - previousStopTime;
         break;
      }
   }
   // send solver terminationstate in this point or deconstructor
   return;
}


///
///
///
void CMapLapParaSolver::solve(
      )
{
}


///
/// constructor
///
CMapLapParaSolver::CMapLapParaSolver(
      int argc,
      char **argv,
      UG::ParaComm     *comm,
      CMapLapParaSolverLocalComm *inLocalComm,
      UG::ParaParamSet *inParaParamSet,
      UG::ParaInstance *inParaInstance,
      UG::ParaDeterministicTimer *inDetTimer,
      double timeOffset
      )
      :
      UG::ParaSolver(argc, argv, N_CMAP_LAP_TAGS, comm, inParaParamSet, inParaInstance, inDetTimer),
      localComm(inLocalComm),
      rank(-1),
      nThreadsPerRank(1),
      quiet(false),
      logname(0),
      problemName(0),
      terminationState(0),
      sendBasisRequest(false),
      previousIReceiveTime(0.0),
      previousShareVectorsTime(0.0),
      previousSendBasisTime(0.0),
      receivedBasis(nullptr),
      receiveVectorsIsProcessed(false),
      receivedPackedVector(nullptr),
      currentMinObjValue(DBL_MAX),
      currentMinEnumCost(DBL_MAX),
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
      nSolverStateSent(0),
      osLog(0)
{

   DEF_CMAP_LAP_PARA_COMM(cmapLapParaComm, comm);
   cmapLapParaComm->setLocalComm(inLocalComm);

   notificationInterval = paraParams->getRealParamValue(UG::NotificationInterval);
   messageQueueSizeHistory = MessageQueueSizeHistory();

   receiveVectorsData.resize(2);

   ///
   /// create timer for this BbParaSolver
   ///
   localMessageHandler = new CMapLapMessageHandlerFunctionPointer[N_CMAP_LAP_TAGS];
   for( int i = 0; i < N_CMAP_LAP_TAGS; i++ )
   {
      localMessageHandler[i] = 0;
   }

   localMessageHandler[UG::TagTask] = &ParaCMapLAP::CMapLapParaSolver::processLocalTagTask;
   localMessageHandler[UG::TagSolution] = &ParaCMapLAP::CMapLapParaSolver::processLocalTagSolution;
   localMessageHandler[UG::TagNotificationId] = &ParaCMapLAP::CMapLapParaSolver::processLocalTagNotificationId;
   localMessageHandler[UG::TagTerminateRequest] = &ParaCMapLAP::CMapLapParaSolver::processLocalTagTerminateRequest;
   localMessageHandler[UG::TagInterruptRequest] = &ParaCMapLAP::CMapLapParaSolver::processLocalTagInterruptRequest;
   localMessageHandler[TagTimeLimitRequest] = &ParaCMapLAP::CMapLapParaSolver::processLocalTagTimeLimitRequest;
   localMessageHandler[TagCMapLapPackedVector] = &ParaCMapLAP::CMapLapParaSolver::processLocalTagCMapLapPackedVector;
   localMessageHandler[TagBasisEnumCost] = &ParaCMapLAP::CMapLapParaSolver::processLocalTagBasisEnumCost;
   localMessageHandler[TagBasisRequest] = &ParaCMapLAP::CMapLapParaSolver::processLocalTagBasisRequest;
   localMessageHandler[TagBasis] = &ParaCMapLAP::CMapLapParaSolver::processLocalTagBasis;
   localMessageHandler[TagUpdateNotificationInterval] = &ParaCMapLAP::CMapLapParaSolver::processLocalTagUpdateNotificationInterval;

   if( comm )
   {
      CMapLapMessageHandlerFunctionPointer *cmapLapMessageHandler = reinterpret_cast<CMapLapMessageHandlerFunctionPointer *>(messageHandler);

      cmapLapMessageHandler[UG::TagTask] = &ParaCMapLAP::CMapLapParaSolver::processTagTask;
      cmapLapMessageHandler[UG::TagTaskReceived] = 0;
      cmapLapMessageHandler[UG::TagSolution] = &ParaCMapLAP::CMapLapParaSolver::processTagSolution;
      cmapLapMessageHandler[UG::TagNotificationId] = &ParaCMapLAP::CMapLapParaSolver::processTagNotificationId;
      cmapLapMessageHandler[UG::TagTerminateRequest] = &ParaCMapLAP::CMapLapParaSolver::processTagTerminateRequest;
      cmapLapMessageHandler[UG::TagInterruptRequest] = &ParaCMapLAP::CMapLapParaSolver::processTagInterruptRequest;
      cmapLapMessageHandler[TagTimeLimitRequest] = &ParaCMapLAP::CMapLapParaSolver::processTagTimeLimitRequest;
      cmapLapMessageHandler[UG::TagToken] = &ParaCMapLAP::CMapLapParaSolver::processTagToken;
      cmapLapMessageHandler[TagCMapLapPackedVector] = &ParaCMapLAP::CMapLapParaSolver::processTagCMapLapPackedVector;
      cmapLapMessageHandler[TagBasisEnumCost] = &ParaCMapLAP::CMapLapParaSolver::processTagBasisEnumCost;
      cmapLapMessageHandler[TagBasisRequest] = &ParaCMapLAP::CMapLapParaSolver::processTagBasisRequest;
      cmapLapMessageHandler[TagBasis] = &ParaCMapLAP::CMapLapParaSolver::processTagBasis;
      cmapLapMessageHandler[TagUpdateNotificationInterval] = &ParaCMapLAP::CMapLapParaSolver::processTagUpdateNotificationInterval;

      rank = comm->getRank();
      CMapLapParaInstance *cmapLapParaInstance = dynamic_cast< CMapLapParaInstance *>(paraInstance);

      paraTimer->setOffset(timeOffset);
      /* Initialize the CMAP_LAP environment */
      /*********
       * Setup *
       *********/
      /* initialize CMAP_LAP */
      if( paraParams->getRealParamValue(UG::TimeLimit) > 0.0 )
      {
         // double timeRemains =  paraParams->getRealParamValue(UG::TimeLimit) - paraTimer->getElapsedTime() + 10.0;  // 10.0: timming issue
      }

      /* include communication point handler */

      /********************
       * Parse parameters *
       ********************/
      for( int i = 3; i < argc; ++i )   /** the first argument is runtime parameter file for ParaCMAP_LAP */
      {
         if( strcmp(argv[i], "-l") == 0 )
         {
            i++;
            if( i < argc )
               logname = argv[i];
            else
            {
               THROW_LOGICAL_ERROR1("missing log filename after parameter '-l'");
            }
         }
         else if ( strcmp(argv[i], "-ntpr") == 0 )
         {
            i++;
            if( i < argc )
            {
               nThreadsPerRank = atoi(argv[i]);
            }
            else
            {
               std::cerr << "missing number after parameter '-ntpr'" << std::endl;
               exit(1);
            }
         }
         else if( strcmp(argv[i], "-q") == 0 )
            quiet = true;
         // other arguments are omitted in Solver
      }

      /***********************************
       * create log file message handler *
       ***********************************/
      if( paraParams->getBoolParamValue(UG::Quiet) )
      {
         ofsLog.open("/dev/null", std::ios::app );
         if( !ofsLog )
         {
            std::cout << "Solver log file cannot open : file name = /dev/null" << std::endl;
            exit(1);
         }
         osLog = &ofsLog;
      }
      else
      {
         osLog = &std::cout;
         if( logname != NULL )
         {
            std::ostringstream os;
            os << logname << getRank();
            ofsLog.open(os.str().c_str(), std::ios::app );
            if( !ofsLog )
            {
               std::cout << "Solver log file cannot open : file name = " << os.str() << std::endl;
               exit(1);
            }
            osLog = &ofsLog;
            std::cout << "*********** LOG " << logname << std::endl;
         }
         if( quiet )
         {
            ofsLog.open("/dev/null", std::ios::app );
            if( !ofsLog )
            {
               std::cout << "Solver log file cannot open : file name = /dev/null" << std::endl;
               exit(1);
            }
            osLog = &ofsLog;
         }
      }

      int tempIsWarmStarted;
      double startTime = paraTimer->getElapsedTime();
      if( getThreadId() == 0 )
      {
         comm->bcast(&tempIsWarmStarted, 1, UG::ParaINT, 0);
         localComm->bcast(&tempIsWarmStarted, 1, UG::ParaINT, 0);
         warmStarted = (tempIsWarmStarted == 1);
      }
      else
      {
         localComm->bcast(&tempIsWarmStarted, 1, UG::ParaINT, 0);
         warmStarted = (tempIsWarmStarted == 1);
      }
      idleTimeToWaitAckCompletion += ( paraTimer->getElapsedTime() - startTime );

      /** set parameters for CMAP_LAP: this values are reseted before solving */

      /** create problem */

      /** set callback routines **/

      // save parameters
      problemName = new char[strlen(cmapLapParaInstance->getProbName())+1];
      strcpy(problemName, cmapLapParaInstance->getProbName());
      delete paraInstance;
      paraInstance = 0;

      #ifdef PARA_CMAP_LAP_DEBUG
      std::ostringstream os;
      os << "debugR" << getRank()<< ".txt";
      debugOfstream.open(os.str().c_str(), std::ios::app);
      #endif
   }
   else
   {
      /// creation of LC local solver
      rank = getThreadId();
      paraTimer = localComm->createParaTimer();
#ifdef _COMM_CPP11
      // In ParaTimerMpi, created time is start time
      paraTimer->init(localComm);
#endif

      // save parameters
      CMapLapParaInstance *cmapLapParaInstance = dynamic_cast< CMapLapParaInstance *>(paraInstance);
      problemName = new char[strlen(cmapLapParaInstance->getProbName())+1];
      strcpy(problemName, cmapLapParaInstance->getProbName());
   }

}


///
/// deconstructor
///
CMapLapParaSolver::~CMapLapParaSolver(
      )
{

   if( paraComm )
   {
      DEF_CMAP_LAP_PARA_COMM(cmapLapParaComm, paraComm);

      std::shared_ptr<CMapLapParaSolverTerminationState> termState{
         dynamic_cast<CMapLapParaSolverTerminationState *>
         (cmapLapParaComm->createParaSolverTerminationState(
            getThreadId(),
            terminationMode,
            getRank(),
            nParaTasksReceived,
            nParaTasksSolved,
            paraTimer->getElapsedTime(),
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
            0.0,  // detTime
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
            ))};

      double startTime = paraTimer->getElapsedTime();
      cmapLapParaComm->lockApp();
      termState->send(cmapLapParaComm, 0, UG::TagTerminated);
      cmapLapParaComm->unlockApp();
      idleTimeToWaitAckCompletion += ( paraTimer->getElapsedTime() - startTime );
   }

   /// Note: To not delete paraComm until all solvers whose threadId
   /// greater than 1, threadId 0 solver waits second teraminate request.
   if( getThreadId() == 0 )
   {
      int source = -1, tag = -1;
      // wait second termination request
      (void)paraComm->probe(&source, &tag);
      assert( tag == UG::TagTerminateRequest );
      paraComm->receive( NULL, 0, UG::ParaBYTE, source, UG::TagTerminateRequest);
   }

   // if( logname ) delete [] logname;  // do not delete, since logname points to argv
   if( problemName ) delete [] problemName;

   #ifdef PARA_CMAP_LAP_DEBUG
   debugOfstream.close();
   #endif

   if( localMessageHandler ) delete[] localMessageHandler;

}

int
CMapLapParaSolver::processTagTask(
      int source,    ///< source rank
      int tag        ///< TagNode
      )
{
   DEF_CMAP_LAP_PARA_COMM(cmapLapParaComm, paraComm);

   CMapLapParaTask *cmapLapParaTask = dynamic_cast<CMapLapParaTask *>(cmapLapParaComm->createParaTask());
   paraComm->lockApp();
   cmapLapParaTask->receive(cmapLapParaComm, source);
   paraComm->unlockApp();

   if( cmapLapParaTask->getThreadId() != getThreadId() )
   {
      std::shared_ptr<CMapLapParaTask> __cmapLapParaTask(cmapLapParaTask);
      __cmapLapParaTask->send(localComm, cmapLapParaTask->getThreadId());
      return 0;
   }

   if( currentTask )
   {
      newTask = cmapLapParaTask;
   }
   else
   {
      currentTask = cmapLapParaTask;
   }
   nParaTasksReceived++;

   return 0;
}

int
CMapLapParaSolver::processLocalTagTask(
      int source,    ///< source rank
      int tag        ///< TagNode
      )
{
   CMapLapParaTask *cmapLapParaTask = dynamic_cast<CMapLapParaTask *>(localComm->createParaTask());
   double startTime = paraTimer->getElapsedTime();
   cmapLapParaTask->receive(localComm, source);
   idleTimeToWaitAckCompletion += ( paraTimer->getElapsedTime() - startTime );

   if( currentTask )
   {
      newTask = cmapLapParaTask;
   }
   else
   {
      currentTask = cmapLapParaTask;
   }
   nParaTasksReceived++;

   return 0;
}

int
CMapLapParaSolver::processTagTaskReceived(
      int source,
      int tag
      )
{
   paraComm->lockApp();
   PARA_COMM_CALL(
         paraComm->receive( NULL, 0, UG::ParaBYTE, source, UG::TagTaskReceived )
         );
   paraComm->unlockApp();
   return 0;
}

///
/// (measure idleTimeToWaitSolution)
///
int
CMapLapParaSolver::processTagSolution(
      int source,
      int tag
      )
{
   assert( getThreadId() == 0 );

   DEF_CMAP_LAP_PARA_COMM(cmapLapParaComm, paraComm);

   CMapLapParaSolution *sol = dynamic_cast<CMapLapParaSolution *>(cmapLapParaComm->createParaSolution());
   double startTime = paraTimer->getElapsedTime();
   paraComm->lockApp();
   sol->receive(cmapLapParaComm, source);
   paraComm->unlockApp();
   if( sol->getThreadId() != 0 )
   {
      std::shared_ptr<CMapLapParaSolution> __sol(sol);
      __sol->send(localComm, sol->getThreadId());
      return 0;
   }
   idleTimeToWaitSolution += paraTimer->getElapsedTime() - startTime;
   if( globalBestIncumbentSolution )
   {
      if( EPSLE( sol->getObjectiveFunctionValue(), globalBestIncumbentSolution->getObjectiveFunctionValue(), DEFAULT_NUM_EPSILON ) )
         // sometime it is necessary to save the solution in solver side
         //   updateGlobalBestIncumbentValue( sol->getObjectiveFunctionValue() ) ) //  DO NOT UPDATE!!
         //   The timing of the update depends on solver used
      {
         delete globalBestIncumbentSolution;
         globalBestIncumbentSolution = sol;
      }
      else
      {
         delete sol;
      }
   }
   else
   {
      globalBestIncumbentSolution = sol;
   }
   return 0;
}

///
/// (measure idleTimeToWaitSolution)
///
int
CMapLapParaSolver::processLocalTagSolution(
      int source,
      int tag
      )
{
   assert( getThreadId() != 0 );

   CMapLapParaSolution *sol = dynamic_cast<CMapLapParaSolution *>(localComm->createParaSolution());
   double startTime = paraTimer->getElapsedTime();
   sol->receive(localComm, source);
   idleTimeToWaitSolution += paraTimer->getElapsedTime() - startTime;
   if( globalBestIncumbentSolution )
   {
      if( EPSLE( sol->getObjectiveFunctionValue(), globalBestIncumbentSolution->getObjectiveFunctionValue(), DEFAULT_NUM_EPSILON ) )
         // sometime it is necessary to save the solution in solver side
         //   updateGlobalBestIncumbentValue( sol->getObjectiveFunctionValue() ) ) //  DO NOT UPDATE!!
         //   The timing of the update depends on solver used
      {
         delete globalBestIncumbentSolution;
         globalBestIncumbentSolution = sol;
      }
      else
      {
         delete sol;
      }
   }
   else
   {
      globalBestIncumbentSolution = sol;
   }
   return 0;
}


int
CMapLapParaSolver::processTagNotificationId(
      int source,
      int tag
      )
{
   DEF_CMAP_LAP_PARA_COMM(cmapLapParaComm, paraComm);
   unsigned int notificationId[2];   /// notificationId[0]: notification ID
                                     /// notificationId[1]: thread ID

   paraComm->lockApp();
   assert(getThreadId() == 0);
   PARA_COMM_CALL(
         cmapLapParaComm->receive( notificationId, 2, UG::ParaUNSIGNED, source, UG::TagNotificationId)
         );
   paraComm->unlockApp();

   if( notificationId[1] != 0 )
   {
      localComm->send(notificationId, 2, UG::ParaUNSIGNED, notificationId[1], UG::TagNotificationId);
      return 1;
   }

   assert( notificationProcessed );
   if( notificationId[0] == notificationIdGenerator )
   {
      notificationProcessed = false;
   }
   else
   {
      THROW_LOGICAL_ERROR9(paraTimer->getElapsedTime(), " Rank = ", getRank(), ", Thread = ", getThreadId(),
            ", notificationId received is ", notificationId[0], ", but generator value is ", notificationIdGenerator);
   }
   return 0;
}

int
CMapLapParaSolver::processLocalTagNotificationId(
      int source,
      int tag
      )
{
   unsigned int notificationId[2];   /// notificationId[0]: notification ID
                                     /// notificationId[1]: thread ID
   PARA_COMM_CALL(
         localComm->receive( notificationId, 2, UG::ParaUNSIGNED, source, UG::TagNotificationId)
         );

   if( !notificationProcessed )
   {
      THROW_LOGICAL_ERROR9(paraTimer->getElapsedTime(), " notificationProcessed must be true, but got false; Rank = ", getRank(), ", Thread = ", getThreadId(),
            ", notificationId = ", notificationId[0],  ", notificationThread = ", notificationId[1]);
   }

   assert( notificationProcessed );
   assert( static_cast<int>(notificationId[1]) == getThreadId() );
   if( notificationId[0] == notificationIdGenerator )
   {
      notificationProcessed = false;
   }
   else
   {
      THROW_LOGICAL_ERROR9(paraTimer->getElapsedTime(), " Rank = ", getRank(), ", Thread = ", getThreadId(),
            ", notificationId received is ", notificationId[0], ", but generator value is ", notificationIdGenerator);
   }
   return 0;
}

int
CMapLapParaSolver::processTagTerminateRequest(
      int source,
      int tag
      )
{
   assert( getThreadId() == 0 );

   paraComm->lockApp();
   PARA_COMM_CALL(
         paraComm->receive( NULL, 0, UG::ParaBYTE, source, UG::TagTerminateRequest)
         );
   paraComm->unlockApp();

   for( int i = 1; i < localComm->getSize(); i++ )
   {
      localComm->send( NULL, 0, UG::ParaBYTE, i, UG::TagTerminateRequest);
   }

   terminationMode = UG::TimeLimitTerminationMode;
   stayAliveAfterInterrupt = false;
   return 0;
}

int
CMapLapParaSolver::processLocalTagTerminateRequest(
      int source,
      int tag
      )
{
   assert( getThreadId() != 0 );

   PARA_COMM_CALL(
         localComm->receive( NULL, 0, UG::ParaBYTE, source, UG::TagTerminateRequest)
         );

   terminationMode = UG::TimeLimitTerminationMode;
   stayAliveAfterInterrupt = false;
   return 0;
}


///
/// When solver receives processTagInterruptRequest, then it interrupts algorithm and
/// wait next task that will be sent from LoadCoordinator
///
int
CMapLapParaSolver::processTagInterruptRequest(
      int source,
      int tag
      )
{
   int threadId = -1;
   paraComm->lockApp();
   PARA_COMM_CALL(
         paraComm->receive( &threadId, 1, UG::ParaINT, source, UG::TagInterruptRequest)
         );
   paraComm->unlockApp();

   if( threadId != getThreadId() )
   {
      localComm->send(&threadId, 1, UG::ParaINT, threadId, UG::TagInterruptRequest);
      return 0;
   }

   terminationMode = UG::InterruptedTerminationMode;
   terminationState = -1;
   stayAliveAfterInterrupt = true;
   return 0;

}

///
/// When solver receives processTagInterruptRequest, then it interrupts algorithm and
/// wait next task that will be sent from LoadCoordinator
///
int
CMapLapParaSolver::processLocalTagInterruptRequest(
      int source,
      int tag
      )
{
   int threadId = -1;
   PARA_COMM_CALL(
         localComm->receive(&threadId, 1, UG::ParaINT, source, UG::TagInterruptRequest)
         );

   terminationMode = UG::InterruptedTerminationMode;
   terminationState = -1;
   stayAliveAfterInterrupt = true;
   return 0;
}


///
/// When solver receives processTagTimeLimitRequest, solver interrupts their algorithm
/// and processes to finalize
///
int
CMapLapParaSolver::processTagTimeLimitRequest(
      int source,
      int tag
      )
{
   int threadId = -1;
   paraComm->lockApp();
   PARA_COMM_CALL(
         paraComm->receive( &threadId, 1, UG::ParaINT, source, TagTimeLimitRequest)
         );
   paraComm->unlockApp();

   if( threadId != getThreadId() )
   {
      localComm->send(&threadId, 1, UG::ParaINT, threadId, TagTimeLimitRequest);
      return 0;
   }

   terminationMode = UG::TimeLimitTerminationMode;
   terminationState = -1;
   stayAliveAfterInterrupt = true;
   return 0;

}

///
/// When solver receives processTagTimeLimitRequest, solver interrupts their algorithm
/// and processes to finalize
///
int
CMapLapParaSolver::processLocalTagTimeLimitRequest(
      int source,
      int tag
      )
{
   int threadId = -1;
   PARA_COMM_CALL(
         localComm->receive(&threadId, 1, UG::ParaINT, source, TagTimeLimitRequest)
         );
   terminationMode = UG::TimeLimitTerminationMode;
   terminationState = -1;
   stayAliveAfterInterrupt = true;
   return 0;
}


int
CMapLapParaSolver::processTagToken(
      int source,
      int tag
      )
{
   assert( paraParams->getBoolParamValue(UG::Deterministic) );
   assert(localComm->getSize() == 1); /// only work for one thread
   int token[2];
   paraComm->lockApp();
   PARA_COMM_CALL(
         paraComm->receive( token, 2, UG::ParaINT, source, UG::TagToken )
         );
   paraComm->unlockApp();

   paraComm->setToken(getRank(), token);
   return 0;
}


///
/// process TagBasis
/// @return always 0 (for extension)
/// @note measure idleTimeToWaitBasis
///
int
CMapLapParaSolver::processTagBasis(
      int source,
      int tag
      )
{
   assert( getThreadId() == 0 );

   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm);

   std::shared_ptr<CMapLapParaBasis> paraBasis(cmapLapParaComm->createCMapLapParaBasis());
   double startTime = paraTimer->getElapsedTime();
   paraComm->lockApp();
   paraBasis->receive(cmapLapParaComm, source);
   paraComm->unlockApp();
   idleTimeToWaitBasis += paraTimer->getElapsedTime() - startTime;

   if( paraBasis->getThreadId() != 0 )
   {
      paraBasis->send( localComm, paraBasis->getThreadId() );
   }
   else
   {
      receivedBasis = paraBasis;
   }
   return 0;
}

///
/// process TagBasis among local threads
/// @return always 0 (for extension)
///
int
CMapLapParaSolver::processLocalTagBasis(
      int source,
      int tag
      )
{
   assert( getThreadId() != 0 );

   std::shared_ptr<CMapLapParaBasis> paraBasis(localComm->createCMapLapParaBasis());

   paraBasis->receive(localComm, source);

   receivedBasis = paraBasis;
   return 0;
}


///
/// process TagUpdateNotificationInterval
/// @return always 0 (for extension)
/// @note measure idleTimeToWaitAckCompletion
///
int
CMapLapParaSolver::processTagUpdateNotificationInterval(
      int source,
      int tag
      )
{
   assert( getThreadId() == 0 );

   double startTime = paraTimer->getElapsedTime();
   paraComm->lockApp();
   PARA_COMM_CALL(
         paraComm->receive( &notificationInterval, 1, UG::ParaDOUBLE, source, TagUpdateNotificationInterval)
         );
   paraComm->unlockApp();
   idleTimeToWaitAckCompletion += paraTimer->getElapsedTime() - startTime;

   for( int i = 1; i < localComm->getSize(); i++ )
   {
      localComm->send(&notificationInterval, 1, UG::ParaDOUBLE, i, TagUpdateNotificationInterval);
   }
   return 0;
}

///
/// process TagUpdateNotificationInterval among local threads
/// @return always 0 (for extension)
///
int
CMapLapParaSolver::processLocalTagUpdateNotificationInterval(
      int source,
      int tag
      )
{
   assert( getThreadId() != 0 );
   PARA_COMM_CALL(
         localComm->receive( &notificationInterval, 1, UG::ParaDOUBLE, source, TagUpdateNotificationInterval)
         );
   return 0;
}

///
/// @note measure idleTimeToWaitPackedVector)
///
int
CMapLapParaSolver::processTagCMapLapPackedVector(
      int source,
      int tag
      )
{
   assert( getThreadId() == 0 );

   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm);

   std::shared_ptr<CMapLapParaPackedVector> pv{cmapLapParaComm->createCMapLapParaPackedVector()};
   double startTime = paraTimer->getElapsedTime();
   paraComm->lockApp();
   pv->receive(cmapLapParaComm, source);
   paraComm->unlockApp();
   idleTimeToWaitPackedVector += paraTimer->getElapsedTime() - startTime;

   if( pv->getThreadId() != 0 )
   {
      startTime = paraTimer->getElapsedTime();
      pv->send( localComm, pv->getThreadId() );
      idleTimeToWaitPackedVector += paraTimer->getElapsedTime() - startTime;
      return 0;
   }
   assert( !receivedPackedVector );
   receivedPackedVector = std::move(pv);
   receiveVectorsIsProcessed = false;

   // counter
   switch( getSolverType() )
   {
      case DeepBkz:
         nVectorsReceivedDeepBkz += receivedPackedVector->getNVectors();
         break;
      case Enum:
         nVectorsReceivedEnum += receivedPackedVector->getNVectors();
         break;
      case Sieve:
         nVectorsReceivedSieve += receivedPackedVector->getNVectors();
         break;
      case Undefined:
         break;
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaSolver::processTagCMapLapPackedVector: Invalid solver type = ",static_cast<int>(getSolverType()));
   }

   return 0;
}

///
/// @note measure idleTimeToWaitPackedVector
///
int
CMapLapParaSolver::processLocalTagCMapLapPackedVector(
      int source,
      int tag
      )
{
   assert( getThreadId() != 0 );

   std::shared_ptr<CMapLapParaPackedVector> pv{localComm->createCMapLapParaPackedVector()};

   double startTime = paraTimer->getElapsedTime();
   pv->receive(localComm, source);
   idleTimeToWaitPackedVector += paraTimer->getElapsedTime() - startTime;

   assert( !receivedPackedVector );
   receivedPackedVector = std::move(pv);
   receiveVectorsIsProcessed = false;

   return 0;
}


int
CMapLapParaSolver::processTagBasisEnumCost(
      int source,
      int tag
      )
{
   assert( getThreadId() == 0 );

   double enumCost = 0.0;
   double startTime = paraTimer->getElapsedTime();
   paraComm->lockApp();
   PARA_COMM_CALL(
         paraComm->receive( &enumCost, 1, UG::ParaDOUBLE, source, TagBasisEnumCost)
         );
   paraComm->unlockApp();
   idleTimeToWaitBasis += paraTimer->getElapsedTime() - startTime;
   for( int i = 1; i < localComm->getSize(); i++ )
   {
      localComm->send(&enumCost, 1, UG::ParaDOUBLE, i, TagBasisEnumCost);
   }

   currentMinEnumCost = std::min(currentMinEnumCost, enumCost);

   return 0;
}

int
CMapLapParaSolver::processLocalTagBasisEnumCost(
      int source,
      int tag
      )
{
   assert( getThreadId() != 0 );

   double enumCost = 0.0;
   PARA_COMM_CALL(
         localComm->receive( &enumCost, 1, UG::ParaDOUBLE, source, TagBasisEnumCost)
         );

   currentMinEnumCost = std::min(currentMinEnumCost, enumCost);

   return 0;
}


int
CMapLapParaSolver::processTagBasisRequest(
      int source,
      int tag
      )
{
   int threadId = -1;
   paraComm->lockApp();
   PARA_COMM_CALL(
         paraComm->receive(&threadId, 1, UG::ParaINT, source, TagBasisRequest)
         );
   paraComm->unlockApp();
   if( threadId != getThreadId() )
   {
      localComm->send(&threadId, 1, UG::ParaINT, threadId, TagBasisRequest);
      return 0;
   }
   sendBasisRequest = 1;
   return 0;
}

int
CMapLapParaSolver::processLocalTagBasisRequest(
      int source,
      int tag
      )
{
   int threadId = -1;
   PARA_COMM_CALL(
         localComm->receive(&threadId, 1, UG::ParaINT, source, TagBasisRequest)
         );

   sendBasisRequest = 1;
   return 0;
}


bool
CMapLapParaSolver::receiveNewTaskAndReactivate()
{
   for(;;)
   {
      int source;
      int tag;
      int status;
      bool receive;
      /*******************************************
       *  waiting for any message from anywhere  *
       *******************************************/
      if( paraParams->getBoolParamValue(UG::Deterministic) )
      {
         do
         {
            iReceiveMessages();
         } while( !waitToken(getRank()) );
         iReceiveMessages();
         paraDetTimer->update(1.0);
         previousCommTime = paraDetTimer->getElapsedTime();
#ifdef _DEBUG_DET
         std::cout << previousCommTime << " receiveNewNodeAndReactivate R." << getRank() << ": token passed" << std::endl;
#endif
         passToken(getRank());
      }
      else
      {
         if( getThreadId() ==  0 )
         {
            paraComm->lockApp();
            receive = paraComm->iProbe(&source, &tag);
            paraComm->unlockApp();
            if( receive )
            {
               if( messageHandler[tag] )
               {
                  status = (this->*messageHandler[tag])(source, tag);
                  if( status != 0 && status != 1 )
                  {
                     std::ostringstream s;
                     s << "[ERROR RETURN from Message Hander]:" <<  __FILE__ <<  "] func = "
                       << __func__ << ", line = " << __LINE__ << " - "
                       << "process tag = " << tag << std::endl;
                     abort();
                  }
               }
               else
               {
                  THROW_LOGICAL_ERROR3( "No message hander for ", tag, " is not registered" );
               }
            }
         }
         else
         {
            (void)localComm->probe(&source, &tag);
            if( localMessageHandler[tag] )
            {
               status = (this->*localMessageHandler[tag])(source, tag);
               if( status )
               {
                  std::ostringstream s;
                  s << "[ERROR RETURN from Message Hander]:" <<  __FILE__ <<  "] func = "
                    << __func__ << ", line = " << __LINE__ << " - "
                    << "process tag = " << tag << std::endl;
                  abort();
               }
            }
            else
            {
               THROW_LOGICAL_ERROR3( "No local message hander for ", tag, " is not registered" );
            }
         }
      }
      if( currentTask )
      {
         if( paraParams->getBoolParamValue(UG::Deterministic)  )
         {
            previousNotificationTime = paraDetTimer->getElapsedTime();
         }
         else
         {
            previousNotificationTime = paraTimer->getElapsedTime();
         }
         return true;
      }
      if( terminationMode && !stayAliveAfterInterrupt ) break;
   }
   return false;
}

void
CMapLapParaSolver::iReceiveMessages(
      )
{

   /************************************************************************
    * This fucntion is called from a callback routine of the target solver *
    * **********************************************************************/
   int source;
   int tag;
   int status;
   /************************************
    * check if there are some messages *
    ************************************/
   if( getThreadId() == 0 )
   {
      paraComm->lockApp();
      while( paraComm->iProbe(&source, &tag) )
      {
         paraComm->unlockApp();
         if( messageHandler[tag] )
         {
            status = (this->*messageHandler[tag])(source, tag);
            if( status != 0 && status != 1 )
            {
               std::ostringstream s;
               s << "[ERROR RETURN from Message Hander]:" <<  __FILE__ <<  "] func = "
                 << __func__ << ", line = " << __LINE__ << " - "
                 << "process tag = " << tag << std::endl;
               abort();
            }
         }
         else
         {
            sleep(2);
            THROW_LOGICAL_ERROR3( "No message hander for ", tag, " is not registered" );
         }
         paraComm->lockApp();
      }
      paraComm->unlockApp();
      previousIReceiveTime = paraTimer->getElapsedTime();
   }
   else
   {
      while( localComm->iProbe(&source, &tag) )
      {
         if( localMessageHandler[tag] )
         {
            status = (this->*localMessageHandler[tag])(source, tag);
            if( status )
            {
               std::ostringstream s;
               s << "[ERROR RETURN from Local Message Hander]:" <<  __FILE__ <<  "] func = "
                 << __func__ << ", line = " << __LINE__ << " - "
                 << "process tag = " << tag << std::endl;
               abort();
            }
         }
         else
         {
            THROW_LOGICAL_ERROR3( "No local message hander for ", tag, " is not registered" );
         }
      }
   }

#ifdef _COMM_MPI_WORLD
   double startTime = paraTimer->getElapsedTime();
   paraComm->lockApp();
   int messageQueueSize = dynamic_cast<CMapLapParaCommMpi *>(paraComm)->testAllIsends();
   paraComm->unlockApp();
   idleTimeToWaitIsend += ( paraTimer->getElapsedTime() - startTime );
   messageQueueSizeHistory.push(messageQueueSize);
#endif
}


void
CMapLapParaSolver::waitNotificationIdMessage(
      )
{
   int tag;
   int status;
   double startTime = paraTimer->getElapsedTime();

   if( getThreadId() == 0 )
   {
      for(;;)
      {
         paraComm->lockApp();
         paraComm->waitSpecTagFromSpecSource(0, UG::TagAny, &tag);
         paraComm->unlockApp();
         while( tag != UG::TagNotificationId )
         {
            idleTimeToWaitNotificationId += ( paraTimer->getElapsedTime() - startTime );
            status = (this->*messageHandler[tag])(0, tag);
            if( status )
            {
               std::ostringstream s;
               s << "[ERROR RETURN from Message Hander]:" <<  __FILE__ <<  "] func = "
                 << __func__ << ", line = " << __LINE__ << " - "
                 << "process tag = " << tag << std::endl;
               abort();
            }
            startTime = paraTimer->getElapsedTime();
            paraComm->lockApp();
            paraComm->waitSpecTagFromSpecSource(0, UG::TagAny, &tag);
            paraComm->unlockApp();
         }
         idleTimeToWaitNotificationId += ( paraTimer->getElapsedTime() - startTime );
         //
         // receive the notification Id message
         //
         status = processTagNotificationId(0, tag);
         if( status == 0 ) break;
         if( status != 1 )
         {
            std::ostringstream s;
            s << "[ERROR RETURN from Message Hander]:" <<  __FILE__ <<  "] func = "
              << __func__ << ", line = " << __LINE__ << " - "
              << "process tag = " << tag << std::endl;
            abort();
         }
      }
      waitingSpecificMessage = false;
   }
   else
   {
      std::ostringstream s;
      s << "waitNotificationIdMessage is called from Thread ID != 0, " <<  __FILE__ <<  "] func = "
        << __func__ << ", line = " << __LINE__ << std::endl;
      abort();
   }
}


///
/// check if a notification message needs to send or not
/// @return true if the notification message needs to send else false
///
bool
CMapLapParaSolver::notificationIsNecessary(
      )
{
   if( paraParams->getBoolParamValue(UG::Deterministic) )
   {
      return ( (paraDetTimer->getElapsedTime() - previousNotificationTime) > notificationInterval );
   }
   else
   {
      return ( (paraTimer->getElapsedTime() - previousNotificationTime) > notificationInterval );
   }
   return false;
}


///
/// check if solver needs to call iReceive()
/// @return true if solver needs to call iReceive() else false
///
bool
CMapLapParaSolver::iReceiveIsNecessary(
      )
{
   if( getThreadId() == 0 ) return true;  // thread = 0 solver should call iReceiveMessages frequency
   if( paraParams->getBoolParamValue(UG::Deterministic) )
   {
      if( ( paraDetTimer->getElapsedTime() - previousIReceiveTime )
            > paraParams->getRealParamValue(IReceiveInterval) )
      {
         return true;
      }
   }
   else
   {
      if( ( paraTimer->getElapsedTime() - previousIReceiveTime )
            > paraParams->getRealParamValue(IReceiveInterval) )
      {
         return true;
      }
   }
   return false;
}


///
/// check if solver needs to send and receive vectors
/// @return true if solver needs to send and receive vectors else false
///
bool
CMapLapParaSolver::shareVectorsIsNecessary(
      )
{
   if( paraParams->getBoolParamValue(UG::Deterministic) )
   {
      if( ( paraDetTimer->getElapsedTime() - previousShareVectorsTime )
            > paraParams->getRealParamValue(ShareVectorsInterval) )
      {
         return true;
      }
   }
   else
   {
      if( ( paraTimer->getElapsedTime() - previousShareVectorsTime )
            > paraParams->getRealParamValue(ShareVectorsInterval) )
      {
         return true;
      }
   }
   return false;
}



///
/// wait a notification id message if it is needed to synchronize with LoadCoordinaor
///
void
CMapLapParaSolver::waitMessageIfNecessary(
      )
{
   if( paraParams->getIntParamValue(UG::NotificationSynchronization) == 0 ||
          !rampUp ||
         ( paraParams->getIntParamValue(UG::NotificationSynchronization) == 1 ) )
   {
      waitNotificationIdMessage();
   }
}

void
CMapLapParaSolver::waitLocalMessageIfNecessary(
      )
{
   if( paraParams->getIntParamValue(UG::NotificationSynchronization) == 0 ||
          !rampUp ||
         ( paraParams->getIntParamValue(UG::NotificationSynchronization) == 1 ) )
   {
      int tag;
      int status;
      double startTime = paraTimer->getElapsedTime();
      localComm->waitSpecTagFromSpecSource(0, UG::TagNotificationId, &tag);
      idleTimeToWaitNotificationId += ( paraTimer->getElapsedTime() - startTime );
      startTime = paraTimer->getElapsedTime();
      //
      // receive the notification Id message
      //
      status = processLocalTagNotificationId(0,tag);
      if( status )
      {
         std::ostringstream s;
         s << "[ERROR RETURN from Local Message Hander]:" <<  __FILE__ <<  "] func = "
           << __func__ << ", line = " << __LINE__ << " - "
           << "process tag = " << tag << std::endl;
         abort();
      }
      waitingSpecificMessage = false;
      idleTimeToWaitNotificationId += ( paraTimer->getElapsedTime() - startTime );
   }
}


///
/// @note measure idleTimeToWaitAckCompletion
///
bool
CMapLapParaSolver::checkIfProjNormIsUpdated(
      double *projNorm,
      LatticeVector<int> *vector
      )
{
   double startTime = paraTimer->getElapsedTime();
   iReceiveMessages();
   idleTimeToWaitAckCompletion += paraTimer->getElapsedTime() - startTime;

   CMapLapParaSolution *cmapLapSolution = dynamic_cast<CMapLapParaSolution *>(globalBestIncumbentSolution);
   if( cmapLapSolution && cmapLapSolution->getObjectiveFunctionValue() < currentMinObjValue )
   {
      currentMinObjValue = cmapLapSolution->getObjectiveFunctionValue();
      *projNorm = currentMinObjValue;
      *vector = cmapLapSolution->getVector();
      return true;
   }
   else
   {
      return false;
   }
}


///
/// @brief send vector-request to Load Cooardinator
/// @param[in] nRequestVectors number of requested lattice vector
///
void
CMapLapParaSolver::requestVectors(
      int nRequestVectors
      )
{
   assert( nRequestVectors > 0 );
   if( receivedPackedVector || receiveVectorsIsProcessed ){ return; }

   receiveVectorsData[0] = nRequestVectors;
   receiveVectorsData[1] = getThreadId();

   double startTime = paraTimer->getElapsedTime();
   paraComm->lockApp();
   paraComm->send(&receiveVectorsData[0], 2, UG::ParaINT, 0, TagVectorRequest);
   paraComm->unlockApp();
   idleTimeToWaitPackedVector += paraTimer->getElapsedTime() - startTime;

   receiveVectorsIsProcessed = true;
   previousShareVectorsTime = paraTimer->getElapsedTime();
}


///
/// getter of received basis
/// @param[out] vectors received vectors as column-based matrix form
///
template<typename BasisFloat>
bool
CMapLapParaSolver::getReceivedVectors(
      LatticeBasis<BasisFloat> &vectors
      )
{
   assert( receivedPackedVector );
   int nReceiveVectors = receivedPackedVector->getNVectors();
   if( nReceiveVectors == 0 )
   {
      vectors.resize(0, 0);
      receiveVectorsIsProcessed = false;
      receivedPackedVector = nullptr;
      return false;
   }
   else
   {
      int dimension = receivedPackedVector->getDimension();
      vectors.resize(nReceiveVectors, dimension);
      for( int k = 0; k < nReceiveVectors; ++k )
      {
         for( int i = 0; i < dimension; ++i )
         {
            vectors(k, i) = receivedPackedVector->getVectorElements()[k]->vector(i);
         }
      }
      receiveVectorsIsProcessed = false;
      receivedPackedVector = nullptr;
      return true;
   }
}


///
/// send SolverState DeepBkz
///
void
CMapLapParaSolver::sendSolverState(
      LatticeBasis<int>& inBasis,    ///< current basis
      int      inCurrentBlockSize,   ///< current DeepBkz block size
      int      inTour,               ///< number of DeepBkz loop
      double   inElapsedTime         ///< elapsed time
      )
{
   CMapLapParaTask *cmapLapParaTask = dynamic_cast<CMapLapParaTask *>(currentTask);
   if( cmapLapParaTask->isLocalTask() ) return;
   if( paraParams->getIntParamValue(MaxSizeOfMessageQueue) > 0 &&
         paraComm->getNumOfMessagesWaitingToSend(0)
         > paraParams->getIntParamValue(MaxSizeOfMessageQueue) )
   {
      return;
   }

   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm );

   if( !notificationProcessed && !waitingSpecificMessage )
   {
      if( !paraParams->getBoolParamValue(NoWaitNotificationId) )
      {
         waitingSpecificMessage = true;
      }

      int *basis = 0;
      int basisRows = inBasis.rows();
      basis = new int [basisRows*inBasis.cols()];
      int index = 0;
      for( int row = 0; row < basisRows; row++ )
      {
         for( int col = 0; col < inBasis.cols(); col++ )
         {
            basis[index++] = inBasis(row, col);
         }
      }

      Lattice<int, double> L(inBasis);

      std::shared_ptr<CMapLapParaSolverState> solverState(
         dynamic_cast<CMapLapParaSolverState *>(cmapLapParaComm->createParaSolverState(
            ++notificationIdGenerator,
            currentTask->getLcId(),
            currentTask->getGlobalSubtaskIdInLc(),
            paraTimer->getElapsedTime(),
            getThreadId(),
            inBasis.cols(),  // dimension
            basisRows,
            messageQueueSizeHistory.getMeanSize(),
            messageQueueSizeHistory.getMaxSize(),
            basis,
            inCurrentBlockSize,
            inTour,
            inElapsedTime,
            L.shortestNorm(),
            L.approxFactor(),
            L.hermiteFactor(),
            L.rootHermiteFactor(),
            std::log(L.enumCost()),
            std::log(L.enumCost(L.GH)),
            L.slopeGSA(),
            L.slopeGSA(std::ceil(L.m/2)),
            L.logOrthogonalityDefect()
            )));
      double startTime = paraTimer->getElapsedTime();
      paraComm->lockApp();
      solverState->send(paraComm, 0, UG::TagSolverState);
      paraComm->unlockApp();
      idleTimeToWaitSolverState += paraTimer->getElapsedTime() - startTime;
      if( paraParams->getBoolParamValue(UG::Deterministic) )
      {
         previousNotificationTime = paraDetTimer->getElapsedTime();
      }
      else
      {
         previousNotificationTime = paraTimer->getElapsedTime();
      }
      if( waitingSpecificMessage )
      {
         notificationProcessed = true;
         startTime = paraTimer->getElapsedTime();
         if( getThreadId() == 0 )
         {
            waitMessageIfNecessary();       // inside of waitMessageIfNecessary(), paraComm->unlockApp() is issued
         }
         else
         {
            waitLocalMessageIfNecessary();  // inside of waitLocalMessageIfNecessary(), paraComm->unlockApp() is issued
         }
         idleTimeToWaitSolverState += paraTimer->getElapsedTime() - startTime;
      }
      nSolverStateSent++;
   }
}

///
/// send SolverState Enum
///
void
CMapLapParaSolver::sendSolverState(
      double   inElapsedTime,          ///< elapsed time
      LatticeVector<int> &inCoeff,     ///< coefficients of current search node
      int      inDepth,                ///< depth of current searched node
      double   inEnumCost,             ///< log of approximted nodes of enumeration tree with incumbent radius
      double   inShortestNorm,         ///< the shortest norm found
      long int inNumSearchedNodes      ///< number of searched nodes in the enumeration tree
      )
{
   CMapLapParaTask *cmapLapParaTask = dynamic_cast<CMapLapParaTask *>(currentTask);
   if( cmapLapParaTask->isLocalTask() ) return;
   if( paraParams->getIntParamValue(MaxSizeOfMessageQueue) > 0 &&
         paraComm->getNumOfMessagesWaitingToSend(0)
         > paraParams->getIntParamValue(MaxSizeOfMessageQueue) )
   {
      return;
   }

   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm);

   int *pcoeff = new int [inCoeff.size()];
   Eigen::Map<LatticeVector<int>>(&pcoeff[0], inCoeff.size()) = inCoeff;

   if( !notificationProcessed && !waitingSpecificMessage )
   {
      if( !paraParams->getBoolParamValue(NoWaitNotificationId) )
      {
         waitingSpecificMessage = true;
      }

      Lattice<int, double> L(dynamic_cast<CMapLapParaTask *>(currentTask)->getBasis());

      std::shared_ptr<CMapLapParaSolverState> solverState(
         dynamic_cast<CMapLapParaSolverState *>(cmapLapParaComm->createParaSolverState(
            ++notificationIdGenerator,
            currentTask->getLcId(),
            currentTask->getGlobalSubtaskIdInLc(),
            paraTimer->getElapsedTime(),
            getThreadId(),
            cmapLapParaTask->getDimension(),
            cmapLapParaTask->getNumBasisVectors(),
            messageQueueSizeHistory.getMeanSize(),
            messageQueueSizeHistory.getMaxSize(),
            inElapsedTime,
            pcoeff,
            inDepth,
            inEnumCost,
            inShortestNorm,
            L.approxFactor(inShortestNorm),
            L.hermiteFactor(inShortestNorm),
            L.rootHermiteFactor(inShortestNorm),
            inNumSearchedNodes
            )));
      double startTime = paraTimer->getElapsedTime();
      paraComm->lockApp();
      solverState->send(paraComm, 0, UG::TagSolverState);
      paraComm->unlockApp();
      idleTimeToWaitSolverState += paraTimer->getElapsedTime() - startTime;
      if( paraParams->getBoolParamValue(UG::Deterministic) )
      {
         previousNotificationTime = paraDetTimer->getElapsedTime();
      }
      else
      {
         previousNotificationTime = paraTimer->getElapsedTime();
      }
      if( waitingSpecificMessage )
      {
         notificationProcessed = true;
         if( getThreadId() == 0 )
         {
            waitMessageIfNecessary();
         }
         else
         {
            waitLocalMessageIfNecessary();
         }
      }
      nSolverStateSent++;
   }
}

///
/// send SolverState Sieve
///
void
CMapLapParaSolver::sendSolverState(
      double       inElapsedTime,       ///< elapsed time
      int          inBlockSize,         ///< block size
      long int     inNLoop,             ///< number of Sieve algorithm loops
      int          inListSize,          ///< size of List L
      int          inStackSize,         ///< size of Stack S
      int          inMaxListSize,       ///< maximum size used by List L
      int          inNCollisions,       ///< number of collisions
      double       inShortestNorm       ///< the shortest norm found
      )
{
   CMapLapParaTask *cmapLapParaTask = dynamic_cast<CMapLapParaTask *>(currentTask);
   if( cmapLapParaTask->isLocalTask() ) return;
   if( paraParams->getIntParamValue(MaxSizeOfMessageQueue) > 0 &&
         paraComm->getNumOfMessagesWaitingToSend(0)
         > paraParams->getIntParamValue(MaxSizeOfMessageQueue) )
   {
      return;
   }

   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm);
   if( !notificationProcessed && !waitingSpecificMessage )
   {

      if( !paraParams->getBoolParamValue(NoWaitNotificationId) )
      {
         waitingSpecificMessage = true;
      }

      Lattice<int, double> L(dynamic_cast<CMapLapParaTask *>(currentTask)->getBasis());

      std::shared_ptr<CMapLapParaSolverState> solverState(
         dynamic_cast<CMapLapParaSolverState *>(cmapLapParaComm->createParaSolverState(
            ++notificationIdGenerator,
            currentTask->getLcId(),
            currentTask->getGlobalSubtaskIdInLc(),
            paraTimer->getElapsedTime(),
            getThreadId(),
            0, // dimension
            0, // basisRows
            messageQueueSizeHistory.getMeanSize(),
            messageQueueSizeHistory.getMaxSize(),
            inElapsedTime,
            inBlockSize,
            inNLoop,
            inListSize,
            inStackSize,
            inMaxListSize,
            inNCollisions,
            inShortestNorm,
            L.approxFactor(inShortestNorm),
            L.hermiteFactor(inShortestNorm),
            L.rootHermiteFactor(inShortestNorm)
            )));
      double startTime = paraTimer->getElapsedTime();
      cmapLapParaComm->lockApp();
      solverState->send(paraComm, 0, UG::TagSolverState);
      cmapLapParaComm->unlockApp();
      idleTimeToWaitSolverState += paraTimer->getElapsedTime() - startTime;
      if( paraParams->getBoolParamValue(UG::Deterministic) )
      {
         previousNotificationTime = paraDetTimer->getElapsedTime();
      }
      else
      {
         previousNotificationTime = paraTimer->getElapsedTime();
      }
      if( waitingSpecificMessage )
      {
         notificationProcessed = true;
         if( getThreadId() == 0 )
         {
            waitMessageIfNecessary();
         }
         else
         {
            waitLocalMessageIfNecessary();
         }
      }
      nSolverStateSent++;
   }
}


///
/// send CalcutationState DeepBkz
///
void
CMapLapParaSolver::sendDeepBkzCalculationState(
      LatticeBasis<int>& inBasis,   ///< current basis
      int      inCurrentBlockSize,  ///< current DeepBkz block size
      int      inTour,              ///< number of DeepBkz loop
      double   inElapsedTime        ///< elapsed time
      )
{
   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm);

   Lattice<int, double> L(inBasis);

   if( cmapLapParaComm )
   {
      std::shared_ptr<CMapLapParaCalculationState> paraCalculationState(
         dynamic_cast<CMapLapParaCalculationState *>(cmapLapParaComm->createParaCalculationState(
            terminationState,
            getThreadId(),
            inCurrentBlockSize,
            inTour,
            inElapsedTime,
            L.shortestNorm(),
            L.approxFactor(),
            L.hermiteFactor(),
            L.rootHermiteFactor(),
            std::log(L.enumCost()),
            std::log(L.enumCost(L.GH)),
            L.slopeGSA(),
            L.slopeGSA(std::ceil(L.m/2)),
            L.logOrthogonalityDefect()
            )));
      double startTime = paraTimer->getElapsedTime();
      paraComm->lockApp();
      paraCalculationState->send(paraComm, 0, UG::TagCompletionOfCalculation);
      paraComm->unlockApp();
      idleTimeToWaitAckCompletion += ( paraTimer->getElapsedTime() - startTime );
   }
   else
   {
      std::shared_ptr<CMapLapParaCalculationState> paraCalculationState(
         dynamic_cast<CMapLapParaCalculationState *>(localComm->createParaCalculationState(
                  terminationState,
                  getThreadId(),
                  inCurrentBlockSize,
                  inTour,
                  inElapsedTime,
                  L.shortestNorm(),
                  L.approxFactor(),
                  L.hermiteFactor(),
                  L.rootHermiteFactor(),
                  std::log(L.enumCost()),
                  std::log(L.enumCost(L.GH)),
                  L.slopeGSA(),
                  L.slopeGSA(std::ceil(L.m/2)),
                  L.logOrthogonalityDefect()
                  )));
      double startTime = paraTimer->getElapsedTime();
      paraComm->lockApp();
      paraCalculationState->send(localComm, 0, UG::TagCompletionOfCalculation);
      paraComm->unlockApp();
      idleTimeToWaitAckCompletion += ( paraTimer->getElapsedTime() - startTime );
      // delete paraCalculationState;
   }
}


///
/// send CalcutationState Enum
///
void
CMapLapParaSolver::sendEnumCalculationState(
      LatticeBasis<int>& inBasis,   ///< current basis
      double   inElapsedTime,       ///< elapsed time
      double   inShortestNorm,      ///< the shortest norm found
      long int inNSearch            ///< number of search nodes
      )
{
   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm);

   Lattice<int, double> L(inBasis);

   if( cmapLapParaComm )
   {
      std::shared_ptr<CMapLapParaCalculationState> paraCalculationState(
         dynamic_cast<CMapLapParaCalculationState *>(cmapLapParaComm->createParaCalculationState(
                  terminationState,
                  getThreadId(),
                  inElapsedTime,
                  inShortestNorm,
                  L.approxFactor(inShortestNorm),
                  L.hermiteFactor(inShortestNorm),
                  L.rootHermiteFactor(inShortestNorm),
                  inNSearch
                  )));
      double startTime = paraTimer->getElapsedTime();
      paraComm->lockApp();
      paraCalculationState->send(paraComm, 0, UG::TagCompletionOfCalculation);
      paraComm->unlockApp();
      idleTimeToWaitAckCompletion += ( paraTimer->getElapsedTime() - startTime );
   }
   else
   {
      std::shared_ptr<CMapLapParaCalculationState> paraCalculationState(
         dynamic_cast<CMapLapParaCalculationState *>(localComm->createParaCalculationState(
                  terminationState,
                  getThreadId(),
                  inElapsedTime,
                  inShortestNorm,
                  L.approxFactor(inShortestNorm),
                  L.hermiteFactor(inShortestNorm),
                  L.rootHermiteFactor(inShortestNorm),
                  inNSearch
                  )));
      double startTime = paraTimer->getElapsedTime();
      paraComm->lockApp();
      paraCalculationState->send(localComm, 0, UG::TagCompletionOfCalculation);
      paraComm->unlockApp();
      idleTimeToWaitAckCompletion += ( paraTimer->getElapsedTime() - startTime );
      // delete paraCalculationState;
   }
}


///
/// send CalcutationState Sieve
///
void
CMapLapParaSolver::sendSieveCalculationState(
      LatticeBasis<int>& inBasis,      ///< current basis
      double       inElapsedTime,      ///< elapsed time
      int          inBlockSize,        ///< block size
      long int     inNLoop,            ///< number of Sieve algorithm loop
      int          inListSize,         ///< size of List L
      int          inStackSize,        ///< size of Stack S
      int          inMaxListSize,      ///< maximum size of List L up to the point
      int          inNCollisions,      ///< number of collision
      double       inShortestNorm      ///< the shortest norm found
      )
{
   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm);

   Lattice<int, double> L(inBasis);

   if( cmapLapParaComm )
   {
      std::shared_ptr<CMapLapParaCalculationState> paraCalculationState(
         dynamic_cast<CMapLapParaCalculationState *>(cmapLapParaComm->createParaCalculationState(
                  terminationState,
                  getThreadId(),
                  inElapsedTime,
                  inBlockSize,
                  inNLoop,
                  inListSize,
                  inStackSize,
                  inMaxListSize,
                  inNCollisions,
                  inShortestNorm,
                  L.approxFactor(inShortestNorm),
                  L.hermiteFactor(inShortestNorm),
                  L.rootHermiteFactor(inShortestNorm)
                  )));
      double startTime = paraTimer->getElapsedTime();
      paraComm->lockApp();
      paraCalculationState->send(paraComm, 0, UG::TagCompletionOfCalculation);
      paraComm->unlockApp();
      idleTimeToWaitAckCompletion += ( paraTimer->getElapsedTime() - startTime );
   }
   else
   {
      std::shared_ptr<CMapLapParaCalculationState> paraCalculationState(
         dynamic_cast<CMapLapParaCalculationState *>(localComm->createParaCalculationState(
                  terminationState,
                  getThreadId(),
                  inElapsedTime,
                  inBlockSize,
                  inNLoop,
                  inListSize,
                  inStackSize,
                  inMaxListSize,
                  inNCollisions,
                  inShortestNorm,
                  L.approxFactor(inShortestNorm),
                  L.hermiteFactor(inShortestNorm),
                  L.rootHermiteFactor(inShortestNorm)
                  )));
      double startTime = paraTimer->getElapsedTime();
      paraComm->lockApp();
      paraCalculationState->send(localComm, 0, UG::TagCompletionOfCalculation);
      paraComm->unlockApp();
      idleTimeToWaitAckCompletion += ( paraTimer->getElapsedTime() - startTime );
      // delete paraCalculationState;
   }
}


///
/// @brief wait UG::TagAckCompletion tag from rank 0
///
void
CMapLapParaSolver::waitAckCompletion(
      )
{
   double inIdleTime = paraTimer->getElapsedTime();
   int tag;
   // int hstatus;
   paraComm->waitSpecTagFromSpecSource(0, UG::TagAckCompletion, &tag);
   double current = paraTimer->getElapsedTime();
   double idleTime = current - inIdleTime;
   idleTimeToWaitAckCompletion += idleTime;
   if( paraParams->getBoolParamValue(UG::DynamicAdjustNotificationInterval) &&
         current < paraParams->getRealParamValue(UG::CheckpointInterval) &&
         idleTime > 1.0 )
   {
      paraParams->setRealParamValue(UG::NotificationInterval, notificationInterval + idleTime );
   }
   //
   // receive the notification Id message
   //
   PARA_COMM_CALL(
         paraComm->receive( NULL, 0, UG::ParaBYTE, 0, UG::TagAckCompletion)
         );
}


///
/// @brief send solution of vector
/// @param[in] v lattice vector
/// @param[in] objValue objective function value
/// @return true if solution is sent else false
///
bool
CMapLapParaSolver::sendSolution(
      LatticeVector<int> &v,  ///< lattice vector
      double objValue         ///< objective function value
      )
{
   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm);

   // create cmapLapParaSolution
   std::shared_ptr<CMapLapParaSolution> cmapLapParaSolution = nullptr;
   if( cmapLapParaComm )
   {
      cmapLapParaSolution = std::shared_ptr<CMapLapParaSolution>(
         cmapLapParaComm->createCMapLapParaSolution(getThreadId(), v, objValue
         ));
   }
   else
   {
      cmapLapParaSolution = std::shared_ptr<CMapLapParaSolution>(
         localComm->createCMapLapParaSolution(getThreadId(), v, objValue
         ));
   }

   // delete cmapLapParaSolution or (save it to globalBestIncumbentSolution and send it to Load Coordinator)
   if( !globalBestIncumbentSolution ||
       ( globalBestIncumbentSolution &&
         EPSLT( cmapLapParaSolution->getObjectiveFunctionValue(), globalBestIncumbentSolution->getObjectiveFunctionValue(), eps ) ) )
   {
      if( globalBestIncumbentSolution ) delete globalBestIncumbentSolution;
      globalBestIncumbentSolution = cmapLapParaSolution->clone(paraComm);
      globalBestIncumbentValue = cmapLapParaSolution->getObjectiveFunctionValue();
      assert( localIncumbentSolution == 0 );
      double startTime = paraTimer->getElapsedTime();
      if( paraComm )
      {
         paraComm->lockApp();
         cmapLapParaSolution->send(paraComm, 0);
         paraComm->unlockApp();
      }
      else
      {
         iReceiveMessages();
         cmapLapParaSolution->send(localComm, 0);
         // delete cmapLapParaSolution;
      }
      idleTimeToWaitSolution += paraTimer->getElapsedTime() - startTime;
      nImprovedIncumbent++;
      switch ( getSolverType() )
      {
         case DeepBkz:
            nVectorsSentDeepBkz++;
            break;
         case Enum:
            nVectorsSentEnum++;
            break;
         case Sieve:
            nVectorsSentSieve++;
            break;
         case Undefined:
            break;
         default:
            THROW_LOGICAL_ERROR2("CMapLapParaSolver::sendSolution: Invalid solver type = ",static_cast<int>(getSolverType()));
      }
      return true;
   }
   else
   {
      // delete cmapLapParaSolution;
      return false;
   }
}


///
/// @brief send inK vectors to Load Coordinator from solver
/// @param[in] vectors  sent vectors as marix form with column-major
/// @note measure idleTimeToWaitPackedVector
///
void
CMapLapParaSolver::sendVectors(
      LatticeBasis<int> &vectors
   )
{

   if( paraParams->getIntParamValue(MaxSizeOfMessageQueue) > 0 &&
         paraComm->getNumOfMessagesWaitingToSend(0)
         > paraParams->getIntParamValue(MaxSizeOfMessageQueue) )
   {
      return;
   }

   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm);
   std::shared_ptr<CMapLapParaPackedVector> pv(
         cmapLapParaComm->createCMapLapParaPackedVector(
            getThreadId(), vectors, rank, nThreadsPerRank
            )
         );
   double startTime = paraTimer->getElapsedTime();
   paraComm->lockApp();
   pv->send(cmapLapParaComm, 0);
   paraComm->unlockApp();
   idleTimeToWaitPackedVector += paraTimer->getElapsedTime() - startTime;

   previousShareVectorsTime = paraTimer->getElapsedTime();
}


///
/// @brief send basis to Load Coordinator
/// @param[in] basis lattice basis
/// @param[in] enumCost enumeration cost. It value is negative, enumeration cost is calculated in this function
/// @note measure idleTimeToWaitBasis
///
void
CMapLapParaSolver::sendBasis(
      LatticeBasis<int> &basis,
      double enumCost
      )
{
   if( paraParams->getIntParamValue(MaxSizeOfMessageQueue) > 0 &&
         paraComm->getNumOfMessagesWaitingToSend(0)
         > paraParams->getIntParamValue(MaxSizeOfMessageQueue) )
   {
      return;
   }
   if( enumCost < 0 )
   {
      enumCost = Lattice<int, double>(basis).enumCost();
   }

   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm);
   std::shared_ptr<CMapLapParaBasis> paraBasis(
         cmapLapParaComm->createCMapLapParaBasis(getThreadId(), enumCost, basis)
         );
   double startTime = paraTimer->getElapsedTime();
   paraComm->lockApp();
   paraBasis->send(cmapLapParaComm, 0);
   paraComm->unlockApp();
   idleTimeToWaitBasis += paraTimer->getElapsedTime() - startTime;
   if( getSolverType() == DeepBkz ) nBasesSentDeepBkz++;
   sendBasisRequest = false;
}


///
/// instantiation
///
template LatticeBasis<int> CMapLapParaSolver::getReceivedBasis<int>();
template LatticeBasis<long int> CMapLapParaSolver::getReceivedBasis<long int>();
template bool CMapLapParaSolver::getReceivedVectors<int>(LatticeBasis<int> &vectors);
template bool CMapLapParaSolver::getReceivedVectors<long int>(LatticeBasis<long int> &vectors);

} // namespace ParaCMapLAP
