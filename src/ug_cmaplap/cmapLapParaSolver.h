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

/**@file    cmapLapParaSolver.h
 * @brief   ParaSolver extension for CMAP-LAP
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_SOLVER_H__
#define __CMAP_LAP_PARA_SOLVER_H__


#include <vector>
#include <deque>
#include "ug/paraSolver.h"
#include "ug/paraPthLock.h"
#include "cmapLapParaDef.h"
#include "cmapLapParaBasis.h"
#include "cmapLapParaSolverLocalComm.h"


#define ENFORCED_THRESHOLD 5


namespace ParaCMapLAP
{

// class CMapLapParaCommPointHdlr;



///
/// @class MessageQueueSizeHistory
/// @brief Recored history of message queue
///
class MessageQueueSizeHistory
{

private:
   std::deque<int>   queueSizes;
   int               maxHistorySize;
   int               maxQueueSize;
   int               totalSizeInQueue;

public:
   MessageQueueSizeHistory(
         int inMaxHistorySize=10
         ) : maxHistorySize(inMaxHistorySize),
             maxQueueSize(0),
             totalSizeInQueue(0)
   {}
   virtual ~MessageQueueSizeHistory(){}

   virtual void push(int queueSize)
   {
      queueSizes.push_back(queueSize);
      maxQueueSize = std::max(maxQueueSize, queueSize);
      totalSizeInQueue += queueSize;
      if( static_cast<int>(queueSizes.size()) > maxHistorySize )
      {
         totalSizeInQueue -= queueSizes.front();
         queueSizes.pop_front();
      }
   }
   virtual int getMaxSize()
   {
      return maxQueueSize;
   }
   virtual double getMeanSize()
   {
      if( queueSizes.empty() ) return 0;
      return totalSizeInQueue / queueSizes.size();
   }
};



///
/// @class CMapLapParaSolver
///
class CMapLapParaSolver : public UG::ParaSolver
{

protected:

   typedef int(CMapLapParaSolver::*CMapLapMessageHandlerFunctionPointer)(int, int);

   CMapLapMessageHandlerFunctionPointer *localMessageHandler;     ///< table for message handlers

   CMapLapParaSolverLocalComm   *localComm;

   int                         rank;

   UG::Lock                    callbackLock;

   int                         nThreadsPerRank;             ///< number of threads in a solver
   bool                        quiet;
   char                        *logname;                    ///< log name for each solver
   char                        *problemName;                ///< keep the name for restart, in case of file read
   int                         terminationState;            ///< termination state, 0: normal, -1: interrupted
   double                      notificationInterval;        ///< interval time of notification to LoadCoordinator
   bool                        sendBasisRequest;            ///< send basis request, 1: has request form LC, 0: donot has request

   double                      previousIReceiveTime;        ///< previous iReceive time
   double                      previousShareVectorsTime;    ///< previous send and receive vectors time
   double                      previousSendBasisTime;       ///< previous send basis time

   std::shared_ptr<CMapLapParaBasis> receivedBasis;         ///< basis received from LC
   bool                        receiveVectorsIsProcessed;   ///< request state, 1: requesting, 0: otherwise
   std::vector<int>            receiveVectorsData;          ///< buffer for TagVectorRequest

   std::shared_ptr<CMapLapParaPackedVector> receivedPackedVector; ///< packed vector received from LC

   double                      currentMinObjValue;
   double                      currentMinEnumCost;

   double                      idleTimeToWaitSolverState;   ///< idle time to wait solver state
   double                      idleTimeToWaitPackedVector;  ///< idle time to wait packed vector
   double                      idleTimeToWaitSolution;      ///< idle time to wait solution
   double                      idleTimeToWaitBasis;         ///< idle time to wait basis
   double                      idleTimeToWaitIsend;         ///< idle time to wait Isend

   int                         nParaTasksDeepBkzReceived;   ///< number of DeepBkz ParaTasks received in this solver
   int                         nParaTasksEnumReceived;      ///< number of Enum ParaTasks received in this solver
   int                         nParaTasksSieveReceived;     ///< number of Sieve ParaTasks received in this solver
   double                      runningTimeDeepBkz;          ///< this solver running time of DeepBkz
   double                      runningTimeEnum;             ///< this solver running time of Enum
   double                      runningTimeSieve;            ///< this solver running time of Sieve
   int                         nVectorsReceivedDeepBkz;     ///< number of vectors received in DeepBkz
   int                         nVectorsReceivedEnum;        ///< number of vectors received in Enum
   int                         nVectorsReceivedSieve;       ///< number of vectors received in Sieve
   int                         nVectorsSentDeepBkz;         ///< number of vectors sent in DeepBkz
   int                         nVectorsSentEnum;            ///< number of vectors sent in Enum
   int                         nVectorsSentSieve;           ///< number of vectors sent in Sieve
   int                         nBasesSentDeepBkz;           ///< number of vectors sent in DeepBkz
   int                         nSolverStateSent;            ///< number of solver states sent

   MessageQueueSizeHistory     messageQueueSizeHistory;     ///< history of the message queue size

   std::ofstream               ofsLog;
   std::ostream                *osLog;

   #ifdef PARA_CMAP_LAP_DEBUG
   std::ofstream debugOfstream;  ///< ofstream for debugging
   #endif

   virtual void solve();

   ///-------------------
   /// Message handlers
   ///-------------------

   ///
   /// process TagNode
   /// @return always 0 (for extension)
   ///
   virtual int processTagTask(
         int source,    ///< source rank
         int tag        ///< TagNode
         );

   ///
   /// process TagNode among local threads
   /// @return always 0 (for extension)
   ///
   virtual int processLocalTagTask(
         int source,    ///< source rank
         int tag        ///< TagNode
         );

   ///
   /// process TagTaskReceived
   /// @return always 0 (for extension)
   ///
   virtual int processTagTaskReceived(
         int source,    ///< source rank
         int tag        ///< UG::TagTaskReceived
         );

   ///
   /// process UG::TagSolution
   /// @return always 0 (for extension)
   ///
   virtual int processTagSolution(
         int source,    ///< source rank
         int tag        ///< UG::TagSolution
         );

   ///
   /// process UG::TagSolution among local threads
   /// @return always 0 (for extension)
   ///
   virtual int processLocalTagSolution(
         int source,    ///< source rank
         int tag        ///< UG::TagSolution
         );

   ///
   /// process TagNotificationId
   /// @return always 0 (for extension)
   ///
   virtual int processTagNotificationId(
         int source,    ///< source rank
         int tag        ///< UG::TagNotificationId
         );

   ///
   /// process TagNotificationId among local threads
   /// @return always 0 (for extension)
   ///
   virtual int processLocalTagNotificationId(
         int source,    ///< source rank
         int tag        ///< UG::TagNotificationId
         );

   ///
   /// process TagTerminateRequest
   /// @return always 0 (for extension)
   ///
   virtual int processTagTerminateRequest(
         int source,    ///< source rank
         int tag        ///< UG::TagTerminateRequest
         );

   ///
   /// process TagTerminateRequest among local threads
   /// @return always 0 (for extension)
   ///
   virtual int processLocalTagTerminateRequest(
         int source,    ///< source rank
         int tag        ///< UG::TagTerminateRequest
         );


   ///
   /// process TagInterruptRequest
   /// @return always 0 (for extension)
   ///
   virtual int processTagInterruptRequest(
         int source,    ///< source rank
         int tag        ///< UG::TagInterruptRequest
         );

   ///
   /// process TagInterruptRequest among local threads
   /// @return always 0 (for extension)
   ///
   virtual int processLocalTagInterruptRequest(
         int source,    ///< source rank
         int tag        ///< UG::TagInterruptRequest
         );


   ///
   /// process TagTimeLimitRequest
   /// @return always 0 (for extension)
   ///
   virtual int processTagTimeLimitRequest(
         int source,    ///< source rank
         int tag        ///< TagTimeLimitRequest
         );

   ///
   /// process TagTimeLimitRequest among local threads
   /// @return always 0 (for extension)
   ///
   virtual int processLocalTagTimeLimitRequest(
         int source,    ///< source rank
         int tag        ///< TagTimeLimitRequest
         );


   ///
   /// process TagToken
   /// @return always 0 (for extension)
   ///
   virtual int processTagToken(
         int source,    ///< source rank
         int tag        ///< TagToken
         );


   ///
   /// process TagBasis
   /// @return always 0 (for extension)
   /// @note measure idleTimeToWaitBasis
   ///
   virtual int processTagBasis(
         int source,    ///< source rank
         int tag        ///< TagBasis
         );

   ///
   /// process TagBasis among local threads
   /// @return always 0 (for extension)
   ///
   virtual int processLocalTagBasis(
         int source,    ///< source rank
         int tag        ///< TagBasis
         );


   ///
   /// process TagUpdateNotificationInterval
   /// @return always 0 (for extension)
   /// @note measure idleTimeToWaitAckCompletion
   ///
   virtual int processTagUpdateNotificationInterval(
         int source,    ///< source rank
         int tag        ///< TagBasis
         );

   ///
   /// process TagUpdateNotificationInterval among local threads
   /// @return always 0 (for extension)
   ///
   virtual int processLocalTagUpdateNotificationInterval(
         int source,    ///< source rank
         int tag        ///< TagBasis
         );


   ///
   /// process TagCMapLapPackedVector
   /// @return always 0 (for extension)
   ///
   virtual int processTagCMapLapPackedVector(
         int source,    ///< source rank
         int tag        ///< TagPackedVector
         );

   ///
   /// process TagCMapLapPackedVector among local threads
   /// @return always 0 (for extension)
   ///
   virtual int processLocalTagCMapLapPackedVector(
         int source,    ///< source rank
         int tag        ///< TagPackedVector
         );

   ///
   /// process TagBasisEnumCost
   /// @return always 0 (for extension)
   /// @note measure idleTimeToWaitBasis
   ///
   virtual int processTagBasisEnumCost(
         int source,    ///< source rank
         int tag        ///< TagBasisEnumCost
         );

   ///
   /// process TagBasisEnumCost among local threads
   /// @return always 0 (for extension)
   ///
   virtual int processLocalTagBasisEnumCost(
         int source,    ///< source rank
         int tag        ///< TagBasisEnumCost
         );


   ///
   /// process TagBasisRequest
   /// @return always 0 (for extension)
   ///
   virtual int processTagBasisRequest(
         int source,    ///< source rank
         int tag        ///< TagBasisRequest
         );

   ///
   /// process TagBasisRequest among local threads
   /// @return always 0 (for extension)
   ///
   virtual int processLocalTagBasisRequest(
         int source,    ///< source rank
         int tag        ///< TagBasisRequest
         );


   ///
   /// wait for receiving a new node and reactivate solver
   /// @return true if a new node is received, false esle
   ///
   virtual bool receiveNewTaskAndReactivate(
         );

   ///
   /// wait notification id message to synchronized with LoadCoordinator
   ///
   virtual void waitNotificationIdMessage(
         );

   ///
   /// wait notification id message to synchronized with LoadCoordinator
   ///
   virtual void waitLocalMessageIfNecessary(
         );

   ///
   /// send completion of calculation
   ///
   virtual void sendCompletionOfCalculation(
         double stopTime       ///< stopping time
         )
   {
   }

   ///
   /// wait ack completion to synchronized with LoadCoordinator
   ///
   virtual void waitAckCompletion(
         );

public:
   CMapLapParaSolver(
         int argc,
         char **argv,
         UG::ParaComm *comm,
         CMapLapParaSolverLocalComm *localComm,
         UG::ParaParamSet *paraParamSet,
         UG::ParaInstance *paraInstance,
         UG::ParaDeterministicTimer *detTimer,
         double timeOffset);

   virtual ~CMapLapParaSolver(
         );

   ///
   /// @return rank of this solver
   ///
   virtual int getRank(
         )
   {
      return rank;
   }

   ///
   /// @return thread id of this solver
   ///
   virtual int getThreadId(
         )
   {
      return localComm->getThreadId();
   }

   ///
   /// @return parameter set
   ///
   virtual UG::ParaParamSet *getParaParamSet(
         )
   {
      return paraParams;
   }

   ///
   /// run this Solver
   ///
   using ParaSolver::run;
   virtual void run();

   ///
   /// non-blocking receive messages
   ///
   virtual void iReceiveMessages(
         );

   ///
   /// check if a notification message needs to send or not
   /// @return true if the notification message needs to send else false
   ///
   virtual bool notificationIsNecessary(
         );

   ///
   /// check if solver needs to call iReceive()
   /// @return true if solver needs to call iReceive() else false
   ///
   virtual bool iReceiveIsNecessary(
         );

   ///
   /// check if solver needs to send and receive vectors
   /// @return true if solver needs to send and receive vectors else false
   ///
   virtual bool shareVectorsIsNecessary(
         );

   ///
   /// wait a notification id message if it is needed to synchronize with LoadCoordinaor
   ///
   virtual void waitMessageIfNecessary(
         );


   ///
   /// @return true if this solver has received basis else false
   ///
   virtual bool hasReceivedBasis(
         )
   {
      return ( receivedBasis != nullptr );
   }

   ///
   /// getter of received basis
   /// @return Eigen::LatticeBasis<int>
   ///
   template<typename BasisFloat> LatticeBasis<BasisFloat> getReceivedBasis(
         )
   {
      assert( receivedBasis != nullptr );
      LatticeBasis<BasisFloat> basis = receivedBasis->getBasis().cast<BasisFloat>();
      receivedBasis = nullptr;
      return basis;
   }


   ///
   /// @brief send vector-request to Load Cooardinator
   /// @param[in] nRequestVectors number of requested lattice vector
   ///
   virtual void requestVectors(
         int nRequestVectors
         );

   ///
   /// @return true if this solver has received vectors else false
   ///
   virtual bool hasReceivedVectors(
         )
   {
      return ( receivedPackedVector != nullptr );
   }

   ///
   /// getter of received basis
   /// @param[out] vectors received vectors as column-based matrix form
   /// @param[out] nReceiveVectors number of received vectors
   ///
   template<typename BasisFloat> bool getReceivedVectors(
         LatticeBasis<BasisFloat> &vectors
         );

   ///
   /// @return solverType of current task
   ///
   virtual SolverType getSolverType(
         )
   {
      assert( currentTask );
      return (dynamic_cast<CMapLapParaTask *>(currentTask))->getSolverType();
   }

   virtual void writeCurrentTaskProblem(const std::string& filename){ throw "** writeCurrentNodeProblem is called **"; };
   virtual void tryNewSolution(UG::ParaSolution *sol){}
   virtual void writeSubproblem(){}
   virtual bool wasTerminatedNormally(){ return true; }

   ///
   /// send Solver State DeepBkz
   ///
   virtual void sendSolverState(
      LatticeBasis<int>& inBasis,   ///< current basis
      int      inCurrentBlockSize,  ///< current DeepBkz block size
      int      inTour,              ///< number of DeepBkz loop
      double   inElapsedTime        ///< elapsed time
      );

   ///
   /// send Solver State Enum
   ///
   virtual void sendSolverState(
      double   inElapsedTime,       ///< elapsed time
      LatticeVector<int> &inCoeff,  ///< coefficients of current search node
      int      inStart,             ///< the index of current search node in the enumeration tree
      double   inEnumCost,          ///< approximated time of full Enum
      double   inShortestNorm,      ///< the shortest norm found
      long int inNumSearchedNodes   ///< number of searched nodes in the enumeration tree
      );

   ///
   /// send Solver State Sieve
   ///
   virtual void sendSolverState(
      double       inElapsedTime,       ///< elapsed time
      int          inBlockSize,         ///< block size
      long int     inNLoop,             ///< number of Sieve algorithm loops
      int          inListSize,          ///< size of List L
      int          inStackSize,         ///< size of Stack S
      int          inMaxListSize,       ///< maximum size used by List L
      int          inNCollisions,       ///< number of collisions
      double       inShortestNorm       ///< the shortest norm found
      );


   ///
   /// send CalcutationState DeepBkz
   ///
   virtual void sendDeepBkzCalculationState(
         LatticeBasis<int>& inBasis,    ///< current basis
         int      inCurrentBlockSize,   ///< current DeepBkz block size
         int      inTour,               ///< number of DeepBkz loop
         double   inElapsedTime         ///< elapsed time
         );


   ///
   /// send CalcutationState Enum
   ///
   virtual void sendEnumCalculationState(
         LatticeBasis<int>& inBasis,   ///< current basis
         double   inElapsedTime,       ///< elapsed time
         double   inShortestNorm,      ///< the shortest norm found
         long int inNSearch            ///< number of search nodes
         );


   ///
   /// send CalcutationState Sieve
   ///
   virtual void sendSieveCalculationState(
         LatticeBasis<int>& inBasis,      ///< current basis
         double       inElapsedTime,      ///< elapsed time
         int          inBlockSize,        ///< block size
         long int     inNLoop,            ///< number of Sieve algorithm loop
         int          inListSize,         ///< size of List L
         int          inStackSize,        ///< size of Stack S
         int          inMaxListSize,      ///< maximum size of List L up to the point
         int          inNCollisions,      ///< number of collision
         double       inShortestNorm      ///< the shortest norm found
         );


   ///
   /// check whether solver should interrupt or not
   /// @return true if solver should interrupt else false
   ///
   virtual bool isInterrupted(
      )
   {
      return ( terminationState == -1 );
   }

   ///
   /// check whether solver should send basis or not
   /// @return true if solver should interrupt else false
   ///
   virtual bool hasSendBasisRequest(
      )
   {
      return ( sendBasisRequest == true );
   }

   ///
   /// @note measure idleTimeToWaitAckCompletion
   ///
   virtual bool checkIfProjNormIsUpdated(
         double *projNorm,
         LatticeVector<int> *vector
         );

   ///
   /// @brief send solution of vector
   /// @brief send solution of vector
   /// @param[in] v lattice vector
   /// @param[in] objValue objective function value
   /// @return true if solution is sent else false
   ///
   virtual bool sendSolution(
         LatticeVector<int> &v,  ///< lattice vector
         double objValue         ///< objective function value
         );

   ///
   /// @brief send inK vectors to Load Coordinator from solver
   /// @param[in] vectors  sent vectors as marix form with column-major
   /// @note measure idleTimeToWaitPackedVector
   ///
   virtual void sendVectors(
         LatticeBasis<int> &vectors
         );

   ///
   /// @brief send basis to Load Coordinator
   /// @param[in] basis lattice basis
   /// @param[in] enumCost enumeration cost
   /// @note measure idleTimeToWaitBasis
   ///
   virtual void sendBasis(
         LatticeBasis<int> &basis,
         double enumCost=-1
         );

};

} // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_SOLVER_H__
