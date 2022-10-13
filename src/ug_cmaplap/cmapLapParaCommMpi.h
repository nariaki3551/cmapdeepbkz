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

/**@file    cmapLapParaCommMpi.h
 * @brief   CMapLapParaComm extension for MPI communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_COMM_MPI_H__
#define __CMAP_LAP_PARA_COMM_MPI_H__


#include <utility>
#include <memory>
#include <deque>
#include "ug/paraCommMpi.h"
#include "cmapLapParaParamSet.h"
#include "cmapLapParaInstanceMpi.h"
#include "cmapLapParaSolutionMpi.h"
namespace ParaCMapLAP { class VectorElement; }
namespace ParaCMapLAP { class CMapLapParaBasis; }
namespace ParaCMapLAP { class CMapLapParaSolver; }
namespace ParaCMapLAP { class CMapLapParaIsendRequest; }
namespace ParaCMapLAP { class CMapLapParaPackedVector; }
namespace ParaCMapLAP { class CMapLapParaSolverLocalComm; }


namespace ParaCMapLAP
{


class CMapLapParaCommMpi : public UG::ParaCommMpi
{

   static const char          *tagStringTable[];                  ///< tag name string table
   static thread_local CMapLapParaSolverLocalComm *localComm;     ///< for get rank and threadId when lock / unlock App
   std::deque<std::shared_ptr<CMapLapParaIsendRequest>> cmapLapISendRequestQueue;

   ///
   /// check if tag string table (for debugging) set up correctly
   /// @return true if tag string table is set up correctly, false otherwise
   ///
   virtual bool tagStringTableIsSetUpCorrectly(
         );

public:

   CMapLapParaCommMpi(
         )
   {
   }

   virtual ~CMapLapParaCommMpi(
         );

   ///
   /// setter of cmapLapParaSolverLocalComm
   ///
   virtual void setLocalComm(
         CMapLapParaSolverLocalComm *inLocalComm
      );

   ///
   /// lock UG application to synchronize with other threads
   ///
   // virtual void lockApp(
   //       );
   virtual void lockApp(
         )
   {
      applicationLockMutex.lock();
   }

   ///
   /// unlock UG application to synchronize with other threads
   ///
   // virtual void unlockApp(
   //       );
   virtual void unlockApp(
         )
   {
      applicationLockMutex.unlock();
   }

   virtual void init(int argc, char** argv)
   {
      UG::ParaCommMpi::init(argc,argv);
      myComm = MPI_COMM_WORLD;
      MPI_CALL(
         MPI_Comm_rank(myComm, &myRank)
      );
      MPI_CALL(
         MPI_Comm_size(myComm, &myCommSize)
      );
   }

   /////////////////////////////////
   ///
   /// transfer object factory
   ///
   /// /////////////////////////////

   ///
   /// create ParaCalculationState object by default constructor
   /// @return pointer to ParaCalculationState object
   ///
   virtual UG::ParaCalculationState *createParaCalculationState(
         );

   ///
   /// create ParaRacingRampUpParamSet object (Currently not used)
   /// @return pointer to ParaRacingRampUpParamSet object
   ///
   virtual UG::ParaRacingRampUpParamSet* createParaRacingRampUpParamSet(
         )
   {
      throw "** createParaRacingRampUpParamSet is called **";
   }

   ///
   /// create ParaTask object by default constructor
   /// @return pointer to ParaTask object
   ///
   virtual UG::ParaTask *createParaTask(
         );

   ///
   /// create ParaSolverState object by default constructor
   /// @return pointer to ParaSolverState object
   ///
   virtual UG::ParaSolverState *createParaSolverState(
         );

   virtual UG::ParaInstance *createParaInstance(
         )
   {
      return new CMapLapParaInstanceMpi();
   }

   virtual UG::ParaSolution *createParaSolution(
         )
   {
      return new CMapLapParaSolutionMpi();
   }

   virtual CMapLapParaInstance *createCMapLapParaInstance(
         char      *inProbFileName,
         char      *inProbName
         )
   {
      return new CMapLapParaInstanceMpi(inProbFileName, inProbName);
   }

   virtual CMapLapParaInstance *createCMapLapParaInstance(
         char      *inProbName
         )
   {
      return new CMapLapParaInstanceMpi(inProbName);
   }

   virtual CMapLapParaSolution *createCMapLapParaSolution(
         int inThreadId,         ///< thread Id
         LatticeVector<int> inV, ///< lattice vector
         double inObjValue       ///< objective function value
         );

   virtual CMapLapParaBasis *createCMapLapParaBasis(
         int inThreadId,            ///< thread Id
         double inEnumCost,         ///< enumeration cost
         LatticeBasis<int> inBasis  ///< lattice basis
         );

   virtual CMapLapParaBasis *createCMapLapParaBasis(
         );

   virtual CMapLapParaPackedVector *createCMapLapParaPackedVector(
         int      inThreadId,
         LatticeBasis<int> &inVectors,
         int      inRank,
         int      nThreadsPerRank
         );

   virtual CMapLapParaPackedVector *createCMapLapParaPackedVector(
         int inThreadId,
         std::vector<std::shared_ptr<VectorElement>> &inVectorElements
         );

   virtual CMapLapParaPackedVector *createCMapLapParaPackedVector(
         int inThreadId,
         int inNVectors
         );

   virtual CMapLapParaPackedVector *createCMapLapParaPackedVector(
         );

   ///
   /// create ParaTask DeepBkz object by constructor
   /// @return pointer to ParaTask object
   ///
   virtual UG::ParaTask *createParaTask(
         UG::TaskId inTaskId,                ///< node id
         UG::TaskId inGeneratorTaskId,       ///< generator node id
         double inEstimatedValue,            ///< estimated value
         int inThreadId,                     ///< thread ID
         int inBegin,                        ///< 1st index of block matrix
         int inEnd,                          ///< last index of block matrix
         int inBlockSize,                    ///< the current blocksize
         int inU,                            ///< the unimodular matrix size
         int inSeed,                         ///< the seed of randomize
         std::shared_ptr<LatticeBasis<int>> inBasis   ///< lattice basis
         );

   ///
   /// create ParaTask Enum object by constructor
   /// @return pointer to ParaTask object
   ///
   virtual UG::ParaTask *createParaTask(
         UG::TaskId inTaskId,                ///< node id
         UG::TaskId inGeneratorTaskId,       ///< generator task id
         double   inEstimatedValue,          ///< estimated value
         int      inThreadId,                ///< thread ID
         int      inBegin,                   ///< 1st index of block matrix
         int      inEnd,                     ///< last index of block matrix
         int      inStart,                   ///< the index of current search node in the enumeration tree
         int      inLast,                    ///< search will be terminated when node reaches last-index depth node
         double   inProb,                    ///< the probability of the extreme pruning
         std::shared_ptr<LatticeVector<int>> inCoeffs,   ///< coefficients of basis vectors
         std::shared_ptr<LatticeBasis<int>> inBasis      ///< lattice basis
         );

   ///
   /// create ParaTask Sieve object by constructor
   /// @return pointer to ParaTask object
   ///
   virtual UG::ParaTask *createParaTask(
         UG::TaskId inTaskId,                ///< node id
         UG::TaskId inGeneratorTaskId,       ///< generator task id
         double inEstimatedValue,            ///< estimated value
         int inThreadId,                     ///< thread ID
         int inBegin,                        ///< 1st index of block matrix
         int inEnd,                          ///< last index of block matrix
         std::shared_ptr<LatticeBasis<int>> inBasis   ///< lattice basis
         );


   ///
   /// create ParaSolverState DeepBkz object by default constructor
   /// @return pointer to ParaSolverState object
   ///
   virtual UG::ParaSolverState *createParaSolverState(
         unsigned int  inNotificationId,        ///< id for this notification
         int           inLcId,                  ///< lc id of current ParaTask
         int           inGlobalSubtreeId,       ///< global subtree id of current ParaTask
         double        inDetTime,               ///< deterministic time, -1: should be non-deterministic
         int           inThreadId,              ///< thread id
         int           inDimension,             ///< dimension of SVP instance
         int           inBasisRows,             ///< the number of basis vectors of the lattice basis
         int           inMeanMessageQueueSize,  ///< mean of the message queue size
         int           inMaxMessageQueueSize,   ///< max of the message queue size
         int*          inBasis,                 ///< 1-dimension array of lattice basis row major
         int           inCurrentBlockSize,      ///< current DeepBkz block size
         int           inTour,                  ///< number of DeepBkz loop
         double        inElapsedTime,           ///< elapsed time
         double        inShortestNorm,          ///< the shortest norm found
         double        inApproxFactor,          ///< approximated factor
         double        inHermiteFactor,         ///< hermite factor
         double        inRootHermiteFactor,     ///< (hermite factor)^(1/dim)
         double        inEnumCost,              ///< log of approximted nodes of enumeration tree with incumbent radius
         double        inEnumCostGH,            ///< log of approximted nodes of enumeration tree with GH radius
         double        inSlopeGSA,              ///< slope of GSA line
         double        inTopHalfSlopeGSA,       ///< slope of top-half GSA line
         double        inOrthogonalFactor       ///< orthogonal factor
         );

   ///
   /// create ParaSolverState Enum object by default constructor
   /// @return pointer to ParaSolverState object
   ///
   virtual UG::ParaSolverState *createParaSolverState(
         unsigned int  inNotificationId,        ///< id for this notification
         int           inLcId,                  ///< lc id of current ParaTask
         int           inGlobalSubtreeId,       ///< global subtree id of current ParaTask
         double        inDetTime,               ///< deterministic time, -1: should be non-deterministic
         int           inThreadId,              ///< thread id
         int           inDimension,             ///< dimension of SVP instance
         int           inBasisRows,             ///< the number of basis vectors of the lattice basis
         int           inMeanMessageQueueSize,  ///< mean of the message queue size
         int           inMaxMessageQueueSize,   ///< max of the message queue size
         double        inElapsedTime,           ///< elapsed time
         int*          inCoeffs,                ///< size is dimension, coefficient of a current search node
         int           inDepth,                 ///< current depth of searched node in enumeration tree
         double        inEnumCost,              ///< log of approximted nodes of enumeration tree with incumbent radius
         double        inShortestNorm,          ///< the shortest norm found
         double        inApproxFactor,          ///< approximate factor
         double        inHermiteFactor,         ///< hermite factor
         double        inRootHermiteFactor,     ///< (hermite false)^(1/dim)
         long int      inNumSearchedNodes       ///< number of searched nodes in the enumeration tree
         );

   ///
   /// create ParaSolverState Sieve object by default constructor
   /// @return pointer to ParaSolverState object
   ///
   virtual UG::ParaSolverState *createParaSolverState(
         unsigned int  inNotificationId,        ///< id for this notification
         int           inLcId,                  ///< lc id of current ParaTask
         int           inGlobalSubtreeId,       ///< global subtree id of current ParaTask
         double        inDetTime,               ///< deterministic time, -1: should be non-deterministic
         int           inThreadId,              ///< thread id
         int           inDimension,             ///< dimension of SVP instance
         int           inBasisRows,             ///< the number of basis vectors of the lattice basis
         int           inMeanMessageQueueSize,  ///< mean of the message queue size
         int           inMaxMessageQueueSize,   ///< max of the message queue size
         double        inElapsedTime,           ///< elapsed time
         int           inBlockSize,             ///< block size
         long int      inNLoop,                 ///< number of Sieve algorithm loop
         int           inListSize,              ///< size of List L
         int           inStackSize,             ///< size of Stack S
         int           inMaxListSize,           ///< maximum size of List L up to the point
         int           inNCollisions,           ///< number of collision
         double        inShortestNorm,          ///< the shortest norm found
         double        inApproxFactor,          ///< approximated factor
         double        inHermiteFactor,         ///< hermite factor
         double        inRootHermiteFactor      ///< (hermite factor)^(1/dim)
         );

   ///
   /// create CMapLapParaCalculationState DeepBkz object
   /// @return pointer to CMapLapParaCalculationState object
   ///
   virtual UG::ParaCalculationState *createParaCalculationState(
         int          inTermState,          ///< termination status, 0: normal, -1: interrupted
         int          inThreadId,           ///< thread id
         int          inCurrentBlockSize,   ///< current DeepBkz block size
         int          inTour,               ///< number of DeepBkz loop
         double       inElapsedTime,        ///< elapsed time
         double       inShortestNorm,       ///< the shortest norm found
         double       inApproxFactor,       ///< approximated factor
         double       inHermiteFactor,      ///< hermite factor
         double       inRootHermiteFactor,  ///< (hermite factor)^(1/dim)
         double       inEnumCost,           ///< log of approximted nodes of enumeration tree with incumbent radius
         double       inEnumCostGH,         ///< log of approximted nodes of enumeration tree with GH radius
         double       inSlopeGSA,           ///< slope of GSA line
         double       inTopHalfSlopeGSA,    ///< slope of top-half GSA line
         double       inOrthogonalFactor    ///< orthogonal factor
         );

   ///
   /// create CMapLapParaCalculationState Enum object
   /// @return pointer to CMapLapParaCalculationState object
   ///
   virtual UG::ParaCalculationState *createParaCalculationState(
         int          inTermState,          ///< termination status, 0: normal, -1: interrupted
         int          inThreadId,           ///< thread id
         double       inElapsedTime,        ///< elapsed time
         double       inShortestNorm,       ///< the shortest norm found
         double       inApproxFactor,       ///< approximated factor
         double       inHermiteFactor,      ///< hermite factor
         double       inRootHermiteFactor,  ///< (hermite factor)^(1/dim)
         long int     inNumSearchedNodes    ///< number of searched nodes in the enumeration tree
         );

   ///
   /// create CMapLapParaCalculationState Sieve object
   /// @return pointer to CMapLapParaCalculationState object
   ///
   virtual UG::ParaCalculationState *createParaCalculationState(
         int          inTermState,          ///< termination status, 0: normal, -1: interrupted
         int          inThreadId,           ///< thread id
         double       inElapsedTime,        ///< elapsed time
         int          inBlockSize,          ///< block size
         long int     inNLoop,              ///< number of Sieve algorithm loop
         int          inListSize,           ///< size of List L
         int          inStackSize,          ///< size of Stack S
         int          inMaxListSize,        ///< maximum size of List L up to the point
         int          inNCollisions,        ///< number of collision
         double       inShortestNorm,       ///< the shortest norm found
         double       inApproxFactor,       ///< approximated factor
         double       inHermiteFactor,      ///< hermite factor
         double       inRootHermiteFactor   ///< (hermite factor)^(1/dim)
         );


   ///
   /// create ParaSolverTerminationState object by default constructor
   /// @return pointer to ParaSolverTerminationState object
   ///
   virtual UG::ParaSolverTerminationState *createParaSolverTerminationState(
         );

   ///
   /// create ParaSolverTerminationState object by default constructor
   /// @return pointer to ParaSolverTerminationState object
   ///
   virtual UG::ParaSolverTerminationState *createParaSolverTerminationState(
         int          inThreadId,                     ///< thread id
         int          inInterrupted,                  ///< indicate that this solver is interrupted or not.
                                                      ///< 0: not interrupted,
                                                      ///< 1: interrupted
                                                      ///< 2: checkpoint,
                                                      ///< 3: racing-ramp up
         int          inRank,                         ///< rank of this solver
         int          inNParaTasksReceived,           ///< number of ParaTasks received in this ParaSolver
         int          inNParaTasksSolved,             ///< number of ParaTasks solved ( received ) in this ParaSolver
         double       inRunningTime,                  ///< this solver running time
         double       inIdleTimeToFirstParaTask,      ///< idle time to start solving the first ParaTask
         double       inIdleTimeBetweenParaTasks,     ///< idle time between ParaTasks processing
         double       inIdleTimeAfterLastParaTask,    ///< idle time after the last ParaTask was solved
         double       inIdleTimeToWaitNotificationId, ///< idle time to wait notification Id messages
         double       inIdleTimeToWaitAckCompletion,  ///< idle time to wait ack completion message
         double       inIdleTimeToWaitToken,          ///< idle time to wait token
         double       inIdleTimeToSolverState,        ///< idle time to wait solver state
         double       inIdleTimeToWaitPackedVector,   ///< idle time to wait packed vector
         double       inIdleTimeToWaitSolution,       ///< idle time to wait solution
         double       inIdleTimeToWaitBasis,          ///< idle time to wait basis
         double       inidleTimeToWaitIsend,          ///< idle time to wait Isend
         double       inDetTime,                      ///< deterministic time, -1: should be non-deterministic
         int          inNParaTasksDeepBkzReceived,    ///< number of DeepBkz ParaTasks received in this solver
         int          inNParaTasksEnumReceived,       ///< number of Enum ParaTasks received in this solver
         int          inNParaTasksSieveReceived,      ///< number of Sieve ParaTasks received in this solver
         double       inRunningTimeDeepBkz,           ///< this solver running time of DeepBkz
         double       inRunningTimeEnum,              ///< this solver running time of Enum
         double       inRunningTimeSieve,             ///< this solver running time of Sieve
         int          inNVectorsReceivedDeepBkz,      ///< number of vectors received in DeepBkz
         int          inNVectorsReceivedEnum,         ///< number of vectors received in Enum
         int          inNVectorsReceivedSieve,        ///< number of vectors received in Sieve
         int          inNVectorsSentDeepBkz,          ///< number of vectors sent in DeepBkz
         int          inNVectorsSentEnum,             ///< number of vectors sent in Enum
         int          inNVectorsSentSieve,            ///< number of vectors sent in Sieve
         int          inNBasesSentDeepBkz,            ///< number of vectors sent in DeepBkz
         int          inNSolverStateSent              ///< number of solver states sent
         );

   ///
   /// Crrently not used virtual functions
   ///
   virtual UG::ParaInitialStat* createParaInitialStat(){ throw "** createParaInitialStat is called **"; }

   /// send function for standard ParaData types
   /// @return always 0 (for future extensions)
   ///
   virtual int send(
        void* buffer,
        int count,
        const int datatypeId,
        int dest,
        const int tag
        );


   ///
   /// push CMapLapParaIsendRequest into iSendRequestQueue
   /// @param[in] iSendReq
   ///
   virtual void pushISendRequest(
         std::shared_ptr<CMapLapParaIsendRequest>& iSendReq
         );


   ///
   /// get Tag string for debugging
   /// @return string which shows Tag
   ///
   virtual const char *getTagString(
         int tag                 /// tag to be converted to string
         );


   ///
   /// test isend message
   ///
   virtual int testAllIsends(
         );


   ///
   /// test isend message
   ///
   virtual void waitAllIsends(
         );

};

#define DEF_CMAP_LAP_PARA_COMM( cmapLap_para_comm, comm ) CMapLapParaCommMpi *cmapLap_para_comm = dynamic_cast< CMapLapParaCommMpi* >(comm)

} // namespace ParaCMapLAP


#endif // __CMAP_LAP_PARA_COMM_MPI_H__
