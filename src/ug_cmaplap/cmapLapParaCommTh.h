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

/**@file    cmapLapParaCommTh.h
 * @brief   CMapLapParaComm extension for threads communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_COMM_TH_H__
#define __CMAP_LAP_PARA_COMM_TH_H__


#include <vector>
#include "cmapLapParaCalculationStateTh.h"
#include "cmapLapParaSolverStateTh.h"
#include "cmapLapParaSolverTerminationStateTh.h"
#include "ug/paraCalculationState.h"
#include "ug/paraCommCPP11.h"
#include "ug/paraInstance.h"
#include "ug/paraParamSet.h"
#include "ug/paraRacingRampUpParamSet.h"
#include "ug/paraSolution.h"
#include "ug/paraSolverState.h"
#include "ug/paraSolverTerminationState.h"
#include "ug/paraTask.h"
namespace ParaCMapLAP { class CMapLapParaBasis; }
namespace ParaCMapLAP { class CMapLapParaInstance; }
namespace ParaCMapLAP { class CMapLapParaSolution; }
namespace ParaCMapLAP { class CMapLapParaPackedVector; }
namespace ParaCMapLAP { class VectorElement; }


namespace ParaCMapLAP
{

class CMapLapParaSolverLocalComm;

///
/// cmapLap user defined transfer data types
///
#ifndef __CMAP_LAP_SOLVER_LOCAL_COMM_H_
static const int CMAP_LAP_USER_TYPE_FIRST       = UG::UG_USER_TYPE_LAST + 1;
static const int ParaBasisType                  = CMAP_LAP_USER_TYPE_FIRST + 0;
static const int CMapLapParaPackedVectorType    = CMAP_LAP_USER_TYPE_FIRST + 1;
static const int ParaCheckpointElementType      = CMAP_LAP_USER_TYPE_FIRST + 2;
static const int CMAP_LAP_USER_TYPE_LAST        = CMAP_LAP_USER_TYPE_FIRST + 2;
#endif

class CMapLapParaCommTh : public UG::ParaCommCPP11
{

   static const char          *tagStringTable[];               ///< tag name string table
   static thread_local CMapLapParaSolverLocalComm *localComm;   ///< for get rank and threadId when lock / unlock App
   std::mutex                 *applicationLockMutex;           ///< mutex for applications


   ///
   /// check if tag string table (for debugging) set up correctly
   /// @return true if tag string table is set up correctly, false otherwise
   ///
   virtual bool tagStringTableIsSetUpCorrectly(
         );

public:

   CMapLapParaCommTh(
         )
         : UG::ParaCommCPP11(), applicationLockMutex(0)
   {
   }

   virtual ~CMapLapParaCommTh(
         )
   {
      if( applicationLockMutex ) delete[] applicationLockMutex;
   }

   ///
   /// initializer of this communicator
   ///
   virtual void init(
         int argc,             ///< the number of arguments
         char **argv           ///< pointers to the arguments
         )
   {
      ParaCommCPP11::init(argc, argv);
      applicationLockMutex = new std::mutex[comSize];
   }

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
      applicationLockMutex[localRank].lock();
   }

   ///
   /// unlock UG application to synchronize with other threads
   ///
   // virtual void unlockApp(
   //       );
   virtual void unlockApp(
         )
   {
      applicationLockMutex[localRank].unlock();
   }


   virtual void solverInit(
         UG::ParaParamSet *paraParamSet   ///< UG parameter set
         )
   {
   }

   ///
   /// initializer for a specific slave thread Solver
   ///
   virtual void solverInit(
         int rank,                        ///< rank of the Solver
         UG::ParaParamSet *paraParamSet   ///< UG parameter set
         )
   {
      localRank = rank;
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
         )
   {
      return new CMapLapParaCalculationStateTh();
   }

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
         )
   {
      return new CMapLapParaSolverStateTh();
   }

   ///
   /// create ParaSolverTerminationState object by default constructor
   /// @return pointer to ParaSolverTerminationState object
   ///
   virtual UG::ParaSolverTerminationState *createParaSolverTerminationState(
         )
   {
      return new CMapLapParaSolverTerminationStateTh();
   }

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
         double       inIdleTimeToWaitIsend,          ///< idle time to wait Isend
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


   virtual UG::ParaInstance *createParaInstance(
         );


   virtual UG::ParaSolution *createParaSolution(
         );

   ///
   /// create ParaTask DeepBkz object by constructor
   /// @return pointer to ParaTask object
   ///
   virtual UG::ParaTask *createParaTask(
         UG::TaskId inTaskId,                ///< task id
         UG::TaskId inGeneratorTaskId,       ///< generator task id
         double   inEstimatedValue,          ///< estimated value
         int      inThreadId,                ///< thread ID
         int      inBegin,                   ///< 1st index of block matrix
         int      inEnd,                     ///< last index of block matrix
         int      inBlockSize,               ///< the current blocksize
         int      inU,                       ///< the unimodular matrix size
         int      inSeed,                    ///< the seed of randomize
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
         UG::TaskId inTaskId,                ///< task id
         UG::TaskId inGeneratorTaskId,       ///< generator task id
         double   inEstimatedValue,          ///< estimated value
         int      inThreadId,                ///< thread ID
         int      inBegin,                   ///< 1st index of block matrix
         int      inEnd,                     ///< last index of block matrix
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
         int           inNumBasisVectors,       ///< the number of basis vectors of the lattice basis
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
         )
   {
      return new CMapLapParaSolverStateTh(
            inNotificationId,
            inLcId,
            inGlobalSubtreeId,
            inDetTime,
            inThreadId,
            inDimension,
            inNumBasisVectors,
            inMeanMessageQueueSize,
            inMaxMessageQueueSize,
            inBasis,
            inCurrentBlockSize,
            inTour,
            inElapsedTime,
            inShortestNorm,
            inApproxFactor,
            inHermiteFactor,
            inRootHermiteFactor,
            inEnumCost,
            inEnumCostGH,
            inSlopeGSA,
            inTopHalfSlopeGSA,
            inOrthogonalFactor
            );
   }

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
         int           inNumBasisVectors,       ///< the number of basis vectors of the lattice basis
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
         )
   {
      return new CMapLapParaSolverStateTh(
            inNotificationId,
            inLcId,
            inGlobalSubtreeId,
            inDetTime,
            inThreadId,
            inDimension,
            inNumBasisVectors,
            inMeanMessageQueueSize,
            inMaxMessageQueueSize,
            inElapsedTime,
            inCoeffs,
            inDepth,
            inEnumCost,
            inShortestNorm,
            inApproxFactor,
            inHermiteFactor,
            inRootHermiteFactor,
            inNumSearchedNodes
            );
   }

   ///
   /// create ParaSolverState Sieve object by default constructor
   /// @return pointer to ParaSolverState object
   ///
   virtual UG::ParaSolverState *createParaSolverState(
         unsigned int  inNotificationId,        ///< id for this notification
         int           inLcId,                  ///< lc id of current ParaNode
         int           inGlobalSubtreeId,       ///< global subtree id of current ParaNode
         double        inDetTime,               ///< deterministic time, -1: should be non-deterministic
         int           inThreadId,              ///< thread id
         int           inDimension,             ///< dimension of SVP instance
         int           inNumBasisVectors,       ///< the number of basis vectors of the lattice basis
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
         )
   {
      return new CMapLapParaSolverStateTh(
            inNotificationId,
            inLcId,
            inGlobalSubtreeId,
            inDetTime,
            inThreadId,
            inDimension,
            inNumBasisVectors,
            inMeanMessageQueueSize,
            inMaxMessageQueueSize,
            inElapsedTime,
            inBlockSize,
            inNLoop,
            inListSize,
            inStackSize,
            inMaxListSize,
            inNCollisions,
            inShortestNorm,
            inApproxFactor,
            inHermiteFactor,
            inRootHermiteFactor
            );
   }


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
         )
   {
      return new CMapLapParaCalculationStateTh(
            inTermState,
            inThreadId,
            inCurrentBlockSize,
            inTour,
            inElapsedTime,
            inShortestNorm,
            inApproxFactor,
            inHermiteFactor,
            inRootHermiteFactor,
            inEnumCost,
            inEnumCostGH,
            inSlopeGSA,
            inTopHalfSlopeGSA,
            inOrthogonalFactor
            );
   }

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
         )
   {
      return new CMapLapParaCalculationStateTh(
            inTermState,
            inThreadId,
            inElapsedTime,
            inShortestNorm,
            inApproxFactor,
            inHermiteFactor,
            inRootHermiteFactor,
            inNumSearchedNodes
            );
   }

   ///
   /// create CMapLapParaCalculationState Sieve object
   /// @return pointer to CMapLapParaCalculationState object
   ///
   virtual UG::ParaCalculationState *createParaCalculationState(
         int          inTermState,          ///< termination status
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
         )
   {
      return new CMapLapParaCalculationStateTh(
            inTermState,
            inThreadId,
            inElapsedTime,
            inBlockSize,
            inNLoop,
            inListSize,
            inStackSize,
            inMaxListSize,
            inNCollisions,
            inShortestNorm,
            inApproxFactor,
            inHermiteFactor,
            inRootHermiteFactor
            );
   }


   virtual CMapLapParaInstance *createCMapLapParaInstance(
         char      *inProbFileName,
         char      *probName
         );

   virtual CMapLapParaInstance *createCMapLapParaInstance(
         char      *probName
         );

   virtual CMapLapParaSolution *createCMapLapParaSolution(
         int      inThreadId,          ///< thread Id
         LatticeVector<int> inV,       ///< updated vector
         double inObjValue             ///< objective function value
         );

   ///
   /// @brief create ParaBasis object
   ///
   virtual CMapLapParaBasis *createCMapLapParaBasis(
         int inThreadId,            ///< thread Id
         double inEnumCost,         ///< enumeration cost
         LatticeBasis<int> &inBasis ///< lattice basis
         );

   ///
   /// @brief create ParaBasis object
   ///
   virtual CMapLapParaBasis *createCMapLapParaBasis(
         );


   ///
   /// @brief create ParaPackedVector object
   ///
   virtual CMapLapParaPackedVector *createCMapLapParaPackedVector(
         int      inThreadId,
         LatticeBasis<int> &inVectors,
         int      inRank,
         int      nThreadsPerRank
         );

   ///
   /// @brief create ParaPackedVector object
   ///
   virtual CMapLapParaPackedVector *createCMapLapParaPackedVector(
         int inThreadId,
         std::vector<std::shared_ptr<VectorElement>> &inVectorElements
         );

   ///
   /// @brief create ParaPackedVector object
   ///
   virtual CMapLapParaPackedVector *createCMapLapParaPackedVector(
         int inThreadId,
         int inNVectors
         );

   ///
   /// @brief create ParaPackedVector object
   ///
   virtual CMapLapParaPackedVector *createCMapLapParaPackedVector(
         );

   ///
   /// get Tag string for debugging
   /// @return string which shows Tag
   ///
   virtual const char *getTagString(
         int tag                 /// tag to be converted to string
         );

};

#define DEF_CMAP_LAP_PARA_COMM( cmapLap_para_comm, comm ) CMapLapParaCommTh *cmapLap_para_comm = dynamic_cast< CMapLapParaCommTh* >(comm)

} // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_COMM_TH_H__
