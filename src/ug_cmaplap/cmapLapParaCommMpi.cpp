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

/**@file    cmapLapParaCommMpi.cpp
 * @brief   CMapLapParaComm extension for MPI communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "cmapLapParaCommMpi.h"
#include "ug/paraIsendRequest.h"
#include "cmapLapParaTagDef.h"
#include "cmapLapParaTaskMpi.h"
#include "cmapLapParaBasisMpi.h"
#include "cmapLapParaSolutionMpi.h"
#include "cmapLapParaIsendRequest.h"
#include "cmapLapParaSolverStateMpi.h"
#include "cmapLapParaPackedVectorMpi.h"
#include "cmapLapParaCalculationStateMpi.h"
#include "cmapLapParaSolverTerminationStateMpi.h"
namespace ParaCMapLAP { class CMapLapParaBasis; }
namespace ParaCMapLAP { class CMapLapParaPackedVector; }
namespace ParaCMapLAP { class CMapLapParaSolution; }
namespace ParaCMapLAP { class VectorElement; }


namespace ParaCMapLAP
{

const char *
CMapLapParaCommMpi::tagStringTable[] = {
  TAG_STR(TagCMapLapPackedVector),
  TAG_STR(TagBasisEnumCost),
  TAG_STR(TagTimeLimitRequest),
  TAG_STR(TagBasisRequest),
  TAG_STR(TagBasis),
  TAG_STR(TagVectorRequest),
  TAG_STR(TagUpdateNotificationInterval),
  TAG_STR(TagParaInstance),
  TAG_STR(TagTask1),
  TAG_STR(TagTask2),
  TAG_STR(TagTask3),
  TAG_STR(TagSolution1),
  TAG_STR(TagSolverState1),
  TAG_STR(TagCompletionOfCalculation1),
  TAG_STR(TagCMapLapPackedVector1),
  TAG_STR(TagBasis1)
};

thread_local CMapLapParaSolverLocalComm *CMapLapParaCommMpi::localComm;

CMapLapParaCommMpi::~CMapLapParaCommMpi(
      )
{
}

///
/// setter of cmapLapParaSolverLocalComm
///
void
CMapLapParaCommMpi::setLocalComm(
      CMapLapParaSolverLocalComm *inLocalComm
   )
{
   localComm = inLocalComm;
}

///
/// lock UG application to synchronize with other threads
/// this is for debug
///
// void
// CMapLapParaCommMpi::lockApp(
//       )
// {
//    if( localComm )
//    {
//       std::cout << "[Rank = " << getRank() << ", Thread = " << localComm->getThreadId() << "] going to lockApp" << std::endl;
//    }
//    applicationLockMutex.lock();
//    if( localComm )
//    {
//       std::cout << "[Rank = " << getRank() << ", Thread = " << localComm->getThreadId() << "] done lockApp" << std::endl;
//    }
// }

///
/// unlock UG application to synchronize with other threads
/// this is for debug
///
// void
// CMapLapParaCommMpi::unlockApp(
//       )
// {
//    applicationLockMutex.unlock();
//    if( localComm )
//    {
//       std::cout << "[Rank = " << getRank() << ", Thread = " << localComm->getThreadId() << "] done unlockApp" << std::endl;
//    }
// }

bool
CMapLapParaCommMpi::tagStringTableIsSetUpCorrectly(
      )
{
   return ( sizeof(tagStringTable)/sizeof(char*) == (N_CMAP_LAP_MPI_TAGS - UG::N_MPI_TAGS) );
}

const char *
CMapLapParaCommMpi::getTagString(
      int tag                 /// tag to be converted to string
      )
{
   assert( tag >= 0 && tag < N_CMAP_LAP_MPI_TAGS );
   if( tag >= 0 && tag < TAG_CMAP_LAP_FIRST )
   {
      return ParaCommMpi::getTagString(tag);
   }
   else
   {
      return tagStringTable[(tag - TAG_CMAP_LAP_FIRST)];
   }
}

/*******************************************************************************
* transfer object factory
*******************************************************************************/

UG::ParaTask *
CMapLapParaCommMpi::createParaTask(
      )
{
   return new CMapLapParaTaskMpi();
}

///
/// create ParaTask DeepBkz
///
UG::ParaTask *
CMapLapParaCommMpi::createParaTask(
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
      )
{
   return new CMapLapParaTaskMpi(inTaskId, inGeneratorTaskId, inEstimatedValue, inThreadId,
         inBegin, inEnd, inBlockSize, inU, inSeed, inBasis);
}

///
/// create ParaTask Enum
///
UG::ParaTask *
CMapLapParaCommMpi::createParaTask(
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
      )
{
   return new CMapLapParaTaskMpi(inTaskId, inGeneratorTaskId, inEstimatedValue, inThreadId,
         inBegin, inEnd, inStart, inLast, inProb, inCoeffs, inBasis);
}

///
/// create ParaTask Sieve
///
UG::ParaTask *
CMapLapParaCommMpi::createParaTask(
      UG::TaskId inTaskId,                ///< task id
      UG::TaskId inGeneratorTaskId,       ///< generator task id
      double inEstimatedValue,            ///< estimated value
      int inThreadId,                     ///< thread ID
      int inBegin,                        ///< 1st index of block matrix
      int inEnd,                          ///< last index of block matrix
      std::shared_ptr<LatticeBasis<int>> inBasis   ///< lattice basis
      )
{
   return new CMapLapParaTaskMpi(inTaskId, inGeneratorTaskId, inEstimatedValue, inThreadId,
         inBegin, inEnd, inBasis);
}



UG::ParaSolverState *
CMapLapParaCommMpi::createParaSolverState(
      )
{
   return new CMapLapParaSolverStateMpi();
}

///
/// create ParaSolverState DeepBkz
///
UG::ParaSolverState *
CMapLapParaCommMpi::createParaSolverState(
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
      )
{
   return new CMapLapParaSolverStateMpi(
         inNotificationId,
         inLcId,
         inGlobalSubtreeId,
         inDetTime,
         inThreadId,
         inDimension,
         inBasisRows,
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
/// create ParaSolverState Enum
///
UG::ParaSolverState *
CMapLapParaCommMpi::createParaSolverState(
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
      double        inEnumCost,              ///< approximated time of full Enum
      double        inShortestNorm,          ///< the shortest norm found
      double        inApproxFactor,          ///< approximate factor
      double        inHermiteFactor,         ///< hermite factor
      double        inRootHermiteFactor,     ///< (hermite false)^(1/dim)
      long int      inNumSearchedNodes       ///< number of searched nodes in the enumeration tree
      )
{
   return new CMapLapParaSolverStateMpi(
         inNotificationId,
         inLcId,
         inGlobalSubtreeId,
         inDetTime,
         inThreadId,
         inDimension,
         inBasisRows,
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
/// create ParaSolverState Sieve
///
UG::ParaSolverState *
CMapLapParaCommMpi::createParaSolverState(
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
      )
{
   return new CMapLapParaSolverStateMpi(
         inNotificationId,
         inLcId,
         inGlobalSubtreeId,
         inDetTime,
         inThreadId,
         inDimension,
         inBasisRows,
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



UG::ParaCalculationState *
CMapLapParaCommMpi::createParaCalculationState(
      )
{
   return new CMapLapParaCalculationStateMpi();
}

///
/// create CalculationState DeepBkz
///
UG::ParaCalculationState *
CMapLapParaCommMpi::createParaCalculationState(
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
   return new CMapLapParaCalculationStateMpi(
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
/// create CalculationState Enum
///
UG::ParaCalculationState *
CMapLapParaCommMpi::createParaCalculationState(
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
   return new CMapLapParaCalculationStateMpi(
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
/// create CalculationState Sieve
///
UG::ParaCalculationState *
CMapLapParaCommMpi::createParaCalculationState(
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
      )
{
   return new CMapLapParaCalculationStateMpi(
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



///
/// create ParaSolution
///
CMapLapParaSolution*
CMapLapParaCommMpi::createCMapLapParaSolution(
      int inThreadId,               ///< thread Id
      LatticeVector<int> inV,       ///< lattice vector
      double inObjValue             ///< objective function value
      )
{
   return ( new CMapLapParaSolutionMpi(inThreadId, inV, inObjValue) );
}


///
/// @brief create ParaBasis
///
CMapLapParaBasis *
CMapLapParaCommMpi::createCMapLapParaBasis(
      int inThreadId,               ///< thread Id
      double inEnumCost,            ///< enumeration cost
      LatticeBasis<int> inBasis     ///< lattice basis
      )
{
   return new CMapLapParaBasisMpi(inThreadId, inEnumCost, inBasis);
}

///
/// @brief create ParaBasis
///
CMapLapParaBasis *
CMapLapParaCommMpi::createCMapLapParaBasis(
      )
{
   return new CMapLapParaBasisMpi();
}


CMapLapParaPackedVector *
CMapLapParaCommMpi::createCMapLapParaPackedVector(
      int      inThreadId,
      LatticeBasis<int> &inVectors,
      int      inRank,
      int      nThreadsPerRank
      )
{
   return new CMapLapParaPackedVectorMpi(inThreadId, inVectors, inRank, nThreadsPerRank);
}

CMapLapParaPackedVector *
CMapLapParaCommMpi::createCMapLapParaPackedVector(
      int inThreadId,
      std::vector<std::shared_ptr<VectorElement>> &inVectorElements
      )
{
   return new CMapLapParaPackedVectorMpi(inThreadId, inVectorElements);
}

CMapLapParaPackedVector *
CMapLapParaCommMpi::createCMapLapParaPackedVector(
      int inThreadId,
      int inNVectors
      )
{
   return new CMapLapParaPackedVectorMpi(inThreadId, inNVectors);
}

CMapLapParaPackedVector *
CMapLapParaCommMpi::createCMapLapParaPackedVector()
{
  return new CMapLapParaPackedVectorMpi();
}



UG::ParaSolverTerminationState *
CMapLapParaCommMpi::createParaSolverTerminationState(
      )
{
   return new CMapLapParaSolverTerminationStateMpi();
}

UG::ParaSolverTerminationState *
CMapLapParaCommMpi::createParaSolverTerminationState(
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
      )
{
      return new CMapLapParaSolverTerminationStateMpi(
            inThreadId,
            inInterrupted,
            inRank,
            inNParaTasksReceived,
            inNParaTasksSolved,
            inRunningTime,
            inIdleTimeToFirstParaTask,
            inIdleTimeBetweenParaTasks,
            inIdleTimeAfterLastParaTask,
            inIdleTimeToWaitNotificationId,
            inIdleTimeToWaitAckCompletion,
            inIdleTimeToWaitToken,
            inIdleTimeToSolverState,
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
            );
}


int
CMapLapParaCommMpi::send(
   void* buffer,
   int count,
   const int datatypeId,
   int dest,
   const int tag
   )
{
   UG::ParaIsendRequest *request = nullptr;
   switch ( datatypeId )
   {
   case UG::ParaCHAR:
   {
      char *pChar = nullptr;
      if( buffer )
      {
         pChar = new char[count];
         memcpy(pChar, buffer, sizeof(char)*count);
      }
      request = new UG::ParaIsendRequest(UG::ParaIsendRequest::ParaCHAR, new MPI_Request, pChar);
      break;
   }
   case UG::ParaSHORT:
   {
      short *pShort = nullptr;
      if( buffer )
      {
         pShort = new short[count];
         memcpy(pShort, buffer, sizeof(short)*count);
      }
      request = new UG::ParaIsendRequest(UG::ParaIsendRequest::ParaSHORT, new MPI_Request, pShort);
      break;
   }
   case UG::ParaINT:
   {
      int *pInt = nullptr;
      if( buffer )
      {
         pInt = new int[count];
         memcpy(pInt, buffer, sizeof(int)*count);
      }
      request = new UG::ParaIsendRequest(UG::ParaIsendRequest::ParaINT, new MPI_Request, pInt);
      break;
   }
   case UG::ParaLONG:
   {
      long *pLong = nullptr;
      if( buffer )
      {
         pLong = new long[count];
         memcpy(pLong, buffer, sizeof(long)*count);
      }
      request = new UG::ParaIsendRequest(UG::ParaIsendRequest::ParaLONG, new MPI_Request, pLong);
      break;
   }
   case UG::ParaLONG_LONG:
   {
      long long *pLongLong = nullptr;
      if( buffer )
      {
         pLongLong = new long long[count];
         memcpy(pLongLong, buffer, sizeof(long long)*count);
      }
      request = new UG::ParaIsendRequest(UG::ParaIsendRequest::ParaLONG_LONG, new MPI_Request, pLongLong);
      break;
   }
   case UG::ParaSIGNED_CHAR:
   {
      signed char *pSignedChar = nullptr;
      if( buffer )
      {
         pSignedChar = new signed char[count];
         memcpy(pSignedChar, buffer, sizeof(signed char)*count);
      }
      request = new UG::ParaIsendRequest(UG::ParaIsendRequest::ParaSIGNED_CHAR, new MPI_Request, pSignedChar);
      break;
   }
   case UG::ParaUNSIGNED_CHAR:
   {
      unsigned char *pUnsignedChar = nullptr;
      if( buffer )
      {
         pUnsignedChar = new unsigned char[count];
         memcpy(pUnsignedChar, buffer, sizeof(unsigned char)*count);
      }
      request = new UG::ParaIsendRequest(UG::ParaIsendRequest::ParaUNSIGNED_CHAR, new MPI_Request, pUnsignedChar);
      break;
   }
   case UG::ParaUNSIGNED_SHORT:
   {
      unsigned short *pUnsignedShort = nullptr;
      if( buffer )
      {
         pUnsignedShort = new unsigned short[count];
         memcpy(pUnsignedShort, buffer, sizeof(unsigned short)*count);
      }
      request = new UG::ParaIsendRequest(UG::ParaIsendRequest::ParaUNSIGNED_SHORT, new MPI_Request, pUnsignedShort);
      break;
   }
   case UG::ParaUNSIGNED:
   {
      unsigned *pUnsigned = nullptr;
      if( buffer )
      {
         pUnsigned = new unsigned[count];
         memcpy(pUnsigned, buffer, sizeof(unsigned)*count);
      }
      request = new UG::ParaIsendRequest(UG::ParaIsendRequest::ParaUNSIGNED, new MPI_Request, pUnsigned);
      break;
   }
   case UG::ParaUNSIGNED_LONG:
   {
      unsigned long *pUnsignedLong = nullptr;
      if( buffer )
      {
         pUnsignedLong = new unsigned long[count];
         memcpy(pUnsignedLong, buffer, sizeof(unsigned long)*count);
      }
      request = new UG::ParaIsendRequest(UG::ParaIsendRequest::ParaUNSIGNED_LONG, new MPI_Request, pUnsignedLong);
      break;
   }
   case UG::ParaUNSIGNED_LONG_LONG:
   {
      unsigned long long *pUnsignedLongLong = nullptr;
      if( buffer )
      {
         pUnsignedLongLong = new unsigned long long[count];
         memcpy(pUnsignedLongLong, buffer, sizeof(unsigned long long)*count);
      }
      request = new UG::ParaIsendRequest(UG::ParaIsendRequest::ParaUNSIGNED_LONG_LONG, new MPI_Request, pUnsignedLongLong);
      break;
   }
   case UG::ParaFLOAT:
   {
      float *pFloat = nullptr;
      if( buffer )
      {
         pFloat = new float[count];
         memcpy(pFloat, buffer, sizeof(float)*count);
      }
      request = new UG::ParaIsendRequest(UG::ParaIsendRequest::ParaFLOAT, new MPI_Request, pFloat);
      break;
   }
   case UG::ParaDOUBLE:
   {
      double *pDouble = nullptr;
      if( buffer )
      {
         pDouble = new double[count];
         memcpy(pDouble, buffer, sizeof(double)*count);
      }
      request = new UG::ParaIsendRequest(UG::ParaIsendRequest::ParaDOUBLE, new MPI_Request, pDouble);
      break;
   }
   case UG::ParaLONG_DOUBLE:
   {
      long double *pLongDouble = nullptr;
      if( buffer )
      {
         pLongDouble = new long double[count];
         memcpy(pLongDouble, buffer, sizeof(long double)*count);
      }
      request = new UG::ParaIsendRequest(UG::ParaIsendRequest::ParaLONG_DOUBLE, new MPI_Request, pLongDouble);
      break;
   }
   case UG::ParaBOOL:
   {
      bool *pBool = nullptr;
      if( buffer )
      {
         pBool = new bool[count];
         memcpy(pBool, buffer, sizeof(bool)*count);
      }
      request = new UG::ParaIsendRequest(UG::ParaIsendRequest::ParaBOOL, new MPI_Request, pBool);
      break;
   }
   case UG::ParaBYTE:
   {
      assert(!buffer);
      char *pByte = nullptr;
      if( buffer )
      {
         pByte = new char[count];
         memcpy(pByte, buffer, sizeof(char)*count);
      }
      request = new UG::ParaIsendRequest(UG::ParaIsendRequest::ParaBYTE, new MPI_Request, pByte);
      break;
   }
   default:
   {
      std::cerr << "Invalid dataType = " << datatypeId << std::endl;
      abort();
   }
   }

   MPI_CALL(
      MPI_Isend( request->buffer(), count, datatypes[datatypeId], dest, tag, myComm, request->req )
   );
   TAG_TRACE (SendIsend, To, dest, tag);
   iSendRequestDeque.push_back(request);
   return 0;
}


///
/// push CMapLapParaIsendRequest into iSendRequestQueue
/// @param[in] iSendReq
///
void
CMapLapParaCommMpi::pushISendRequest(
      std::shared_ptr<CMapLapParaIsendRequest>& iSendReq
      )
{
   cmapLapISendRequestQueue.push_back(iSendReq);
}

///
/// test isend message
///
int
CMapLapParaCommMpi::testAllIsends(
      )
{
   UG::ParaCommMpi::testAllIsends();
   if( !cmapLapISendRequestQueue.empty() )
   {
      auto it = cmapLapISendRequestQueue.begin();
      while( it != cmapLapISendRequestQueue.end() )
      {
         if( (*it)->test() ) // it is std::shared_ptr<CMapLapParaIsendRequest> object
         {
            it = cmapLapISendRequestQueue.erase(it);
         }
         else
         {
            it++;
         }
      }
   }
   return cmapLapISendRequestQueue.size() + iSendRequestDeque.size();
}


///
/// test isend message
///
void
CMapLapParaCommMpi::waitAllIsends(
      )
{
   if( !cmapLapISendRequestQueue.empty() )
   {
      for( auto const &it : cmapLapISendRequestQueue )
      {
         it->wait(); // it is std::shared_ptr<CMapLapParaIsendRequest> object
      }
      cmapLapISendRequestQueue.clear();
   }
}


} // namespace ParaCMapLAP
