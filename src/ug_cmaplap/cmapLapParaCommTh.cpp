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

/**@file    cmapLapParaCommTh.cpp
 * @brief   CMapLapParaComm extension for threads communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <assert.h>
#include "cmapLapParaInstanceTh.h"
#include "cmapLapParaPackedVectorTh.h"
#include "cmapLapParaParamSet.h"
#include "cmapLapParaSolutionTh.h"
#include "cmapLapParaBasisTh.h"
#include "cmapLapParaCommTh.h"
#include "cmapLapParaDef.h"
#include "cmapLapParaTagDef.h"
#include "cmapLapParaTaskTh.h"
namespace ParaCMapLAP { class CMapLapParaBasis; }
namespace ParaCMapLAP { class CMapLapParaInstance; }
namespace ParaCMapLAP { class CMapLapParaSolution; }
namespace ParaCMapLAP { class CMapLapParaSolverLocalComm; }
namespace ParaCMapLAP { class CMapLapParaPackedVector; }
namespace ParaCMapLAP { class VectorElement; }


namespace ParaCMapLAP
{

const char *
CMapLapParaCommTh::tagStringTable[] = {
   TAG_STR(TagCMapLapPackedVector),
   TAG_STR(TagBasisEnumCost),
   TAG_STR(TagTimeLimitRequest),
   TAG_STR(TagBasisRequest),
   TAG_STR(TagBasis),
   TAG_STR(TagVectorRequest),
   TAG_STR(TagUpdateNotificationInterval)
};

thread_local CMapLapParaSolverLocalComm *CMapLapParaCommTh::localComm;

///
/// setter of cmapLapParaSolverLocalComm
///
void
CMapLapParaCommTh::setLocalComm(
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
// CMapLapParaCommTh::lockApp(
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
// CMapLapParaCommTh::unlockApp(
//       )
// {
//    applicationLockMutex.unlock();
//    if( localComm )
//    {
//       std::cout << "[Rank = " << getRank() << ", Thread = " << localComm->getThreadId() << "] done unlockApp" << std::endl;
//    }
// }

bool
CMapLapParaCommTh::tagStringTableIsSetUpCorrectly(
      )
{
   // if( !UG::ParaCommCPP11::tagStringTableIsSetUpCorrectly() ) return false;
   // std::cout << "size = " << sizeof(tagStringTable)/sizeof(char*)
   //       << ", (N_CMAP_LAP_TH_TAGS - UG::N_TH_TAGS) = " <<  (N_CMAP_LAP_TH_TAGS - UG::N_TH_TAGS) << std::endl;
   return ( sizeof(tagStringTable)/sizeof(char*) == (N_CMAP_LAP_TH_TAGS - UG::N_TH_TAGS) );
}

const char *
CMapLapParaCommTh::getTagString(
      int tag                 ///< tag to be converted to string
      )
{
   assert( tag >= 0 && tag < N_CMAP_LAP_TH_TAGS );
   if( tag >= 0 && tag < TAG_CMAP_LAP_FIRST )
   {
      return ParaCommCPP11::getTagString(tag);
   }
   else
   {
      return tagStringTable[(tag - TAG_CMAP_LAP_FIRST)];
   }
}


//thread_local int
//CMapLapParaCommTh::localThreadId = -1;     /*< local thread id

/*******************************************************************************
* transfer object factory
*******************************************************************************/
UG::ParaInstance*
CMapLapParaCommTh::createParaInstance(
      )
{
   return new CMapLapParaInstanceTh();
}

UG::ParaSolution*
CMapLapParaCommTh::createParaSolution(
      )
{
   return new CMapLapParaSolutionTh();
}

CMapLapParaInstance*
CMapLapParaCommTh::createCMapLapParaInstance(
      char      *inProbFileName,
      char      *inProbName
      )
{
   return new CMapLapParaInstanceTh(inProbFileName, inProbName);
}

CMapLapParaInstance*
CMapLapParaCommTh::createCMapLapParaInstance(
      char      *inProbName
      )
{
   return new CMapLapParaInstanceTh(inProbName);
}

CMapLapParaSolution*
CMapLapParaCommTh::createCMapLapParaSolution(
      int inThreadId,               ///< thread Id
      LatticeVector<int> inV,       ///< updated vector
      double inObjValue             ///< objective function value
      )
{
   return new CMapLapParaSolutionTh(inThreadId, inV, inObjValue);
}

///
/// create ParaSolverTerminationState object by default constructor
/// @return pointer to ParaSolverTerminationState object
///
UG::ParaSolverTerminationState *
CMapLapParaCommTh::createParaSolverTerminationState(
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
      return new CMapLapParaSolverTerminationStateTh(
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


UG::ParaTask *
CMapLapParaCommTh::createParaTask(
      )
{
   return new CMapLapParaTaskTh();
}


///
/// create ParaTask DeepBkz object by constructor
/// @return pointer to ParaTask object
///
UG::ParaTask *
CMapLapParaCommTh::createParaTask(
      UG::TaskId inTaskId,                ///< task id
      UG::TaskId inGeneratorTaskId,       ///< generator task id
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
   return new CMapLapParaTaskTh(inTaskId, inGeneratorTaskId, inEstimatedValue, inThreadId,
         inBegin, inEnd, inBlockSize, inU, inSeed, inBasis);
}

///
/// create ParaTask Enum object by constructor
/// @return pointer to ParaTask object
///
UG::ParaTask *
CMapLapParaCommTh::createParaTask(
      UG::TaskId inTaskId,                ///< task id
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
   return new CMapLapParaTaskTh(inTaskId, inGeneratorTaskId, inEstimatedValue, inThreadId,
         inBegin, inEnd, inStart, inLast, inProb, inCoeffs, inBasis);
}

///
/// create ParaTask Sieve object by constructor
/// @return pointer to ParaTask object
///
UG::ParaTask *
CMapLapParaCommTh::createParaTask(
      UG::TaskId inTaskId,                ///< task id
      UG::TaskId inGeneratorTaskId,       ///< generator task id
      double inEstimatedValue,            ///< estimated value
      int inThreadId,                     ///< thread ID
      int inBegin,                        ///< 1st index of block matrix
      int inEnd,                          ///< last index of block matrix
      std::shared_ptr<LatticeBasis<int>> inBasis   ///< lattice basis
      )
{
   return new CMapLapParaTaskTh(inTaskId, inGeneratorTaskId, inEstimatedValue, inThreadId,
         inBegin, inEnd, inBasis);
}


///
/// @brief create ParaBasis object
///
CMapLapParaBasis *
CMapLapParaCommTh::createCMapLapParaBasis(
      int inThreadId,            ///< thread Id
      double inEnumCost,         ///< enumeration cost
      LatticeBasis<int> &inBasis ///< lattice basis
      )
{
   return new CMapLapParaBasisTh(inThreadId, inEnumCost, inBasis);
}

///
/// @brief create ParaBasis object
///
CMapLapParaBasis *
CMapLapParaCommTh::createCMapLapParaBasis(
      )
{
   return new CMapLapParaBasisTh();
}


///
/// @brief create ParaPackedVector object
///
CMapLapParaPackedVector *
CMapLapParaCommTh::createCMapLapParaPackedVector(
      int      inThreadId,
      LatticeBasis<int> &inVectors,
      int      inRank,
      int      nThreadsPerRank
      )
{
   return new CMapLapParaPackedVectorTh(inThreadId, inVectors, inRank, nThreadsPerRank);
}

///
/// @brief create ParaPackedVector object
///
CMapLapParaPackedVector *
CMapLapParaCommTh::createCMapLapParaPackedVector(
      int inThreadId,
      std::vector<std::shared_ptr<VectorElement>> &inVectorElements
      )
{
   return new CMapLapParaPackedVectorTh(inThreadId, inVectorElements);
}

///
/// @brief create ParaPackedVector object
///
CMapLapParaPackedVector *
CMapLapParaCommTh::createCMapLapParaPackedVector(
      int inThreadId,
      int inNVectors
      )
{
   return new CMapLapParaPackedVectorTh(inThreadId, inNVectors);
}

///
/// @brief create ParaPackedVector object
///
CMapLapParaPackedVector *
CMapLapParaCommTh::createCMapLapParaPackedVector(
      )
{
   return new CMapLapParaPackedVectorTh();
}

} // namespace ParaCMapLAP
