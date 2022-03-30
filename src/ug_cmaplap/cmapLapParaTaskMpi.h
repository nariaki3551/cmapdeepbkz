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

/**@file    cmapLapParaTaskMpi.h
 * @brief   CMapLapParaTask extension for MIP communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_NODE_MPI_H__
#define __CMAP_LAP_PARA_NODE_MPI_H__

#include <memory>
#include <vector>
#include <mpi.h>
#include "ug/paraDef.h"
#include "ug/paraTask.h"
#include "cmapLapParaDef.h"
#include "cmapLapParaTask.h"
namespace UG { class ParaComm; }


namespace ParaCMapLAP
{

///
/// class CMapLapParaNodeMpi
///
class CMapLapParaTaskMpi : public CMapLapParaTask
{

   int iSolverType;     ///< int solver type: 0:DeepBkz, 1:Enum, 2:Sieve, -1: not set
   int nRows;           ///< number of rows of basis
   int nCols;           ///< number of cols of basis

   std::vector<int> vec;
   std::vector<int> box;

   ///
   /// create DeepBkz CMapLapParaSolverStateMpi datatype
   /// @return MPI datatype for CMapLapParaSolverStateMpi
   ///
   virtual MPI_Datatype createDatatypeDeepBkz(
         );

   ///
   /// create Enum CMapLapParaSolverStateMpi datatype
   /// @return MPI datatype for CMapLapParaSolverStateMpi
   ///
   virtual MPI_Datatype createDatatypeEnum(
         );

   ///
   /// create Sieve CMapLapParaSolverStateMpi datatype
   /// @return MPI datatype for CMapLapParaSolverStateMpi
   ///
   virtual MPI_Datatype createDatatypeSieve(
         );

   ///
   /// create CMapLapParaNode datatype
   /// @return pointer to CMapLapParaNodeTh object
   ///
   virtual CMapLapParaTaskMpi *createThreadDatatype(
         UG::ParaComm *comm      ///< communicator used
         )
   {
      return dynamic_cast<CMapLapParaTaskMpi *>(clone(comm));
   }

public :

   ///
   /// default constructor
   ///
   CMapLapParaTaskMpi(
         )
         : iSolverType(static_cast<int>(Undefined))
   {
   }

   ///
   /// constructor of DeepBkz
   ///
   CMapLapParaTaskMpi(
         UG::TaskId intaskId,                ///< task id
         UG::TaskId inGeneratorTaskId,       ///< generator task id
         double   inEstimatedValue,          ///< estimated value
         int      inThreadId,                ///< thread ID
         int      inBegin,                   ///< 1st index of block matrix
         int      inEnd,                     ///< last index of block matrix
         int      inBlockSize,               ///< the current blocksize
         int      inU,                       ///< the unimodular matrix size
         int      inSeed,                    ///< the seed of randomize
         std::shared_ptr<LatticeBasis<int>> inBasis   ///< lattice basis
         )
         : CMapLapParaTask(intaskId, inGeneratorTaskId, inEstimatedValue, inThreadId,
                          inBegin, inEnd, inBlockSize, inU, inSeed, inBasis),
           iSolverType(static_cast<int>(DeepBkz)),
           nRows(inBasis->rows()),
           nCols(inBasis->cols())
   {
   }

   ///
   /// constructor of Enum
   ///
   CMapLapParaTaskMpi(
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
         : CMapLapParaTask(inTaskId, inGeneratorTaskId, inEstimatedValue, inThreadId,
               inBegin, inEnd, inStart, inLast, inProb, inCoeffs, inBasis),
           iSolverType(static_cast<int>(Enum)),
           nRows(inBasis->rows()),
           nCols(inBasis->cols())
   {
   }

   ///
   /// constructor of Sieve
   ///
   CMapLapParaTaskMpi(
         UG::TaskId intaskId,                ///< task id
         UG::TaskId inGeneratorTaskId,       ///< generator task id
         double   inEstimatedValue,          ///< estimated value
         int      inThreadId,                ///< thread ID
         int      inBegin,                   ///< 1st index of block matrix
         int      inEnd,                     ///< last index of block matrix
         std::shared_ptr<LatticeBasis<int>> inBasis ///< lattice basis
         )
         : CMapLapParaTask(intaskId, inGeneratorTaskId, inEstimatedValue, inThreadId,
                          inBegin, inEnd, inBasis),
           iSolverType(static_cast<int>(Sieve)),
           nRows(inBasis->rows()),
           nCols(inBasis->cols())
   {
   }


   ///
   /// destructor
   ///
   virtual ~CMapLapParaTaskMpi(
         )
   {
   }


   ///
   /// clone this CMapLapParaNodeMpi
   /// @return pointer to cloned CMapLapParaNodeMpi object
   ///
   virtual UG::ParaTask *clone(
         UG::ParaComm *comm      ///< communicator used
         )
   {
      switch( solverType )
      {
      case DeepBkz:
      {
         return ( new
            CMapLapParaTaskMpi(taskId, generatorTaskId, estimatedValue, threadId,
                  cmapLapParaTaskDeepBkz.begin,
                  cmapLapParaTaskDeepBkz.end,
                  cmapLapParaTaskDeepBkz.blocksize,
                  cmapLapParaTaskDeepBkz.u,
                  cmapLapParaTaskDeepBkz.seed,
                  std::make_shared<LatticeBasis<int>>(*(cmapLapParaTaskDeepBkz.basis))
                  )
            );
      }
      case Enum:
      {
         return ( new
            CMapLapParaTaskMpi(taskId, generatorTaskId, estimatedValue, threadId,
                  cmapLapParaTaskEnum.begin,
                  cmapLapParaTaskEnum.end,
                  cmapLapParaTaskEnum.start,
                  cmapLapParaTaskEnum.last,
                  cmapLapParaTaskEnum.prob,
                  std::make_shared<LatticeVector<int>>(*(cmapLapParaTaskEnum.coeffs)),
                  std::make_shared<LatticeBasis<int>>(*(cmapLapParaTaskEnum.basis))
                  )
            );
      }
      case Sieve:
      {
         return ( new
            CMapLapParaTaskMpi(taskId, generatorTaskId, estimatedValue, threadId,
                  cmapLapParaTaskSieve.begin,
                  cmapLapParaTaskSieve.end,
                  std::make_shared<LatticeBasis<int>>(*(cmapLapParaTaskSieve.basis))
                  ));
      }
      default:
      {
         THROW_LOGICAL_ERROR2("CMapLapParaTaskMpi::clone: Invalid solver type = ", static_cast<int>(solverType));
      }
      }
      return 0;
   }

   ///
   /// broadcast this object
   /// @return always 0 (for future extensions)
   ///
   virtual int bcast(
         UG::ParaComm *comm,     ///< communicator used
         int root                ///< root rank of broadcast
         );

   ///
   /// send this object
   /// @return always 0 (for future extensions)
   ///
   virtual int send(
         UG::ParaComm *comm,     ///< communicator used
         int destination         ///< destination rank
         );

   ///
   /// receive this object
   /// @return always 0 (for future extensions)
   ///
   virtual int receive(
         UG::ParaComm *comm,     ///< communicator used
         int source              ///< source rank
         );
};

} // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_NODE_MPI_H__
