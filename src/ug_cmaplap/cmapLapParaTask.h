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

/**@file    cmapLapParaTask.h
 * @brief   Base class for CMapLapParaTask.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_TASK_H__
#define __CMAP_LAP_PARA_TASK_H__

#include <cassert>
#include <iostream>
#include <memory>
#include "ug/paraDef.h"
#include "ug/paraTask.h"
#include "cmapLapParaDef.h"
#include "cmapLapParaLog.h"
#ifdef UG_WITH_ZLIB
#include "ug/gzstream.h"
#endif
namespace UG { class ParaComm; }
namespace ParaCMapLAP { class CMapLapParaTaskLC; }
namespace ParaCMapLAP { class CMapLapParaTaskMpi; }
namespace ParaCMapLAP { class CMapLapParaTaskTh; }


namespace ParaCMapLAP
{

struct CMapLapParaTaskDeepBkz {
   int         begin;            ///< 1st index of block matrix
   int           end;            ///< last index of block matrix
   int     blocksize;            ///< the current blocksize
   int             u;            ///< the size of unimodular matrix
   int          seed;            ///< the seed of randomize
   std::shared_ptr<LatticeBasis<int>> basis; ///< lattice basis
};

struct CMapLapParaTaskEnum {
   int         begin;            ///< 1st index of block matrix
   int           end;            ///< last index of block matrix
   int         start;            ///< index of current search node in the enumeration tree
   int          last;            ///< search will be terminated when node reaches last-index depth node
   double       prob;            ///< the probability of the extreme pruning
   std::shared_ptr<LatticeVector<int>> coeffs;  ///< the coefficients of basis vectors
   std::shared_ptr<LatticeBasis<int>> basis;    ///< lattice basis
};

struct CMapLapParaTaskSieve {
   int         begin;            ///< 1st index of block matrix
   int           end;            ///< last index of block matrix
   std::shared_ptr<LatticeBasis<int>> basis; ///< lattice basis
};


///
/// class CMapLapParaTask
///
class CMapLapParaTask : public UG::ParaTask, public std::enable_shared_from_this<CMapLapParaTask>
{

protected:

   int threadId;  ///< threadId = -1 when solverType is Sieve
   SolverType solverType;
   CMapLapParaTaskDeepBkz cmapLapParaTaskDeepBkz;
   CMapLapParaTaskEnum cmapLapParaTaskEnum;
   CMapLapParaTaskSieve cmapLapParaTaskSieve;

public:

   friend CMapLapParaTaskTh;
   friend CMapLapParaTaskMpi;
   friend CMapLapParaTaskLC;
   ///
   /// default constructor
   ///
   CMapLapParaTask(
         )
         :
         threadId(-1),
         solverType(Undefined)
   {
   }

   ///
   ///  constructor of DeepBkz
   ///
   CMapLapParaTask(
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
         )
         : UG::ParaTask(inTaskId, inGeneratorTaskId, inEstimatedValue, 0),
           threadId(inThreadId),
           solverType(DeepBkz)
   {
      cmapLapParaTaskDeepBkz.begin     = inBegin;
      cmapLapParaTaskDeepBkz.end       = inEnd;
      cmapLapParaTaskDeepBkz.blocksize = inBlockSize;
      cmapLapParaTaskDeepBkz.u         = inU;
      cmapLapParaTaskDeepBkz.seed      = inSeed;
      cmapLapParaTaskDeepBkz.basis     = inBasis;
   }

   ///
   ///  constructor of Enum
   ///
   CMapLapParaTask(
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
         : UG::ParaTask(inTaskId, inGeneratorTaskId, inEstimatedValue, 0),
           threadId(inThreadId),
           solverType(Enum)
   {
      cmapLapParaTaskEnum.begin   = inBegin;
      cmapLapParaTaskEnum.end     = inEnd;
      cmapLapParaTaskEnum.start   = inStart;
      cmapLapParaTaskEnum.last    = inLast;
      cmapLapParaTaskEnum.prob    = inProb;
      cmapLapParaTaskEnum.coeffs  = inCoeffs;
      cmapLapParaTaskEnum.basis   = inBasis;
   }

   ///
   ///  constructor of Sieve
   ///
   CMapLapParaTask(
         UG::TaskId inTaskId,                ///< task id
         UG::TaskId inGeneratorTaskId,       ///< generator task id
         double   inEstimatedValue,          ///< estimated value
         int      inThreadId,                ///< thread ID
         int      inBegin,                   ///< 1st index of block matrix
         int      inEnd,                     ///< last index of block matrix
         std::shared_ptr<LatticeBasis<int>> inBasis   ///< lattice basis
         )
         : UG::ParaTask(inTaskId, inGeneratorTaskId, inEstimatedValue, 0),
           threadId(inThreadId),
           solverType(Sieve)
   {
      cmapLapParaTaskSieve.begin = inBegin;
      cmapLapParaTaskSieve.end   = inEnd;
      cmapLapParaTaskSieve.basis = inBasis;
   }

   ///
   ///  deconstructor
   ///
   virtual ~CMapLapParaTask(
         )
   {
   }

   ///
   /// upudate paraSolverTask DeepBkz
   ///
   virtual void updateDeepBkz(
         int inBlockSize,
         int *inBasis
         )
   {
      assert( solverType == DeepBkz );
      cmapLapParaTaskDeepBkz.blocksize = inBlockSize;
      if( inBasis )
      {
         int nRow = cmapLapParaTaskDeepBkz.basis->rows();
         int nCol = cmapLapParaTaskDeepBkz.basis->cols();
         *cmapLapParaTaskDeepBkz.basis = Eigen::Map<LatticeBasis<int>>(inBasis, nRow, nCol);
      }
   }

   ///
   /// upudate paraSolverTask DeepBkz
   ///
   virtual void updateEnum(
         int  inStart,
         int* inCoeffs
         )
   {
      assert( solverType == Enum );
      cmapLapParaTaskEnum.start = inStart;
      int n = cmapLapParaTaskEnum.coeffs->size();
      *cmapLapParaTaskEnum.coeffs = Eigen::Map<LatticeVector<int>>(inCoeffs, n);
   }

   ///
   /// update paraSolver Sieve
   ///
   virtual void updateSieve(
         int *inBasis
         )
   {
      int nRow = cmapLapParaTaskSieve.basis->rows();
      int nCol = cmapLapParaTaskSieve.basis->cols();
      *cmapLapParaTaskSieve.basis = Eigen::Map<LatticeBasis<int>>(inBasis, nRow, nCol);
   }

   ///
   /// get therad ID
   /// @return thread ID
   ///
   virtual int getThreadId(
         )
   {
      return threadId;
   }

   ///
   /// set therad ID
   ///
   virtual void setThreadId(
         int inThreadId      ///< thread ID
         )
   {
      threadId = inThreadId;
   }

   ///
   /// get solver type
   /// @return solver type
   ///
   virtual SolverType getSolverType(
         )
   {
      return solverType;
   }

   ///
   /// getter for CMapLapParaTaskLC
   ///
   virtual bool isLocalTask(
         )
   {
      return false;
   }

#ifdef UG_WITH_ZLIB

   ///
   /// write to checkpoint file
   ///
   virtual void write(
         gzstream::ogzstream &out   ///< gzstream for output
         )
   {
      THROW_LOGICAL_ERROR1("***** we use 2 argments version defined by ourselves *****");
   }

   ///
   /// write to checkpoint file
   ///
   virtual void write(
         gzstream::ogzstream &out,   ///< gzstream for output
         int activeFlag
         );

   ///
   /// read from checkpoint file
   ///
   virtual bool read(
         UG::ParaComm *comm,            ///< communicator used
         gzstream::igzstream &in        ///< gzstream for input
         );

#endif

   ///
   /// stringfy CMapLapParaNode
   /// @return string to show inside of this object
   ///
   virtual const std::string toString(
         )
   {
      std::ostringstream s;
      s << "CMapLapParaTaskId = " << (taskId.toString()) << ", GeneratorTaskId = " << (generatorTaskId.toString())
      << ", estimated value = " << estimatedValue << std::endl;
      s << toSimpleString();
      return s.str();
   }

   ///
   /// stringfy CMapLapParaNode as simple string
   /// @return string to show inside of this object
   ///
   virtual const std::string toSimpleString(
         )
   {
      std::ostringstream s;
      s << taskId.toString()
            << ", "
            << generatorTaskId.toString()
            << ", "
            << estimatedValue;
      switch( solverType )
      {
      case DeepBkz:
         s << " DeepBkz"
           << ", begin = "       << cmapLapParaTaskDeepBkz.begin
           << ", end = "         << cmapLapParaTaskDeepBkz.end
           << ", u = "           << cmapLapParaTaskDeepBkz.u
           << ", seed = "        << cmapLapParaTaskDeepBkz.seed
           << ", basis.size = (" << cmapLapParaTaskDeepBkz.basis->rows()
           << ", "               << cmapLapParaTaskDeepBkz.basis->cols()
           << ")"
           << std::endl;
         break;
      case Enum:
         s << " Enum"
           << ", begin = "       << cmapLapParaTaskEnum.begin
           << ", end = "         << cmapLapParaTaskEnum.end
           << ", start = "       << cmapLapParaTaskEnum.start
           << ", last = "        << cmapLapParaTaskEnum.last
           << ", prob = "        << cmapLapParaTaskEnum.prob
           << ", coeffs = "      << cmapLapParaTaskEnum.coeffs->transpose()
           << ", basis.size = (" << cmapLapParaTaskEnum.basis->rows()
           << ", "               << cmapLapParaTaskEnum.basis->cols()
           << ")"
           << std::endl;
         break;
      case Sieve:
         s << " Sieve"
           << ", begin = "       << cmapLapParaTaskSieve.begin
           << ", end = "         << cmapLapParaTaskSieve.end
           << ", basis.size = (" << cmapLapParaTaskSieve.basis->rows()
           << ", "               << cmapLapParaTaskSieve.basis->cols()
           << ")"
           << std::endl;
         break;
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaTask::toSimpleString: Invalid solver type = ", static_cast<int>(solverType));
      }

      return s.str();
   }


   virtual const std::string toStartString(
         std::string delimiter=""
         )
   {
      switch( solverType )
      {
      case DeepBkz:
         return toStringLogDeepBkz(delimiter);
      case Enum:
         return toStringLogEnum(delimiter);
      case Sieve:
         return toStringLogSieve(delimiter);
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaTask::toStartString: Invalid solver type = ",static_cast<int>(solverType));
      }
   }

   ///
   /// stringfy CMapLapParaSolverStateData
   /// @return string to show inside of CMapLapParaSolverStateDeepBkz
   ///
   virtual const std::string toStringLogDeepBkz(
         std::string delimiter=""
         )
   {
      std::string taskName             = "DeepBkz";
      double      elapsedTime          = 0.0;
      int         size                 = cmapLapParaTaskDeepBkz.blocksize;
      long int    iter                 = 0;
      double      progress             = 0.0;
      double      leftTime             = -1.0;
      double      logCost              = -1.0;
      double      shortestNorm         = -1.0;
      double      approxFactor         = -1.0;
      double      hermiteFactor        = -1.0;
      double      rootHermiteFactor    = -1.0;
      double      orthogonalFactor     = -1.0;
      int         meanMessageQueueSize = -1;
      int         maxMessageQueueSize  = -1;
      std::string appendix             = "";

      return Logging::toStringLogBase(taskName, elapsedTime, size, iter, progress, leftTime, logCost,
            shortestNorm, approxFactor, hermiteFactor, rootHermiteFactor, orthogonalFactor,
            meanMessageQueueSize, maxMessageQueueSize,
            appendix, delimiter);
   }

   ///
   /// stringfy CMapLapParaSolverStateData
   /// @return string to show inside of CMapLapParaSolverStateEnum
   ///
   virtual const std::string toStringLogEnum(
         std::string delimiter=""
         )
   {
      std::string taskName             = "Enum";
      double      elapsedTime          = 0.0;
      int         size                 = cmapLapParaTaskEnum.end
                                       - cmapLapParaTaskEnum.begin + 1;
      long int    iter                 = 0;
      double      progress             = 0.0;
      double      leftTime             = -1.0;
      double      logCost              = -1.0;
      double      shortestNorm         = -1.0;
      double      approxFactor         = -1.0;
      double      hermiteFactor        = -1.0;
      double      rootHermiteFactor    = -1.0;
      double      orthogonalFactor     = -1.0;
      int         meanMessageQueueSize = -1;
      int         maxMessageQueueSize  = -1;
      std::string appendix             = "";

      return Logging::toStringLogBase(taskName, elapsedTime, size, iter, progress, leftTime, logCost,
            shortestNorm, approxFactor, hermiteFactor, rootHermiteFactor, orthogonalFactor,
            meanMessageQueueSize, maxMessageQueueSize,
            appendix, delimiter);
   }

   ///
   /// stringfy CMapLapParaSolverStateData
   /// @return string to show inside of CMapLapParaSolverStateSieve
   ///
   virtual const std::string toStringLogSieve(
         std::string delimiter=""
         )
   {
      std::string taskName             = "Sieve";
      double      elapsedTime          = 0.0;
      int         size                 = cmapLapParaTaskSieve.end
                                       - cmapLapParaTaskSieve.begin + 1;
      long int    iter                 = 0;
      double      progress             = 0.0;
      double      leftTime             = -1.0;
      double      logCost              = -1.0;
      double      shortestNorm         = -1.0;
      double      approxFactor         = -1.0;
      double      hermiteFactor        = -1.0;
      double      rootHermiteFactor    = -1.0;
      double      orthogonalFactor     = -1.0;
      int         meanMessageQueueSize = -1;
      int         maxMessageQueueSize  = -1;
      std::string appendix             = "";

      return Logging::toStringLogBase(taskName, elapsedTime, size, iter, progress, leftTime, logCost,
            shortestNorm, approxFactor, hermiteFactor, rootHermiteFactor, orthogonalFactor,
            meanMessageQueueSize, maxMessageQueueSize,
            appendix, delimiter);
   }


   ///
   /// getter of dimension of basis
   ///
   virtual int getDimension(
         )
   {
      switch( solverType )
      {
      case DeepBkz:
         return cmapLapParaTaskDeepBkz.basis->cols();
      case Enum:
         return cmapLapParaTaskEnum.basis->cols();
      case Sieve:
         return cmapLapParaTaskSieve.basis->cols();
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaTask::getNumBasisVectors: Invalid solver type = ",static_cast<int>(solverType));
      }
   }


   ///
   /// getter of number of vectors in basis = number of rows of basis
   ///
   virtual int getNumBasisVectors(
         )
   {
      switch( solverType )
      {
      case DeepBkz:
         return cmapLapParaTaskDeepBkz.basis->rows();
      case Enum:
         return cmapLapParaTaskEnum.basis->rows();
      case Sieve:
         return cmapLapParaTaskSieve.basis->rows();
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaTask::getNumBasisVectors: Invalid solver type = ", static_cast<int>(solverType));
      }
   }

   ///
   /// getter of the basis of lattice
   /// @return the basis of lattice
   ///
   virtual std::shared_ptr<LatticeBasis<int>> getBasis()
   {
      switch( solverType )
      {
      case DeepBkz:
         return cmapLapParaTaskDeepBkz.basis;
      case Enum:
         return cmapLapParaTaskEnum.basis;
      case Sieve:
         return cmapLapParaTaskSieve.basis;
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaTask::getBasis: Invalid solver type = ", static_cast<int>(solverType));
      }
   }

   ///
   /// getter of the coefficients of lattice basis
   /// @return the coefficients of lattice basis
   ///
   virtual std::shared_ptr<LatticeVector<int>> getCoeffs()
   {
      assert( solverType == Enum );
      return cmapLapParaTaskEnum.coeffs;
   }

   ///
   /// getter of 1st index of block matrix
   /// @return 1st index of block matrix
   ///
   virtual int getBegin()
   {
      switch( solverType )
      {
      case DeepBkz:
         return cmapLapParaTaskDeepBkz.begin;
      case Enum:
         return cmapLapParaTaskEnum.begin;
      case Sieve:
         return cmapLapParaTaskSieve.begin;
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaTask::getBegin: Invalid solver type = ",static_cast<int>(solverType));
      }
   }

   ///
   /// getter of last index of block matrix
   /// @return the last index of block matrix
   ///
   virtual int getEnd()
   {
      switch( solverType )
      {
      case DeepBkz:
         return cmapLapParaTaskDeepBkz.end;
      case Enum:
         return cmapLapParaTaskEnum.end;
      case Sieve:
         return cmapLapParaTaskSieve.end;
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaTask::getEnd: Invalid solver type = ",static_cast<int>(solverType));
      }
   }

   ///
   /// getter of index of current search node in the enumeration tree
   /// @return start
   ///
   virtual int getStart()
   {
      assert( solverType == Enum );
      return cmapLapParaTaskEnum.start;
   }

   ///
   /// getter of last attributes in CMapLapParaTaskEnum
   /// @return last attributes
   ///
   virtual int getLast()
   {
      assert( solverType == Enum );
      return cmapLapParaTaskEnum.last;
   }

   ///
   /// getter of the size of unimodular matrix
   /// @return the size of unimodular matrix
   ///
   virtual int getU()
   {
      assert( solverType == DeepBkz );
      return cmapLapParaTaskDeepBkz.u;
   }

   ///
   /// setter of the size of unimodular matrix
   ///
   virtual void setU(
         int inU
         )
   {
      assert( solverType == DeepBkz );
      cmapLapParaTaskDeepBkz.u = inU;
   }

   ///
   /// getter of the seed which is a randomize parameter
   /// @return the seed which is a randomize parameter
   ///
   virtual int getSeed()
   {
      assert( solverType == DeepBkz );
      return cmapLapParaTaskDeepBkz.seed;
   }

   ///
   /// getter of the probability of the extreme pruning
   /// @return prob
   ///
   virtual double getProb()
   {
      assert( solverType == Enum );
      return cmapLapParaTaskEnum.prob;
   }

};



///
/// class CMapLapParaTask for Local Solver
///
class CMapLapParaTaskLC : public CMapLapParaTask
{
   CMapLapParaTaskLC *createDatatype(
         UG::ParaComm *comm
         );

public:
   ///
   /// constructor DeepBkz
   ///
   CMapLapParaTaskLC(
         UG::TaskId inTaskId,                ///< task id
         UG::TaskId inGeneratorTaskId,       ///< generator task id
         double inEstimatedValue,            ///< estimated value
         int inThreadId,                     ///< thread ID
         int inBegin,                        ///< 1st index of block matrix
         int inEnd,                          ///< last index of block matrix
         int inBlockSize,                    ///< the current blocksize
         int inU,                            ///< the unimodular matrix size
         int inSeed,                         ///< the seed of randomize
         std::shared_ptr<LatticeBasis<int>> inBasis   ///< sub lattice basis (row <= col)
         )
         :
         CMapLapParaTask(
               inTaskId,
               inGeneratorTaskId,
               inEstimatedValue,
               inThreadId,
               inBegin,
               inEnd,
               inBlockSize,
               inU,
               inSeed,
               inBasis)
   {}

   virtual ~CMapLapParaTaskLC(
      )
   {
      assert( solverType == DeepBkz );
   }

   ///
   /// getter for CMapLapParaTaskLC
   ///
   virtual bool isLocalTask(
         )
   {
      return true;
   }

   ///
   /// clone this ParaTask
   /// @return pointer to cloned ParaTask object
   ///
   virtual ParaTask* clone(
         UG::ParaComm *comm   ///< communicator used
         )
   {
      assert( solverType == DeepBkz );
      return ( new
         CMapLapParaTaskLC(taskId, generatorTaskId, estimatedValue, threadId,
               cmapLapParaTaskDeepBkz.begin,
               cmapLapParaTaskDeepBkz.end,
               cmapLapParaTaskDeepBkz.blocksize,
               cmapLapParaTaskDeepBkz.u,
               cmapLapParaTaskDeepBkz.seed,
               std::make_shared<LatticeBasis<int>>(*(cmapLapParaTaskDeepBkz.basis))
            )
         );
   }

   ///
   /// broadcast this object
   /// @return always 0 (for future extensions)
   ///
   virtual int bcast(
         UG::ParaComm *comm,  ///< communicator used
         int root             ///< root rank of broadcast
         );

   ///
   /// send this object
   /// @return always 0 (for future extensions)
   ///
   virtual int send(
         UG::ParaComm *comm,  ///< communicator used
         int destination      ///< destination rank
         );

   ///
   /// receive this object
   /// @return always 0 (for future extensions)
   ///
   virtual int receive(
         UG::ParaComm *comm,  ///< communicator used
         int source           ///< source rank
         );

};

} // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_TASK_H__
