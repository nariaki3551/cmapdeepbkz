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

/**@file    cmapLapParaTask.cpp
 * @brief   Base class for CMapLapParaTask.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "cmapLapParaTask.h"
#include <vector>
#include <float.h>
#include "ug/paraComm.h"
#include "ug/paraTagDef.h"
#include "cmapLapParaDef.h"
#include "cmapLapParaSolverLocalComm.h"
#ifdef UG_WITH_ZLIB
#include "ug/gzstream.h"
#endif


namespace ParaCMapLAP
{

#ifdef UG_WITH_ZLIB
void
CMapLapParaTask::write(
      gzstream::ogzstream &out,
      int activeFlag
      )
{
   out.write(reinterpret_cast<char *>(&taskId.subtaskId.lcId), sizeof(int));
   out.write(reinterpret_cast<char *>(&taskId.subtaskId.globalSubtaskIdInLc), sizeof(int));
   out.write(reinterpret_cast<char *>(&taskId.subtaskId.solverId), sizeof(int));
   out.write(reinterpret_cast<char *>(&taskId.seqNum), sizeof(long long));
   out.write(reinterpret_cast<char *>(&generatorTaskId.subtaskId.lcId), sizeof(int));
   out.write(reinterpret_cast<char *>(&generatorTaskId.subtaskId.globalSubtaskIdInLc), sizeof(int));
   out.write(reinterpret_cast<char *>(&generatorTaskId.subtaskId.solverId), sizeof(int));
   out.write(reinterpret_cast<char *>(&generatorTaskId.seqNum), sizeof(long long));
   if( activeFlag == 1)
   {
      // active is high priority
      double highestPriorityValue = -DBL_MAX;
      out.write(reinterpret_cast<char *>(&highestPriorityValue), sizeof(double));
   }
   else
   {
      out.write(reinterpret_cast<char *>(&estimatedValue), sizeof(double));
   }
   out.write(reinterpret_cast<char *>(&diffSubproblemInfo), sizeof(int));
   assert( diffSubproblemInfo == 0 );
   switch( solverType )
   {
   case DeepBkz:
   {
      int iSolverType = static_cast<int>(DeepBkz);
      out.write(reinterpret_cast<char *>(&iSolverType), sizeof(int));
      out.write(reinterpret_cast<char *>(&cmapLapParaTaskDeepBkz.begin), sizeof(int));
      out.write(reinterpret_cast<char *>(&cmapLapParaTaskDeepBkz.end), sizeof(int));
      out.write(reinterpret_cast<char *>(&cmapLapParaTaskDeepBkz.blocksize), sizeof(int));
      out.write(reinterpret_cast<char *>(&cmapLapParaTaskDeepBkz.u), sizeof(int));
      out.write(reinterpret_cast<char *>(&cmapLapParaTaskDeepBkz.seed), sizeof(int));

      int m = cmapLapParaTaskDeepBkz.basis->rows();
      int n = cmapLapParaTaskDeepBkz.basis->cols();
      out.write(reinterpret_cast<char *>(&m), sizeof(int));
      out.write(reinterpret_cast<char *>(&n), sizeof(int));

      std::vector<int> basis(m*n);
      Eigen::Map<LatticeBasis<int>>(&basis[0], m, n) = *cmapLapParaTaskDeepBkz.basis;
      out.write(reinterpret_cast<char *>(&basis[0]), sizeof(int)*(m*n));
      break;
   }
   case Enum:
   {
      int iSolverType = static_cast<int>(Enum);
      out.write(reinterpret_cast<char *>(&iSolverType), sizeof(int));
      out.write(reinterpret_cast<char *>(&cmapLapParaTaskEnum.begin), sizeof(int));
      out.write(reinterpret_cast<char *>(&cmapLapParaTaskEnum.end), sizeof(int));
      out.write(reinterpret_cast<char *>(&cmapLapParaTaskEnum.start), sizeof(int));
      out.write(reinterpret_cast<char *>(&cmapLapParaTaskEnum.last), sizeof(int));
      out.write(reinterpret_cast<char *>(&cmapLapParaTaskEnum.prob), sizeof(double));

      int m = cmapLapParaTaskEnum.basis->rows();
      int n = cmapLapParaTaskEnum.basis->cols();
      out.write(reinterpret_cast<char *>(&m), sizeof(int));
      out.write(reinterpret_cast<char *>(&n), sizeof(int));

      std::vector<int> coeffs(m);
      Eigen::Map<LatticeVector<int>>(&coeffs[0], m) = *cmapLapParaTaskEnum.coeffs;
      out.write(reinterpret_cast<char *>(&coeffs[0]), sizeof(int)*(m));

      std::vector<int> basis(m*n);
      Eigen::Map<LatticeBasis<int>>(&basis[0], m, n) = *cmapLapParaTaskEnum.basis;
      out.write(reinterpret_cast<char *>(&basis[0]), sizeof(int)*(m*n));
      break;
   }
   case Sieve:
   {
      int iSolverType = static_cast<int>(Sieve);
      out.write(reinterpret_cast<char *>(&iSolverType), sizeof(int));
      out.write(reinterpret_cast<char *>(&cmapLapParaTaskSieve.begin), sizeof(int));
      out.write(reinterpret_cast<char *>(&cmapLapParaTaskSieve.end), sizeof(int));

      int m = cmapLapParaTaskSieve.basis->rows();
      int n = cmapLapParaTaskSieve.basis->cols();
      out.write(reinterpret_cast<char *>(&m), sizeof(int));
      out.write(reinterpret_cast<char *>(&n), sizeof(int));

      std::vector<int> basis(m*n);
      Eigen::Map<LatticeBasis<int>>(&basis[0], m, n) = *cmapLapParaTaskSieve.basis;
      out.write(reinterpret_cast<char *>(&basis[0]), sizeof(int)*(m*n));
      break;
   }
   default:
      THROW_LOGICAL_ERROR2("MCMapLapParaTask::write: Invalid solver type = ", static_cast<int>(solverType));
   }
}

bool
CMapLapParaTask::read(
      UG::ParaComm *comm,
      gzstream::igzstream &in
      )
{
   in.read(reinterpret_cast<char *>(&taskId.subtaskId.lcId), sizeof(int));
   if( in.eof() ) return false;
   in.read(reinterpret_cast<char *>(&taskId.subtaskId.globalSubtaskIdInLc), sizeof(int));
   in.read(reinterpret_cast<char *>(&taskId.subtaskId.solverId), sizeof(int));
   in.read(reinterpret_cast<char *>(&taskId.seqNum), sizeof(long long));
   in.read(reinterpret_cast<char *>(&generatorTaskId.subtaskId.lcId), sizeof(int));
   in.read(reinterpret_cast<char *>(&generatorTaskId.subtaskId.globalSubtaskIdInLc), sizeof(int));
   in.read(reinterpret_cast<char *>(&generatorTaskId.subtaskId.solverId), sizeof(int));
   in.read(reinterpret_cast<char *>(&generatorTaskId.seqNum), sizeof(long long));
   in.read(reinterpret_cast<char *>(&estimatedValue), sizeof(double));
   in.read(reinterpret_cast<char *>(&diffSubproblemInfo), sizeof(int));
   assert( diffSubproblemInfo == 0);

   int iSolverType = 0;
   in.read(reinterpret_cast<char *>(&iSolverType), sizeof(int));
   switch( iSolverType )
   {
      case static_cast<int>(DeepBkz):
         solverType = DeepBkz;
         break;
      case static_cast<int>(Enum):
         solverType = Enum;
         break;
      case static_cast<int>(Sieve):
         solverType = Sieve;
         break;
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaTask::read: Invalid solver type = ", static_cast<int>(iSolverType));
   }

   switch( solverType )
   {
   case DeepBkz:
   {
      in.read(reinterpret_cast<char *>(&cmapLapParaTaskDeepBkz.begin), sizeof(int));
      in.read(reinterpret_cast<char *>(&cmapLapParaTaskDeepBkz.end), sizeof(int));
      in.read(reinterpret_cast<char *>(&cmapLapParaTaskDeepBkz.blocksize), sizeof(int));
      in.read(reinterpret_cast<char *>(&cmapLapParaTaskDeepBkz.u), sizeof(int));
      in.read(reinterpret_cast<char *>(&cmapLapParaTaskDeepBkz.seed), sizeof(int));

      int m = 0, n = 0;
      in.read(reinterpret_cast<char *>(&m), sizeof(int));
      in.read(reinterpret_cast<char *>(&n), sizeof(int));

      std::vector<int> basis(m*n);
      in.read(reinterpret_cast<char *>(&basis[0]), sizeof(int)*(m*n));
      cmapLapParaTaskDeepBkz.basis = std::make_shared<LatticeBasis<int>>(Eigen::Map<LatticeBasis<int>>(&basis[0], m, n));
      break;
   }
   case Enum:
   {
      in.read(reinterpret_cast<char *>(&cmapLapParaTaskEnum.begin), sizeof(int));
      in.read(reinterpret_cast<char *>(&cmapLapParaTaskEnum.end), sizeof(int));
      in.read(reinterpret_cast<char *>(&cmapLapParaTaskEnum.start), sizeof(int));
      in.read(reinterpret_cast<char *>(&cmapLapParaTaskEnum.last), sizeof(int));
      in.read(reinterpret_cast<char *>(&cmapLapParaTaskEnum.prob), sizeof(double));

      int m = 0, n = 0;
      in.read(reinterpret_cast<char *>(&m), sizeof(int));
      in.read(reinterpret_cast<char *>(&n), sizeof(int));

      std::vector<int> coeffs(m);
      in.read(reinterpret_cast<char *>(&coeffs[0]), sizeof(int)*(m));
      cmapLapParaTaskEnum.coeffs = std::make_shared<LatticeVector<int>>(Eigen::Map<LatticeVector<int>>(&coeffs[0], m));

      std::vector<int> basis(m*n);
      in.read(reinterpret_cast<char *>(&basis[0]), sizeof(int)*(m*n));
      cmapLapParaTaskEnum.basis = std::make_shared<LatticeBasis<int>>(Eigen::Map<LatticeBasis<int>>(&basis[0], m, n));
      break;
   }
   case Sieve:
   {
      in.read(reinterpret_cast<char *>(&cmapLapParaTaskSieve.begin), sizeof(int));
      in.read(reinterpret_cast<char *>(&cmapLapParaTaskSieve.end), sizeof(int));

      int m = 0, n = 0;
      in.read(reinterpret_cast<char *>(&m), sizeof(int));
      in.read(reinterpret_cast<char *>(&n), sizeof(int));

      std::vector<int> basis(m*n);
      in.read(reinterpret_cast<char *>(&basis[0]), sizeof(int)*(m*n));
      cmapLapParaTaskSieve.basis = std::make_shared<LatticeBasis<int>>(Eigen::Map<LatticeBasis<int>>(&basis[0], m, n));
      break;
   }
   default:
      THROW_LOGICAL_ERROR2("CMapLapParaTask::read: Invalid solver type = ", static_cast<int>(solverType));
   }

   return true;
}

#endif

CMapLapParaTaskLC *
CMapLapParaTaskLC::createDatatype(
      UG::ParaComm *comm
      )
{
   return dynamic_cast<CMapLapParaTaskLC *>(clone(comm));
}

int
CMapLapParaTaskLC::bcast(
      UG::ParaComm *comm,
      int root
      )
{
   CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
   if( localComm->getRank() == root )
   {
      for( int i = 0; i < localComm->getSize(); i++ )
      {
         if( i != root )
         {
            CMapLapParaTaskLC *sent;
            sent = createDatatype(localComm);
            PARA_COMM_CALL(
                  localComm->uTypeSend((void *)sent, ParaCMapLAP::ParaTaskType, i, UG::TagTask)
                  );
         }
      }
   }
   else
   {
      CMapLapParaTaskLC *received;
      PARA_COMM_CALL(
            localComm->uTypeReceive((void **)&received, ParaCMapLAP::ParaTaskType, root, UG::TagTask)
            );
      taskId = received->taskId;
      generatorTaskId = received->generatorTaskId;
      estimatedValue = received->estimatedValue;
      diffSubproblemInfo = received->diffSubproblemInfo;
      delete received;
   }

   return 0;
}

int
CMapLapParaTaskLC::send(
      UG::ParaComm *comm,
      int destination
      )
{

   CMapLapParaTaskLC *sent;

   CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
   sent = createDatatype(localComm);
   PARA_COMM_CALL(
         localComm->uTypeSend((void *)sent, ParaCMapLAP::ParaTaskType, destination, UG::TagTask)
         );

   return 0;
}

int
CMapLapParaTaskLC::receive(
      UG::ParaComm *comm,
      int source
      )
{

   std::unique_ptr<CMapLapParaTaskLC> received;

   CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
   PARA_COMM_CALL(
         localComm->uTypeReceive((void **)&received, ParaCMapLAP::ParaTaskType, source, UG::TagTask)
   );

   taskId            = received->taskId;
   generatorTaskId   = received->generatorTaskId;
   estimatedValue    = received->estimatedValue;
   threadId          = received->threadId;
   solverType        = received->solverType;

   assert( solverType == DeepBkz );

   cmapLapParaTaskDeepBkz.begin     = received->cmapLapParaTaskDeepBkz.begin;
   cmapLapParaTaskDeepBkz.end       = received->cmapLapParaTaskDeepBkz.end;
   cmapLapParaTaskDeepBkz.blocksize = received->cmapLapParaTaskDeepBkz.blocksize;
   cmapLapParaTaskDeepBkz.u         = received->cmapLapParaTaskDeepBkz.u;
   cmapLapParaTaskDeepBkz.seed      = received->cmapLapParaTaskDeepBkz.seed;
   cmapLapParaTaskDeepBkz.basis     = std::make_shared<LatticeBasis<int>>(*(received->cmapLapParaTaskDeepBkz.basis));

   return 0;
}

} // namespace ParaCMapLAP
