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

/**@file    cmapLapParaTaskTh.cpp
 * @brief   CMapLapParaTask extension for threads communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "cmapLapParaTaskTh.h"
#include <memory>
#include "ug/paraDef.h"
#include "ug/paraComm.h"
#include "ug/paraTask.h"
#include "ug/paraTagDef.h"
#include "cmapLapParaDef.h"
#include "cmapLapParaTask.h"
#include "cmapLapParaSolverLocalComm.h"


namespace ParaCMapLAP
{

CMapLapParaTaskTh *
CMapLapParaTaskTh::createDatatype(
      UG::ParaComm *comm
      )
{
   return dynamic_cast<CMapLapParaTaskTh *>(clone(comm));
}

int
CMapLapParaTaskTh::bcast(
      UG::ParaComm *comm,
      int root
      )
{

#ifdef _COMM_CPP11
   CMapLapParaCommTh *commTh = dynamic_cast<CMapLapParaCommTh*>(comm);
   if( commTh )
   {
      if( commTh->getRank() == root )
      {
         for( int i = 0; i < commTh->getSize(); i++ )
         {
            if( i != root )
            {
               CMapLapParaTaskTh *sent;
               sent = createDatatype(comm);
               PARA_COMM_CALL(
                     commTh->uTypeSend((void *)sent, ParaTaskType, i, UG::TagTask)
                     );
            }
         }
      }
      else
      {
         std::unique_ptr<CMapLapParaTaskTh> received;
         PARA_COMM_CALL(
               commTh->uTypeReceive((void **)&received, ParaTaskType, root, UG::TagTask)
               );
         taskId = received->taskId;
         generatorTaskId = received->generatorTaskId;
         estimatedValue = received->estimatedValue;
         diffSubproblemInfo = received->diffSubproblemInfo;
      }
   }
   else
   {
      CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
      if( localComm->getRank() == root )
      {
         for( int i = 0; i < localComm->getSize(); i++ )
         {
            if( i != root )
            {
               CMapLapParaTaskTh *sent;
               sent = createDatatype(localComm);
               PARA_COMM_CALL(
                     localComm->uTypeSend((void *)sent, ParaTaskType, i, UG::TagTask)
                     );
            }
         }
      }
      else
      {
         std::unique_ptr<CMapLapParaTaskTh> received;
         PARA_COMM_CALL(
               localComm->uTypeReceive((void **)&received, ParaTaskType, root, UG::TagTask)
               );
         taskId = received->taskId;
         generatorTaskId = received->generatorTaskId;
         estimatedValue = received->estimatedValue;
         diffSubproblemInfo = received->diffSubproblemInfo;
      }
   }
#else
   CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
   if( localComm->getRank() == root )
   {
      for( int i = 0; i < localComm->getSize(); i++ )
      {
         if( i != root )
         {
            CMapLapParaTask *sent;
            sent = createDatatype(localComm);
            PARA_COMM_CALL(
                  localComm->uTypeSend((void *)sent, ParaCMapLAP::ParaTaskType, i, UG::TagTask)
                  );
         }
      }
   }
   else
   {
      std::unique_ptr<CMapLapParaTask> received;
      PARA_COMM_CALL(
            localComm->uTypeReceive((void **)&received, ParaCMapLAP::ParaTaskType, root, UG::TagTask)
            );
      taskId = received->taskId;
      generatorTaskId = received->generatorTaskId;
      estimatedValue = received->estimatedValue;
      diffSubproblemInfo = received->diffSubproblemInfo;
   }
#endif

   return 0;
}

int
CMapLapParaTaskTh::send(
      UG::ParaComm *comm,
      int destination
      )
{
   CMapLapParaTask *sent;

#ifdef _COMM_CPP11
   CMapLapParaCommTh *commTh = dynamic_cast<CMapLapParaCommTh*>(comm);
   if( commTh )
   {
      sent = createDatatype(commTh);
      PARA_COMM_CALL(
            commTh->uTypeSend((void *)sent, ParaTaskType, destination, UG::TagTask)
            );
   }
   else
   {
      CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
      sent = createDatatype(localComm);
      PARA_COMM_CALL(
            localComm->uTypeSend((void *)sent, ParaTaskType, destination, UG::TagTask)
            );
   }
#else
   CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
   sent = createDatatype(localComm);
   PARA_COMM_CALL(
         localComm->uTypeSend((void *)sent, ParaCMapLAP::ParaTaskType, destination, UG::TagTask)
         );
#endif

   return 0;
}

int
CMapLapParaTaskTh::receive(
      UG::ParaComm *comm,
      int source
      )
{

   std::unique_ptr<CMapLapParaTask> received;

#ifdef _COMM_CPP11
   CMapLapParaCommTh *commTh = dynamic_cast<CMapLapParaCommTh*>(comm);
   if( commTh )
   {
      PARA_COMM_CALL(
         commTh->uTypeReceive((void **)&received, ParaTaskType, source, UG::TagTask)
      );
   }
   else
   {
      CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
      PARA_COMM_CALL(
            localComm->uTypeReceive((void **)&received, ParaTaskType, source, UG::TagTask)
      );
   }
#else
   CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
   PARA_COMM_CALL(
         localComm->uTypeReceive((void **)&received, ParaCMapLAP::ParaTaskType, source, UG::TagTask)
   );
#endif

   taskId          = received->taskId;
   generatorTaskId = received->generatorTaskId;
   estimatedValue  = received->estimatedValue;
   threadId        = received->threadId;
   solverType      = received->solverType;
   switch( solverType )
   {
   case DeepBkz:
      cmapLapParaTaskDeepBkz.begin     = received->cmapLapParaTaskDeepBkz.begin;
      cmapLapParaTaskDeepBkz.end       = received->cmapLapParaTaskDeepBkz.end;
      cmapLapParaTaskDeepBkz.blocksize = received->cmapLapParaTaskDeepBkz.blocksize;
      cmapLapParaTaskDeepBkz.u         = received->cmapLapParaTaskDeepBkz.u;
      cmapLapParaTaskDeepBkz.seed      = received->cmapLapParaTaskDeepBkz.seed;
      cmapLapParaTaskDeepBkz.basis     = std::make_shared<LatticeBasis<int>>(*(received->cmapLapParaTaskDeepBkz.basis));
      break;
   case Enum:
      cmapLapParaTaskEnum.begin   = received->cmapLapParaTaskEnum.begin;
      cmapLapParaTaskEnum.end     = received->cmapLapParaTaskEnum.end;
      cmapLapParaTaskEnum.start   = received->cmapLapParaTaskEnum.start;
      cmapLapParaTaskEnum.last    = received->cmapLapParaTaskEnum.last;
      cmapLapParaTaskEnum.prob    = received->cmapLapParaTaskEnum.prob;
      cmapLapParaTaskEnum.coeffs  = std::make_shared<LatticeVector<int>>(*(received->cmapLapParaTaskEnum.coeffs));
      cmapLapParaTaskEnum.basis   = std::make_shared<LatticeBasis<int>>(*(received->cmapLapParaTaskEnum.basis));
      break;
   case Sieve:
      cmapLapParaTaskSieve.begin  = received->cmapLapParaTaskSieve.begin;
      cmapLapParaTaskSieve.end    = received->cmapLapParaTaskSieve.end;
      cmapLapParaTaskSieve.basis  = std::make_shared<LatticeBasis<int>>(*(received->cmapLapParaTaskSieve.basis));
      break;
   default:
      THROW_LOGICAL_ERROR2("CMapLapParaTaskTh::clone: Invalid solver type = ", static_cast<int>(solverType));
   }

   return 0;
}

} // namespace ParaCMapLAP
