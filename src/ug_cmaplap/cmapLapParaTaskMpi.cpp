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

/**@file    cmapLapParaTaskMpi.cpp
 * @brief   CMapLapParaTask extension for MIP communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "cmapLapParaTaskMpi.h"
#include <memory>
#include <vector>
#include <mpi.h>
#include "ug/paraComm.h"
#include "ug/paraTagDef.h"
#include "cmapLapParaComm.h"
#include "cmapLapParaTask.h"
#include "cmapLapParaTagDef.h"
#include "cmapLapParaIsendRequest.h"
#include "cmapLapParaSolverLocalComm.h"


namespace ParaCMapLAP
{

MPI_Datatype
CMapLapParaTaskMpi::createDatatypeDeepBkz(
      )
{
   const int nBlocks = 17;

   MPI_Datatype datatype;

   MPI_Aint startAddress = 0;
   MPI_Aint address = 0;

   int blockLengths[nBlocks];
   MPI_Aint displacements[nBlocks];
   MPI_Datatype types[nBlocks];

   for( int i = 0; i < nBlocks; i++ ){
      blockLengths[i] = 1;
      types[i] = MPI_INT;
   }

   MPI_CALL(
      MPI_Get_address( &taskId.subtaskId.lcId, &startAddress )
   );
   displacements[0] = 0;

   MPI_CALL(
      MPI_Get_address( &taskId.subtaskId.globalSubtaskIdInLc, &address )
   );
   displacements[1] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &taskId.subtaskId.solverId, &address )
   );
   displacements[2] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &taskId.seqNum, &address )
   );
   displacements[3] = address - startAddress;
   types[3] = MPI_LONG_LONG;

   MPI_CALL(
      MPI_Get_address( &generatorTaskId.subtaskId.lcId, &address )
   );
   displacements[4] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &generatorTaskId.subtaskId.globalSubtaskIdInLc, &address )
   );
   displacements[5] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &generatorTaskId.subtaskId.solverId, &address )
   );
   displacements[6] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &generatorTaskId.seqNum, &address )
   );
   displacements[7] = address - startAddress;
   types[7] = MPI_LONG_LONG;

   MPI_CALL(
      MPI_Get_address( &estimatedValue, &address )
   );
   displacements[8] = address - startAddress;
   types[8] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &threadId, &address )
   );
   displacements[9] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaTaskDeepBkz.begin, &address )
   );
   displacements[10] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaTaskDeepBkz.end, &address )
   );
   displacements[11] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaTaskDeepBkz.blocksize, &address )
   );
   displacements[12] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaTaskDeepBkz.u, &address )
   );
   displacements[13] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaTaskDeepBkz.seed, &address )
   );
   displacements[14] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &nRows, &address )
   );
   displacements[15] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &nCols, &address )
   );
   displacements[16] = address - startAddress;

   MPI_CALL(
         MPI_Type_create_struct(nBlocks, blockLengths, displacements, types, &datatype)
         );

   return datatype;
}

MPI_Datatype
CMapLapParaTaskMpi::createDatatypeEnum(
      )
{
   const int nBlocks = 17;

   MPI_Datatype datatype;

   MPI_Aint startAddress = 0;
   MPI_Aint address = 0;

   int blockLengths[nBlocks];
   MPI_Aint displacements[nBlocks];
   MPI_Datatype types[nBlocks];

   for( int i = 0; i < nBlocks; i++ ){
      blockLengths[i] = 1;
      types[i] = MPI_INT;
   }

   MPI_CALL(
      MPI_Get_address( &taskId.subtaskId.lcId, &startAddress )
   );
   displacements[0] = 0;

   MPI_CALL(
      MPI_Get_address( &taskId.subtaskId.globalSubtaskIdInLc, &address )
   );
   displacements[1] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &taskId.subtaskId.solverId, &address )
   );
   displacements[2] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &taskId.seqNum, &address )
   );
   displacements[3] = address - startAddress;
   types[3] = MPI_LONG_LONG;

   MPI_CALL(
      MPI_Get_address( &generatorTaskId.subtaskId.lcId, &address )
   );
   displacements[4] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &generatorTaskId.subtaskId.globalSubtaskIdInLc, &address )
   );
   displacements[5] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &generatorTaskId.subtaskId.solverId, &address )
   );
   displacements[6] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &generatorTaskId.seqNum, &address )
   );
   displacements[7] = address - startAddress;
   types[7] = MPI_LONG_LONG;

   MPI_CALL(
      MPI_Get_address( &estimatedValue, &address )
   );
   displacements[8] = address - startAddress;
   types[8] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &threadId, &address )
   );
   displacements[9] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaTaskEnum.begin, &address )
   );
   displacements[10] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaTaskEnum.end, &address )
   );
   displacements[11] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaTaskEnum.start, &address )
   );
   displacements[12] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaTaskEnum.last, &address )
   );
   displacements[13] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaTaskEnum.prob, &address)
   );
   displacements[14] = address - startAddress;
   types[14] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &nRows, &address )
   );
   displacements[15] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &nCols, &address )
   );
   displacements[16] = address - startAddress;

   MPI_CALL(
         MPI_Type_create_struct(nBlocks, blockLengths, displacements, types, &datatype)
         );

   return datatype;
}

MPI_Datatype
CMapLapParaTaskMpi::createDatatypeSieve(
      )
{
   const int nBlocks = 14;

   MPI_Datatype datatype;

   MPI_Aint startAddress = 0;
   MPI_Aint address = 0;

   int blockLengths[nBlocks];
   MPI_Aint displacements[nBlocks];
   MPI_Datatype types[nBlocks];

   for( int i = 0; i < nBlocks; i++ ){
      blockLengths[i] = 1;
      types[i] = MPI_INT;
   }

   MPI_CALL(
      MPI_Get_address( &taskId.subtaskId.lcId, &startAddress )
   );
   displacements[0] = 0;

   MPI_CALL(
      MPI_Get_address( &taskId.subtaskId.globalSubtaskIdInLc, &address )
   );
   displacements[1] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &taskId.subtaskId.solverId, &address )
   );
   displacements[2] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &taskId.seqNum, &address )
   );
   displacements[3] = address - startAddress;
   types[3] = MPI_LONG_LONG;

   MPI_CALL(
      MPI_Get_address( &generatorTaskId.subtaskId.lcId, &address )
   );
   displacements[4] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &generatorTaskId.subtaskId.globalSubtaskIdInLc, &address )
   );
   displacements[5] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &generatorTaskId.subtaskId.solverId, &address )
   );
   displacements[6] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &generatorTaskId.seqNum, &address )
   );
   displacements[7] = address - startAddress;
   types[7] = MPI_LONG_LONG;

   MPI_CALL(
      MPI_Get_address( &estimatedValue, &address )
   );
   displacements[8] = address - startAddress;
   types[8] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &threadId, &address )
   );
   displacements[9] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaTaskSieve.begin, &address )
   );
   displacements[10] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaTaskSieve.end, &address )
   );
   displacements[11] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &nRows, &address )
   );
   displacements[12] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &nCols, &address )
   );
   displacements[13] = address - startAddress;

   MPI_CALL(
         MPI_Type_create_struct(nBlocks, blockLengths, displacements, types, &datatype)
         );

   return datatype;
}


int
CMapLapParaTaskMpi::bcast(
      UG::ParaComm *comm,
      int root
      )
{
   UG::ParaCommMpi *commMpi = dynamic_cast< UG::ParaCommMpi* >(comm);

   if( iSolverType == static_cast<int>(Undefined) )
   {
      iSolverType = static_cast<int>(solverType);
   }

   if( commMpi )
   {
      MPI_CALL(
            commMpi->bcast(&iSolverType, 1, UG::ParaINT, root)
            );

      switch( solverType )
      {
      case DeepBkz:
      {
         solverType = DeepBkz;
         if( comm->getRank() == root )
         {
            MPI_Datatype datatypeDeepBkz = createDatatypeDeepBkz();
            MPI_CALL(
                  MPI_Type_commit( &datatypeDeepBkz )
                  );
            PARA_COMM_CALL(
                  commMpi->ubcast(&taskId.subtaskId.lcId, 1, datatypeDeepBkz, root)
                  );
            MPI_CALL(
                  MPI_Type_free( &datatypeDeepBkz )
                  );

            box.resize(nRows*nCols);
            Eigen::Map<LatticeBasis<int>>(&box[0], nRows, nCols) = *cmapLapParaTaskDeepBkz.basis;
            PARA_COMM_CALL(
                  commMpi->bcast(&box[0], nRows*nCols, UG::ParaINT, root)
                  );
         }
         else
         {
            MPI_Datatype datatypeDeepBkz = createDatatypeDeepBkz();
            MPI_CALL(
                  MPI_Type_commit( &datatypeDeepBkz )
                  );
            PARA_COMM_CALL(
                  commMpi->ubcast(&taskId.subtaskId.lcId, 1, datatypeDeepBkz, root)
                  );
            MPI_CALL(
                  MPI_Type_free( &datatypeDeepBkz )
                  );

            box.resize(nRows*nCols);
            PARA_COMM_CALL(
                  commMpi->bcast(&box[0], nRows*nCols, UG::ParaINT, root)
                  );
            cmapLapParaTaskDeepBkz.basis = std::make_shared<LatticeBasis<int>>(
                  Eigen::Map<LatticeBasis<int>>(&box[0], nRows, nCols)
                  );
         }
         break;
      }
      case Enum:
      {
         solverType = Enum;
         if( comm->getRank() == root )
         {
            MPI_Datatype datatypeEnum = createDatatypeEnum();
            MPI_CALL(
                  MPI_Type_commit( &datatypeEnum )
                  );
            PARA_COMM_CALL(
                  commMpi->ubcast(&taskId.subtaskId.lcId, 1, datatypeEnum, root)
                  );
            MPI_CALL(
                  MPI_Type_free( &datatypeEnum )
                  );

            vec.resize(nCols);
            Eigen::Map<LatticeVector<int>>(&vec[0], nCols) = *cmapLapParaTaskEnum.coeffs;
            PARA_COMM_CALL(
                  commMpi->bcast(&vec[0], nCols, UG::ParaINT, root)
                  );

            box.resize(nRows*nCols);
            Eigen::Map<LatticeBasis<int>>(&box[0], nRows, nCols) = *cmapLapParaTaskEnum.basis;
            PARA_COMM_CALL(
                  commMpi->bcast(&box[0], nRows*nCols, UG::ParaINT, root)
                  );
         }
         else
         {
            MPI_Datatype datatypeDeepBkz = createDatatypeDeepBkz();
            MPI_CALL(
                  MPI_Type_commit( &datatypeDeepBkz )
                  );
            PARA_COMM_CALL(
                  commMpi->ubcast(&taskId.subtaskId.lcId, 1, datatypeDeepBkz, root)
                  );
            MPI_CALL(
                  MPI_Type_free( &datatypeDeepBkz )
                  );

            vec.resize(nCols);
            PARA_COMM_CALL(
                  commMpi->bcast(&vec[0], nCols, UG::ParaINT, root)
                  );
            cmapLapParaTaskEnum.coeffs = std::make_shared<LatticeVector<int>>(
                  Eigen::Map<LatticeVector<int>>(&vec[0], nCols)
                  );

            box.resize(nRows*nCols);
            PARA_COMM_CALL(
                  commMpi->bcast(&box[0], nRows*nCols, UG::ParaINT, root)
                  );
            cmapLapParaTaskEnum.basis = std::make_shared<LatticeBasis<int>>(
                  Eigen::Map<LatticeBasis<int>>(&box[0], nRows, nCols)
                  );
        }
         break;
      }
      case Sieve:
      {
         solverType = Sieve;
         if( comm->getRank() == root )
         {
            MPI_Datatype datatypeSieve = createDatatypeSieve();
            MPI_CALL(
                  MPI_Type_commit( &datatypeSieve )
                  );
            PARA_COMM_CALL(
                  commMpi->ubcast(&taskId.subtaskId.lcId, 1, datatypeSieve, root)
                  );
            MPI_CALL(
                  MPI_Type_free( &datatypeSieve )
                  );

            box.resize(nRows*nCols);
            Eigen::Map<LatticeBasis<int>>(&box[0], nRows, nCols) = *cmapLapParaTaskSieve.basis;
            PARA_COMM_CALL(
                  commMpi->bcast(&box[0], nRows*nCols, UG::ParaINT, root)
                  );
         }
         else
         {
            MPI_Datatype datatypeSieve = createDatatypeSieve();
            MPI_CALL(
                  MPI_Type_commit( &datatypeSieve )
                  );
            PARA_COMM_CALL(
                  commMpi->ubcast(&taskId.subtaskId.lcId, 1, datatypeSieve, root)
                  );
            MPI_CALL(
                  MPI_Type_free( &datatypeSieve )
                  );

            box.resize(nRows*nCols);
            PARA_COMM_CALL(
                  commMpi->bcast(&box[0], nRows*nCols, UG::ParaINT, root)
                  );
            cmapLapParaTaskSieve.basis = std::make_shared<LatticeBasis<int>>(
                  Eigen::Map<LatticeBasis<int>>(&box[0], nRows, nCols)
                  );
         }
         break;
      }
      default:
      {
         THROW_LOGICAL_ERROR2("CMapLapParaTaskMpi::bcast: Invalid iSolver type = ", iSolverType);
      }
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
               CMapLapParaTaskMpi *sent;
               sent = createThreadDatatype(localComm);
               PARA_COMM_CALL(
                     localComm->uTypeSend((void *)sent, ParaCMapLAP::ParaTaskType, i, UG::TagTask)
                     );
            }
         }
      }
      else
      {
         std::unique_ptr<CMapLapParaTaskMpi> received;
         PARA_COMM_CALL(
               localComm->uTypeReceive((void **)&received, ParaCMapLAP::ParaTaskType, root, UG::TagTask)
               );
         taskId = received->taskId;
         generatorTaskId = received->generatorTaskId;
         estimatedValue = received->estimatedValue;
         diffSubproblemInfo = received->diffSubproblemInfo;
      }
   }
   return 0;
}

int
CMapLapParaTaskMpi::send(
      UG::ParaComm *comm,
      int destination
      )
{
   CMapLapParaCommMpi *commMpi = dynamic_cast< CMapLapParaCommMpi* >(comm);

   iSolverType = static_cast<int>(solverType);

   if( commMpi )
   {
	   auto req = std::make_shared<MPI_Request>();
      MPI_CALL(
            commMpi->iSend(&iSolverType, 1, UG::ParaINT, destination, UG::TagTask, req.get())
            );
      auto iSendReq = std::make_shared<CMapLapParaIsendRequest>(req, shared_from_this());
      commMpi->pushISendRequest(iSendReq);

      switch( solverType )
      {
      case DeepBkz:
      {
         MPI_Datatype datatypeDeepBkz = createDatatypeDeepBkz();
         MPI_CALL(
               MPI_Type_commit( &datatypeDeepBkz )
               );
         auto req1 = std::make_shared<MPI_Request>();
         PARA_COMM_CALL(
               commMpi->iUsend(&taskId.subtaskId.lcId, 1, datatypeDeepBkz, destination, TagTask1, req1.get())
               );
         auto iSendReq1 = std::make_shared<CMapLapParaIsendRequest>(req1, shared_from_this());
         commMpi->pushISendRequest(iSendReq1);
         MPI_CALL(
               MPI_Type_free( &datatypeDeepBkz )
               );

         box.resize(nRows*nCols);
         Eigen::Map<LatticeBasis<int>>(&box[0], nRows, nCols) = *cmapLapParaTaskDeepBkz.basis;
         auto req2 = std::make_shared<MPI_Request>();
         PARA_COMM_CALL(
               commMpi->iSend(&box[0], nRows*nCols, UG::ParaINT, destination, TagTask2, req2.get())
               );
         auto iSendReq2 = std::make_shared<CMapLapParaIsendRequest>(req2, shared_from_this());
         commMpi->pushISendRequest(iSendReq2);
         break;
      }
      case Enum:
      {
         MPI_Datatype datatypeEnum = createDatatypeEnum();
         MPI_CALL(
               MPI_Type_commit( &datatypeEnum )
               );
         auto req1 = std::make_shared<MPI_Request>();
         PARA_COMM_CALL(
               commMpi->iUsend(&taskId.subtaskId.lcId, 1, datatypeEnum, destination, TagTask1, req1.get())
               );
         auto iSendReq1 = std::make_shared<CMapLapParaIsendRequest>(req1, shared_from_this());
         commMpi->pushISendRequest(iSendReq1);
         MPI_CALL(
               MPI_Type_free( &datatypeEnum )
               );

         vec.resize(nCols);
         Eigen::Map<LatticeVector<int>>(&vec[0], nCols) = *cmapLapParaTaskEnum.coeffs;
         auto req2 = std::make_shared<MPI_Request>();
         PARA_COMM_CALL(
               commMpi->iSend(&vec[0], nCols, UG::ParaINT, destination, TagTask2, req2.get())
               );
         auto iSendReq2 = std::make_shared<CMapLapParaIsendRequest>(req2, shared_from_this());
         commMpi->pushISendRequest(iSendReq2);

         box.resize(nRows*nCols);
         Eigen::Map<LatticeBasis<int>>(&box[0], nRows, nCols) = *cmapLapParaTaskEnum.basis;
         auto req3 = std::make_shared<MPI_Request>();
         PARA_COMM_CALL(
               commMpi->iSend(&box[0], nRows*nCols, UG::ParaINT, destination, TagTask3, req3.get())
               );
         auto iSendReq3 = std::make_shared<CMapLapParaIsendRequest>(req3, shared_from_this());
         commMpi->pushISendRequest(iSendReq3);
         break;
      }
      case Sieve:
      {
         MPI_Datatype datatypeSieve = createDatatypeSieve();
         MPI_CALL(
               MPI_Type_commit( &datatypeSieve )
               );
         auto req1 = std::make_shared<MPI_Request>();
         PARA_COMM_CALL(
               commMpi->iUsend(&taskId.subtaskId.lcId, 1, datatypeSieve, destination, TagTask1, req1.get())
               );
         auto iSendReq1 = std::make_shared<CMapLapParaIsendRequest>(req1, shared_from_this());
         commMpi->pushISendRequest(iSendReq1);
         MPI_CALL(
               MPI_Type_free( &datatypeSieve )
               );

         box.resize(nRows*nCols);
         Eigen::Map<LatticeBasis<int>>(&box[0], nRows, nCols) = *cmapLapParaTaskSieve.basis;
         auto req2 = std::make_shared<MPI_Request>();
         PARA_COMM_CALL(
               commMpi->iSend(&box[0], nRows*nCols, UG::ParaINT, destination, TagTask2, req2.get())
               );
         auto iSendReq2 = std::make_shared<CMapLapParaIsendRequest>(req2, shared_from_this());
         commMpi->pushISendRequest(iSendReq2);
         break;
      }
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaSolverStateMpi::send, invalid solver type = ", static_cast<int>(solverType) );
      }
   }
   else
   {
      CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
      CMapLapParaTaskMpi *sent = createThreadDatatype(localComm);
      PARA_COMM_CALL(
            localComm->uTypeSend((void *)sent, ParaCMapLAP::ParaTaskType, destination, UG::TagTask)
            );
   }

   return 0;
}

int
CMapLapParaTaskMpi::receive(
      UG::ParaComm *comm,
      int source
      )
{
   UG::ParaCommMpi *commMpi = dynamic_cast< UG::ParaCommMpi* >(comm);

   if( commMpi )
   {
      MPI_CALL(
            commMpi->receive(&iSolverType, 1, UG::ParaINT, source, UG::TagTask)
            );

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
      {
         THROW_LOGICAL_ERROR2("CMapLapParaTaskMpi::receive: Invalid iSolver type = ", iSolverType);
      }
      }

      int tempTag;
      commMpi->waitSpecTagFromSpecSource(source, TagTask1, &tempTag );

      switch(solverType)
      {
      case DeepBkz:
      {
         solverType = DeepBkz;
         MPI_Datatype datatypeDeepBkz = createDatatypeDeepBkz();
         MPI_CALL(
            MPI_Type_commit( &datatypeDeepBkz )
         );
         PARA_COMM_CALL(
            commMpi->ureceive(&taskId.subtaskId.lcId, 1, datatypeDeepBkz, source, TagTask1)
         );
         MPI_CALL(
            MPI_Type_free( &datatypeDeepBkz )
         );

         int tempTag2;
         commMpi->waitSpecTagFromSpecSource(source, TagTask2, &tempTag2 );

         box.resize(nRows*nCols);
         PARA_COMM_CALL(
            commMpi->receive(&box[0], nRows*nCols, UG::ParaINT, source, TagTask2)
         );
         cmapLapParaTaskDeepBkz.basis = std::make_shared<LatticeBasis<int>>(
               Eigen::Map<LatticeBasis<int>>(&box[0], nRows, nCols)
               );
         break;
      }
      case Enum:
      {
         solverType = Enum;
         MPI_Datatype datatypeEnum = createDatatypeEnum();
         MPI_CALL(
            MPI_Type_commit( &datatypeEnum )
         );
         PARA_COMM_CALL(
            commMpi->ureceive(&taskId.subtaskId.lcId, 1, datatypeEnum, source, TagTask1)
         );
         MPI_CALL(
            MPI_Type_free( &datatypeEnum )
         );

         int tempTag2;
         commMpi->waitSpecTagFromSpecSource(source, TagTask2, &tempTag2 );

         vec.resize(nCols);
         PARA_COMM_CALL(
            commMpi->receive(&vec[0], nCols, UG::ParaINT, source, TagTask2)
         );
         cmapLapParaTaskEnum.coeffs = std::make_shared<LatticeVector<int>>(
               Eigen::Map<LatticeVector<int>>(&vec[0], nCols)
               );

         int tempTag3;
         commMpi->waitSpecTagFromSpecSource(source, TagTask3, &tempTag3 );

         box.resize(nRows*nCols);
         PARA_COMM_CALL(
            commMpi->receive(&box[0], nRows*nCols, UG::ParaINT, source, TagTask3)
         );
         cmapLapParaTaskEnum.basis = std::make_shared<LatticeBasis<int>>(
               Eigen::Map<LatticeBasis<int>>(&box[0], nRows, nCols)
               );
         break;
      }
      case Sieve:
      {
         solverType = Sieve;
         MPI_Datatype datatypeSieve = createDatatypeSieve();
         MPI_CALL(
            MPI_Type_commit( &datatypeSieve )
         );
         PARA_COMM_CALL(
            commMpi->ureceive(&taskId.subtaskId.lcId, 1, datatypeSieve, source, TagTask1)
         );
         MPI_CALL(
            MPI_Type_free( &datatypeSieve )
         );

         int tempTag2;
         commMpi->waitSpecTagFromSpecSource(source, TagTask2, &tempTag2 );

         box.resize(nRows*nCols);
         PARA_COMM_CALL(
            commMpi->receive(&box[0], nRows*nCols, UG::ParaINT, source, TagTask2)
         );
         cmapLapParaTaskSieve.basis = std::make_shared<LatticeBasis<int>>(
               Eigen::Map<LatticeBasis<int>>(&box[0], nRows, nCols)
               );
         break;
      }
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaTaskMpi::receive: Invalid iSolver type = ", static_cast<int>(iSolverType));
      }
   }
   else
   {
      std::unique_ptr<CMapLapParaTaskMpi> received;
      CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
      PARA_COMM_CALL(
            localComm->uTypeReceive((void **)&received, ParaCMapLAP::ParaTaskType, source, UG::TagTask)
      );
      taskId = received->taskId;
      generatorTaskId = received->generatorTaskId;
      estimatedValue = received->estimatedValue;
      threadId = received->threadId;
      solverType = received->solverType;
      switch( solverType )
      {
      case DeepBkz:
      {
         cmapLapParaTaskDeepBkz.begin      = received->cmapLapParaTaskDeepBkz.begin;
         cmapLapParaTaskDeepBkz.end        = received->cmapLapParaTaskDeepBkz.end;
         cmapLapParaTaskDeepBkz.blocksize  = received->cmapLapParaTaskDeepBkz.blocksize;
         cmapLapParaTaskDeepBkz.u          = received->cmapLapParaTaskDeepBkz.u;
         cmapLapParaTaskDeepBkz.seed       = received->cmapLapParaTaskDeepBkz.seed;
         cmapLapParaTaskDeepBkz.basis      = std::make_shared<LatticeBasis<int>>(*(received->cmapLapParaTaskDeepBkz.basis));
         break;
      }
      case Enum:
      {
         cmapLapParaTaskEnum.begin   = received->cmapLapParaTaskEnum.begin;
         cmapLapParaTaskEnum.end     = received->cmapLapParaTaskEnum.end;
         cmapLapParaTaskEnum.start   = received->cmapLapParaTaskEnum.start;
         cmapLapParaTaskEnum.last    = received->cmapLapParaTaskEnum.last;
         cmapLapParaTaskEnum.prob    = received->cmapLapParaTaskEnum.prob;
         cmapLapParaTaskEnum.coeffs  = std::make_shared<LatticeVector<int>>(*(received->cmapLapParaTaskEnum.coeffs));
         cmapLapParaTaskEnum.basis   = std::make_shared<LatticeBasis<int>>(*(received->cmapLapParaTaskEnum.basis));
         break;
      }
      case Sieve:
      {
         cmapLapParaTaskSieve.begin  = received->cmapLapParaTaskSieve.begin;
         cmapLapParaTaskSieve.end    = received->cmapLapParaTaskSieve.end;
         cmapLapParaTaskSieve.basis  = std::make_shared<LatticeBasis<int>>(*(received->cmapLapParaTaskSieve.basis));
         break;
      }
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaTaskMpi::receive: Invalid solver type = ", static_cast<int>(solverType));
      }
   }

   return 0;
}

} // namespace ParaCMapLAP
