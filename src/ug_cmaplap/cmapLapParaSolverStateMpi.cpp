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

/**@file    cmapLapParaSolverStateMpi.cpp
 * @brief   CMapLapParaSolverState extension for MPI communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "cmapLapParaSolverStateMpi.h"
#include "ug/paraComm.h"
#include "ug/paraDef.h"
#include "ug/paraTagDef.h"
#include "cmapLapParaCommMpi.h"
#include "cmapLapParaIsendRequest.h"
#include "cmapLapParaTagDef.h"
#include "cmapLapParaSolverState.h"


namespace ParaCMapLAP
{

MPI_Datatype
CMapLapParaSolverStateMpi::createDatatypeDeepBkz(
      )
{

   const int nBlocks = 20;

   MPI_Datatype datatype;

   MPI_Aint startAddress = 0;
   MPI_Aint address = 0;

   int blockLengths[nBlocks];
   MPI_Aint displacements[nBlocks];
   MPI_Datatype types[nBlocks];

   for( int i = 0; i < nBlocks; i++ )
   {
      blockLengths[i] = 1;
      types[i] = MPI_INT;
   }

   MPI_CALL(
      MPI_Get_address( &racingStage, &startAddress )
   );
   displacements[0] = 0;

   MPI_CALL(
      MPI_Get_address( &notificationId, &address )
   );
   displacements[1] = address - startAddress;
   types[1] = MPI_UNSIGNED;

   MPI_CALL(
      MPI_Get_address( &lcId, &address )
   );
   displacements[2] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &globalSubtreeIdInLc, &address )
   );
   displacements[3] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &detTime, &address )
   );
   displacements[4] = address - startAddress;
   types[4] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &threadId, &address )
   );
   displacements[5] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &meanMessageQueueSize, &address )
   );
   displacements[6] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &maxMessageQueueSize, &address )
   );
   displacements[7] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( cmapLapParaSolverStateData.cmapLapParaSolverStateDeepBkz.basis, &address )
   );
   displacements[8] = address - startAddress;
   blockLengths[8] = dimension * basisRows;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateDeepBkz.currentBlockSize, &address )
   );
   displacements[9] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateDeepBkz.tour, &address )
   );
   displacements[10] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateDeepBkz.elapsedTime, &address )
   );
   displacements[11] = address - startAddress;
   types[11] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateDeepBkz.shortestNorm, &address )
   );
   displacements[12] = address - startAddress;
   types[12] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateDeepBkz.approxFactor, &address )
   );
   displacements[13] = address - startAddress;
   types[13] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateDeepBkz.hermiteFactor, &address )
   );
   displacements[14] = address - startAddress;
   types[14] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateDeepBkz.rootHermiteFactor, &address )
   );
   displacements[15] = address - startAddress;
   types[15] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateDeepBkz.enumCost, &address )
   );
   displacements[16] = address - startAddress;
   types[16] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateDeepBkz.slopeGSA, &address )
   );
   displacements[17] = address - startAddress;
   types[17] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateDeepBkz.topHalfSlopeGSA, &address )
   );
   displacements[18] = address - startAddress;
   types[18] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateDeepBkz.orthogonalFactor, &address )
   );
   displacements[19] = address - startAddress;
   types[19] = MPI_DOUBLE;

   MPI_CALL(
         MPI_Type_create_struct(nBlocks, blockLengths, displacements, types, &datatype)
         );

   return datatype;
}

MPI_Datatype
CMapLapParaSolverStateMpi::createDatatypeEnum(
      )
{
   const int nBlocks = 17;

   MPI_Datatype datatype;

   MPI_Aint startAddress = 0;
   MPI_Aint address = 0;

   int blockLengths[nBlocks];
   MPI_Aint displacements[nBlocks];
   MPI_Datatype types[nBlocks];

   for( int i = 0; i < nBlocks; i++ )
   {
      blockLengths[i] = 1;
      types[i] = MPI_INT;
   }

   MPI_CALL(
      MPI_Get_address( &racingStage, &startAddress )
   );
   displacements[0] = 0;

   MPI_CALL(
      MPI_Get_address( &notificationId, &address )
   );
   displacements[1] = address - startAddress;
   types[1] = MPI_UNSIGNED;

   MPI_CALL(
      MPI_Get_address( &lcId, &address )
   );
   displacements[2] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &globalSubtreeIdInLc, &address )
   );
   displacements[3] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &detTime, &address )
   );
   displacements[4] = address - startAddress;
   types[4] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &threadId, &address )
   );
   displacements[5] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &meanMessageQueueSize, &address )
   );
   displacements[6] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &maxMessageQueueSize, &address )
   );
   displacements[7] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.elapsedTime, &address )
   );
   displacements[8] = address - startAddress;
   types[8] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.coeffs, &address )
   );
   displacements[9] = address - startAddress;
   blockLengths[9] = dimension;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.searchNodeIndex, &address )
   );
   displacements[10] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.enumCost, &address )
   );
   displacements[11] = address - startAddress;
   types[11] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.shortestNorm, &address )
   );
   displacements[12] = address - startAddress;
   types[12] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.approxFactor, &address )
   );
   displacements[13] = address - startAddress;
   types[13] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.hermiteFactor, &address )
   );
   displacements[14] = address - startAddress;
   types[14] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.rootHermiteFactor, &address )
   );
   displacements[15] = address - startAddress;
   types[15] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.numSearchedNodes, &address )
   );
   displacements[16] = address - startAddress;
   types[16] = MPI_LONG_INT;

   MPI_CALL(
         MPI_Type_create_struct(nBlocks, blockLengths, displacements, types, &datatype)
         );

   return datatype;
}

MPI_Datatype
CMapLapParaSolverStateMpi::createDatatypeSieve(
      )
{

   const int nBlocks = 19;

   MPI_Datatype datatype;

   MPI_Aint startAddress = 0;
   MPI_Aint address = 0;

   int blockLengths[nBlocks];
   MPI_Aint displacements[nBlocks];
   MPI_Datatype types[nBlocks];

   for( int i = 0; i < nBlocks; i++ )
   {
      blockLengths[i] = 1;
      types[i] = MPI_INT;
   }

   MPI_CALL(
      MPI_Get_address( &racingStage, &startAddress )
   );
   displacements[0] = 0;

   MPI_CALL(
      MPI_Get_address( &notificationId, &address )
   );
   displacements[1] = address - startAddress;
   types[1] = MPI_UNSIGNED;

   MPI_CALL(
      MPI_Get_address( &lcId, &address )
   );
   displacements[2] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &globalSubtreeIdInLc, &address )
   );
   displacements[3] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &detTime, &address )
   );
   displacements[4] = address - startAddress;
   types[4] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &threadId, &address )
   );
   displacements[5] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &meanMessageQueueSize, &address )
   );
   displacements[6] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &maxMessageQueueSize, &address )
   );
   displacements[7] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.elapsedTime, &address )
   );
   displacements[8] = address - startAddress;
   types[8] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.blockSize, &address )
   );
   displacements[9] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.nLoop, &address )
   );
   displacements[10] = address - startAddress;
   types[10] = MPI_LONG_INT;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.listSize, &address )
   );
   displacements[11] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.stackSize, &address )
   );
   displacements[12] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.maxListSize, &address )
   );
   displacements[13] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.nCollisions, &address )
   );
   displacements[14] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.shortestNorm, &address )
   );
   displacements[15] = address - startAddress;
   types[15] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.approxFactor, &address )
   );
   displacements[16] = address - startAddress;
   types[16] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.hermiteFactor, &address )
   );
   displacements[17] = address - startAddress;
   types[17] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaSolverStateData.cmapLapParaSolverStateSieve.rootHermiteFactor, &address )
   );
   displacements[18] = address - startAddress;
   types[18] = MPI_DOUBLE;

   MPI_CALL(
         MPI_Type_create_struct(nBlocks, blockLengths, displacements, types, &datatype)
         );

   return datatype;
}

void
CMapLapParaSolverStateMpi::send(
      UG::ParaComm *comm,
      int destination,
      int tag
      )
{
   CMapLapParaCommMpi *commMpi = dynamic_cast< CMapLapParaCommMpi* >(comm);

   typeAndDim[0] = iSolverType;
   typeAndDim[1] = dimension;
   typeAndDim[2] = basisRows;


   auto req = std::make_shared<MPI_Request>();
   MPI_CALL(
         commMpi->iSend(typeAndDim, 3, UG::ParaINT, destination, UG::TagSolverState, req.get())
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
         commMpi->iUsend(&racingStage, 1, datatypeDeepBkz, destination, TagSolverState1, req1.get())
      );
      auto iSendReq1 = std::make_shared<CMapLapParaIsendRequest>(req1, shared_from_this());
      commMpi->pushISendRequest(iSendReq1);
      MPI_CALL(
         MPI_Type_free( &datatypeDeepBkz )
      );
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
         commMpi->iUsend(&racingStage, 1, datatypeEnum, destination, TagSolverState1, req1.get())
      );
      auto iSendReq1 = std::make_shared<CMapLapParaIsendRequest>(req1, shared_from_this());
      commMpi->pushISendRequest(iSendReq1);
      MPI_CALL(
         MPI_Type_free( &datatypeEnum )
      );
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
         commMpi->iUsend(&racingStage, 1, datatypeSieve, destination, TagSolverState1, req1.get())
      );
      auto iSendReq1 = std::make_shared<CMapLapParaIsendRequest>(req1, shared_from_this());
      commMpi->pushISendRequest(iSendReq1);
      MPI_CALL(
         MPI_Type_free( &datatypeSieve )
      );
      break;
   }
   default:
      THROW_LOGICAL_ERROR2("CMapLapParaSolverStateMpi::send, invalid solver type = ", static_cast<int>(solverType) );
   }
}

void
CMapLapParaSolverStateMpi::receive(
      UG::ParaComm *comm,
      int source,
      int tag
      )
{
   UG::ParaCommMpi *commMpi = dynamic_cast< UG::ParaCommMpi* >(comm);

   /*
   std::chrono::system_clock::time_point startTime, endTime;
   startTime = std::chrono::system_clock::now();
   */

   MPI_CALL(
         commMpi->receive(typeAndDim, 3, UG::ParaINT, source, UG::TagSolverState)
         );

   iSolverType    = typeAndDim[0];
   dimension      = typeAndDim[1];
   basisRows      = typeAndDim[2];

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
      THROW_LOGICAL_ERROR2("CMapLapParaSolverStateMpi::receive, invalid iSolver type = ", iSolverType );
   }

   int tempTag;
   commMpi->waitSpecTagFromSpecSource(source, TagSolverState1, &tempTag );

   switch( solverType )
   {
   case DeepBkz:
   {
      if( basisRows > 0)
      {
         cmapLapParaSolverStateData.cmapLapParaSolverStateDeepBkz.basis = new int [dimension * basisRows];
      }
      else
      {
         cmapLapParaSolverStateData.cmapLapParaSolverStateDeepBkz.basis = nullptr;
      }

      MPI_Datatype datatypeDeepBkz = createDatatypeDeepBkz();
      MPI_CALL(
         MPI_Type_commit( &datatypeDeepBkz )
      );
      PARA_COMM_CALL(
         commMpi->ureceive(&racingStage, 1, datatypeDeepBkz, source, TagSolverState1)
      );
      MPI_CALL(
         MPI_Type_free( &datatypeDeepBkz )
      );
      break;
   }
   case Enum:
   {
      cmapLapParaSolverStateData.cmapLapParaSolverStateEnum.coeffs = new int [dimension];

      MPI_Datatype datatypeEnum = createDatatypeEnum();
      MPI_CALL(
         MPI_Type_commit( &datatypeEnum )
      );
      PARA_COMM_CALL(
         commMpi->ureceive(&racingStage, 1, datatypeEnum, source, TagSolverState1)
      );
      MPI_CALL(
         MPI_Type_free( &datatypeEnum )
      );
      break;
   }
   case Sieve:
   {
      MPI_Datatype datatypeSieve = createDatatypeSieve();
      MPI_CALL(
         MPI_Type_commit( &datatypeSieve )
      );
      PARA_COMM_CALL(
         commMpi->ureceive(&racingStage, 1, datatypeSieve, source, TagSolverState1)
      );
      MPI_CALL(
         MPI_Type_free( &datatypeSieve )
      );
      break;
   }
   default:
      THROW_LOGICAL_ERROR2("CMapLapParaSolverStateMpi::receive, invalid iSolver type = ", iSolverType );
   }

   /*
   endTime = std::chrono::system_clock::now();
   std::cout << "CMapLapParaSolverStateMpi::receive from " << source << " time "
      << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(endTime-startTime).count())
         / 100000.0
      << std::endl;
   */
}

} // namespace ParaCMapLAP
