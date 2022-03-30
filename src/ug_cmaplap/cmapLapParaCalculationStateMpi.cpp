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

/**@file    cmapLapParaCalculationStateMpi.cpp
 * @brief   CalcutationState object extension for MPI communication
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "cmapLapParaCalculationStateMpi.h"
#include <memory>
#include "ug/paraComm.h"
#include "ug/paraDef.h"
#include "ug/paraTagDef.h"
#include "cmapLapParaCommMpi.h"
#include "cmapLapParaIsendRequest.h"
#include "cmapLapParaTagDef.h"
#include "cmapLapParaCalculationState.h"

namespace ParaCMapLAP
{

///
/// create MPI datatype of this object
/// @return
///
MPI_Datatype
CMapLapParaCalculationStateMpi::createDatatypeDeepBkz(
      )
{

   const int nBlocks = 15;

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
      MPI_Get_address( &compTime, &startAddress )
   );
   displacements[0] = 0;
   types[0] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &nSolved, &address )
   );
   displacements[1] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &terminationState, &address )
   );
   displacements[2] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &threadId, &address )
   );
   displacements[3] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.currentBlockSize, &address )
   );
   displacements[4] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.tour, &address )
   );
   displacements[5] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.shortestNorm, &address )
   );
   displacements[6] = address - startAddress;
   types[6] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.approxFactor, &address )
   );
   displacements[7] = address - startAddress;
   types[7] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.hermiteFactor, &address )
   );
   displacements[8] = address - startAddress;
   types[8] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.rootHermiteFactor, &address )
   );
   displacements[9] = address - startAddress;
   types[9] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.enumCost, &address )
   );
   displacements[10] = address - startAddress;
   types[10] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.enumCostGH, &address )
   );
   displacements[11] = address - startAddress;
   types[11] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.slopeGSA, &address )
   );
   displacements[12] = address - startAddress;
   types[12] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.topHalfSlopeGSA, &address )
   );
   displacements[13] = address - startAddress;
   types[13] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.orthogonalFactor, &address )
   );
   displacements[14] = address - startAddress;
   types[14] = MPI_DOUBLE;

   MPI_CALL(
         MPI_Type_create_struct(nBlocks, blockLengths, displacements, types, &datatype)
         );

   return datatype;

}

///
/// create MPI datatype of this object
/// @return
///
MPI_Datatype
CMapLapParaCalculationStateMpi::createDatatypeEnum(
      )
{

   const int nBlocks = 9;

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
      MPI_Get_address( &compTime, &startAddress )
   );
   displacements[0] = 0;
   types[0] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &nSolved, &address )
   );
   displacements[1] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &terminationState, &address )
   );
   displacements[2] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &threadId, &address )
   );
   displacements[3] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateEnum.shortestNorm, &address )
   );
   displacements[4] = address - startAddress;
   types[4] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateEnum.approxFactor, &address )
   );
   displacements[5] = address - startAddress;
   types[5] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateEnum.hermiteFactor, &address )
   );
   displacements[6] = address - startAddress;
   types[6] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateEnum.rootHermiteFactor, &address )
   );
   displacements[7] = address - startAddress;
   types[7] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateEnum.numSearchedNodes, &address )
   );
   displacements[8] = address - startAddress;
   types[8] = MPI_LONG_INT;

   MPI_CALL(
         MPI_Type_create_struct(nBlocks, blockLengths, displacements, types, &datatype)
         );

   return datatype;

}

///
/// create MPI datatype of this object
/// @return
///
MPI_Datatype
CMapLapParaCalculationStateMpi::createDatatypeSieve(
      )
{

   const int nBlocks = 14;

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
      MPI_Get_address( &compTime, &startAddress )
   );
   displacements[0] = 0;
   types[0] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &nSolved, &address )
   );
   displacements[1] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &terminationState, &address )
   );
   displacements[2] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &threadId, &address )
   );
   displacements[3] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.blockSize, &address )
   );
   displacements[4] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.nLoop, &address )
   );
   displacements[5] = address - startAddress;
   types[5] = MPI_LONG_INT;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.listSize, &address )
   );
   displacements[6] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.stackSize, &address )
   );
   displacements[7] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.maxListSize, &address )
   );
   displacements[8] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.nCollisions, &address )
   );
   displacements[9] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.shortestNorm, &address )
   );
   displacements[10] = address - startAddress;
   types[10] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.approxFactor, &address )
   );
   displacements[11] = address - startAddress;
   types[11] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.hermiteFactor, &address )
   );
   displacements[12] = address - startAddress;
   types[12] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.rootHermiteFactor, &address )
   );
   displacements[13] = address - startAddress;
   types[13] = MPI_DOUBLE;

   MPI_CALL(
         MPI_Type_create_struct(nBlocks, blockLengths, displacements, types, &datatype)
         );

   return datatype;

}

///
/// send this object to destination
///
void
CMapLapParaCalculationStateMpi::send(
      UG::ParaComm *comm,     ///< communicator used to send this object
      int destination,        ///< destination rank to send
      int tag                 ///< tag to show this object
      )
{
   CMapLapParaCommMpi *commMpi = dynamic_cast< CMapLapParaCommMpi* >(comm);

   auto req = std::make_shared<MPI_Request>();
   MPI_CALL(
         commMpi->iSend(&iSolverType, 1, UG::ParaINT, destination, UG::TagCompletionOfCalculation, req.get())
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
         commMpi->iUsend(&compTime, 1, datatypeDeepBkz, destination, TagCompletionOfCalculation1, req1.get())
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
         commMpi->iUsend(&compTime, 1, datatypeEnum, destination, TagCompletionOfCalculation1, req1.get())
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
         commMpi->iUsend(&compTime, 1, datatypeSieve, destination, TagCompletionOfCalculation1, req1.get())
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

///
/// receive this object from source
///
void
CMapLapParaCalculationStateMpi::receive(
      UG::ParaComm *comm,     ///< communicator used to receive this object
      int source,             ///< source rank to receive this object
      int tag                 ///< tag to show this object
      )
{
   UG::ParaCommMpi *commMpi = dynamic_cast< UG::ParaCommMpi* >(comm);

   MPI_CALL(
         commMpi->receive(&iSolverType, 1, UG::ParaINT, source, UG::TagCompletionOfCalculation)
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
      THROW_LOGICAL_ERROR2("CMapLapParaSolverStateMpi::receive, invalid solver type = ", iSolverType );
   }

   int tempTag;
   commMpi->waitSpecTagFromSpecSource(source, TagCompletionOfCalculation1, &tempTag );

   switch( solverType )
   {
      case DeepBkz:
   {
      MPI_Datatype datatypeDeepBkz = createDatatypeDeepBkz();
      MPI_CALL(
         MPI_Type_commit( &datatypeDeepBkz )
      );
      PARA_COMM_CALL(
         commMpi->ureceive(&compTime, 1, datatypeDeepBkz, source, TagCompletionOfCalculation1)
      );
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
      PARA_COMM_CALL(
         commMpi->ureceive(&compTime, 1, datatypeEnum, source, TagCompletionOfCalculation1)
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
         commMpi->ureceive(&compTime, 1, datatypeSieve, source, TagCompletionOfCalculation1)
      );
      MPI_CALL(
         MPI_Type_free( &datatypeSieve )
      );
      break;
   }
   default:
      THROW_LOGICAL_ERROR2("CMapLapParaSolverStateMpi::send, invalid solver type = ", static_cast<int>(solverType) );
   }
}

} // namespace ParaCMapLAP
