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

/**@file    cmapLapParaPackedVectorMpi.cpp
 * @brief   CMapLapParaPackedVector extension for MPI communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "cmapLapParaPackedVectorMpi.h"
#include <assert.h>
#include <memory>
#include <mpi.h>
#include "ug/paraComm.h"
#include "cmapLapParaTagDef.h"
#include "cmapLapParaCommMpi.h"
#include "cmapLapParaSolverLocalComm.h"
#include "cmapLapParaIsendRequest.h"
#include "cmapLapParaShareData.h"


namespace ParaCMapLAP
{

///
/// @brief create clone of this object
/// @note copy pointer of vector element object because vector element will not be change
///
CMapLapParaPackedVectorMpi *
CMapLapParaPackedVectorMpi::clone(
      UG::ParaComm *comm
      )
{
   std::vector<std::shared_ptr<VectorElement>> inVectorElements(nVectors);
   for (int i = 0; i < nVectors; i++ )
   {
      inVectorElements[i] = vectorElements[i];
   }
   return( new CMapLapParaPackedVectorMpi(threadId, inVectorElements) );
}

///
/// create CMapLapParaSolutioDatatype
///
MPI_Datatype
CMapLapParaPackedVectorMpi::createDatatype(
      int inNVectors,
      int inDimension
      )
{
   int nBlocks = 0;

   MPI_Datatype datatype;

   MPI_Aint startAddress = 0;
   MPI_Aint address = 0;

   int blockLengths[4];
   MPI_Aint displacements[4];
   MPI_Datatype types[4];

   MPI_CALL(
      MPI_Get_address( &nVectors, &startAddress )
   );
   displacements[nBlocks] = 0;
   blockLengths[nBlocks] = 1;
   types[nBlocks++] = MPI_INT;

   MPI_CALL(
      MPI_Get_address( &createdRanks[0], &address )
   );
   displacements[nBlocks] = address - startAddress;
   blockLengths[nBlocks] = inNVectors;
   types[nBlocks++] = MPI_INT;

   MPI_CALL(
      MPI_Get_address( &squaredNorms[0], &address )
   );
   displacements[nBlocks] = address - startAddress;;
   blockLengths[nBlocks] = inNVectors;
   types[nBlocks++] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &allVectorElements[0], &address )
   );
   displacements[nBlocks] = address - startAddress;;
   blockLengths[nBlocks] = inNVectors * inDimension;
   types[nBlocks++] = MPI_INT;

   assert(nBlocks == 4);

   MPI_CALL(
      MPI_Type_create_struct(nBlocks, blockLengths, displacements, types, &datatype)
   );

   return datatype;
}

///
/// send solution data to the rank
///
void
CMapLapParaPackedVectorMpi::send(
      UG::ParaComm *comm,
      int destination
      )
{
   CMapLapParaCommMpi *commMpi = dynamic_cast< CMapLapParaCommMpi* >(comm);

   if( commMpi )
   {
      if( vectorElements.size() > 0 )
      {
         dimension = vectorElements[0]->size();
      }
      else
      {
         assert( nVectors == 0 );
         dimension = 0;
      }

      len[0] = threadId;
      len[1] = nVectors;
      len[2] = dimension;
      auto req = std::make_shared<MPI_Request>();
      PARA_COMM_CALL(
            commMpi->iSend(len, 3, UG::ParaINT, destination, TagCMapLapPackedVector, req.get())
      );
      auto iSendReq = std::make_shared<CMapLapParaIsendRequest>(req, shared_from_this());
      commMpi->pushISendRequest(iSendReq);

      if( nVectors > 0 )
      {
         createdRanks.resize(nVectors);
         squaredNorms.resize(nVectors);
         allVectorElements.resize(nVectors * dimension);

         for(int i = 0; i < nVectors; i++ )
         {
            createdRanks[i] = vectorElements[i]->solverThreadId;
            squaredNorms[i] = vectorElements[i]->squaredNorm();
            for( int j = 0; j < dimension; j++ )
            {
               allVectorElements[i*dimension+j] = vectorElements[i]->vector[j];
            }
         }

         MPI_Datatype datatype = createDatatype(nVectors, dimension);
         MPI_CALL(
            MPI_Type_commit( &datatype )
         );
         auto req1 = std::make_shared<MPI_Request>();
         PARA_COMM_CALL(
               commMpi->iUsend(&nVectors, 1, datatype, destination, TagCMapLapPackedVector1, req1.get())
         );
         auto iSendReq1 = std::make_shared<CMapLapParaIsendRequest>(req1, shared_from_this());
         commMpi->pushISendRequest(iSendReq1);
         MPI_CALL(
            MPI_Type_free( &datatype )
         );

      }
   }
   else
   {
      CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
      PARA_COMM_CALL(
         localComm->uTypeSend((void *)createThreadDatatype(comm), ParaCMapLAP::CMapLapParaPackedVectorType, destination, TagCMapLapPackedVector)
      );
   }


}

///
/// receive solution data from the source rank
///
void
CMapLapParaPackedVectorMpi::receive(
      UG::ParaComm *comm,
      int source
      )
{
   UG::ParaCommMpi *commMpi = dynamic_cast< UG::ParaCommMpi* >(comm);

   if( commMpi )
   {

      PARA_COMM_CALL(
            commMpi->receive(len, 3, UG::ParaINT, source, TagCMapLapPackedVector)
      );

      threadId = len[0];
      nVectors = len[1];
      dimension = len[2];

      if( nVectors > 0 )
      {
         createdRanks.resize(nVectors);
         squaredNorms.resize(nVectors);
         allVectorElements.resize(nVectors * dimension);
         vectorElements.resize(nVectors);

         int tempTag;
         commMpi->waitSpecTagFromSpecSource(source, TagCMapLapPackedVector1, &tempTag);

         MPI_Datatype datatype = createDatatype(nVectors, dimension);
         MPI_CALL(
            MPI_Type_commit( &datatype )
         );
         PARA_COMM_CALL(
               commMpi->ureceive(&nVectors, 1, datatype, source, TagCMapLapPackedVector1)
         );
         MPI_CALL(
            MPI_Type_free( &datatype )
         );

         for(int i = 0; i < nVectors; i++ )
         {
            int rRank = createdRanks[i];
            LatticeVector<int> rVector(dimension);
            for( int j = 0; j < dimension; j++ )
            {
               rVector[j] = allVectorElements[i*dimension+j];
            }
            vectorElements[i] = std::shared_ptr<VectorElement>(new VectorElement(rVector, rRank));
         }
      }
   }
   else
   {
      std::unique_ptr<CMapLapParaPackedVectorMpi> received;
      CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
      PARA_COMM_CALL(
            localComm->uTypeReceive((void **)&received, ParaCMapLAP::CMapLapParaPackedVectorType, source, TagCMapLapPackedVector)
      );
      threadId = received->threadId;
      nVectors = received->nVectors;
      vectorElements.resize(nVectors);
      for (int i = 0; i < nVectors; i++ )
      {
         vectorElements[i] = received->getVectorElements()[i];
      }
   }
}

} // namespace ParaCMapLAP
