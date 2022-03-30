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

/**@file    cmapLapParaBasisMpi.cpp
 * @brief   CMapLapParaBasis extension for MPI communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "cmapLapParaBasisMpi.h"
#include <mpi.h>
#include <memory>
#include <assert.h>
#include "ug/paraComm.h"
#include "cmapLapParaTagDef.h"
#include "cmapLapParaCommMpi.h"
#include "cmapLapParaSolverLocalComm.h"
#include "cmapLapParaIsendRequest.h"


namespace ParaCMapLAP
{

///
/// create clone of this object
///
ParaCMapLAP::CMapLapParaBasisMpi*
CMapLapParaBasisMpi::clone(
      UG::ParaComm *comm
      )
{
   LatticeBasis<int> clonedBasis(basis);
   return( new CMapLapParaBasisMpi(threadId, enumCost, clonedBasis));
}

///
/// create CMapLapParaSolutioDatatype
///
MPI_Datatype
CMapLapParaBasisMpi::createDatatype(
      )
{

   const int nBlocks = 4;

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
      MPI_Get_address( &threadId, &startAddress )
   );
   displacements[0] = 0;

   MPI_CALL(
      MPI_Get_address( &enumCost, &address )
   );
   displacements[1] = address - startAddress;
   types[1] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &nRows, &address )
   );
   displacements[2] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &nCols, &address )
   );
   displacements[3] = address - startAddress;

   MPI_CALL(
      MPI_Type_create_struct(nBlocks, blockLengths, displacements, types, &datatype)
   );

   return datatype;
}

///
/// send solution data to the rank
///
void
CMapLapParaBasisMpi::bcast(
      UG::ParaComm *comm,
      int root
      )
{
   UG::ParaCommMpi *commMpi = dynamic_cast< UG::ParaCommMpi* >(comm);

   if( commMpi )
   {
      if( comm->getRank() == root )
      {
         nRows = basis.rows();
         nCols = basis.cols();

         MPI_Datatype datatype = createDatatype();
         MPI_CALL(
            MPI_Type_commit( &datatype )
         );
         PARA_COMM_CALL(
            commMpi->ubcast(&threadId, 1, datatype, root)
         );
         MPI_CALL(
            MPI_Type_free( &datatype )
         );

         assert( (nRows > 0) && (nCols > 0) );
         box.resize(nRows*nCols);
         Eigen::Map<LatticeBasis<int>>(&box[0], nRows, nCols) = basis;
         PARA_COMM_CALL(
            commMpi->bcast(&box[0], nRows*nCols, UG::ParaINT, root)
         );
      }
      else
      {
         MPI_Datatype datatype = createDatatype();
         MPI_CALL(
            MPI_Type_commit( &datatype )
         );
         PARA_COMM_CALL(
               commMpi->ubcast(&threadId, 1, datatype, root)
         );
         MPI_CALL(
            MPI_Type_free( &datatype )
         );

	      box.resize(nRows*nCols);
         PARA_COMM_CALL(
            commMpi->bcast(&box[0], nRows*nCols, UG::ParaINT, root)
         );
         basis = Eigen::Map<LatticeBasis<int>>(&box[0], nRows, nCols);
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
               PARA_COMM_CALL(
                     localComm->uTypeSend((void *)createThreadDatatype(comm), ParaCMapLAP::ParaBasisType, i, TagBasis)
               );
            }
         }
      }
      else
      {
         std::unique_ptr<CMapLapParaBasisMpi> received;
         PARA_COMM_CALL(
               localComm->uTypeReceive((void **)&received, ParaCMapLAP::ParaBasisType, root, TagBasis)
         );

         threadId = received->threadId;
         enumCost = received->enumCost;
         basis    = LatticeBasis<int>(received->basis);
      }
   }
}


///
/// send solution data to the rank
///
void
CMapLapParaBasisMpi::send(
      UG::ParaComm *comm,
      int destination
      )
{
   CMapLapParaCommMpi *commMpi = dynamic_cast< CMapLapParaCommMpi* >(comm);

   if( commMpi )
   {
      nRows = basis.rows();
      nCols = basis.cols();

      MPI_Datatype datatype = createDatatype();
      MPI_CALL(
         MPI_Type_commit( &datatype )
      );
      auto req = std::make_shared<MPI_Request>();
      PARA_COMM_CALL(
         commMpi->iUsend(&threadId, 1, datatype, destination, TagBasis, req.get())
      );
      auto iSendReq = std::make_shared<CMapLapParaIsendRequest>(req, shared_from_this());
      commMpi->pushISendRequest(iSendReq);
      MPI_CALL(
         MPI_Type_free( &datatype )
      );

      assert( (nRows > 0) && (nCols > 0) );
      box.resize(nRows*nCols);
      Eigen::Map<LatticeBasis<int>>(&box[0], nRows, nCols) = basis;
      auto req1 = std::make_shared<MPI_Request>();
      PARA_COMM_CALL(
         commMpi->iSend(&box[0], nRows*nCols, UG::ParaINT, destination, TagBasis1, req1.get())
      );
      auto iSendReq1 = std::make_shared<CMapLapParaIsendRequest>(req1, shared_from_this());
      commMpi->pushISendRequest(iSendReq1);
   }
   else
   {
      CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
      PARA_COMM_CALL(
         localComm->uTypeSend((void *)createThreadDatatype(comm), ParaCMapLAP::ParaBasisType, destination, TagBasis)
      );
   }

}


///
/// receive solution data from the source rank
///
void
CMapLapParaBasisMpi::receive(
      UG::ParaComm *comm,
      int source
      )
{
   UG::ParaCommMpi *commMpi = dynamic_cast< UG::ParaCommMpi* >(comm);

   if( commMpi )
   {
      MPI_Datatype datatype = createDatatype();
      MPI_CALL(
         MPI_Type_commit( &datatype )
      );
      PARA_COMM_CALL(
            commMpi->ureceive(&threadId, 1, datatype, source, TagBasis)
      );
      MPI_CALL(
         MPI_Type_free( &datatype )
      );

      assert( (nRows > 0) && (nCols > 0) );

      box.resize(nRows*nCols);
      int tempTag;
      commMpi->waitSpecTagFromSpecSource(source, TagBasis1, &tempTag);
      PARA_COMM_CALL(
         commMpi->receive(&box[0], nRows*nCols, UG::ParaINT, source, TagBasis1)
      );
      basis = Eigen::Map<LatticeBasis<int>>(&box[0], nRows, nCols);
   }
   else
   {
      std::unique_ptr<CMapLapParaBasisMpi> received;
      CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
      PARA_COMM_CALL(
            localComm->uTypeReceive((void **)&received, ParaCMapLAP::ParaBasisType, source, TagBasis)
      );

      threadId = received->threadId;
      enumCost = received->enumCost;
      basis    = LatticeBasis<int>(received->basis);
   }
}


} // namespace ParaCMapLAP
