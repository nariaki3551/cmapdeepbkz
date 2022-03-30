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

/**@file    cmapLapParaSolutionMpi.cpp
 * @brief   CMapLapParaSolution extension for MPI communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <mpi.h>
#include "ug/paraTagDef.h"
#include "cmapLapParaCommMpi.h"
#include "cmapLapParaSolverLocalComm.h"
#include "cmapLapParaSolutionMpi.h"
#include "cmapLapParaIsendRequest.h"
#include "cmapLapParaTagDef.h"


namespace ParaCMapLAP
{

///
/// create clone of this object
///
UG::ParaSolution *
CMapLapParaSolutionMpi::clone(
      UG::ParaComm *comm
      )
{
   LatticeVector<int> pv(v);
   return( new CMapLapParaSolutionMpi(threadId, pv, projectedSquaredNorm));
}

///
/// create CMapLapParaSolutioDatatype
///
MPI_Datatype
CMapLapParaSolutionMpi::createDatatype(
      )
{

   const int nBlocks = 3;

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
      MPI_Get_address( &projectedSquaredNorm, &address )
   );
   displacements[1] = address - startAddress;
   types[1] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &lenv, &address )
   );
   displacements[2] = address - startAddress;

   MPI_CALL(
      MPI_Type_create_struct(nBlocks, blockLengths, displacements, types, &datatype)
   );

   return datatype;
}

///
/// send solution data to the rank
///
void
CMapLapParaSolutionMpi::bcast(
      UG::ParaComm *comm,
      int root
      )
{
   UG::ParaCommMpi *commMpi = dynamic_cast< UG::ParaCommMpi* >(comm);

   if( commMpi )
   {
      if( comm->getRank() == root )
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

         vec.resize(lenv);
         Eigen::Map<LatticeVector<int>>(&vec[0], lenv) = v;
         PARA_COMM_CALL(
            commMpi->bcast(&vec[0], lenv, UG::ParaINT, root)
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

         vec.resize(lenv);
         PARA_COMM_CALL(
            commMpi->bcast(&vec[0], lenv, UG::ParaINT, root)
         );
         v = Eigen::Map<LatticeVector<int>>(&vec[0], lenv);
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
                     localComm->uTypeSend((void *)createThreadDatatype(comm), ParaCMapLAP::ParaSolutionType, i, UG::TagSolution)
               );
            }
         }
      }
      else
      {
         std::unique_ptr<CMapLapParaSolutionMpi> received;
         PARA_COMM_CALL(
               localComm->uTypeReceive((void **)&received, ParaCMapLAP::ParaSolutionType, root, UG::TagSolution)
         );

         threadId = received->threadId;
         projectedSquaredNorm = received->projectedSquaredNorm;
         v = LatticeVector<int>(received->v);
      }
   }
}


///
/// send solution data to the rank
///
void
CMapLapParaSolutionMpi::send(
      UG::ParaComm *comm,
      int destination
      )
{
   CMapLapParaCommMpi *commMpi = dynamic_cast< CMapLapParaCommMpi* >(comm);

   if( commMpi )
   {
      MPI_Datatype datatype = createDatatype();
      MPI_CALL(
         MPI_Type_commit( &datatype )
      );
      auto req = std::make_shared<MPI_Request>();
      PARA_COMM_CALL(
         commMpi->iUsend(&threadId, 1, datatype, destination, UG::TagSolution, req.get())
      );
      auto iSendReq = std::make_shared<CMapLapParaIsendRequest>(req, shared_from_this());
      commMpi->pushISendRequest(iSendReq);
      MPI_CALL(
         MPI_Type_free( &datatype )
      );

      vec.resize(lenv);
      Eigen::Map<LatticeVector<int>>(&vec[0], lenv) = v;
      auto req1 = std::make_shared<MPI_Request>();
      PARA_COMM_CALL(
         commMpi->iSend(&vec[0], lenv, UG::ParaINT, destination, TagSolution1, req1.get())
      );
      auto iSendReq1 = std::make_shared<CMapLapParaIsendRequest>(req1, shared_from_this());
      commMpi->pushISendRequest(iSendReq1);
   }
   else
   {
      CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
      PARA_COMM_CALL(
         localComm->uTypeSend((void *)createThreadDatatype(comm), ParaCMapLAP::ParaSolutionType, destination, UG::TagSolution)
      );
   }

}


///
/// receive solution data from the source rank
///
void
CMapLapParaSolutionMpi::receive(
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
            commMpi->ureceive(&threadId, 1, datatype, source, UG::TagSolution)
      );
      MPI_CALL(
         MPI_Type_free( &datatype )
      );

      int tempTag;
      commMpi->waitSpecTagFromSpecSource(source, TagSolution1, &tempTag );

      vec.resize(lenv);
      PARA_COMM_CALL(
         commMpi->receive(&vec[0], lenv, UG::ParaINT, source, TagSolution1)
      );
      v = Eigen::Map<LatticeVector<int>>(&vec[0], lenv);
   }
   else
   {
      std::unique_ptr<CMapLapParaSolutionMpi> received;
      CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
      PARA_COMM_CALL(
            localComm->uTypeReceive((void **)&received, ParaCMapLAP::ParaSolutionType, source, UG::TagSolution)
      );

      threadId = received->threadId;
      projectedSquaredNorm = received->projectedSquaredNorm;
      v = LatticeVector<int>(received->v);
   }
}

} // namespace ParaCMapLAP
