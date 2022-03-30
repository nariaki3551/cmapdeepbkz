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

/**@file    cmapLapParaInstanceMpi.cpp
 * @brief   CMapLapParaInstance extension for MPI communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cstring>
#include <assert.h>
#include "cmapLapParaTagDef.h"
#include "cmapLapParaSolverLocalComm.h"
#include "cmapLapParaInstanceMpi.h"


namespace ParaCMapLAP
{

///
/// create CMapLapInstancePrePreDatatype
///
MPI_Datatype
CMapLapParaInstanceMpi::createDatatype1(
      )
{
   const int nBlocks = 1;
   MPI_Datatype prePreDatatype;

   MPI_Aint startAddress = 0;
   // MPI_Aint address = 0;

   int blockLengths[nBlocks];
   MPI_Aint displacements[nBlocks];
   MPI_Datatype types[nBlocks];

   MPI_CALL(
      MPI_Get_address( &lProbName, &startAddress )
   );
   displacements[0] = 0;
   blockLengths[0] = 1;
   types[0] = MPI_INT;

   MPI_CALL(
      MPI_Type_create_struct(nBlocks, blockLengths, displacements, types, &prePreDatatype)
         );

   return prePreDatatype;
}

void
CMapLapParaInstanceMpi::allocateMemoryForDatatype2(
      )
{
   assert( lProbName > 0 );
   probName = new char[lProbName+1];
}

///
/// create CMapLapInstancePreDatatype
///
MPI_Datatype
CMapLapParaInstanceMpi::createDatatype2(
      bool memAllocNecessary
      )
{
   const int nBlocks = 1;
   MPI_Datatype datatype;

   MPI_Aint startAddress = 0;
   // MPI_Aint address = 0;

   int blockLengths[nBlocks];
   MPI_Aint displacements[nBlocks];
   MPI_Datatype types[nBlocks];

   if( memAllocNecessary )
   {
      allocateMemoryForDatatype2();
   }

   MPI_CALL(
      MPI_Get_address( probName, &startAddress )
   );
   displacements[0] = 0;
   blockLengths[0] = lProbName + 1;
   types[0] = MPI_CHAR;

   MPI_CALL(
         MPI_Type_create_struct(nBlocks, blockLengths, displacements, types, &datatype)
         );

   return datatype;
}

int
CMapLapParaInstanceMpi::bcast(
      UG::ParaComm *comm,
      int root,
      int method
      )
{
   UG::ParaCommMpi *commMpi = dynamic_cast< UG::ParaCommMpi* >(comm);

   if( commMpi )
   {
      MPI_Datatype datatype = createDatatype1();
      MPI_CALL(
         MPI_Type_commit( &datatype )
      );
      PARA_COMM_CALL(
         commMpi->ubcast(&lProbName, 1, datatype, root)
      );
      MPI_CALL(
         MPI_Type_free( &datatype )
      );

      if( commMpi->getRank() == root )
      {
         datatype = createDatatype2(false);
      }
      else
      {
         datatype = createDatatype2(true);
      }
      MPI_CALL(
         MPI_Type_commit( &datatype )
      );

      PARA_COMM_CALL(
         commMpi->ubcast(probName, 1, datatype, root)
      );
      MPI_CALL(
         MPI_Type_free( &datatype )
      );
   }
   else
   {
      CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
      if( localComm->getThreadId() == root )
      {
         for( int i = 0; i < localComm->getSize(); i++ )
         {
            if( i != root )
            {
               PARA_COMM_CALL(
                     localComm->uTypeSend((void *)createThreadDatatype(comm), ParaCMapLAP::ParaInstanceType, i, ParaCMapLAP::TagParaInstance)
               );
            }
         }
      }
      else
      {
         std::unique_ptr<CMapLapParaInstanceMpi> received;
         PARA_COMM_CALL(
               localComm->uTypeReceive((void **)&received, ParaCMapLAP::ParaInstanceType, root, ParaCMapLAP::TagParaInstance)
         );
         lProbName = received->lProbName;
         probName = new char[lProbName + 1];
         strcpy(probName, received->probName);
      }
   }

   return 0;

}

} // namespace ParaCMapLAP