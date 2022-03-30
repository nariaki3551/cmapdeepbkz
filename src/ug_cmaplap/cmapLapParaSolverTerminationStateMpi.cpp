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

/**@file    cmapLapParaSolverTerminationStateMpi.cpp
 * @brief   CMapLapParaSolverTerminationState extension for MIP communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "cmapLapParaSolverTerminationStateMpi.h"
#include "ug/paraComm.h"
#include "cmapLapParaComm.h"


namespace ParaCMapLAP
{

MPI_Datatype
CMapLapParaSolverTerminationStateMpi::createDatatype(){

   const int nBlocks = 32;

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
      MPI_Get_address( &interrupted, &startAddress )
   );
   displacements[0] = 0;

   MPI_CALL(
      MPI_Get_address( &rank, &address )
   );
   displacements[1] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &nParaTasksReceived, &address )
   );
   displacements[2] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &nParaTasksSolved, &address )
   );
   displacements[3] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &runningTime, &address )
   );
   displacements[4] = address - startAddress;
   types[4] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &idleTimeToFirstParaTask, &address )
   );
   displacements[5] = address - startAddress;
   types[5] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &idleTimeBetweenParaTasks, &address )
   );
   displacements[6] = address - startAddress;
   types[6] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &idleTimeAfterLastParaTask, &address )
   );
   displacements[7] = address - startAddress;
   types[7] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &idleTimeToWaitNotificationId, &address )
   );
   displacements[8] = address - startAddress;
   types[8] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &idleTimeToWaitAckCompletion, &address )
   );
   displacements[9] = address - startAddress;
   types[9] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &idleTimeToWaitToken, &address )
   );
   displacements[10] = address - startAddress;
   types[10] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &idleTimeToWaitSolverState, &address )
   );
   displacements[11] = address - startAddress;
   types[11] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &idleTimeToWaitPackedVector, &address )
   );
   displacements[12] = address - startAddress;
   types[12] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &idleTimeToWaitSolution, &address )
   );
   displacements[13] = address - startAddress;
   types[13] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &idleTimeToWaitBasis, &address )
   );
   displacements[14] = address - startAddress;
   types[14] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &idleTimeToWaitIsend, &address )
   );
   displacements[15] = address - startAddress;
   types[15] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &detTime, &address )
   );
   displacements[16] = address - startAddress;
   types[16] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &threadId, &address )
   );
   displacements[17] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &nParaTasksDeepBkzReceived, &address )
   );
   displacements[18] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &nParaTasksEnumReceived, &address )
   );
   displacements[19] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &nParaTasksSieveReceived, &address )
   );
   displacements[20] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &runningTimeDeepBkz, &address )
   );
   displacements[21] = address - startAddress;
   types[21] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &runningTimeEnum, &address )
   );
   displacements[22] = address - startAddress;
   types[22] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &runningTimeSieve, &address )
   );
   displacements[23] = address - startAddress;
   types[23] = MPI_DOUBLE;

   MPI_CALL(
      MPI_Get_address( &nVectorsReceivedDeepBkz, &address )
   );
   displacements[24] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &nVectorsReceivedEnum, &address )
   );
   displacements[25] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &nVectorsReceivedSieve, &address )
   );
   displacements[26] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &nVectorsSentDeepBkz, &address )
   );
   displacements[27] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &nVectorsSentEnum, &address )
   );
   displacements[28] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &nVectorsSentSieve, &address )
   );
   displacements[29] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &nBasesSentDeepBkz, &address )
   );
   displacements[30] = address - startAddress;

   MPI_CALL(
      MPI_Get_address( &nSolverStateSent, &address )
   );
   displacements[31] = address - startAddress;

   MPI_CALL(
         MPI_Type_create_struct(nBlocks, blockLengths, displacements, types, &datatype)
         );

   return datatype;

}

void
CMapLapParaSolverTerminationStateMpi::send(
      UG::ParaComm *comm,
      int destination,
      int tag
      )
{
   UG::ParaCommMpi *commMpi = dynamic_cast< UG::ParaCommMpi* >(comm);

   MPI_Datatype datatype;
   datatype = createDatatype();
   MPI_CALL(
      MPI_Type_commit( &datatype )
   );
   PARA_COMM_CALL(
      commMpi->usend(&interrupted, 1, datatype, destination, tag)
   );
   MPI_CALL(
      MPI_Type_free( &datatype )
   );
}

void
CMapLapParaSolverTerminationStateMpi::receive(
      UG::ParaComm *comm,
      int source,
      int tag
      )
{
   UG::ParaCommMpi *commMpi = dynamic_cast< UG::ParaCommMpi* >(comm);

   MPI_Datatype datatype;
   datatype = createDatatype();
   MPI_CALL(
      MPI_Type_commit( &datatype )
   );
   PARA_COMM_CALL(
      commMpi->ureceive(&interrupted, 1, datatype, source, tag)
   );
   MPI_CALL(
      MPI_Type_free( &datatype )
   );
}

} // namespace ParaCMapLAP
