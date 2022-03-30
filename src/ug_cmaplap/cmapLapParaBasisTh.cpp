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

/**@file    cmapLapParaBasisTh.cpp
 * @brief   CMapLapParaBasis extension for threads communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "cmapLapParaBasisTh.h"
#include <memory>
#include <new>
#include "ug/paraComm.h"
#include "cmapLapParaSolverLocalComm.h"
#include "cmapLapParaTagDef.h"
#include "cmapLapParaBasis.h"


namespace ParaCMapLAP
{

///
/// create clone of this object
///
ParaCMapLAP::CMapLapParaBasisTh*
CMapLapParaBasisTh::clone(
      UG::ParaComm *comm
      )
{
   LatticeBasis<int> clonedBasis(basis);
   return( new CMapLapParaBasisTh(threadId, enumCost, clonedBasis));
}

///
/// create CMapLap Basis Datatype
///
CMapLapParaBasisTh *
CMapLapParaBasisTh::createDatatype(
      UG::ParaComm *comm
      )
{
   return dynamic_cast<CMapLapParaBasisTh *>(clone(comm));
}

///
/// send solution data to the rank
///
void
CMapLapParaBasisTh::bcast(
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
               PARA_COMM_CALL(
                  commTh->uTypeSend((void *)createDatatype(comm), ParaBasisType, i, TagBasis)
               );
            }
         }
      }
      else
      {
         std::unique_ptr<CMapLapParaBasisTh> received;
         PARA_COMM_CALL(
            commTh->uTypeReceive((void **)&received, ParaBasisType, root, TagBasis)
         );
         threadId = received->threadId;
         enumCost = received->enumCost;
         basis    = LatticeBasis<int>(received->basis);
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
                     localComm->uTypeSend((void *)createDatatype(comm), ParaBasisType, i, TagBasis)
               );
            }
         }
      }
      else
      {
         std::unique_ptr<CMapLapParaBasisTh> received;
         PARA_COMM_CALL(
               localComm->uTypeReceive((void **)&received, ParaBasisType, root, TagBasis)
         );
         threadId = received->threadId;
         enumCost = received->enumCost;
         basis    = LatticeBasis<int>(received->basis);
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
            PARA_COMM_CALL(
                  localComm->uTypeSend((void *)createDatatype(comm), ParaCMapLAP::ParaBasisType, i, TagBasis)
            );
         }
      }
   }
   else
   {
      std::unique_ptr<CMapLapParaBasis> received;
      PARA_COMM_CALL(
            localComm->uTypeReceive((void **)&received, ParaCMapLAP::ParaBasisType, root, TagBasis)
      );
      threadId = received->threadId;
      enumCost = received->enumCost;
      basis    = LatticeBasis<int>(received->basis);
   }
#endif

}

///
/// send solution data to the rank
///
void
CMapLapParaBasisTh::send(
      UG::ParaComm *comm,
      int destination
      )
{

#ifdef _COMM_CPP11
   CMapLapParaCommTh *commTh = dynamic_cast<CMapLapParaCommTh*>(comm);
   if( commTh )
   {
      PARA_COMM_CALL(
         commTh->uTypeSend((void *)createDatatype(comm), ParaBasisType, destination, TagBasis)
      );
   }
   else
   {
      CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
      PARA_COMM_CALL(
         localComm->uTypeSend((void *)createDatatype(comm), ParaBasisType, destination, TagBasis)
      );
   }
#else
   CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
   PARA_COMM_CALL(
      localComm->uTypeSend((void *)createDatatype(comm), ParaCMapLAP::ParaBasisType, destination, TagBasis)
   );
#endif
}

///
/// receive solution data from the source rank
///
void
CMapLapParaBasisTh::receive(
      UG::ParaComm *comm,
      int source
      )
{

   std::unique_ptr<CMapLapParaBasis> received;

#ifdef _COMM_CPP11
   CMapLapParaCommTh *commTh = dynamic_cast<CMapLapParaCommTh*>(comm);
   if( commTh )
   {
      PARA_COMM_CALL(
         commTh->uTypeReceive((void **)&received, ParaBasisType, source, TagBasis)
      );
   }
   else
   {
      CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
      PARA_COMM_CALL(
            localComm->uTypeReceive((void **)&received, ParaBasisType, source, TagBasis)
      );
   }
#else
   CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
   PARA_COMM_CALL(
         localComm->uTypeReceive((void **)&received, ParaCMapLAP::ParaBasisType, source, TagBasis)
   );
#endif

   threadId = received->threadId;
   enumCost = received->enumCost;
   basis    = LatticeBasis<int>(received->basis);
}

} // namespace ParaCMapLAP