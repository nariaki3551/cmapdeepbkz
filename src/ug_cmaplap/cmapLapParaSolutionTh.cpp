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

/**@file    cmapLapParaSolutionTh.cpp
 * @brief   CMapLapParaSolution extension for threads communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "cmapLapParaSolutionTh.h"
#include <memory>
#include <new>
#include <utility>
#include "cmapLapParaSolverLocalComm.h"
#include "ug/paraComm.h"
#include "ug/paraTagDef.h"
#include "cmapLapParaSolution.h"


namespace ParaCMapLAP
{

///
/// create clone of this object
///
UG::ParaSolution *
CMapLapParaSolutionTh::clone(UG::ParaComm *comm)
{
   LatticeVector<int> pv(v);
   return( new CMapLapParaSolutionTh(threadId, pv, projectedSquaredNorm));
}

///
/// create CMapLap Solution Datatype
///
CMapLapParaSolutionTh *
CMapLapParaSolutionTh::createDatatype(
      UG::ParaComm *comm
      )
{
   return dynamic_cast<CMapLapParaSolutionTh *>(clone(comm));
}

///
/// send solution data to the rank
///
void
CMapLapParaSolutionTh::bcast(
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
                  commTh->uTypeSend((void *)createDatatype(comm), ParaSolutionType, i, UG::TagSolution)
               );
            }
         }
      }
      else
      {
         std::unique_ptr<CMapLapParaSolutionTh> received;
         PARA_COMM_CALL(
            commTh->uTypeReceive((void **)&received, ParaSolutionType, root, UG::TagSolution)
         );
         threadId = received->threadId;
         projectedSquaredNorm = received->projectedSquaredNorm;
         v = LatticeVector<int>(received->v);
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
                     localComm->uTypeSend((void *)createDatatype(comm), ParaSolutionType, i, UG::TagSolution)
               );
            }
         }
      }
      else
      {
         std::unique_ptr<CMapLapParaSolutionTh> received;
         PARA_COMM_CALL(
               localComm->uTypeReceive((void **)&received, ParaSolutionType, root, UG::TagSolution)
         );
         threadId = received->threadId;
         projectedSquaredNorm = received->projectedSquaredNorm;
         v        = LatticeVector<int>(received->v);
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
                  localComm->uTypeSend((void *)createDatatype(comm), ParaCMapLAP::ParaSolutionType, i, UG::TagSolution)
            );
         }
      }
   }
   else
   {
      std::unique_ptr<CMapLapParaSolution> received;
      PARA_COMM_CALL(
            localComm->uTypeReceive((void **)&received, ParaCMapLAP::ParaSolutionType, root, UG::TagSolution)
      );
      threadId = received->threadId;
      projectedSquaredNorm = received->projectedSquaredNorm;
      v = LatticeVector<int>(received->v);
   }
#endif

}

///
/// send solution data to the rank
///
void
CMapLapParaSolutionTh::send(
      UG::ParaComm *comm,
      int destination
      )
{

#ifdef _COMM_CPP11
   CMapLapParaCommTh *commTh = dynamic_cast<CMapLapParaCommTh*>(comm);
   if( commTh )
   {
      PARA_COMM_CALL(
         commTh->uTypeSend((void *)createDatatype(comm), ParaSolutionType, destination, UG::TagSolution)
      );
   }
   else
   {
      CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
      PARA_COMM_CALL(
         localComm->uTypeSend((void *)createDatatype(comm), ParaSolutionType, destination, UG::TagSolution)
      );
   }
#else
   CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
   PARA_COMM_CALL(
      localComm->uTypeSend((void *)createDatatype(comm), ParaCMapLAP::ParaSolutionType, destination, UG::TagSolution)
   );
#endif
}

///
/// receive solution data from the source rank
///
void
CMapLapParaSolutionTh::receive(
      UG::ParaComm *comm,
      int source
      )
{

   std::unique_ptr<CMapLapParaSolution> received;

#ifdef _COMM_CPP11
   CMapLapParaCommTh *commTh = dynamic_cast<CMapLapParaCommTh*>(comm);
   if( commTh )
   {
      PARA_COMM_CALL(
         commTh->uTypeReceive((void **)&received, ParaSolutionType, source, UG::TagSolution)
      );
   }
   else
   {
      CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
      PARA_COMM_CALL(
            localComm->uTypeReceive((void **)&received, ParaSolutionType, source, UG::TagSolution)
      );
   }
#else
   CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
   PARA_COMM_CALL(
         localComm->uTypeReceive((void **)&received, ParaCMapLAP::ParaSolutionType, source, UG::TagSolution)
   );
#endif

   threadId = received->threadId;
   projectedSquaredNorm = received->projectedSquaredNorm;
   v = LatticeVector<int>(received->v);
}

} // namespace ParaCMapLAP