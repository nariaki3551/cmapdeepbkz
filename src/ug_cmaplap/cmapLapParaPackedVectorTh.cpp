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

/**@file    cmapLapParaPackedVectorTh.cpp
 * @brief   CMapLapParaPackedVector extension for threads communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "cmapLapParaPackedVectorTh.h"
#include <memory>
#include "ug/paraComm.h"
#include "cmapLapParaTagDef.h"
#include "cmapLapParaCommTh.h"
#include "cmapLapParaSolverLocalComm.h"
#include "cmapLapParaPackedVector.h"
namespace ParaCMapLAP { class VectorElement; }


namespace ParaCMapLAP
{

///
/// @brief create clone of this object
/// @note copy pointer of vector element object because vector element will not be change
///
CMapLapParaPackedVectorTh *
CMapLapParaPackedVectorTh::clone(UG::ParaComm *comm)
{
   std::vector<std::shared_ptr<VectorElement>> inVectorElements(nVectors);
   for (int i = 0; i < nVectors; i++ )
   {
      inVectorElements[i] = vectorElements[i];
   }
   return( new CMapLapParaPackedVectorTh(threadId, inVectorElements) );
}

///
/// create CMapLapDiffSubproblemPreDatatype
///
CMapLapParaPackedVectorTh *
CMapLapParaPackedVectorTh::createDatatype(
      UG::ParaComm *comm
      )
{
   return clone(comm);
}

///
/// send solution data to the rank
///
void
CMapLapParaPackedVectorTh::send(
      UG::ParaComm *comm,
      int destination
      )
{

#ifdef _COMM_CPP11
   CMapLapParaCommTh *commTh = dynamic_cast<CMapLapParaCommTh*>(comm);
   if( commTh )
   {
      PARA_COMM_CALL(
         commTh->uTypeSend((void *)createDatatype(comm), CMapLapParaPackedVectorType, destination, TagCMapLapPackedVector)
      );
   }
   else
   {
      CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
      PARA_COMM_CALL(
         localComm->uTypeSend((void *)createDatatype(comm), CMapLapParaPackedVectorType, destination, TagCMapLapPackedVector)
      );
   }
#else
   CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
   PARA_COMM_CALL(
      localComm->uTypeSend((void *)createDatatype(comm), CMapLapParaPackedVectorType, destination, TagCMapLapPackedVector)
   );
#endif

}

///
/// receive solution data from the source rank
///
void
CMapLapParaPackedVectorTh::receive(
      UG::ParaComm *comm,
      int source
      )
{

   std::unique_ptr<CMapLapParaPackedVector> received;

#ifdef _COMM_CPP11
   CMapLapParaCommTh *commTh = dynamic_cast<CMapLapParaCommTh*>(comm);
   if( commTh )
   {
      PARA_COMM_CALL(
         commTh->uTypeReceive((void **)&received, CMapLapParaPackedVectorType, source, TagCMapLapPackedVector)
      );
   }
   else
   {
      CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
      PARA_COMM_CALL(
            localComm->uTypeReceive((void **)&received, CMapLapParaPackedVectorType, source, TagCMapLapPackedVector)
      );
   }
#else
   CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
   PARA_COMM_CALL(
         localComm->uTypeReceive((void **)&received, CMapLapParaPackedVectorType, source, TagCMapLapPackedVector)
   );
#endif

   threadId = received->threadId;
   nVectors = received->nVectors;

   if( nVectors > 0 )
   {
      vectorElements.resize(nVectors);
      for (int i = 0; i < nVectors; i++ )
      {
         vectorElements[i] = received->getVectorElements()[i];
      }
   }
}

} // namespace ParaCMapLAP