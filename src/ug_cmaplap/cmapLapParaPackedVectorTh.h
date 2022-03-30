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

/**@file    cmapLapParaPackedVectorTh.h
 * @brief   CMapLapParaPackedVector extension for threads communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __DEEP_BKZ_PARA_PACKED_VECTOR_TH_H
#define __DEEP_BKZ_PARA_PACKED_VECTOR_TH_H

#ifdef _COMM_CPP11
#include "ug/paraCommCPP11.h"
#endif
#include "cmapLapParaPackedVector.h"
#include "cmapLapParaDef.h"
namespace ParaCMapLAP { class VectorElement; }
namespace UG { class ParaComm; }


namespace ParaCMapLAP
{


class CMapLapParaPackedVectorTh : public CMapLapParaPackedVector
{


public:

   ///
   /// constrctor
   ///
   CMapLapParaPackedVectorTh(
         )
   {
   }

   CMapLapParaPackedVectorTh(
         int      inThradId,
         LatticeBasis<int> &inVectors,
         int      inRank,
         int      nThreadsPerRank
         )
         : CMapLapParaPackedVector(inThradId, inVectors, inRank, nThreadsPerRank)
   {
   }

   CMapLapParaPackedVectorTh(
         int inThreadId,
         std::vector<std::shared_ptr<VectorElement>> &inVectorElements
    		)
         : CMapLapParaPackedVector(inThreadId, inVectorElements)
   {
   }

   CMapLapParaPackedVectorTh(
         int inThreadId,
         int inNVectors
    		)
         : CMapLapParaPackedVector(inThreadId, inNVectors)
   {
   }

   virtual ~CMapLapParaPackedVectorTh(
         )
   {
   }

   virtual CMapLapParaPackedVectorTh *clone(
         UG::ParaComm *comm);

   virtual CMapLapParaPackedVectorTh *createDatatype(
         UG::ParaComm *comm);

   virtual void send(
         UG::ParaComm *comm,
         int destination
         );

   virtual void receive(
         UG::ParaComm *comm,
         int source
         );

};

} // namespace ParaCMapLAP

#endif  // __DEEP_BKZ_PARA_PACKED_VECTOR_TH_H
