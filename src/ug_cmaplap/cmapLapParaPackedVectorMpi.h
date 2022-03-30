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

/**@file    cmapLapParaPackedVectorMpi.h
 * @brief   Base class for communicating vectors.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_PACKED_VECTOR_MPI_H
#define __CMAP_LAP_PARA_PACKED_VECTOR_MPI_H

#include <vector>
#include <mpi.h>
#include "ug/paraComm.h"
#include "cmapLapParaDef.h"
#include "cmapLapParaPackedVector.h"
namespace UG { class ParaComm; }
namespace ParaCMapLAP { class VectorElement; }


namespace ParaCMapLAP
{


class CMapLapParaPackedVectorMpi : public CMapLapParaPackedVector
{
   int len[3];
   std::vector<int> createdRanks;
   std::vector<double> squaredNorms;
   std::vector<int> allVectorElements;
   int dimension;

   virtual CMapLapParaPackedVectorMpi *createThreadDatatype(
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

  public:

   ///
   /// constrctor
   ///
   CMapLapParaPackedVectorMpi(
         ) : dimension(0)
   {
   }

   ///
   /// constructor
   ///
   CMapLapParaPackedVectorMpi(
         int inThreadId,
         LatticeBasis<int> &inVectors,
         int inRank,
         int nThreadsPerRank
         ) : CMapLapParaPackedVector(inThreadId, inVectors, inRank, nThreadsPerRank),
             dimension(0)
   {
   }

   ///
   /// constructor
   ///
   CMapLapParaPackedVectorMpi(
         int inThreadId,
         std::vector<std::shared_ptr<VectorElement>> inVectorElements
         ) : CMapLapParaPackedVector(inThreadId, inVectorElements),
             dimension(0)
   {
   }

   ///
   /// constructor
   ///
   CMapLapParaPackedVectorMpi(
         int inThreadId,
         int inNVectors
         ): CMapLapParaPackedVector(inThreadId, inNVectors),
            dimension(0)
   {
   }

   ///
   /// destructor
   ///
   virtual ~CMapLapParaPackedVectorMpi(
         )
   {
   }

   virtual CMapLapParaPackedVectorMpi *clone(
         UG::ParaComm *comm
         );
   virtual MPI_Datatype createDatatype(
         int inNVectors,
         int inDimension
         );
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

#endif  // __CMAP_LAP_PARA_PACKED_VECTOR_MPI_H
