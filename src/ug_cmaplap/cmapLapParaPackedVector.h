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

/**@file    cmapLapParaPackedVector.h
 * @brief   CMapLapParaPackedVector extension for MPI communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/



#ifndef __CMAP_LAP_PARA_PACKED_VECTOR_H__
#define __CMAP_LAP_PARA_PACKED_VECTOR_H__

#include <memory>
#include <vector>
#include "ug/paraComm.h"
#include "cmapLapParaDef.h"
#include "cmapLapParaShareData.h"


namespace ParaCMapLAP
{

class CMapLapParaPackedVectorTh;
class CMapLapParaPackedVectorMpi;


///
/// @class CMapLapParaPackedVector
/// @brief Container of lattice vectors
///
class CMapLapParaPackedVector : public std::enable_shared_from_this<CMapLapParaPackedVector>
{

protected:

   int threadId;  ///< thread Id
   int nVectors;  ///< number of vectors
   std::vector<std::shared_ptr<VectorElement>> vectorElements;

public:

   friend CMapLapParaPackedVectorTh;
   friend CMapLapParaPackedVectorMpi;

   ///
   /// constrctor
   ///
   CMapLapParaPackedVector(
         )
         :
         threadId(-1),
         nVectors(0),
         vectorElements(0)
   {
   }
   CMapLapParaPackedVector(
         int inThreadId,
         LatticeBasis<int> &vectors,
         int rank,
         int nThreadsPerRank
         )
         :
         threadId(inThreadId),
         nVectors(vectors.rows()),
         vectorElements(0)
   {
      assert( inThreadId >= 0 && rank >= 1 );
      unsigned int solverThreadId = ((rank - 1)*nThreadsPerRank + threadId);
      vectorElements.resize(nVectors);
      for( int i = 0; i < nVectors; i++ )
      {
         LatticeVector<int> vector(vectors.row(i));
         vectorElements[i] = std::make_shared<VectorElement>(vector, solverThreadId);
      }
   }
   CMapLapParaPackedVector(
         int inThreadId,
         std::vector<std::shared_ptr<VectorElement>> inVectorElements
         )
      :
         threadId(inThreadId),
         nVectors(inVectorElements.size()),
         vectorElements(inVectorElements)
   {
   }
   CMapLapParaPackedVector(
         int inThreadId,
         int inNVectors
         )
      :
         threadId(inThreadId),
         nVectors(inNVectors)
   {
      vectorElements.resize(nVectors);
   }


   virtual ~CMapLapParaPackedVector(
         )
   {
   }


   ///
   /// @brief getter threadID
   /// @return thread ID
   ///
   virtual int getThreadId(
         )
   {
      return threadId;
   }


   ///
   /// @brief getter of number of lattice vectors
   /// @return number of lattice vectors
   ///
   virtual int getNVectors(
         )
   {
      return nVectors;
   }


   ///
   /// @brief getter of vector elements
   /// @return vector elements
   ///
   virtual std::vector<std::shared_ptr<VectorElement>> & getVectorElements(
         )
   {
      return vectorElements;
   }

   ///
   /// @brief getter of dimension of lattice vectors
   /// @return dimension of lattice vectors
   ///
   virtual int getDimension(
         )
   {
      assert( nVectors > 0 );
      return vectorElements[0]->size();
   }


   ///
   /// @brief extract a lattice vector
   /// @return a lattice vector
   ///
   virtual std::shared_ptr<VectorElement> extract(
         )
   {
      if( nVectors == 0 ) return nullptr;
      assert( vectorElements[nVectors-1] );
      return vectorElements[--nVectors];
   }


   virtual void send(
         UG::ParaComm *comm,
         int destination) = 0;


   virtual void receive(
         UG::ParaComm *comm,
         int source) = 0;


   virtual const std::string toString(
         )
   {
      std::ostringstream s;
      s << std::endl << std::endl;
      s << "===== Packed Vector =====" << std::endl;
      s << "nVectors " << nVectors << std::endl;
      for( int i = 0; i < nVectors; i++ )
      {
         s << "i " << i << std::endl;
         s << vectorElements[i]->toString();
         s << std::endl;
      }
      s << "===========================" << std::endl;
      return s.str();
   }
};

} // namespace ParaCMapLAP

#endif  // __CMAP_LAP_PARA_PACKED_VECTOR_H__
