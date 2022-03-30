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

/**@file    cmapLapParaSolutionMpi.h
 * @brief   CMapLapParaSolution extension for MPI communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_SOLUTION_MPI_H__
#define __CMAP_LAP_PARA_SOLUTION_MPI_H__

#include <mpi.h>
#include <vector>
#include "cmapLapParaDef.h"
#include "cmapLapParaSolution.h"
namespace UG { class ParaComm; }

namespace ParaCMapLAP
{

///
/// CMapLapSolution class
///
class CMapLapParaSolutionMpi : public CMapLapParaSolution
{
   int lenv;
   std::vector<int> vec;   // allocated memory for vector sending

   ///
   /// create CMapLapParaSolutionDatatype
   ///
   virtual MPI_Datatype createDatatype();

   ///
   /// create cmapLapSolutionDatatype
   ///
   virtual CMapLapParaSolutionMpi *createThreadDatatype(
         UG::ParaComm *comm
         )
   {
      return dynamic_cast<CMapLapParaSolutionMpi *>(clone(comm));
   }

public:

   ///
   /// default constructor
   ///
   CMapLapParaSolutionMpi(
	      ) : lenv(0)
   {
   }

   ///
   /// constructor
   ///
   CMapLapParaSolutionMpi(
         int inThreadId,               ///< thread Id
         LatticeVector<int> &inV,      ///< updated vector
         double inProjectedSquaredNorm ///< objective function value
         )
        : CMapLapParaSolution(inThreadId, inV, inProjectedSquaredNorm)
   {
      lenv = inV.size();
   }

   ///
   /// destructor
   ///
   virtual ~CMapLapParaSolutionMpi(
         )
   {
   }

   ///
   /// create clone of this object
   ///
   virtual UG::ParaSolution *clone(
         UG::ParaComm *comm
         );

   ///
   /// broadcast solution data to from the root rank
   ///
   virtual void bcast(
         UG::ParaComm *comm,
         int root
         );

   ///
   /// send solution data to the rank
   ///
   virtual void send(
         UG::ParaComm *comm,
         int destination
         );

   ///
   /// receive solution data from the source rank
   ///
   virtual void receive(
         UG::ParaComm *comm,
         int source
         );

};

} // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_SOLUTION_MPI_H__
