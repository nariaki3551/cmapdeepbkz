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

/**@file    cmapLapParaInstanceMpi.h
 * @brief   CMapLapParaInstance extension for MPI communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_INSTANCE_MPI_H__
#define __CMAP_LAP_PARA_INSTANCE_MPI_H__

#include <cstring>
#include <mpi.h>
#include "cmapLapParaInstance.h"
namespace UG { class ParaComm; }

namespace ParaCMapLAP
{

///
/// CMapLapInstanceMpi
///
class CMapLapParaInstanceMpi : public CMapLapParaInstance
{

   ///
   /// create deepParaInstancePreType
   ///
   virtual MPI_Datatype createDatatype1();

   ///
   /// create deepParaInstanceType
   ///
   virtual MPI_Datatype createDatatype2(bool memAllocNecessary);

   virtual void allocateMemoryForDatatype2();

   ///
   /// create clone of this object
   ///
   virtual CMapLapParaInstanceMpi *clone(UG::ParaComm *comm)
   {
      return( new CMapLapParaInstanceMpi(probName) );
   }

   ///
   /// create instance object
   ///
   virtual CMapLapParaInstanceMpi *createThreadDatatype(
         UG::ParaComm *comm
         )
   {
      return clone(comm);
   }

public:
   ///
   /// constructor
   ///
   CMapLapParaInstanceMpi(
         )
   {
   }

   /// constructor : only called from CMapLapInitiator

   ///
   /// constructor for cloning
   ///
   CMapLapParaInstanceMpi(
         char      *inProbFileName,
         char      *inProbName
         ) : CMapLapParaInstance(inProbFileName, inProbName)
   {
   }

   ///
   /// constructor for cloning
   ///
   CMapLapParaInstanceMpi(
         char      *inProbName
         ) : CMapLapParaInstance(inProbName)
   {
   }

   ///
   /// destractor
   ///
   virtual ~CMapLapParaInstanceMpi(
	        )
   {
   }

   ///
   /// broadcasts instance to all solvers
   ///
   virtual int bcast(
         UG::ParaComm *comm,
         int rank,
         int method
         );

};

} // namespace ParaCMapLAP

#endif  // __CMAP_LAP_PARA_INSTANCE_MPI_H__
