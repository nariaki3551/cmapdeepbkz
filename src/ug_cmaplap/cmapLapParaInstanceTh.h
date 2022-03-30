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

/**@file    cmapLapParaInstanceTh.h
 * @brief   CMapLapParaInstance extension for threads communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_INSTANCE_TH_H__
#define __CMAP_LAP_PARA_INSTANCE_TH_H__

#ifdef _COMM_PTH
#include "ug/paraCommPth.h"
#endif
#ifdef _COMM_CPP11
#include "ug/paraCommCPP11.h"
#endif

#include "cmapLapParaInstance.h"  // for CMapLapParaInstance


namespace ParaCMapLAP
{

///
/// CMapLapInstanceTh
///
class CMapLapParaInstanceTh : public CMapLapParaInstance
{
   ///
   /// create clone of this object
   ///
   virtual CMapLapParaInstanceTh *clone(UG::ParaComm *comm)
   {
      return( new CMapLapParaInstanceTh(probName) );
   }

   ///
   /// create instance object
   ///
   virtual CMapLapParaInstanceTh *createDatatype(
         UG::ParaComm *comm
         )
   {
      return clone(comm);
   }

public:
   ///
   /// constructor
   ///
   CMapLapParaInstanceTh(
         )
   {
   }
   CMapLapParaInstanceTh(
         const CMapLapParaInstance& cmapLapParaInstance
         ) : CMapLapParaInstance(cmapLapParaInstance)
   {
   }

   /// constructor : only called from CMapLapInitiator

   ///
   /// constructor for cloning
   ///
   CMapLapParaInstanceTh(
         char      *inProbFileName,
         char      *inProbName
         ) : CMapLapParaInstance(inProbFileName, inProbName)
   {
   }

   ///
   /// constructor for cloning
   ///
   CMapLapParaInstanceTh(
         char      *inProbName
         ) : CMapLapParaInstance(inProbName)
   {
   }

   ///
   /// destractor
   ///
   virtual ~CMapLapParaInstanceTh(
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

#endif  // __CMAP_LAP_PARA_INSTANCE_TH_H__

