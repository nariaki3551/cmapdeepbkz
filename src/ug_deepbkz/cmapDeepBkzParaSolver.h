/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*          This file is part of the program and software framework          */
/*  CMAP-DeepBKZ --- Configurable Massively Parallel Solver for DeepBKZ      */
/*                                                                           */
/*  Copyright Written by Nariaki Tateiwa <n-tateiwa@kyudai.jp>,              */
/*                       Yuji Shinano <shinano@zib.de>,                      */
/*            Copyright (C) 2021 by Zuse Institute Berlin,                   */
/*            licensed under LGPL version 3 or later.                        */
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

/**@file    cmapDeepBkzParaSolver.h
 * @brief   ParaSolver extension for CMAP-DeepBKZ
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_DEEPBKZ_PARA_SOLVER_H__
#define __CMAP_DEEPBKZ_PARA_SOLVER_H__

#include "cmapLapParaSolver.h"


///
/// @class CMapDeepBkzParaSolver
///
class CMapDeepBkzParaSolver : public ParaCMapLAP::CMapLapParaSolver
{

protected:

   void solve();
   void runDeepBkz();
   void runExDeepBkz();


public:

   CMapDeepBkzParaSolver(
         int argc,
         char **argv,
         UG::ParaComm *comm,
         ParaCMapLAP::CMapLapParaSolverLocalComm *localComm,
         UG::ParaParamSet *paraParamSet,
         UG::ParaInstance *paraInstance,
         UG::ParaDeterministicTimer *detTimer,
         double timeOffset);

   virtual ~CMapDeepBkzParaSolver(
         );

};


#endif // __CMAP_DEEPBKZ_PARA_SOLVER_H__
