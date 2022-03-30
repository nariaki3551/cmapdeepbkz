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

/**@file    cmapDeepBkzParaLoadCoordinator.h
 * @brief   Load coordinator extension for CMAP-DeepBKZ.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_DEEPBKZ_PARA_LOADCOORDINATOR_H__
#define __CMAP_DEEPBKZ_PARA_LOADCOORDINATOR_H__

#include "cmapLapParaLoadCoordinator.h"



class CMapDeepBkzParaLoadCoordinator : public ParaCMapLAP::CMapLapParaLoadCoordinator
{

   std::unique_ptr<ParaCMapLAP::Lattice<int, double>> sharedLattice; ///< shared lattice amoung solvers
   int dimensionOfSharedLattice;    ///< dimension of shring lattice
   int numOfUnupdateSolverState;    ///< number of times the nSolver has been received since the last sharedLattice update.


   ///
   /// Objects to output statistics and progress state
   ///
   std::ofstream      ofsLogSharedLattice;               ///< ofstream for shared lattice log
   std::ostream       *osLogSharedLattice;               ///< ostream for shared lattice log


   ///////////////////////
   ///
   /// Message handlers
   ///
   ///////////////////////

   ///
   /// function to process TagSolverState message
   /// @return always 0 (for extension)
   ///
   int processTagSolverState(
         int source,                                      ///< source solver rank
         int tag                                          ///< TagSolverState
         ) override;

public:

   CMapDeepBkzParaLoadCoordinator(
      UG::ParaComm *inComm,
      UG::ParaParamSet *inParaParamSet,
      UG::ParaInitiator *inParaInitiator,
      bool *inRacingSolversExist,
      UG::ParaTimer *inParaTimer,
      UG::ParaDeterministicTimer *detTimer,
      int inNThreadsInSolver
      );

   ~CMapDeepBkzParaLoadCoordinator(
         );

};


#endif // __CMAP_DEEPBKZ_PARA_LOADCOORDINATOR_H__

