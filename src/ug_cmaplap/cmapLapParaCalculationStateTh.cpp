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

/**@file    cmapLapParaCalculationStateTh.cpp
 * @brief   CalcutationStte object extension for threads communication
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "cmapLapParaSolverLocalComm.h"
#include "cmapLapParaCalculationStateTh.h"

namespace ParaCMapLAP
{

CMapLapParaCalculationStateTh*
CMapLapParaCalculationStateTh::createDatatype(
      )
{
   switch( solverType )
   {
   case DeepBkz:
      return new CMapLapParaCalculationStateTh(
            terminationState,
            threadId,
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.currentBlockSize,
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.tour,
            compTime,  // elapsedTime
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.shortestNorm,
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.approxFactor,
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.hermiteFactor,
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.rootHermiteFactor,
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.enumCost,
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.enumCostGH,
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.slopeGSA,
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.topHalfSlopeGSA,
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.orthogonalFactor
            );
   case Enum:
      return new CMapLapParaCalculationStateTh(
            terminationState,
            threadId,
            compTime,  // elapsedTime
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateEnum.shortestNorm,
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateEnum.approxFactor,
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateEnum.hermiteFactor,
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateEnum.rootHermiteFactor,
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateEnum.numSearchedNodes
            );
   case Sieve:
      return new CMapLapParaCalculationStateTh(
            terminationState,
            threadId,
            compTime,  // elapsedTime
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.blockSize,
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.nLoop,
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.listSize,
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.stackSize,
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.maxListSize,
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.nCollisions,
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.shortestNorm,
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.approxFactor,
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.hermiteFactor,
            cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.rootHermiteFactor
            );
   default:
      break;
   }
   return 0;
}

void
CMapLapParaCalculationStateTh::send(
      UG::ParaComm *comm,
      int destination,
      int tag
      )
{

#ifdef _COMM_CPP11
   CMapLapParaCommTh *commTh = dynamic_cast<CMapLapParaCommTh*>(comm);
   if( commTh )
   {
      PARA_COMM_CALL(
         commTh->uTypeSend(createDatatype(), ParaCalculationStateType, destination, tag)
      );
   }
   else
   {
      CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
      PARA_COMM_CALL(
         localComm->uTypeSend(createDatatype(), ParaCalculationStateType, destination, tag)
      );
   }
#else
   CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
   PARA_COMM_CALL(
         localComm->uTypeSend(createDatatype(), ParaCalculationStateType, destination, tag)
   );
#endif
}

void
CMapLapParaCalculationStateTh::receive(
      UG::ParaComm *comm,
      int source,
      int tag
      )
{

   std::unique_ptr<CMapLapParaCalculationStateTh> received;

#ifdef _COMM_CPP11
   CMapLapParaCommTh *commTh = dynamic_cast<CMapLapParaCommTh*>(comm);
   if( commTh )
   {
      PARA_COMM_CALL(
         commTh->uTypeReceive((void **)&received, ParaCalculationStateType, source, tag)
      );
   }
   else
   {
      CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
      PARA_COMM_CALL(
         localComm->uTypeReceive((void **)&received, ParaCalculationStateType, source, tag)
      );
   }
#else
   CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
   PARA_COMM_CALL(
         localComm->uTypeReceive((void **)&received, ParaCalculationStateType, source, tag)
   );
#endif
   compTime = received->compTime;
   nSolved = received->nSolved;
   terminationState = received->terminationState;
   threadId = received->threadId;
   solverType = received->solverType;
   cmapLapParaCalculationStateData = received->cmapLapParaCalculationStateData;
}

} // namespace ParaCMapLAP