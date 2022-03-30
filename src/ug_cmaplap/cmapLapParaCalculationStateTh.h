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

/**@file    cmapLapParaCalculationStateTh.h
 * @brief   CalcutationState object extension for threads communication
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_CALCULATION_STATE_TH_H__
#define __CMAP_LAP_PARA_CALCULATION_STATE_TH_H__


#include "cmapLapParaCalculationState.h"
namespace UG { class ParaComm; }


namespace ParaCMapLAP
{

///
/// Calculation state object for thread communications
///
class CMapLapParaCalculationStateTh : public CMapLapParaCalculationState
{

   ///
   /// create thread object datatype of this object
   /// @return thread object dataypte of this object
   ///
   virtual CMapLapParaCalculationStateTh* createDatatype(
         );

public:

   ///
   /// default constructor of this object
   ///
   CMapLapParaCalculationStateTh(
         )
   {
   }

   ///
   /// constructor Deepbkz of this object
   ///
   CMapLapParaCalculationStateTh(
         int          inTermState,          ///< termination status, 0: normal, -1: interrupted
         int          inThreadId,           ///< thread id
         int          inCurrentBlockSize,   ///< current DeepBkz block size
         int          inTour,               ///< number of DeepBkz loop
         double       inElapsedTime,        ///< elapsed time
         double       inShortestNorm,       ///< the shortest norm found
         double       inApproxFactor,       ///< approximated factor
         double       inHermiteFactor,      ///< hermite factor
         double       inRootHermiteFactor,  ///< (hermite factor)^(1/dim)
         double       inEnumCost,           ///< log of approximted nodes of enumeration tree with incumbent radius
         double       inEnumCostGH,         ///< log of approximted nodes of enumeration tree with GH radius
         double       inSlopeGSA,           ///< slope of GSA line
         double       inTopHalfSlopeGSA,    ///< slope of top-half GSA line
         double       inOrthogonalFactor    ///< orthogonal factor
         )
         : CMapLapParaCalculationState(
               inTermState,
               inThreadId,
               inCurrentBlockSize,
               inTour,
               inElapsedTime,
               inShortestNorm,
               inApproxFactor,
               inHermiteFactor,
               inRootHermiteFactor,
               inEnumCost,
               inEnumCostGH,
               inSlopeGSA,
               inTopHalfSlopeGSA,
               inOrthogonalFactor
               )
   {
   }

   ///
   /// constructor Enum of this object
   ///
   CMapLapParaCalculationStateTh(
         int          inTermState,          ///< termination status, 0: normal, -1: interrupted
         int          inThreadId,           ///< thread id
         double       inElapsedTime,        ///< elapsed time
         double       inShortestNorm,       ///< the shortest norm found
         double       inApproxFactor,       ///< approximated factor
         double       inHermiteFactor,      ///< hermite factor
         double       inRootHermiteFactor,  ///< (hermite factor)^(1/dim)
         long int     inNumSearchedNodes    ///< number of searched nodes in the enumeration tree
         )
         : CMapLapParaCalculationState(
               inTermState,
               inThreadId,
               inElapsedTime,
               inShortestNorm,
               inApproxFactor,
               inHermiteFactor,
               inRootHermiteFactor,
               inNumSearchedNodes
               )
   {
   }

   ///
   /// constructor of Sieve
   ///
   CMapLapParaCalculationStateTh(
         int          inTermState,          ///< termination status, 0: normal, -1: interrupted
         int          inThreadId,           ///< thread id
         double       inElapsedTime,        ///< elapsed time
         int&         inBlockSize,          ///< block size
         long int&    inNLoop,              ///< number of Sieve algorithm loop
         int&         inListSize,           ///< size of List L
         int&         inStackSize,          ///< size of Stack S
         int&         inMaxListSize,        ///< maximum size of List L up to the point
         int&         inNCollisions,        ///< number of collision
         double&      inShortestNorm,       ///< the shortest norm found
         double&      inApproxFactor,       ///< approximated factor
         double&      inHermiteFactor,      ///< hermite factor
         double&      inRootHermiteFactor   ///< (hermite factor)^(1/dim)
         )
         : CMapLapParaCalculationState(
               inTermState,
               inThreadId,
               inElapsedTime,
               inBlockSize,
               inNLoop,
               inListSize,
               inStackSize,
               inMaxListSize,
               inNCollisions,
               inShortestNorm,
               inApproxFactor,
               inHermiteFactor,
               inRootHermiteFactor
               )
   {
   }

   ///
   /// destructor of this object
   ///
   virtual ~CMapLapParaCalculationStateTh(
         )
   {
   }

   ///
   /// send this object to destination
   ///
   virtual void send(
         UG::ParaComm *comm,  ///< communicator used to send this object
         int destination,     ///< destination rank to send
         int tag              ///< tag to show this object
         );

   ///
   /// receive this object from source
   ///
   virtual void receive(
         UG::ParaComm *comm,  ///< communicator used to receive this object
         int source,          ///< source rank to receive this object
         int tag              ///< tag to show this object
         );

};

} // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_CALCULATION_STATE_TH_H__
