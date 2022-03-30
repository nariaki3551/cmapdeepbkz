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

/**@file    cmapLapParaCalculationState.h
 * @brief   Base class for calculation state.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_CALCULATION_STATE_H__
#define __CMAP_LAP_PARA_CALCULATION_STATE_H__

#include <iostream>
#include <iomanip>
#include <climits>
#include <cfloat>
#include <memory>
#include "ug/paraComm.h"
#include "ug/paraCalculationState.h"
#include "cmapLapParaLog.h"

namespace ParaCMapLAP
{

typedef struct CMapLapParaCalculationStateDeepBkz_ {
   int      currentBlockSize;   ///< current DeepBkz block size
   int      tour;               ///< number of DeepBkz loop
   double   shortestNorm;       ///< the shortest norm found
   double   approxFactor;       ///< approximate factor
   double   hermiteFactor;      ///< hermite factor
   double   rootHermiteFactor;  ///< (hermite factor)^(1/dim)
   double   enumCost;           ///< approximated nodes of enumeration tree with incumbent radius
   double   enumCostGH;         ///< approximated nodes of enumeration tree with GH radius
   double   slopeGSA;           ///< slope of GSA line [Pool] better a than b when a > b
   double   topHalfSlopeGSA;    ///< slope of top-half GSA line [Pool] better a than b when a > b
   double   orthogonalFactor;   ///< orthogonal factor [Pool] better a than b when a < b
} CMapLapParaCalculationStateDeepBkz;

typedef struct CMapLapParaCalculationStateEnum_ {
   double   shortestNorm;       ///< the shortest norm found
   double   approxFactor;       ///< approximate factor
   double   hermiteFactor;      ///< hermite factor
   double   rootHermiteFactor;  ///< (hermite false)^(1/dim)
   long int numSearchedNodes;   ///< number of searched nodes in the enumeration tree [Pool] better a than b when b > a
} CMapLapParaCalculationStateEnum;

typedef struct CMapLapParaCalculationStateSieve_ {
   int         blockSize;          ///< block size
   long int    nLoop;              ///< number of Sieve loop
   int         listSize;           ///< size of List L
   int         stackSize;          ///< size of Stack S
   int         maxListSize;        ///< maximum size of List L up to the point
   int         nCollisions;        ///< number of collision
   double      shortestNorm;       ///< the shortest norm found
   double      approxFactor;       ///< approximate factor
   double      hermiteFactor;      ///< hermite factor
   double      rootHermiteFactor;  ///< (hermite false)^(1/dim)
} CMapLapParaCalculationStateSieve;

typedef union CMapLapParaCalculationStateData_ {
   CMapLapParaCalculationStateDeepBkz cmapLapParaCalculationStateDeepBkz;
   CMapLapParaCalculationStateEnum    cmapLapParaCalculationStateEnum;
   CMapLapParaCalculationStateSieve   cmapLapParaCalculationStateSieve;
} CMapLapParaCalculationStateData;


///
/// \class CMapLapParaCalculationState
/// Base class of Calculation state in a ParaSolver
///
class CMapLapParaCalculationState : public UG::ParaCalculationState, public std::enable_shared_from_this<CMapLapParaCalculationState>
{
protected:

   int threadId;  ///< threadId = -1 when solverType is Sieve
   SolverType solverType;
   CMapLapParaCalculationStateData cmapLapParaCalculationStateData;

public:

   ///
   /// Default Constructor
   ///
   CMapLapParaCalculationState(
         )
         : ParaCalculationState(),
           threadId(-1),
           solverType(Undefined)
   {
   }

   ///
   /// Constructor
   ///
   CMapLapParaCalculationState(
         double inCompTime,                   ///< computation time of this ParaNode
         int    inNSolved,                    ///< the number of nodes solved
         int    inTerminationState,           ///< indicate whether if this computation is terminationState or not. 0: no, 1: terminationState
         int    inThreadId,                   ///< thread ID
         SolverType inSolverType              ///< solver type
         )
         : ParaCalculationState(inCompTime, inNSolved, inTerminationState),
           threadId(inThreadId),
           solverType(inSolverType)
   {
   }

   ///
   /// Constructor of DeepBkz
   ///
   CMapLapParaCalculationState(
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
         : ParaCalculationState(inElapsedTime, 1, inTermState),
           threadId(inThreadId),
           solverType(DeepBkz)
   {
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.currentBlockSize = inCurrentBlockSize;
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.tour = inTour;
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.shortestNorm = inShortestNorm;
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.approxFactor = inApproxFactor;
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.hermiteFactor = inHermiteFactor;
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.rootHermiteFactor = inRootHermiteFactor;
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.enumCost = inEnumCost;
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.enumCostGH = inEnumCostGH;
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.slopeGSA = inSlopeGSA;
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.topHalfSlopeGSA = inTopHalfSlopeGSA;
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.orthogonalFactor = inOrthogonalFactor;
   }

   ///
   /// Constructor of Enum
   ///
   CMapLapParaCalculationState(
         int          inTermState,          ///< termination status
         int          inThreadId,           ///< thread id
         double       inElapsedTime,        ///< elapsed time
         double       inShortestNorm,       ///< the shortest norm found
         double       inApproxFactor,       ///< approximate factor
         double       inHermiteFactor,      ///< hermite factor
         double       inRootHermiteFactor,  ///< (hermite false)^(1/dim)
         long int     inNumSearchedNodes    ///< number of searched nodes in the enumeration tree
         )
         : ParaCalculationState(inElapsedTime, 1, inTermState),
           threadId(inThreadId),
           solverType(Enum)
   {
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateEnum.shortestNorm = inShortestNorm;
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateEnum.approxFactor = inApproxFactor;
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateEnum.hermiteFactor = inHermiteFactor;
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateEnum.rootHermiteFactor = inRootHermiteFactor;
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateEnum.numSearchedNodes = inNumSearchedNodes;
   }

   ///
   /// Constructor of Sieve
   ///
   CMapLapParaCalculationState(
         int          inTermState,          ///< termination status
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
         : ParaCalculationState(inElapsedTime, 1, inTermState),
           threadId(inThreadId),
           solverType(Sieve)
   {
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.blockSize = inBlockSize;
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.nLoop = inNLoop;
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.listSize = inListSize;
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.stackSize = inStackSize;
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.maxListSize = inMaxListSize;
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.nCollisions = inNCollisions;
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.shortestNorm = inShortestNorm;
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.approxFactor = inApproxFactor;
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.hermiteFactor = inHermiteFactor;
      cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.rootHermiteFactor = inRootHermiteFactor;
   }


   ///
   /// Destructor
   ///
   virtual
   ~CMapLapParaCalculationState(
         )
   {
   }


   ///
   /// get solverType
   /// @return solverType
   ///
   virtual SolverType getSolverType(
         )
   {
      return solverType;
   }

   ///
   /// get thread ID
   /// @return thread ID
   ///
   virtual int getThreadId(
         )
   {
      return threadId;
   }

   ///
   /// stringfy CMapLapParaCalculationState
   /// @return string to show this object
   ///
   virtual std::string toString(
         )
   {
      std::ostringstream s;
      if( terminationState )
      {
         s << "Termination state of this computation was " << terminationState << " : [ "
         << compTime << " sec. computed ]"
         << nSolved << " nodes were solved. ";
      }
      else
      {
         s << "Computation was normally terminated: [ "
         << compTime << " sec. computed ]"
         << nSolved << " nodes were solved. ";
      }
      return s.str();
   }

   ///
   /// stringfy CMapLapParaCalculationState (simple string version)
   /// @return simple string to show this object
   ///
   virtual std::string toSimpleString(
         )
   {
      std::ostringstream s;

      s << compTime
            << ", "
            << nSolved
            << ", ";

      return s.str();
   }

   ///
   /// stringfy CMapLapParaCalculationState
   /// @return simple string to show CMapLapParaCalculationState
   ///
   virtual std::string toStopString(
         std::string delimiter=""
         )
   {
      switch( solverType )
      {
      case DeepBkz:
         return toStringLogDeepBkz(delimiter);
      case Enum:
         return toStringLogEnum(delimiter);
      case Sieve:
         return toStringLogSieve(delimiter);
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaSolverPoolElement: Invalid solver type = ",static_cast<int>(solverType));
      }
   }

   ///
   /// stringfy CMapLapParaCalculationStateDeepBkz
   /// @return simple string to show CMapLapParaCalculationStateDeepBkz
   ///
   virtual const std::string toStringLogDeepBkz(
         std::string delimiter=""
         )
   {
      std::string taskName       = "DeepBkz";
      double   elapsedTime       = getCompTime();
      int      size              = cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.currentBlockSize;
      long int iter              = cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.tour;
      double   progress          = 0.0;
      double   leftTime          = cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.enumCost;
      double   logCost           = cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.enumCost;
      double   shortestNorm      = cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.shortestNorm;
      double   approxFactor      = cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.approxFactor;
      double   hermiteFactor     = cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.hermiteFactor;
      double   rootHermiteFactor = cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.rootHermiteFactor;
      double   orthogonalFactor  = cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.orthogonalFactor;
      int      meanMessageQueueSize = -1;
      int      maxMessageQueueSize  = -1;
      std::ostringstream s_appendix;
      s_appendix << " -- rho " << std::fixed << std::setprecision(5)
                 << cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.slopeGSA
                 << ", half rho " << std::fixed << std::setprecision(5)
                 << cmapLapParaCalculationStateData.cmapLapParaCalculationStateDeepBkz.topHalfSlopeGSA;

      return Logging::toStringLogBase(taskName, elapsedTime, size, iter, progress, leftTime, logCost,
            shortestNorm, approxFactor, hermiteFactor, rootHermiteFactor, orthogonalFactor,
            meanMessageQueueSize, maxMessageQueueSize,
            s_appendix.str(), delimiter);
   }

   ///
   /// stringfy CMapLapParaCalculationStateEnum
   /// @return simple string to show CMapLapParaCalculationStateEnum
   ///
   virtual const std::string toStringLogEnum(
         std::string delimiter=""
         )
   {
      std::string taskName       = "Enum";
      double   elapsedTime       = getCompTime();
      int      size              = 0;
      long int iter              = cmapLapParaCalculationStateData.cmapLapParaCalculationStateEnum.numSearchedNodes;
      double   logCost           = -1.0;
      double   progress          = 0.0;
      double   leftTime          = 0.0;
      double   shortestNorm      = cmapLapParaCalculationStateData.cmapLapParaCalculationStateEnum.shortestNorm;
      double   approxFactor      = cmapLapParaCalculationStateData.cmapLapParaCalculationStateEnum.approxFactor;
      double   hermiteFactor     = cmapLapParaCalculationStateData.cmapLapParaCalculationStateEnum.hermiteFactor;
      double   rootHermiteFactor = cmapLapParaCalculationStateData.cmapLapParaCalculationStateEnum.rootHermiteFactor;
      double   orthogonalFactor  = -1;
      int      meanMessageQueueSize = -1;
      int      maxMessageQueueSize  = -1;
      std::ostringstream s_appendix;

      return Logging::toStringLogBase(taskName, elapsedTime, size, iter, progress, leftTime, logCost,
            shortestNorm, approxFactor, hermiteFactor, rootHermiteFactor, orthogonalFactor,
            meanMessageQueueSize, maxMessageQueueSize,
            s_appendix.str(), delimiter);
   }


   ///
   /// stringfy CMapLapParaCalculationStateSieve
   /// @return simple string to show CMapLapParaCalculationStateSieve
   ///
   virtual const std::string toStringLogSieve(
         std::string delimiter=""
         )
   {
      std::string taskName       = "Sieve";
      double   elapsedTime       = getCompTime();
      int      size              = cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.blockSize;
      long int iter              = cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.nLoop;
      double   progress          = 0.0;
      double   leftTime          = 0.0;
      double   logCost           = 0.0;
      double   shortestNorm      = cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.shortestNorm;
      double   approxFactor      = cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.approxFactor;
      double   hermiteFactor     = cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.hermiteFactor;
      double   rootHermiteFactor = cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.rootHermiteFactor;
      double   orthogonalFactor  = -1;
      int      meanMessageQueueSize = -1;
      int      maxMessageQueueSize  = -1;
      std::ostringstream s_appendix;
      s_appendix << " -- |L| "
                 << cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.listSize
                 << ", |S| "
                 << cmapLapParaCalculationStateData.cmapLapParaCalculationStateSieve.stackSize;

      return Logging::toStringLogBase(taskName, elapsedTime, size, iter, progress, leftTime, logCost,
            shortestNorm, approxFactor, hermiteFactor, rootHermiteFactor, orthogonalFactor,
            meanMessageQueueSize, maxMessageQueueSize,
            s_appendix.str(), delimiter);
   }

};

} // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_CALCULATION_STATE_H__
