/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*          This file is part of the program and software framework          */
/*                 CMAP-TEST --- Test configure for CMAP-LAP                 */
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

/**@file    cmapTestEnumeration.h
 * @brief   Enumeration extension for CMAP-TEST.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_TEST_ENUM_H__
#define __CMAP_TEST_ENUM_H__

#include <memory>
#include "Enumeration.h"
#include "Config.h"
#include "Def.h"
#include "Lattice.h"
#include "Timer.h"
#include "VectorElement.h"
#include "VectorElementPool.h"
#include "cmapTestParaSolver.h"



template<typename BasisFloat=int, typename GSFloat=double, typename EnumGSFloat=double>
class CmapEnumeration : public LapTools::Enumeration<BasisFloat, GSFloat, EnumGSFloat>
{

using VectorElementType = LapTools::VectorElementBase<BasisFloat>;
using VectorElementPtr = std::shared_ptr<VectorElementType>;
using VectorElementQueue = LapTools::VectorElementBasePool<BasisFloat>;
using LatticePtr = std::shared_ptr<LapTools::Lattice<BasisFloat, GSFloat>>;


private:

   LatticePtr L;              ///< lattice
   LapTools::Config config;   ///< hyper parameters
   int rank;                  ///< solver-id
   int threadId;              ///< solver-id
   int verbose;               ///< verbose <= 0: none, 1: light, 2: medium, 3: heavy


protected:

   CMapTestParaSolver *cmapLapParaSolver; ///< communicator to LC

   int nSendVectors;                                  ///< number of vectors per one-communication
   int samplingDepth;                                 ///< depth of sampling
   VectorElementQueue SendStack;                      ///< stack
   LapTools::LatticeVector<BasisFloat> sendVector;    ///< vector for stock

   double nextIReceiveMessagesTime;                   ///< previous time when it run iReceiveMessages()
   double iReceiveMessagesInterval;                   ///< interval time for iReceiveMessages()
   LapTools::LatticeVector<int> bestVector;           ///< shortest vector


public:

   ///
   /// @brief Constructor
   ///
   CmapEnumeration(){};


   ///
   /// @brief Constructor
   /// @param[in] inL         Lattice
   /// @param[in] inVerbose   <= 0 : not, 1 : light, 2: heavy
   /// @param[in] inRank      solver-id
   /// @param[in] inThreadId  solver-id
   /// @param[in] inIReceiveMessagesInterval
   ///
   CmapEnumeration(
         LatticePtr inL,
         CMapTestParaSolver *inCmapLapParaSolver,
         int inRank=1,
         int inThreadId=0,
         int inVerbose=0,
         double inIReceiveMessagesInterval=1
         )
      :
         LapTools::Enumeration<BasisFloat, GSFloat, EnumGSFloat>(inL, inRank, inThreadId, inVerbose),
         rank(inRank),
         threadId(inThreadId),
         verbose(inVerbose),
         nSendVectors(0),
         samplingDepth(-1),
         iReceiveMessagesInterval(inIReceiveMessagesInterval)
   {
      L = inL;
      config = inL->config;
      cmapLapParaSolver = inCmapLapParaSolver;
      nextIReceiveMessagesTime = LapTools::Timer::getElapsedTime();
      SendStack = VectorElementQueue(nSendVectors * 10);
   }


   ///
   /// @brief post process of search node
   ///
   bool postProcess(
         );


   ///
   /// @brief communicate with LC
   /// @param[out] shouldAbort true if it should abort else false
   ///
   bool communicate(
         bool &shouldAbort
         );


   ///
   /// @brief update bestObjectiveValue
   /// @param[in] sigh line-header character
   /// @param[in] objectiveValue new objectiveValue
   ///
   bool updateBestObjectiveValue(
         char sigh,
         EnumGSFloat objectiveValue
         );


   ///
   /// @brief send SolverState of Enumeration
   ///
   void sendStatus(
         );


   ///
   /// @breif setter of nSendVectors
   ///
   void setNSendVectors(
         int inNSendVectors
         )
   {
      nSendVectors = inNSendVectors;
      SendStack.resetLimit(nSendVectors * 10);
   }


   ///
   /// @brief setter of samplingDepth
   ///
   void setSamplingDepth(
         int inSamplingDepth
         )
   {
      samplingDepth = inSamplingDepth;
   }


}; // class CmapEnumeration


#endif // __CMAP_TEST_ENUM_H__
