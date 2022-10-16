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

/**@file    cmapTestDeepLll.h
 * @brief   DeepLll extension for CMAP-TEST.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_TEST_DEEP_LLL_H__
#define __CMAP_TEST_DEEP_LLL_H__


#include "DeepLll.h"
#include <memory>
#include "Config.h"
#include "Def.h"
#include "Lattice.h"
#include "cmapTestParaSolver.h"


template<typename BasisFloat=int, typename GSFloat=double>
class CmapDeepLll : public LapTools::DeepLll<BasisFloat, GSFloat>
{

using LatticePtr = std::shared_ptr<LapTools::Lattice<BasisFloat, GSFloat>>;


private:

   LatticePtr L;              ///< lattice
   LapTools::Config config;   ///< hyper parameters
   int rank;                  ///< rank
   int threadId;              ///< thread Id
   int verbose;               ///< verbose <= 0: none, 1: light, 2: medium, 3: heavy

   double bestObjectiveValue; ///< bestObjectiveValue

   CMapTestParaSolver *cmapLapParaSolver; ///< communicator to LC

   bool mergeBasisFromLC;           ///< merge basis from LoadCoordinator
   int iReceiveMessagesCounter;     ///< counter of iReceive message
   int iReceiveMessagesInterval;    ///< interval of iReceive message
   int nReceivedVectors;            ///< number of receiving vectors in one communication
   LapTools::LatticeVector<int> v;  ///< allocated memory for bestVector


public:

   ///
   /// @brief Constructor
   ///
   CmapDeepLll(){};


   ///
   /// @brief Constructor
   /// @param[in] inL         Lattice
   /// @param[in] inVerbose   <= 0 : not, 1 : light, 2: heavy
   /// @param[in] inRank      solver-id
   /// @param[in] inThreadId  solver-id
   /// @param[in] inMergeIncumbentFromLC  if it is true, then merge incumbent vector from LoadCoordinator
   /// @param[in] inMergeBasisFromLC if it is true, then merge basis from LoadCoordinator
   /// @param[in] inIReceiveMessagesInterval number of interval loop of communicaion with LoadCoordinator
   ///
   CmapDeepLll(
         LatticePtr inL,
         CMapTestParaSolver *inCmapLapParaSolver,
         int inRank=1,
         int inThreadId=0,
         int inVerbose=0,
         bool inMergeIncumbentFromLC=false,
         int inMergeBasisFromLC=false,
         int inIReceiveMessagesInterval=10000
         )
      :
         LapTools::DeepLll<BasisFloat, GSFloat>(inL, inRank, inThreadId, inVerbose),
         rank(inRank),
         threadId(inThreadId),
         verbose(inVerbose),
         mergeBasisFromLC(inMergeBasisFromLC),
         iReceiveMessagesCounter(0),
         iReceiveMessagesInterval(inIReceiveMessagesInterval),
         nReceivedVectors(0)
   {
      L = inL;
      config = inL->config;
      bestObjectiveValue = inL->B(0);
      cmapLapParaSolver = inCmapLapParaSolver;
   }


   ///
   /// @brief communicate with LC
   /// @param[out] shouldAbort true if it should abort else false
   /// @return true if it communicated with LC else false
   ///
   bool communicate(
         bool &shouldAbort
         );


   ///
   /// @brief update bestObjectiveValue
   /// @param[in] sigh line-header character
   ///
   bool updateBestObjectiveValue(
         char sigh
         );


   ///
   /// @breif send SolverState of DeepLll
   ///
   void sendStatus(
         );


   ///
   /// @breif setter of nReceivedVectors
   /// @param[in] inNReceivedVectors number of receiving vectors in one communication
   ///
   void setNReceivedVectors(
         int inNReceivedVectors
         )
   {
      nReceivedVectors = inNReceivedVectors;
   }


}; // class CmapDeepLll


#endif // __CMAP_TEST_DEEP_LLL_H__
