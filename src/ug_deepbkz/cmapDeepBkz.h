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

/**@file    cmapDeepBkz.h
 * @brief   DeepBkz extension for CMAP-DeepBKZ.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_DEEPBKZ_DEEP_BKZ_H__
#define __CMAP_DEEPBKZ_DEEP_BKZ_H__

#include <memory>

#include "Lattice.h"
#include "DeepBkz.h"
#include "Enumeration.h"
#include "VectorElement.h"
#include "VectorElementPool.h"

#include "cmapDeepLll.h"
#include "cmapDeepBkzParaSolver.h"


template<typename BasisFloat=int, typename GSFloat=double, typename EnumGSFloat=double>
class CmapDeepBkz : public LapTools::DeepBkz<BasisFloat, GSFloat, EnumGSFloat>
{

using VectorElementType = LapTools::VectorElementBase<BasisFloat>;
using VectorElementPtr = std::shared_ptr<VectorElementType>;
using VectorElementQueue = LapTools::VectorElementBasePool<BasisFloat>;
using LatticePtr = std::shared_ptr<LapTools::Lattice<BasisFloat, GSFloat>>;


private:

   LatticePtr  L;             ///< lattice
   Config      config;        ///< hyper parameters
   int         rank;          ///< rank
   int         threadId;      ///< thread Id
   int         verbose;       ///< verbose <= 0: none, 1: light, 2: medium, 3: heavy

   ParaCMapLAP::CMapLapParaSolver *cmapLapParaSolver; ///< communicator to LC
   bool mergeBasisFromLC;                             ///< merge basis from LoadCoordinator
   LapTools::LatticeBasis<BasisFloat> prevTourBasis;  ///< basis when previous tour
   int nSendVectors;                                  ///< number of vectors per one-communication
   VectorElementQueue SendStack;                      ///< stack

   CmapDeepLll<BasisFloat, GSFloat> lllObj;                          ///< lll tool
   LapTools::Enumeration<BasisFloat, GSFloat, EnumGSFloat> enumObj;  ///< enumeration tool

public:

   ///
   /// @brief Constructor
   ///
   CmapDeepBkz(){};


   ///
   /// @brief Constructor
   /// @param[in] inL         Lattice
   /// @param[in] inVerbose   <= 0 : not, 1 : light, 2: heavy
   /// @param[in] inRank      solver-id
   /// @param[in] inThreadId  solver-id
   ///
   CmapDeepBkz(
         LatticePtr inL,
         ParaCMapLAP::CMapLapParaSolver *inCmapLapParaSolver,
         int inRank=1,
         int inThreadId=0,
         int inVerbose=0,
         bool inMergeBasisFromLC=false,
         int inNSendVectors=0
         )
      :
         LapTools::DeepBkz<BasisFloat, GSFloat, EnumGSFloat>(inL, inRank, inThreadId, inVerbose),
         rank(inRank),
         threadId(inThreadId),
         verbose(inVerbose),
         mergeBasisFromLC(inMergeBasisFromLC),
         nSendVectors(inNSendVectors)
   {
      L = inL;
      config = inL->config;
      int lllVerbose = 0;
      bool lllMergeBasisFromLC = false;
      lllObj = CmapDeepLll<BasisFloat, GSFloat>{inL, inCmapLapParaSolver, inRank, inThreadId, lllVerbose, lllMergeBasisFromLC};
      cmapLapParaSolver = inCmapLapParaSolver;
      SendStack = VectorElementQueue(10 * nSendVectors);
      prevTourBasis = L->basis;
      // SendStack = VectorElementQueue(-1);
   }


   ///
   /// @brief preprocess of DeepBkz main loop
   /// @param[out] shouldAbort true if it should abort else false
   ///
   bool preprocess(
         bool &shouldAbort
         );


   ///
   /// @brief communicate with LC
   /// @param[out] shouldAbort true if it should abort else false
   /// @return true if it communicated with LC else false
   ///
   bool communicate(
         bool &shouldAbort
         );


   ///
   /// @brief communicate with LC
   /// @param[out] shouldAbort true if it should abort else false
   ///
   bool communicateInTour(
         bool &shouldAbort
         );


   ///
   /// @brief DeepLll reduction algorithm
   /// @return bool true: normal terminate, false: abnormal one
   ///
   bool deeplll(
         )
   {
      int lllVerbose = lllObj.getVerbose();
      lllObj.setVerbose(verbose);
      bool status = lllObj.deeplll();
      lllObj.setVerbose(lllVerbose);
      return status;
   }


   ///
   /// @brief update bestObjectiveValue
   /// @param[in] sigh line-header character
   ///
   bool updateBestObjectiveValue(
         char sigh
         );


   ///
   /// @breif send SolverState of DeepBkz
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
   }


}; // class CmapDeepBkz


#endif // __CMAP_DEEPBKZ_DEEP_BKZ_H__
