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

/**@file    cmapTestGaussSieve.h
 * @brief   GaussSieve extension for CMAP-TEST.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_TEST_GAUSS_SIEVE_H__
#define __CMAP_TEST_GAUSS_SIEVE_H__

#include <memory>

#include "GaussSieve.h"
#include "Config.h"
#include "Lattice.h"
#include "VectorElement.h"
#include "VectorElementPool.h"
#include "cmapTestParaSolver.h"


template<typename BasisFloat=int, typename GSFloat=double>
class CmapGaussSieve : public LapTools::GaussSieve<BasisFloat, GSFloat>
{

using VectorElementType = LapTools::VectorElementBase<BasisFloat>;
using VectorElementPtr = std::shared_ptr<LapTools::VectorElementBase<BasisFloat>>;
using VectorElementQueue = LapTools::VectorElementBasePool<BasisFloat>;
using LatticePtr = std::shared_ptr<LapTools::Lattice<BasisFloat, GSFloat>>;

private:

   LatticePtr L;              ///< lattice
   LapTools::Config config;   ///< hyper parameters
   int rank;                  ///< solver-id
   int threadId;              ///< solver-id
   int verbose;               ///< verbose level


protected:

   CMapTestParaSolver *cmapLapParaSolver; ///< communicator to LC

   int   nSendVectors;           ///< number of sending vectors to LC in one communication
   int   nReceiveVectors;        ///< number of receiving vectors to LC in one communication
   VectorElementQueue SendStack; ///< stack of vectors that will be send


public:

   ///
   /// @brief Constructor
   ///
   CmapGaussSieve(){};


   ///
   /// @brief Constructor
   /// @param[in] inL         Lattice
   /// @param[in] inVerbose   <= 0 : not, 1 : light, 2: heavy
   /// @param[in] inRank      solver-id
   /// @param[in] inThreadId  solver-id
   /// @param[in] inNSendVectors number of sent vectors to LC per one communication
   ///
   CmapGaussSieve(
         LatticePtr inL,
         CMapTestParaSolver *inCmapLapParaSolver,
         int inRank=1,
         int inThreadId=0,
         int inVerbose=0
         )
      :
         LapTools::GaussSieve<BasisFloat, GSFloat>(inL, inRank, inThreadId, inVerbose),
         rank(inRank),
         threadId(inThreadId),
         verbose(inVerbose),
         nSendVectors(0),
         nReceiveVectors(0)
   {
      L = inL;
      config = inL->config;
      cmapLapParaSolver = inCmapLapParaSolver;
   }


   ///
   /// @brief set size of list and statck
   /// @param[in] listsize size of List
   /// @param[in] stacksize size of Stack
   /// @param[in] inMaxCollision
   /// @param[in] sendStacksize size of SendStack
   ///
   void init(
         int listsize,
         int stacksize,
         int inMaxCollision=-1,
         int inNSendVectors=1,
         int sendStacksize=-1
       )
   {
      LapTools::GaussSieve<BasisFloat, GSFloat>::init(listsize, stacksize, inMaxCollision);
      if( inNSendVectors > 0 ) nSendVectors = inNSendVectors;
      if( sendStacksize < 0 ) sendStacksize = nSendVectors * 10;
      // SendStack = VectorElementQueue{sendStacksize};
      SendStack = VectorElementQueue(-1);
   }


   ///
   /// @brief update bestObjectiveValue
   /// @param[in] sigh line-header character
   /// @param[in] shortest vector element
   ///
   bool updateBestObjectiveValue(
         char sigh,
         VectorElementPtr v
         );


   ///
   /// @brief communicate with LC
   /// @param[out] shouldAbort true if it should abort else false
   ///
   bool communicate(
         bool &shouldAbort
         );


   ///
   /// @brief post process after reduce
   /// @param[in] v lattice vector has been reduced
   ///
   void postProcessOfReduce(
         VectorElementPtr v
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


   ///
   /// @breif setter of nReceiveVectors
   ///
   void setNReceiveVectors(
         int inNReceiveVectors
         )
   {
      nReceiveVectors = inNReceiveVectors;
   }


}; // class CmapGaussSieve


#endif // __CMAP_TEST_GAUSS_SIEVE_H__
