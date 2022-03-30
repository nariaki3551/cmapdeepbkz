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

/**@file    cmapTestDeepBkz.cpp
 * @brief   DeepBkz extension for CMAP-TEST.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <iostream>
#include "cmapTestDeepBkz.h"
#include "cmapLapParaDef.h"
#include "cmapTestDeepLll.h"


///
/// @brief preprocess of DeepBkz main loop
/// @param[out] shouldAbort true if it should abort else false
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
bool
CmapDeepBkz<BasisFloat, GSFloat, EnumGSFloat>::preprocess(
      bool &shouldAbort
      )
{
   double startTime = LapTools::Timer::getElapsedTime();
   lllObj.deeplll();
   this->runningTime += LapTools::Timer::getElapsedTime() - startTime;
   updateBestObjectiveValue('*');

   if( nSendVectors > 0 )
   {
      bool skip = false;
      for( int i = 0; i < L->m; ++i )
      {
         skip = false;
         for( int j = 0; j < L->m; ++j )
         {
            if( L->basis.row(i) == prevTourBasis.row(j) ){ skip = true; break; }
         }
         if( skip ) continue;
         LapTools::LatticeVector<BasisFloat> vector(L->basis.row(i));
         SendStack.insert(std::make_shared<VectorElementType>(vector));
         if( verbose > 3 )
         {
            std::cout
               << "[Rank=" << rank << ", Thread=" << threadId << "] "
               << "sampling(DeepBKZ) "
               << "[ "
               << "i " << i
               << "; norm " << std::sqrt(vector.squaredNorm())
               << "; approx " << L->approxFactor(std::sqrt(vector.squaredNorm())) << " ] "
               // << sendVector.transpose()
               << std::endl;
         }
      }
      prevTourBasis = L->basis;
      if( verbose > 3 )
      {
         std::cout
            << "SENDSTACK," << LapTools::Timer::getElapsedTime()
            << "," << this->nTour
            << SendStack.toStatString(L->GH)
            << std::endl;
      }
   }
   return true;
}


///
/// @brief communicate with LC
/// @param[out] shouldAbort true if it should abort else false
/// @return true if it communicated with LC else false
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
bool
CmapDeepBkz<BasisFloat, GSFloat, EnumGSFloat>::communicate(
      bool &shouldAbort
      )
{
   double startTime = LapTools::Timer::getElapsedTime();

   // fetch messages from Load Coordinator
   if( cmapLapParaSolver->iReceiveIsNecessary() )
   {
      cmapLapParaSolver->iReceiveMessages();
   }

   // check if it should interrupt
   if( cmapLapParaSolver->isInterrupted() )
   {
      shouldAbort = true;
      return true;
   }

   // check it has received a send-basis request
   if( cmapLapParaSolver->hasSendBasisRequest() )
   {
      ParaCMapLAP::LatticeBasis<int> basis = L->basis.template cast<int>();
      cmapLapParaSolver->sendBasis(basis, L->enumCost());
   }

   // check it has received basis
   if( mergeBasisFromLC && cmapLapParaSolver->hasReceivedBasis() )
   {
      LapTools::LatticeBasis<BasisFloat> receivedSubBasis(cmapLapParaSolver->getReceivedBasis<BasisFloat>());
      LapTools::Lattice<BasisFloat, GSFloat> receivedSubLattice{receivedSubBasis};
      int index = 0;
      if( receivedSubLattice.isMoreReducedThan(*L, index) )
      {
         double startMergeTime = LapTools::Timer::getElapsedTime();
         L->merge(receivedSubLattice);
         this->mergeTime += LapTools::Timer::getElapsedTime() - startMergeTime;
         updateBestObjectiveValue('U');
      }
   }

   // check it should send vectors
   if( cmapLapParaSolver->shareVectorsIsNecessary() )
   {
      // check it has received vectors
      if( cmapLapParaSolver->hasReceivedVectors() )
      {
         ParaCMapLAP::LatticeBasis<BasisFloat> receivedVectors;
         bool received = cmapLapParaSolver->getReceivedVectors<BasisFloat>(receivedVectors);
         if( received )
         {
            int nRows = receivedVectors.rows();
            if( verbose > 3 )
            {
               std::cout << "[Rank=" << rank << ", Thread=" << threadId << "] receive vector; ";
               for( int i = 0; i < nRows; ++i )
                  std::cout << " " << std::sqrt(receivedVectors.row(i).squaredNorm());
               std::cout << std::endl;
            }
            for( int i = 0; i < nRows; ++i )
            {
               LapTools::LatticeVector<BasisFloat> receivedVector = receivedVectors.row(i);
               L->insertMlll(receivedVector, 0, 0, L->m-1);
               updateBestObjectiveValue('U');
            }
         }
      }

      // send some short lattice vectors
      int nRows = std::min(nSendVectors, static_cast<int>(SendStack.size()));
      if( nRows > 0 )
      {
         ParaCMapLAP::LatticeBasis<int> sendVectorsBasis(nRows, L->n);
         for( int i = 0; i < nRows; ++i )
         {
            sendVectorsBasis.row(i) = (SendStack.extract()->getVector()).template cast<int>();
         }
         if( verbose > 3 )
         {
            std::cout
               << "[Rank=" << rank << ", Thread=" << threadId << "] "
               << "send vectors";
            for( int i = 0; i < nRows; ++i )
            {
               std::cout << " " << std::sqrt(sendVectorsBasis.row(i).squaredNorm());
            }
            std::cout << std::endl;
         }
         cmapLapParaSolver->sendVectors(sendVectorsBasis);
      }
   }

   // send lattice vectors request to Load Coordinator
   if( nReceivedVectors > 0)
   {
      cmapLapParaSolver->requestVectors(nReceivedVectors);
   }

   // check to be necessary to send the status
   if( cmapLapParaSolver->notificationIsNecessary() )
   {
      sendStatus();
   }
   cmapLapParaSolver->iReceiveMessages();

   this->runningTime += LapTools::Timer::getElapsedTime() - startTime;
   return true;
}


///
/// @brief communicate with LC
/// @param[out] shouldAbort true if it should abort else false
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
bool
CmapDeepBkz<BasisFloat, GSFloat, EnumGSFloat>::communicateInTour(
      bool &shouldAbort
      )
{
   cmapLapParaSolver->iReceiveMessages();
   // check interrupting
   if( cmapLapParaSolver->isInterrupted() )
   {
      shouldAbort = true;
   }
   return true;
}


///
/// @brief update bestObjectiveValue
/// @param[in] sigh line-header character
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
bool
CmapDeepBkz<BasisFloat, GSFloat, EnumGSFloat>::updateBestObjectiveValue(
      char sigh
      )
{
   if( LapTools::DeepBkz<BasisFloat, GSFloat, EnumGSFloat>::updateBestObjectiveValue(sigh) )
   {
      ParaCMapLAP::LatticeVector<int> u = L->basis.row(0).template cast<int>();
      cmapLapParaSolver->sendSolution(u, this->bestObjectiveValue);
      return true;
   }
   return false;
}


///
/// @breif send SolverState of DeepBkz
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
void
CmapDeepBkz<BasisFloat, GSFloat, EnumGSFloat>::sendStatus(
      )
{
   ParaCMapLAP::LatticeBasis<int> basis = L->basis.template cast<int>();
   cmapLapParaSolver->sendSolverState(
         basis,
         this->blocksize,
         this->nTour,
         this->runningTime
         );
}


///
/// instantiation
///
template class CmapDeepBkz<int, double>;
template class CmapDeepBkz<int, long double>;
template class CmapDeepBkz<long int, long double>;
