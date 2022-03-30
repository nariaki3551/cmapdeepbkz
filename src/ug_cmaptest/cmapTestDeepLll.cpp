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

/**@file    cmapTestDeepLll.cpp
 * @brief   DeepLll extension for CMAP-TEST.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <iostream>
#include "cmapTestDeepLll.h"
#include "cmapLapParaDef.h"
#include "cmapLapParaSolver.h"


///
/// @brief communicate with LC
/// @param[out] shouldAbort true if it should abort else false
/// @return true if it communicated with LC else false
///
template<typename BasisFloat, typename GSFloat>
bool
CmapDeepLll<BasisFloat, GSFloat>::communicate(
      bool &shouldAbort
      )
{
   // communicate with LC
   ++iReceiveMessagesCounter;
   if( iReceiveMessagesCounter < iReceiveMessagesInterval )
   {
      return false;
   }

   // fetch messages from Load Coordinator
   if( cmapLapParaSolver->iReceiveIsNecessary() )
   {
      cmapLapParaSolver->iReceiveMessages();
   }
   iReceiveMessagesCounter = 0;

   // check if it should interrupt
   if( cmapLapParaSolver->isInterrupted() )
   {
      shouldAbort = true;
      return true;
   }

   // check it has received a send-basis request
   if( mergeBasisFromLC && cmapLapParaSolver->hasReceivedBasis() )
   {
      LapTools::LatticeBasis<BasisFloat> receivedSubBasis = cmapLapParaSolver->getReceivedBasis<BasisFloat>();
      LapTools::Lattice<BasisFloat, GSFloat> receivedSubLattice(receivedSubBasis);
      int index = 0;
      if( receivedSubLattice.isMoreReducedThan(*L, index) )
      {
         L->merge(receivedSubLattice);
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
   return true;
}


///
/// @brief update bestObjectiveValue
/// @param[in] sigh line-header character
///
template<typename BasisFloat, typename GSFloat>
bool
CmapDeepLll<BasisFloat, GSFloat>::updateBestObjectiveValue(
      char sigh
      )
{
   if( LapTools::DeepLll<BasisFloat, GSFloat>::updateBestObjectiveValue(sigh) )
   {
      ParaCMapLAP::LatticeVector<int> receivedVector = L->basis.row(0).template cast<int>();
      cmapLapParaSolver->sendSolution(receivedVector, bestObjectiveValue);
      return true;
   }
   return false;
}


///
/// @breif send SolverState of DeepBkz
///
template<typename BasisFloat, typename GSFloat>
void
CmapDeepLll<BasisFloat, GSFloat>::sendStatus(
      )
{
   ParaCMapLAP::LatticeBasis<int> basis = L->basis.template cast<int>();
   cmapLapParaSolver->sendSolverState(
         basis,
         1, // blocksize
         0, // nTour
         0  ///< taskElapsedTiime
         );
}


///
/// instantiation
///
template class CmapDeepLll<int, double>;
template class CmapDeepLll<int, long double>;
template class CmapDeepLll<long int, long double>;
