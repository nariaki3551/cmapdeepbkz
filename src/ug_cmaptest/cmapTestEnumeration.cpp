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

/**@file    cmapTestEnumeration.cpp
 * @brief   Enumeration extension for CMAP-TEST.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "cmapTestEnumeration.h"
#include <iostream>
#include "cmapLapParaDef.h"


///
/// @brief post process of search node
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
bool
CmapEnumeration<BasisFloat, GSFloat, EnumGSFloat>::postProcess(
      )
{
   if( nSendVectors > 0 && samplingDepth >= this->begin && this->k <= samplingDepth )
   {
      // sampling
      for( int j = this->k - 1; j > -1; --j )
      {
         for( int i = L->m-1; i > j; --i )
         {
            this->sigma.coeffRef(i,j) = this->sigma.coeff(i+1,j) + this->coeffs.coeff(i) * this->mu.coeff(i,j);
         }
         this->c.coeffRef(j) = -this->sigma.coeff(j+1,j);
         this->coeffs.coeffRef(j) = std::round( this->c.coeff(j) );
      }
      sendVector.setZero(L->n);
      for( int i = this->begin; i <= this->end; ++i )
      {
         for( int j = 0; j < L->n; ++j )
         {
            sendVector.coeffRef(j) += this->coeffs(i) * L->basis.coeff(i, j);
         }
      }
      if( verbose > 2 )
      {
         std::cout
            << "[Rank=" << rank << ", Thread=" << threadId << "] "
            << "sampling(ENUM) "
            << "[ "
            << "k " << this->k
            << "; norm " << std::sqrt(sendVector.squaredNorm())
            << "; approx " << L->approxFactor(std::sqrt(sendVector.squaredNorm())) << " ] "
            // << sendVector.transpose()
            << std::endl;
      }
      if( sendVector.squaredNorm() > MINEPSILON )
      {
         SendStack.insert(VectorElementPtr(new VectorElementType(sendVector)));
         if( verbose > 2 && SendStack.size() % 10 == 1 )
         {
            std::cout
               << "SENDSTACK," << LapTools::Timer::getElapsedTime()
               << "," << this->nSearch
               << SendStack.toStatString(L->GH)
               << std::endl;
         }
      }
      return true;
   }
   return false;
}


///
/// @brief communicate with LC
/// @param[out] shouldAbort true if it should abort else false
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
bool
CmapEnumeration<BasisFloat, GSFloat, EnumGSFloat>::communicate(
      bool& shouldAbort
      )
{
   // check iReceiveMessage
   double currentTime = LapTools::Timer::getElapsedTime();
   if( currentTime < nextIReceiveMessagesTime )
      return true;
   nextIReceiveMessagesTime = currentTime + iReceiveMessagesInterval;

   if( cmapLapParaSolver->iReceiveIsNecessary() )
   {
      cmapLapParaSolver->iReceiveMessages();
   }

   // check whether receive the interrupt message from LC
   if( cmapLapParaSolver->isInterrupted() )
   {
      shouldAbort = true;
      return true;
   }

   // check incumbent solution
   double globalIncumbent;
   ParaCMapLAP::LatticeVector<int> tmpV;
   bool updated = cmapLapParaSolver->checkIfProjNormIsUpdated(
         &globalIncumbent,
         &tmpV);
   if( updated )
   {
      if( globalIncumbent < this->bestObjectiveValue - MINEPSILON )
      {
         LapTools::Enumeration<BasisFloat, GSFloat, EnumGSFloat>::updateBestObjectiveValue('U', globalIncumbent);
         bestVector = tmpV;
      }
   }

   // send status
   if( cmapLapParaSolver->notificationIsNecessary() )
   {
      sendStatus();
   }

   return true;
}


///
/// @brief update bestObjectiveValue
/// @param[in] sigh line-header character
/// @param[in] objectiveValue new objectiveValue
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
bool
CmapEnumeration<BasisFloat, GSFloat, EnumGSFloat>::updateBestObjectiveValue(
      char sigh,
      EnumGSFloat objectiveValue
      )
{
   if( LapTools::Enumeration<BasisFloat, GSFloat, EnumGSFloat>::updateBestObjectiveValue(sigh, objectiveValue) )
   {
      bestVector = LapTools::Enumeration<BasisFloat, GSFloat, EnumGSFloat>::getBestVector().template cast<int>();
      cmapLapParaSolver->sendSolution(bestVector, objectiveValue);
      return true;
   }
   return false;
}


///
/// @breif send SolverState of DeepBkz
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
void
CmapEnumeration<BasisFloat, GSFloat, EnumGSFloat>::sendStatus(
      )
{
   ParaCMapLAP::LatticeVector<int> inCoeffs = this->getCoeffs().template cast<int>();
   cmapLapParaSolver->sendSolverState(
      this->runningTime,
      inCoeffs,
      this->getDepth(),
      std::log(this->enumCost),
      static_cast<double>(std::sqrt(this->bestObjectiveValue)),
      this->nSearch  // numSearchedNodes
      );
}


///
/// instantiation
///
template class CmapEnumeration<int, double, double>;
template class CmapEnumeration<int, long double, double>;
template class CmapEnumeration<int, long double, long double>;
template class CmapEnumeration<long int, double, double>;
template class CmapEnumeration<long int, long double, double>;
template class CmapEnumeration<long int, long double, long double>;
