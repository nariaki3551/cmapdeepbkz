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

/**@file    cmapLapParaShareDataPool.cpp
 * @brief   Vector Pool.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "cmapLapParaShareDataPool.h"
#include <algorithm>
#include <iterator>
#include <cfloat>
#include <type_traits>
#include "cmapLapParaShareData.h"
#include "cmapLapParaSolverLocalComm.h"


namespace ParaCMapLAP
{


///
/// @brief insert lattice vector into share data pool
/// @param[in] vectorElement lattice vector
///
bool
ShareDataPool::insert(
      std::shared_ptr<VectorElement> vectorElement
      )
{
   startTime = std::chrono::system_clock::now();

   bool inserted = true;
   if( ascendingPool.find(vectorElement) != ascendingPool.end() )
   {
      // there exists a vector in the pool whose norm is equal to the new vector
      return false;
   }
   else if( limitUsageOfPool > 0 )
   {
      // the new vector is not in the pool
      ascendingPool.insert(vectorElement);
      solverTypeConter[static_cast<int>(vectorElement->getSolverType())]++;
      if( limitUsageOfPool > 0 && ascendingPool.size() > limitUsageOfPool )
      {
         auto i = std::prev(ascendingPool.end());
         solverTypeConter[static_cast<int>((*i)->getSolverType())]--;
         ascendingPool.erase(i);
      }
   }

   endTime = std::chrono::system_clock::now();
   insertingMilliTime += static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(endTime-startTime).count())
                         / 1000000.0;
   nCallInsert++;

   return inserted;
}


///
/// get copied k-element from pool
///
int
ShareDataPool::getCopiedElement(
      int k,
      int solverThreadId,
      std::vector<std::shared_ptr<VectorElement>> &vectorElements,   ///< memory should be reserved from caller for k pointers array
      CMapLapParaSolverLocalComm *comm                      ///< communicator used
      )
{
   startTime = std::chrono::system_clock::now();

   int i = 0;
   for( const auto &p : ascendingPool )
   {
      if( i == k ) break;
      if( !(p->isInSentSolverThreadId(solverThreadId)) )
      {
         vectorElements[i++] = p;
      }
   }
   vectorElements.resize(i);

   // To avoid data conflicts with the writeVectorElementsToCheckpointFile()
   if( comm ) comm->lockThread();
   for( int j = 0; j < i; ++j )
   {
      vectorElements[j]->addSentSolverThreadId(solverThreadId);
   }
   if( comm ) comm->unlockThread();

   endTime = std::chrono::system_clock::now();
   extractingMilliTime += static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(endTime-startTime).count())
                          / 100000.0;
   nCallExtract++;

   return i;
}


///
/// erase a solver thread id from set of sent thread ids of all elements
///
void
ShareDataPool::eraseSolverThreadId(
      unsigned int inSolverThreadId,
      CMapLapParaSolverLocalComm *comm    ///< communicator used
      )
{
   // To avoid data conflicts with the writeVectorElementsToCheckpointFile()
   if( comm ) comm->lockThread();
   for( const auto &p : ascendingPool )
   {
      p->eraseSolverThreadId(inSolverThreadId);
   }
   if( comm ) comm->unlockThread();
}


///
/// @brief get quantile sqnorm value
/// @param[in] GH estimated shortest non-zero vector's norm of lattice
///
ShareDataPoolData
ShareDataPool::getPoolData(
      double GH
      )
{
   ShareDataPoolData data;
   data.size = size();
   data.max_value = 0.0;
   data.min_value = DBL_MAX;
   data.nLessThanAFOnePointFive = 0;
   data.nLessThanAFOnePointFour = 0;
   data.nLessThanAFOnePointThree = 0;
   data.nLessThanAFOnePointTwo = 0;
   data.nLessThanAFOnePointOne = 0;
   data.nLessThanAFOnePointZeroSeven = 0;
   data.nLessThanAFOnePointZeroFive = 0;
   data.nLessThanAFOnePointZeroTwo = 0;
   data.nLessThanAFOnePointZeroOne = 0;

   std::vector<double> sqnorms(size(), 0);
   int i = -1;
   for( const auto &p : ascendingPool )
   {
      sqnorms[++i] = p->squaredNorm();
      if( sqnorms[i] > data.max_value){ data.max_value = sqnorms[i]; }
      if( sqnorms[i] < data.min_value){ data.min_value = sqnorms[i]; }
      if( sqnorms[i] < GH * GH * 1.5 * 1.5 )
         data.nLessThanAFOnePointFive++;
      if( sqnorms[i] < GH * GH * 1.4 * 1.4 )
         data.nLessThanAFOnePointFour++;
      if( sqnorms[i] < GH * GH * 1.3 * 1.3 )
         data.nLessThanAFOnePointThree++;
      if( sqnorms[i] < GH * GH * 1.2 * 1.2 )
         data.nLessThanAFOnePointTwo++;
      if( sqnorms[i] < GH * GH * 1.1 * 1.1 )
         data.nLessThanAFOnePointOne++;
      if( sqnorms[i] < GH * GH * 1.07 * 1.07 )
         data.nLessThanAFOnePointZeroSeven++;
      if( sqnorms[i] < GH * GH * 1.05 * 1.05 )
         data.nLessThanAFOnePointZeroFive++;
      if( sqnorms[i] < GH * GH * 1.02 * 1.02 )
         data.nLessThanAFOnePointZeroTwo++;
      if( sqnorms[i] < GH * GH * 1.01 * 1.01 )
         data.nLessThanAFOnePointZeroOne++;
   }
   auto Q1 = sqnorms.size() / 4;
   auto Q2 = sqnorms.size() / 2;
   auto Q3 = Q1 + Q2;
   std::nth_element(sqnorms.begin(),          sqnorms.begin() + Q1, sqnorms.end());
   std::nth_element(sqnorms.begin() + Q1 + 1, sqnorms.begin() + Q2, sqnorms.end());
   std::nth_element(sqnorms.begin() + Q2 + 1, sqnorms.begin() + Q3, sqnorms.end());
   data.Q1_value     = sqnorms[Q1];
   data.median_value = sqnorms[Q2];
   data.Q3_value     = sqnorms[Q3];
   return data;
}


///
/// output vectorElements in ascendingPool
///
const std::string
ShareDataPool::toString(
      )
{
   std::ostringstream s;
   s << std::endl << std::endl;
   s << "===== Share Data Pool (Size " << ascendingPool.size() << ") =====" << std::endl;
   int index = 0;
   for( const auto &p : ascendingPool )
   {
      s << "index " << index << std::endl;
      s << p->toSimpleString();
      s << std::endl;
      index++;
   }
   s << "=====================================" << std::endl;
   return s.str();
}


///
/// output statistics information of pool
///
std::string
ShareDataPool::toStatStringHeader(
      )
{
   return "ShareDataPoolStat"
      ",time"
      ",insertingTime"
      ",nCallInsert"
      ",extractingTime"
      ",nCallExtract"
      ",size"
      ",min_value"
      ",Q1_value"
      ",median_value"
      ",Q3_value"
      ",max_value"
      ",nLessThanAFOnePointFive"
      ",nLessThanAFOnePointFour"
      ",nLessThanAFOnePointThree"
      ",nLessThanAFOnePointTwo"
      ",nLessThanAFOnePointOne"
      ",nLessThanAFOnePointZeroSeven"
      ",nLessThanAFOnePointZeroFive"
      ",nLessThanAFOnePointZeroTwo"
      ",nLessThanAFOnePointZeroOne";
}


///
/// output statistics information of pool
///
std::string
ShareDataPool::toStatString(
      double GH
      )
{
   std::ostringstream sHead, s;
   s << insertingMilliTime / 1000.0
     << "," << nCallInsert
     << "," << extractingMilliTime / 1000.0
     << "," << nCallExtract;
   auto poolData = getPoolData(GH);
   s << "," << poolData.size
     << "," << poolData.min_value
     << "," << poolData.Q1_value
     << "," << poolData.median_value
     << "," << poolData.Q3_value
     << "," << poolData.max_value
     << "," << poolData.nLessThanAFOnePointFive
     << "," << poolData.nLessThanAFOnePointFour
     << "," << poolData.nLessThanAFOnePointThree
     << "," << poolData.nLessThanAFOnePointTwo
     << "," << poolData.nLessThanAFOnePointOne
     << "," << poolData.nLessThanAFOnePointZeroSeven
     << "," << poolData.nLessThanAFOnePointZeroFive
     << "," << poolData.nLessThanAFOnePointZeroTwo
     << "," << poolData.nLessThanAFOnePointZeroOne;
   return s.str();
}

} // namespace ParaCMapLAP
