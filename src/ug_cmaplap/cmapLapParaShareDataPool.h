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

/**@file    cmapLapParaShareDataPool.h
 * @brief   Base class for Share Data Pool.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_SHARE_DATA_POOL_H__
#define __CMAP_LAP_PARA_SHARE_DATA_POOL_H__


#include <queue>
#include <vector>
#include <set>
#include <chrono>
#include <memory>
#include <assert.h>
#include "cmapLapParaShareData.h"
#include "cmapLapParaSolverLocalComm.h"
#include "cmapLapParaDef.h"
namespace UG { class ParaComm; }


namespace ParaCMapLAP
{


#ifdef UG_WITH_ZLIB


///
/// @class VectorElementQueue
/// @note this container is used when checkpointing
///
class VectorElementQueue
{

private:

   std::queue< std::shared_ptr<VectorElement> > vectorElementQueue;


public:
   VectorElementQueue(){};
   virtual ~VectorElementQueue(){};

   virtual int size()
   {
      return vectorElementQueue.size();
   }

   ///
   /// add vectorElement in the queue
   ///
   virtual void push(
         std::shared_ptr<VectorElement> vectorElement
         )
   {
      vectorElementQueue.push(vectorElement);
   }

   ///
   /// write vectorElements to checkpoint file
   /// @return number of vectorElements written
   ///
   virtual int writeVectorElementsToCheckpointFile(
         gzstream::ogzstream &out,  ///< gzstream for output
         int nWrite                 ///< -1: write all elements
         )
   {
      int n = 0;
      while( !vectorElementQueue.empty() )
      {
         if( nWrite > 0 && n >= nWrite ) break;
         vectorElementQueue.front()->write(out);
         vectorElementQueue.pop();
         n++;
      }
      return n;
   }
};
#endif



class VectorElementSortCriterion
{
public:
   virtual bool operator()(const std::shared_ptr<VectorElement>& n1, const std::shared_ptr<VectorElement>& n2) const
   {
      //return (n1->squaredNorm() < n2->squaredNorm());
      if( n1->squaredNorm() < n2->squaredNorm() ) return true;
      if( n1->squaredNorm() > n2->squaredNorm() ) return false;
      auto a = n1->getVector();
      auto b = n2->getVector();
      for( auto i = 0; i < a.size(); ++i )
      {
         if( a[i] < b[i] ) return true;
         if( a[i] > b[i] ) return false;
      }
      return false;
   }
};



///
/// @struct ShareDataPoolData
/// @brief Data Structure for statistics information of ShareDataPool
///
struct ShareDataPoolData
{
   size_t size;                           ///< size of pool
   double min_value;                      ///< minimum value
   double Q1_value;                       ///< first quantile value
   double median_value;                   ///< median value
   double Q3_value;                       ///< third quantile value
   double max_value;                      ///< maximum value
   long int nLessThanAFOnePointFive;      ///< number of vectors whose norm < 1.5
   long int nLessThanAFOnePointFour;      ///< number of vectors whose norm < 1.4
   long int nLessThanAFOnePointThree;     ///< number of vectors whose norm < 1.3
   long int nLessThanAFOnePointTwo;       ///< number of vectors whose norm < 1.2
   long int nLessThanAFOnePointOne;       ///< number of vectors whose norm < 1.1
   long int nLessThanAFOnePointZeroSeven; ///< number of vectors whose norm < 1.07
   long int nLessThanAFOnePointZeroFive;  ///< number of vectors whose norm < 1.05
   long int nLessThanAFOnePointZeroTwo;   ///< number of vectors whose norm < 1.02
   long int nLessThanAFOnePointZeroOne;   ///< number of vectors whose norm < 1.01
};



///
/// @class ShareDataPool
/// @brief Pool of lattice vector as Share-Data
///
class ShareDataPool
{

protected:

   size_t   limitUsageOfPool; // limit usage of this pool, if negative, this is ignore
   std::set<std::shared_ptr<VectorElement>, VectorElementSortCriterion> ascendingPool;
   std::vector<int> solverTypeConter;  // counter number of vectors in vector pool

   /// for statistics
   double   insertingMilliTime;     // total process time of insert
   double   extractingMilliTime;    // total process time of extract
   int      nCallInsert;            // number of call insert()
   int      nCallExtract;           // number of call getCopiedElement()
   std::chrono::system_clock::time_point startTime, endTime;


public:
   ///
   /// constructor
   ///
   ShareDataPool(
         int limit
         ) : limitUsageOfPool(limit),
             insertingMilliTime(0),
             extractingMilliTime(0),
             nCallInsert(0),
             nCallExtract(0)
   {
      solverTypeConter = std::vector<int>(4);
   };

   ///
   /// deconstructor
   ///
   virtual ~ShareDataPool(
         )
   {
   }

   ///
   /// @brief insert lattice vector into share data pool
   /// @param[in] vectorElement lattice vector
   ///
   virtual bool insert(
         std::shared_ptr<VectorElement> vectorElement
         );

   ///
   /// @return true if this pool is empty else false
   ///
   virtual bool empty(
         )
   {
      return ( ascendingPool.size() == 0 );
   }

   ///
   /// @return size of pool
   ///
   virtual size_t size(
         )
   {
      return ascendingPool.size();
   }

   ///
   /// extract VectorElement
   /// This function should be called after checking pool is not empty.
   ///
   virtual std::shared_ptr<VectorElement> extractVectorElement(
         )
   {
      assert( ascendingPool.size() > 0 );
      auto ve = *ascendingPool.begin();
      ascendingPool.erase(ascendingPool.begin());
      return ve;
   }

   ///
   /// get begin of iterator
   ///
   virtual std::set<std::shared_ptr<VectorElement>, VectorElementSortCriterion >::iterator getIterBegin(
         )
   {
      return ascendingPool.begin();
   }

   ///
   /// get end of iterator
   ///
   virtual std::set<std::shared_ptr<VectorElement>, VectorElementSortCriterion >::iterator getIterEnd(
         )
   {
      return ascendingPool.end();
   }

   ///
   /// erase vectorElement from ascendingPool
   ///
   virtual void erase(
         std::set<std::shared_ptr<VectorElement>, VectorElementSortCriterion >::iterator& p
         )
   {
      ascendingPool.erase(p++);
   }

   ///
   /// get copied k-element from pool
   ///
   virtual int getCopiedElement(
         int k,
         int solverThreadId,
         std::vector<std::shared_ptr<VectorElement>> &vectorElements,   ///< memory should be reserved from caller for k pointers array
         CMapLapParaSolverLocalComm *comm                      ///< communicator used
         );

   ///
   /// erase a solver thread id from set of sent thread ids of all elements
   ///
   virtual void eraseSolverThreadId(
         unsigned int inSolverThreadId,
         CMapLapParaSolverLocalComm *comm    ///< communicator used
         );

   ///
   /// creates subBasis using vectors in pool
   /// @return matrix whose #row = ascendingPool.size() and #col = dimension
   ///
   virtual std::shared_ptr<LatticeBasis<int>> createSubBasis(
         )
   {
      int n = (*ascendingPool.begin())->size();
      int m = ascendingPool.size();
      auto subBasis = std::make_shared<LatticeBasis<int>>(m, n);

      int row = 0;
      for( const auto &p : ascendingPool )
      {
         subBasis->row(row++) = LatticeVector<int>(p->getVector());
      }
      return subBasis;
   }

#ifdef UG_WITH_ZLIB
   ///
   /// get copied vector elements
   /// @return the VectorElementQueue pointer
   ///
   virtual VectorElementQueue *getVectorElementQueue(
         UG::ParaComm *comm,     ///< communicator used
         int nWrite=-1           ///< -1: write all elements
         )
   {
      int n = 0;
      VectorElementQueue *vectorElementQueue = new VectorElementQueue();
      for(const auto &p : ascendingPool )
      {
         if( nWrite > 0 && n >= nWrite ) break;
         vectorElementQueue->push(p);
         n++;
      }
      return vectorElementQueue;
   }

   ///
   /// write VectorElements in pool
   /// @return the number of vector elements witten
   ///
   virtual int writeVectorElementsToCheckpointFile(
         gzstream::ogzstream &out,
         CMapLapParaSolverLocalComm *comm,   ///< communicator used
         int nWrite                          ///< -1: write all elements
         )
   {
      int n = 0;
      if( comm ) comm->lockThread();
      for( const auto &p : ascendingPool )
      {
         if( nWrite > 0 && n >= nWrite ) break;
         p->write(out);
         n++;
      }
      if( comm ) comm->unlockThread();
      return n;
   }

   ///
   /// read VectorElements from Checkpoint file
   /// return the number of vector elements inserted
   ///
   virtual int readVectorElementsFromCheckpointFile(
         UG::ParaComm *comm,
         gzstream::igzstream &in,
         bool clearSolverThreadIds=false
         )
   {
      int n = 0;
      auto vectorElement = std::make_shared<VectorElement>();
      while( vectorElement->read(comm, in) )
      {
         if( clearSolverThreadIds )
         {
            vectorElement->clearSolverThreadIds();
         }
         insert(vectorElement);
         n++;
         vectorElement = std::make_shared<VectorElement>();
      }
      return n;
   }

#endif // UG_WITH_ZLIB


   ///
   /// @brief get quantile sqnorm value
   /// @param[in] GH estimated shortest non-zero vector's norm of lattice
   ///
   virtual ShareDataPoolData getPoolData(
         double GH
         );

   ///
   /// write Vector info in pool
   ///
   virtual void writeVectorNormsToOstream(
         std::ostream *os
         )
   {
      for( const auto &p : ascendingPool )
      {
         *os << p->squaredNorm() << "(" << static_cast<int>(p->getSolverType()) << ") ";
      }
      // *os << "; ";
      // *os << solverTypeConter[static_cast<int>(DeepBkz)] << " "
      //     << solverTypeConter[static_cast<int>(Enum)]    << " "
      //     << solverTypeConter[static_cast<int>(Sieve)]   << " "
      //     << solverTypeConter[static_cast<int>(Undefined)];
      *os << std::endl;
   }

   ///
   /// output vectorElements in ascendingPool
   ///
   virtual const std::string toString(
         );

   ///
   /// output statistics information of pool
   ///
   virtual std::string toStatStringHeader(
         );

   ///
   /// output statistics information of pool
   ///
   virtual std::string toStatString(
         double GH
         );

   };

} // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_SHARE_DATA_POOL_H__
