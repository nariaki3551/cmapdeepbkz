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

/**@file    cmapLapParaTaskPool.h
 * @brief   CMapLapParaTask Pool.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_TASK_POOL_H__
#define __CMAP_LAP_PARA_TASK_POOL_H__

#include <map>
#include <queue>
#include <cfloat>
#include <cmath>
#include "cmapLapParaTask.h"
#include "cmapLapParaInstance.h"
#include "cmapLapParaComm.h"


namespace ParaCMapLAP
{

static const double eps = 1.0e-12;


///
/// class CMapLapParaTaskSortCriterion
///
class CMapLapParaTaskSortCriterion
{

public:

   ///
   /// ()operator
   /// @return true if n1 < n2
   ///
   virtual bool operator()(
         const std::shared_ptr<CMapLapParaTask>& n1,       ///< CMapLapParaTask pointer n1
         const std::shared_ptr<CMapLapParaTask>& n2        ///< CMapLapParaTask pointer n2
         ) const
   {

      return EPSLT(n1->getEstimatedValue(),n2->getEstimatedValue(),eps);

   }
};


///
/// class CMapLapParaTaskPool
///
class CMapLapParaTaskPool
{

protected:

   size_t        maxUsageOfPool;     ///< maximum usage of this pool

public:

   ///
   /// constructor
   ///
   CMapLapParaTaskPool(
         )
         : maxUsageOfPool(0)
   {
   }

   ///
   /// destructor
   ///
   virtual ~CMapLapParaTaskPool(
         )
   {
   }

   ///
   /// insert CMapLapParaTask to this pool
   ///
   virtual void insert(
         std::shared_ptr<CMapLapParaTask> task  ///< pointer to CMapLapParaTask object
         ) = 0;

   ///
   /// check if this pool is empty or not
   /// @return true if this pool is empty
   ///
   virtual bool empty(
         ) = 0;

   ///
   /// extract a CMapLapParaTask object from this pool
   /// @return pointer to CMapLapParaTask object extracted
   ///
   virtual std::shared_ptr<CMapLapParaTask> extractTask(
         ) = 0;

   ///
   /// get best dual bound value of CMapLapParaTask object in this pool
   /// @return best dual bound value
   ///
   virtual double getBestEstimatedValue(
         ) = 0;

   ///
   /// get number of CMapLapParaTasks in this pool
   /// @return number of CMapLapParaTasks
   ///
   virtual size_t getNumOfTasks(
         ) = 0;

   ///
   /// stringfy this object
   /// @return string which shows inside of this object as string
   ///
   virtual const std::string toString(
         ) = 0;

   ///
   /// get maximum usage of this pool
   /// @return the maximum number of CMapLapParaTasks that were in this pool
   ///
   virtual size_t getMaxUsageOfPool(
         )
   {
      return maxUsageOfPool;
   }

};



#ifdef UG_WITH_ZLIB

///
/// class CMapLapParaTaskQueue
/// @note this container is used when checkpointing
///
class CMapLapParaTaskQueue
{
private:

   std::queue< std::shared_ptr<CMapLapParaTask> > paraTaskQueue;

public:
   CMapLapParaTaskQueue(){}
   virtual ~CMapLapParaTaskQueue()
   {
      while( !paraTaskQueue.empty() )
      {
         paraTaskQueue.pop();
      }
   }

   virtual int size(
         )
   {
      return paraTaskQueue.size();
   }

   ///
   /// add MapsvpParaTask in the queue
   ///
   virtual void push(
         std::shared_ptr<CMapLapParaTask> cmapLapParaTask
         )
   {
      paraTaskQueue.push(cmapLapParaTask);
   }

   ///
   /// write CMapLapParaTasks to checkpoint file
   /// @return number of CMapLapParaTasks written
   ///
   virtual int writeTasksToCheckPointFile(
         gzstream::ogzstream &out,  ///< gzstream for output
         bool active                ///< True: tasks in paraTaskQueue is active task, False: else
         )
   {
      std::shared_ptr<CMapLapParaTask> cmapLapParaTask;
      while( !paraTaskQueue.empty() )
      {
         cmapLapParaTask = paraTaskQueue.front();
         cmapLapParaTask->write(out, active);
         paraTaskQueue.pop();
      }
      return paraTaskQueue.size();
   }
};
#endif


///
/// class CMapLapParaTaskPoolInAscendingOrder
/// @note only minimization pool was written, since all problem is converted to minimization problem inside of UG solvers
///
class CMapLapParaTaskPoolInAscendingOrder : virtual public CMapLapParaTaskPool
{

   std::multimap<
      std::shared_ptr<CMapLapParaTask>,
      std::shared_ptr<CMapLapParaTask>,
      CMapLapParaTaskSortCriterion > ascendingPool;    ///< asnding pool

public:

   ///
   /// constructor
   ///
   CMapLapParaTaskPoolInAscendingOrder(
         )
         : CMapLapParaTaskPool()
   {
   }

   ///
   /// destructor
   ///
   virtual ~CMapLapParaTaskPoolInAscendingOrder(
         )
   {
   }

   ///
   /// insert a CMapLapParaTask object to this pool
   ///
   virtual void insert(
         std::shared_ptr<CMapLapParaTask> paraTask          ///< pointer to CMapLapParaTask object to insert
         )
   {
      ascendingPool.insert(std::make_pair(paraTask, paraTask));
      if( maxUsageOfPool < ascendingPool.size() )
      {
         maxUsageOfPool = ascendingPool.size();
      }
   }

   ///
   /// check if this pool is empty or not
   /// @return true if this pool is empty
   ///
   virtual bool empty(
         )
   {
      return ( ascendingPool.size() == 0 );
   }

   ///
   /// extract a CMapLapParaTask from this pool
   /// @return pointer to CMapLapParaTask object extracted
   ///
   virtual std::shared_ptr<CMapLapParaTask> extractTask(
         )
   {
      std::shared_ptr<CMapLapParaTask> extracted(nullptr);
      if( ascendingPool.size() > 0 )
      {
         auto p = ascendingPool.begin();
         extracted = p->second;
         ascendingPool.erase(p);
      }
      return extracted;
   }

   ///
   /// get best dual bound value of CMapLapParaTask object in this pool
   /// @return best dual bound value
   ///
   virtual double getBestEstimatedValue(
         )
   {
      if( ascendingPool.size() > 0 )
      {
         auto p = ascendingPool.begin();
         return p->second->getEstimatedValue();
      }
      else
      {
         return DBL_MAX;  // no tasks exist
      }
   }

   ///
   /// get number of CMapLapParaTasks in this pool
   /// @return number of CMapLapParaTasks
   ///
   virtual size_t getNumOfTasks(
         )
   {
      return ascendingPool.size();
   }

   ///
   /// stringfy this object
   /// @return string which shows inside of this object as string
   ///
   virtual const std::string toString(
         )
   {
      std::ostringstream s;
      for( const auto &p : ascendingPool )
      {
         s << p.second->toString();
      }
      return s.str();
   }

#ifdef UG_WITH_ZLIB

   ///
   /// get copied task objects
   /// @return the CMapLapParaTaskQueue pointer
   ///
   virtual std::shared_ptr<CMapLapParaTaskQueue> getParaTaskQueue(
         UG::ParaComm *comm      ///< communicator used
         )
   {
      std::shared_ptr<CMapLapParaTaskQueue> cmapLapParaTaskQueue(new CMapLapParaTaskQueue());
      for( const auto &p : ascendingPool )
      {
         std::shared_ptr<CMapLapParaTask> cmapLapParaTask(
               dynamic_cast<CMapLapParaTask *>(p.second->clone(comm))
               );
         cmapLapParaTaskQueue->push(cmapLapParaTask);
      }
      return cmapLapParaTaskQueue;
   }

   ///
   /// write CMapLapParaTasks to checkpoint file
   /// @return number of CMapLapParaTasks written
   ///
   virtual int writeTasksToCheckPointFile(
         gzstream::ogzstream &out       ///< gzstream for output
         )
   {
      int n = 0;
      for( const auto &p : ascendingPool )
      {
         p.second->write(out, 0);  // 0 means inactive
         n++;
      }
      return n;
   }

   ///
   /// read CMapLapParaTasks to checkpoint file
   /// @return number of CMapLapParaTasks read
   ///
   virtual int readTasksFromCheckPointFile(
         UG::ParaComm *comm,
         gzstream::igzstream &in
         )
   {
      DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, comm);

      int n = 0;
      std::shared_ptr<CMapLapParaTask> cmapLapParaTask(
            dynamic_cast<CMapLapParaTask *>(cmapLapParaComm->createParaTask())
            );
      while( cmapLapParaTask->read(comm, in) )
      {
         insert(cmapLapParaTask);
         cmapLapParaTask = std::shared_ptr<CMapLapParaTask>(
               dynamic_cast<CMapLapParaTask *>(cmapLapParaComm->createParaTask())
               );
         n++;
      }
      return n;
   }

#endif

};

} // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_TASK_POOL_H__
