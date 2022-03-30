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

/**@file    cmapLapParaSolverPool.cpp
 * @brief   Solver pool.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <iostream>
#include <cassert>
#include <unordered_map>
#include "ug/paraTagDef.h"
#include "cmapLapParaComm.h"
#include "cmapLapParaTask.h"
#include "cmapLapParaSolverPool.h"
#include "cmapLapParaTagDef.h"
namespace ParaCMapLAP { class CMapLapParaTaskPool; }


namespace ParaCMapLAP
{


///
/// resize selection heap
///
void
SelectionHeap::resize(
      std::size_t size                               ///< new heap size
      )
{
   CMapLapParaSolverPoolElementPtr *tempHeap =
           new CMapLapParaSolverPoolElementPtr[size+1]; /// index 0 is a sentinel
   std::memcpy(tempHeap, heap, (unsigned long int)(sizeof(CMapLapParaSolverPoolElementPtr) * (maxHeapSize+1)));
   delete[] heap;
   heap = tempHeap;
   maxHeapSize = size;
}

///
/// insert BbParaSolverPoolElementPtr to selection heap
///
SelectionHeap::ResultOfInsert
SelectionHeap::insert(
      CMapLapParaSolverPoolElementPtr inSolver      ///< BbParaSolverPoolElementPtr to be inserted
      )
{
   if( heapSize >= maxHeapSize ) return FAILED_BY_FULL;
   heap[++heapSize] = inSolver;
   upHeap(heapSize);
   return SUCCEEDED;
}

/// remove the top priority element form selection heap
CMapLapParaSolverPoolElementPtr
SelectionHeap::remove(
      )
{
   CMapLapParaSolverPoolElementPtr solver = heap[1];

   heap[1] = heap[heapSize];
   heapSize--;
   downHeap(1);
   if( heapSize == 0 )
   {
      heap[1] = 0;
   }
   solver->setSelectionHeapElement(0);
   return solver;
}

///
/// stringfy selection hea
///
const std::string
SelectionHeap::toString(
      )
{
   std::ostringstream os;
   os << "--- selection heap ---" << std::endl;
   os << "maxHeapSize: " << maxHeapSize << std::endl;
   os << "heapSize: " << heapSize << std::endl;
   for( std::size_t i = 1; i <= heapSize; i++ )
   {
      os << "heap[" << i << "]->rank: "
         << heap[i]->getRank() << std::endl
         << "heap[" << i << "]->status: "
         << static_cast<int>(heap[i]->getStatus()) << std::endl
         << "heap[" << i << "]->bestBound: "
         << heap[i]->getInterruptPriorityValue() << std::endl;
   }
   return os.str();
}



///
/// update dual bound value of the solver in heap
///
void
AscendingSelectionHeap::updateInterruptPriorityValue(
      CMapLapParaSolverPoolElementPtr solver,           ///< pointer to solver pool element whose dual bound is updated
      double newValue                                  ///< new interrupt priority value
      )
{
   int pos = solver->getSelectionHeapElement() - heap;
   if( solver->getInterruptPriorityValue() > newValue )
   {
      solver->setInterruptPriorityValue(newValue);
      upHeap(pos);
   }
   else
   {
      solver->setInterruptPriorityValue(newValue);
      downHeap(pos);
   }
}

///
/// delete BbParaSolverPoolElement
///
void
AscendingSelectionHeap::deleteElement(
      CMapLapParaSolverPoolElementPtr solver        ///< BbParaSolverPoolElement to be deleted
      )
{
   std::size_t pos = (solver->getSelectionHeapElement()) - heap;

   if( pos == heapSize )
   {
      /* no need to rearrange heap element */
      heap[heapSize--] = 0;
   }
   else
   {
      if( heap[pos]->getInterruptPriorityValue() > heap[heapSize]->getInterruptPriorityValue() )
      {
         heap[pos] = heap[heapSize];
         heap[heapSize--] = 0;
         upHeap(pos);
      }
      else
      {
         heap[pos] = heap[heapSize];
         heap[heapSize--] = 0;
         downHeap(pos);
      }
   }
   solver->setSelectionHeapElement(0);
}

///
/// up heap
///
void
AscendingSelectionHeap::upHeap(
   std::size_t pos                                        ///< up heap this position element
){
   CMapLapParaSolverPoolElementPtr she;

   she = heap[pos];
   heap[0] = NULL;

   while ( heap[pos/2] != NULL )
   {
      heap[pos] = heap[pos/2];
      heap[pos]->setSelectionHeapElement(&heap[pos]);
      pos = pos/2;
   }
   heap[pos] = she;
   heap[pos]->setSelectionHeapElement(&heap[pos]);
}

///
/// down heap
///
void
AscendingSelectionHeap::downHeap(
   std::size_t pos                                       ///< down heap this position element
){
   std::size_t j;
   CMapLapParaSolverPoolElementPtr she;

   she = heap[pos];
   while ( pos <= (heapSize/2) )
   {
      j = pos + pos;
      if( j < heapSize &&
            ( heap[j]->getInterruptPriorityValue()
                  > heap[j+1]->getInterruptPriorityValue() ) ) j++;
      if( she->getInterruptPriorityValue()
            < heap[j]->getInterruptPriorityValue()  )
         break;
      heap[pos] = heap[j];
      heap[pos]->setSelectionHeapElement(&heap[pos]);
      pos = j;
   }
   heap[pos] = she;
   heap[pos]->setSelectionHeapElement(&heap[pos]);
}


///
/// get an inactive Solver rank
/// @return rank of inactivate Solver
///
int
CMapLapParaSolverPool::getInactiveSolverRankForThreadTask(
      int &threadId           ///< return thread ID
      )
{
   threadId = -1;
   auto p = inactiveSolvers.begin();
   while( p != inactiveSolvers.end() && p->second->isReservedForSieve()  )
   {
      p++;
   }
   if( p == inactiveSolvers.end() )
   {
      return -1;   // no inactive Solver
   }
   threadId = p->second->getThreadId();
   return p->second->getRank();
}


///
/// get an inactive Solver rank
/// @return rank of inactivate Solver
///
int
CMapLapParaSolverPool::getInactiveSolverRankForRankTask(
      int &threadId           ///< return thread ID
      )
{
   threadId = -1;
   auto p = inactiveSolvers.begin();
   while( p != inactiveSolvers.end()
         && !p->second->isReservedForSieve()
         && p->second->getThreadId() != 0 )
   {
      p++;
   }
   if( p == inactiveSolvers.end() )
   {
      return -1;   // no inactive Solver
   }
   threadId = p->second->getThreadId();
   return p->second->getRank();
}


///
/// this function is for restart, asssume that an inactive node exists
/// @return true if the reservation of node for Sieve is succeeded else false
///
bool
CMapLapParaSolverPool::reserveForSieve(
      )
{
   auto p = inactiveSolvers.begin();
   while( p != inactiveSolvers.end()
         && p->second->getThreadId() != 0 )
   {
      p++;
   }
   if( p == inactiveSolvers.end() )
   {
      return false;   // can not reserve for Sieve Solver
   }
   int rank = p->second->getRank();
   assert( p->second->getThreadId() == 0 );

   for( unsigned int threadId = 0; threadId < nThreadsPerRank; threadId++ )
   {
      unsigned int solverThreadId = getSolverThreadId(rank, threadId);
      if( !pool[solverThreadId]->isActive() )
      {
         pool[solverThreadId]->reserveForSieve();
      }
      else
      {
         return false;
      }
   }
   return true;
}


///
/// activate a Solver with specified CMapLapParaTask which is sent within this Solver pool
/// @return rank of Solver which is activated
///
int
CMapLapParaSolverPool::activateSolver(
      SolverType           solverType,             ///< solver type
      std::shared_ptr<CMapLapParaTask> paraTask,   ///< pointer to CMapLapParaTask object to be sent to a Solver
      int                  &threadId               ///< thread Id of solver activated (return value)
      )
{
   threadId = -1;
   int rank = -1;
   unsigned int solverThreadId = 0;

   if( solverType == Sieve )
   {
      rank = getInactiveSolverRankForRankTask(threadId); // threadId is set, too
   }
   else
   {
      rank = getInactiveSolverRankForThreadTask(threadId); // threadId is set, too
   }

   if( rank >= 0 )
   {
      solverThreadId = getSolverThreadId(rank, threadId);

      if( pool[( solverThreadId )]->getStatus() == UG::Inactive &&
            ( ( ( solverType == DeepBkz || solverType == Enum ) &&
               !pool[( solverThreadId )]->isReservedForSieve() ) ||
              ( ( solverType == Sieve ) &&
               ( pool[( solverThreadId )]->isReservedForSieve() || nSolvers == 1) ) ) )
      {
         auto p = inactiveSolvers.find(solverThreadId);
         if( p != inactiveSolvers.end() )
         {
            if( p->second->getRank() != rank ||
                  pool[( solverThreadId )]->getRank() != rank ||
                  pool[( solverThreadId )]->getThreadId() != threadId )
            {
               THROW_LOGICAL_ERROR4("Invalid rank. Rank = ", rank, ", or invalid thread Id = ", threadId);
            }
            inactiveSolvers.erase(p);
         }
         else
         {
            THROW_LOGICAL_ERROR4("Invalid rank. Rank = ", rank, ", or invalid thread Id = ", threadId);
         }
      }
      else
      {
         if( solverType == Sieve )
         {
            // cannot activate since it is not reserved
            return -1;
         }
         else
         {
            THROW_LOGICAL_ERROR4("Invalid rank. Rank = ", rank, ", or invalid thread Id = ", threadId);
         }
      }
   }
   else
   {
      return -1;  // no inactive Solver without reserved for Sieve
   }

   pool[( solverThreadId )]->resetReservedForSieve();

   switch( solverType )
   {
   case DeepBkz:
      activeDeepBkzSolvers.insert(std::make_pair(solverThreadId,pool[( solverThreadId )]));
      deepBkzSelectionHeap->insert(pool[( solverThreadId )]);
      nDeepBkzTasks++;
      break;
   case Enum:
      activeEnumSolvers.insert(std::make_pair(solverThreadId,pool[( solverThreadId )]));
      enumSelectionHeap->insert(pool[( solverThreadId )]);
      nEnumTasks++;
      break;
   case Sieve:
      activeSieveSolvers.insert(std::make_pair(solverThreadId,pool[( solverThreadId )]));
      sieveSelectionHeap->insert(pool[( solverThreadId )]);
      nSieveTasks++;
      break;
   default:
      THROW_LOGICAL_ERROR2("CMapLapParaSolverPool::activateSolver: Invalid solver type = ",static_cast<int>(solverType));
   }
   pool[( solverThreadId )]->activate(solverType, paraTask);
   if( paraTask->isLocalTask() )
   {
      threadId = 1; // in local task, thread Id must be 1
   }
   paraTask->setThreadId(threadId);
   paraTask->send(paraComm, rank);
   return rank;
}


///
/// get interrupt DeepBkz solver
/// @return interrupt DeepBkz solver rank, -1 is reserved for Sieve
///
int
CMapLapParaSolverPool::getInterruptDeepBkzSolver(
      int &threadId           ///< interrupt DeepBkz solver thread id (return value)
      )
{
   assert( deepBkzSelectionHeap->getHeapSize() > 0 );
   if( deepBkzSelectionHeap->top()->isReservedForSieve() ) return -1;
   threadId = deepBkzSelectionHeap->top()->getThreadId();
   return deepBkzSelectionHeap->top()->getRank();
}


///
/// get interrupt Enum solver
/// @return interrupt Enum solver rank
///
int
CMapLapParaSolverPool::getInterruptEnumSolver(
      int &threadId           ///< interrupt Enum solver thread id (return value)
      )
{
   assert( enumSelectionHeap->getHeapSize() > 0 );
   if( enumSelectionHeap->top()->isReservedForSieve() ) return -1;
   threadId = enumSelectionHeap->top()->getThreadId();
   return enumSelectionHeap->top()->getRank();
}


///
/// interrupt DeepBkz Solvers
///
void
CMapLapParaSolverPool::interruptDeepBkzSolvers(
      int nInterrupt     ///< number of interrupt solvers
      )
{
   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm );
   while( static_cast<int>(interruptingDeepBkzSolvers.size()) < nInterrupt && activeDeepBkzSolvers.size() > 0 )
   {
      int threadId = -1;
      int rank = getInterruptDeepBkzSolver(threadId);
      if( rank >= 0 )
      {
         assert(threadId >= 0);
         unsigned int solverThreadId = getSolverThreadId(rank, threadId);
         assert( !pool[solverThreadId]->isReservedForSieve() );
         auto p = activeDeepBkzSolvers.find(solverThreadId);
         if( p != activeDeepBkzSolvers.end() )
         {
            if( p->second->getRank() != rank ||
                  pool[( solverThreadId )]->getRank() != rank ||
                  pool[( solverThreadId )]->getThreadId() != threadId )
            {
               THROW_LOGICAL_ERROR4("Invalid rank. Rank = ", rank, ", or invalid thread Id = ", threadId);
            }
            PARA_COMM_CALL(
                  cmapLapParaComm->send( &threadId, 1, UG::ParaINT, rank, UG::TagInterruptRequest);
            );
            activeDeepBkzSolvers.erase(p);
            deepBkzSelectionHeap->deleteElement(pool[( solverThreadId )]);
            interruptingDeepBkzSolvers.insert(std::make_pair(solverThreadId,pool[( solverThreadId )]));
         }
         else
         {
            THROW_LOGICAL_ERROR2("DeepBkz solver rank = ", rank);
         }
      }
   }
}


///
/// interrupt Enum Solvers
///
void
CMapLapParaSolverPool::interruptEnumSolvers(
      int nInterrupt     ///< number of interrupt solvers
      )
{
   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm );
   while( static_cast<int>(interruptingEnumSolvers.size()) < nInterrupt && activeEnumSolvers.size() > 0 )
   {
      int threadId = -1;
      int rank = getInterruptEnumSolver(threadId);
      if( rank >= 0 )
      {
         assert(threadId >= 0);
         unsigned int solverThreadId = getSolverThreadId(rank, threadId);
         assert( !pool[solverThreadId]->isReservedForSieve() );
         auto p = activeEnumSolvers.find(solverThreadId);
         if( p != activeEnumSolvers.end() )
         {
            if( p->second->getRank() != rank ||
                  pool[( solverThreadId )]->getRank()     != static_cast<int>(rank) ||
                  pool[( solverThreadId )]->getThreadId() != static_cast<int>(threadId) )
            {
               THROW_LOGICAL_ERROR4("Invalid rank. Rank = ", rank, ", or invalid thread Id = ", threadId);
            }
            PARA_COMM_CALL(
                  cmapLapParaComm->send( &threadId, 1, UG::ParaINT, rank, UG::TagInterruptRequest);
            );
            activeEnumSolvers.erase(p);
            enumSelectionHeap->deleteElement(pool[( solverThreadId )]);
            interruptingEnumSolvers.insert(std::make_pair(solverThreadId,pool[( solverThreadId )]));
         }
         else
         {
            THROW_LOGICAL_ERROR2("Enum solver rank = ", rank);
         }
      }
   }
}


///
/// interrupt DeepBkz and Enum Solvers to run a Sieve Solver
///
bool
CMapLapParaSolverPool::interruptDeepBkzAndEnumSolversToRunSieveSolver(
      )
{
   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm );

   CMapLapParaSolverPoolElementPtr mainThreadSolver = 0;

   for(int i = 0; i < static_cast<int>(deepBkzSelectionHeap->getHeapSize()); i++ )
   {
      if( deepBkzSelectionHeap->elem(i+1)->getThreadId() == 0 )
      {
         mainThreadSolver = deepBkzSelectionHeap->elem(i+1);
         break;
      }
   }
   if( !mainThreadSolver )
   {
      for(int i = 0; i < static_cast<int>(enumSelectionHeap->getHeapSize()); i++ )
      {
         if( enumSelectionHeap->elem(i+1)->getThreadId() == 0 )
         {
            mainThreadSolver = enumSelectionHeap->elem(i+1);
            break;
         }
      }
   }

   if( !mainThreadSolver ) return false;

   int rank = mainThreadSolver->getRank();
   for( unsigned int threadId = 0; threadId < nThreadsPerRank; threadId++ )
   {
      unsigned int solverThreadId = getSolverThreadId(rank, threadId);
      if( !pool[solverThreadId]->isActive() )
      {
         pool[solverThreadId]->reserveForSieve();
      }
      else
      {
         switch ( pool[solverThreadId]->getSolverType() )
         {
         case DeepBkz:
         {
            auto p = activeDeepBkzSolvers.find(solverThreadId);
            if( p != activeDeepBkzSolvers.end() )
            {
               if( p->second->getRank() != rank ||
                     pool[( solverThreadId )]->getRank()     != static_cast<int>(rank) ||
                     pool[( solverThreadId )]->getThreadId() != static_cast<int>(threadId) )
               {
                  THROW_LOGICAL_ERROR4("Invalid rank. Rank = ", rank, ", or invalid thread Id = ", threadId);
               }
               PARA_COMM_CALL(
                     cmapLapParaComm->send( &threadId, 1, UG::ParaINT, rank, UG::TagInterruptRequest);
               );
               std::cout << "CMapLapParaSolverPool::interruptDeepBkzAndEnumSolversToRunSieveSolver: send UG::TagInterruptRequest rank " << rank << " threadId " << threadId  << std::endl;
               activeDeepBkzSolvers.erase(p);
               deepBkzSelectionHeap->deleteElement(pool[( solverThreadId )]);
               interruptingDeepBkzSolvers.insert(std::make_pair(solverThreadId,pool[( solverThreadId )]));
            }
            else
            {
               THROW_LOGICAL_ERROR2("DeepBkz solver rank = ", rank);
            }
            pool[solverThreadId]->reserveForSieve();
            break;
         }
         case Enum:
         {
            auto p = activeEnumSolvers.find(solverThreadId);
            if( p != activeEnumSolvers.end() )
            {
               if( p->second->getRank() != rank ||
                     pool[( solverThreadId )]->getRank()     != static_cast<int>(rank) ||
                     pool[( solverThreadId )]->getThreadId() != static_cast<int>(threadId) )
               {
                  THROW_LOGICAL_ERROR4("Invalid rank. Rank = ", rank, ", or invalid thread Id = ", threadId);
               }
               PARA_COMM_CALL(
                     cmapLapParaComm->send( &threadId, 1, UG::ParaINT, rank, UG::TagInterruptRequest);
               );
               activeEnumSolvers.erase(p);
               enumSelectionHeap->deleteElement(pool[( solverThreadId )]);
               interruptingEnumSolvers.insert(std::make_pair(solverThreadId,pool[( solverThreadId )]));
            }
            else
            {
               THROW_LOGICAL_ERROR2("Enum solver rank = ", rank);
            }
            pool[solverThreadId]->reserveForSieve();
            break;
         }
         default:
            THROW_LOGICAL_ERROR2(
                  "Invalid solver type to interrupt. solver type = ",
                  static_cast<int>(pool[solverThreadId]->getSolverType()));
         }
      }
   }
   return true;
}


///
/// interrupt DeepBkz and Enum Solvers to run a Sieve Solver
///
bool
CMapLapParaSolverPool::interruptEnumAndDeepBkzSolversToRunSieveSolver(
      )
{
   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm );

   CMapLapParaSolverPoolElementPtr mainThreadSolver = 0;

   for(int i = 0; i < static_cast<int>(enumSelectionHeap->getMaxHeapSize()); i++ )
   {
      if( enumSelectionHeap->elem(i+1)->getThreadId() == 0 )
      {
         mainThreadSolver = enumSelectionHeap->elem(i+1);
         break;
      }
   }
   if( !mainThreadSolver )
   {
      for(int i = 0; i < static_cast<int>(deepBkzSelectionHeap->getMaxHeapSize()); i++ )
      {
         if( deepBkzSelectionHeap->elem(i+1)->getThreadId() == 0 )
         {
            mainThreadSolver = deepBkzSelectionHeap->elem(i+1);
            break;
         }
      }
   }

   if( !mainThreadSolver ) return false;

   int rank = mainThreadSolver->getRank();
   for( unsigned int threadId = 0; threadId < nThreadsPerRank; threadId++ )
   {
      unsigned int solverThreadId = getSolverThreadId(rank, threadId);
      if( !pool[solverThreadId]->isActive() )
      {
         pool[solverThreadId]->reserveForSieve();
      }
      else
      {
         switch ( pool[solverThreadId]->getSolverType() )
         {
         case DeepBkz:
         {
            auto p = activeDeepBkzSolvers.find(solverThreadId);
            if( p != activeDeepBkzSolvers.end() )
            {
               if( p->second->getRank() != rank ||
                     pool[( solverThreadId )]->getRank()     != static_cast<int>(rank) ||
                     pool[( solverThreadId )]->getThreadId() != static_cast<int>(threadId) )
               {
                  THROW_LOGICAL_ERROR4("Invalid rank. Rank = ", rank, ", or invalid thread Id = ", threadId);
               }
               PARA_COMM_CALL(
                     cmapLapParaComm->send( &threadId, 1, UG::ParaINT, rank, UG::TagInterruptRequest);
               );
               activeDeepBkzSolvers.erase(p);
               deepBkzSelectionHeap->deleteElement(pool[( solverThreadId )]);
               interruptingDeepBkzSolvers.insert(std::make_pair(solverThreadId,pool[( solverThreadId )]));
            }
            else
            {
               THROW_LOGICAL_ERROR2("DeepBkz solver rank = ", rank);
            }
            pool[solverThreadId]->reserveForSieve();
            break;
         }
         case Enum:
         {
            auto p = activeEnumSolvers.find(solverThreadId);
            if( p != activeEnumSolvers.end() )
            {
               if( p->second->getRank() != rank ||
                     pool[( solverThreadId )]->getRank()     != static_cast<int>(rank) ||
                     pool[( solverThreadId )]->getThreadId() != static_cast<int>(threadId) )
               {
                  THROW_LOGICAL_ERROR4("Invalid rank. Rank = ", rank, ", or invalid thread Id = ", threadId);
               }
               PARA_COMM_CALL(
                     cmapLapParaComm->send( &threadId, 1, UG::ParaINT, rank, UG::TagInterruptRequest);
               );
               activeEnumSolvers.erase(p);
               enumSelectionHeap->deleteElement(pool[( solverThreadId )]);
               interruptingEnumSolvers.insert(std::make_pair(solverThreadId,pool[( solverThreadId )]));
            }
            else
            {
               THROW_LOGICAL_ERROR2("Enum solver rank = ", rank);
            }
            pool[solverThreadId]->reserveForSieve();
            break;
         }
         default:
            THROW_LOGICAL_ERROR2(
                  "Invalid solver type to interrupt. solver type = ",
                  static_cast<int>(pool[solverThreadId]->getSolverType()));
         }
      }
   }
   return true;
}



///
/// inactivate the Solver specified by rank
///
void
CMapLapParaSolverPool::inactivateSolver(
      int rank,                                          ///< rank of the Solver to be inactivated
      int threadId,                                      ///< thread ID of the Solver to be inactivated
      std::shared_ptr<CMapLapParaTaskPool> paraNodePool  ///< pointer to CMapLapParaNodePool to change into collecting mode
      )
{
   if( rank < originRank || rank >= (signed)(originRank + nSolvers)  )
   {
      THROW_LOGICAL_ERROR2("Invalid rank. Rank = ", rank);
   }
   if( threadId >=  static_cast<int>(nThreadsPerRank) )
   {
      THROW_LOGICAL_ERROR2("Invalid threadId. ThreadId = ", threadId);
   }

   unsigned int solverThreadId = getSolverThreadId(rank, threadId);

   ActiveSolvers           *activeSolvers;
   AscendingSelectionHeap  *selectionHeap;
   InterruptingSolvers     *interruptingSolvers;
   if ( pool[( solverThreadId )]->getStatus() == UG::Active )
   {
      switch( pool[( solverThreadId )]->getSolverType() )
      {
      case DeepBkz:
         {
            activeSolvers = &activeDeepBkzSolvers;
            selectionHeap = deepBkzSelectionHeap;
            interruptingSolvers = &interruptingDeepBkzSolvers;
            break;
         }
         case Enum:
         {
            activeSolvers = &activeEnumSolvers;
            selectionHeap = enumSelectionHeap;
            interruptingSolvers = &interruptingEnumSolvers;
            break;
         }
         case Sieve:
         {
            activeSolvers = &activeSieveSolvers;
            selectionHeap = sieveSelectionHeap;
            interruptingSolvers = &interruptingSieveSolvers;
            break;
         }
         default:
            THROW_LOGICAL_ERROR2("CMapLapParaSolverPool::interruptAllSolvers: Invalid solver type = ",static_cast<int>(pool[( solverThreadId )]->getSolverType()));
      }

      auto p = activeSolvers->find(solverThreadId);
      if( p != activeSolvers->end() )
      {
         if( p->second->getRank() != rank ||
               pool[( solverThreadId )]->getRank()     != static_cast<int>(rank) ||
               pool[( solverThreadId )]->getThreadId() != static_cast<int>(threadId) )
         {
            THROW_LOGICAL_ERROR4("Invalid rank. Rank = ", rank, ", or invalid thread Id = ", threadId);
         }
         activeSolvers->erase(p);
         selectionHeap->deleteElement(pool[( solverThreadId )]);
      }
      else
      {
         p = interruptingSolvers->find(solverThreadId);
         if( p != interruptingSolvers->end() )
         {
            if( p->second->getRank() != rank ||
                  pool[( solverThreadId )]->getRank() != rank ||
                  pool[( solverThreadId )]->getThreadId() != threadId )
            {
               THROW_LOGICAL_ERROR4("Invalid rank. Rank = ", rank, ", or invalid thread Id = ", threadId);
            }
            interruptingSolvers->erase(p);
         }
         else
         {
            THROW_LOGICAL_ERROR2("Solver rank = ", rank);
         }
      }
   }
   else
   {
      THROW_LOGICAL_ERROR2("Rank = ", rank);
   }
   pool[( solverThreadId )]->inactivate();
   size_t mainSolverThreadId = getSolverThreadId(rank, 0);

   if( pool[( mainSolverThreadId )]->getSolverType() == Sieve )
   {
      assert( pool[( solverThreadId )]->isReservedForSieve() );
      if( mainSolverThreadId == solverThreadId )
      {
         for( unsigned int thId = 0; thId < nThreadsPerRank; thId++ )
         {
            int sieveSolverThreadId = getSolverThreadId(rank, thId);
            pool[( sieveSolverThreadId )]->resetReservedForSieve();
            inactiveSolvers.insert(std::make_pair(sieveSolverThreadId,pool[( sieveSolverThreadId )]) );
         }
      }
   }
   else
   {
      if( pool[( solverThreadId )]->getThreadId() == 0 || !pool[( solverThreadId )]->isReservedForSieve() )
      {
         inactiveSolvers.insert(std::make_pair(solverThreadId,pool[( solverThreadId )]));
      }
   }
}


///
/// getter of basis of current task
///
std::shared_ptr<LatticeBasis<int>>
CMapLapParaSolverPool::getBasis(
      int      rank,                ///< rank of the Solver
      int      threadId             ///< thread ID
      )
{
   if( rank < originRank || rank >= (signed)(originRank + nSolvers)  )
   {
      THROW_LOGICAL_ERROR2("Invalid rank. Rank = ", rank);
   }
   if( threadId >=  static_cast<int>(nThreadsPerRank) )
   {
      THROW_LOGICAL_ERROR2("Invalid threadId. ThreadId = ", threadId);
   }

   unsigned int solverThreadId = getSolverThreadId(rank, threadId);
   assert( pool[(solverThreadId)]->getStatus() == UG::Active );
   assert( pool[(solverThreadId)]->getSolverType() == DeepBkz );

   auto p = activeDeepBkzSolvers.find(solverThreadId);
   if( p == activeDeepBkzSolvers.end() )
   {
      p = interruptingDeepBkzSolvers.find(solverThreadId);
   }
   if( p == interruptingDeepBkzSolvers.end() )
   {
      THROW_LOGICAL_ERROR4("Invalid rank. Rank = ", rank, ", or invalid thread Id = ", threadId);
   }
   return p->second->getBasis();
}


///
/// update DeepBkz Solver status
///
void
CMapLapParaSolverPool::updateSolverStatus(
      int      rank,                ///< rank of the Solver
      int      threadId,            ///< thread ID
      int      inBlockSize,         ///< current blocksize
      int*     inBasis,             //< current basis
      double   inEnumCost,          ///< approximated time of full Enum [Pool] better a than b when a < b
      double   inSlopeGSA,          ///< slope of GSA line [Pool] better a than b when a > b
      double   inTopHalfSlopeGSA,   ///< slope of top-half GSA line [Pool] better a than b when a > b
      double   inOrthogonalFactor   ///< orthogonal factor [Pool] better a than b when a < b
      )
{
   if( rank < originRank || rank >= (signed)(originRank + nSolvers)  )
   {
      THROW_LOGICAL_ERROR2("Invalid rank. Rank = ", rank);
   }
   if( threadId >=  static_cast<int>(nThreadsPerRank) )
   {
      THROW_LOGICAL_ERROR2("Invalid threadId. ThreadId = ", threadId);
   }

   unsigned int solverThreadId = getSolverThreadId(rank, threadId);
   assert( pool[(solverThreadId)]->getStatus() == UG::Active );
   assert( pool[(solverThreadId)]->getSolverType() == DeepBkz );

   auto p = activeDeepBkzSolvers.find(solverThreadId);
   if( p != activeDeepBkzSolvers.end() )
   {
      // assert( (numOfNodesSolved - p->second->getNumOfNodesSolved() ) >= 0 ); // if solver restart, this is not always true
      p->second->updateDeepBkzSolverData(
            inBlockSize,
            inBasis,
            inEnumCost,
            inSlopeGSA,
            inTopHalfSlopeGSA,
            inOrthogonalFactor
            );
      deepBkzSelectionHeap->updateInterruptPriorityValue(*(pool[(solverThreadId)]->getSelectionHeapElement()), -inEnumCost);
   }
   else
   {
      p = interruptingDeepBkzSolvers.find(solverThreadId);
      if( p == interruptingDeepBkzSolvers.end() )
      {
         THROW_LOGICAL_ERROR4("Invalid rank. Rank = ", rank, ", or invalid thread Id = ", threadId);
      }
   }
}


///
/// update Enum Solver status
///
void
CMapLapParaSolverPool::updateSolverStatus(
      int      rank,                ///< rank of the Solver
      int      threadId,            ///< thread ID
      int      inDepth,
      int*     inCoeffs,
      long int inNumSearchedNodes   ///< number of searched nodes in the enumeration tree
      )
{
   if( rank < originRank || rank >= (signed)(originRank + nSolvers)  )
   {
      THROW_LOGICAL_ERROR2("Invalid rank. Rank = ", rank);
   }
   if( threadId >=  static_cast<int>(nThreadsPerRank) )
   {
      THROW_LOGICAL_ERROR2("Invalid threadId. ThreadId = ", threadId);
   }

   unsigned int solverThreadId = getSolverThreadId(rank, threadId);
   assert( pool[(solverThreadId)]->getStatus() == UG::Active );
   assert( pool[(solverThreadId)]->getSolverType() == Enum );

   auto p = activeEnumSolvers.find(solverThreadId);
   if( p != activeEnumSolvers.end() )
   {
      p->second->updateEnumSolverData(
            inDepth,
            inCoeffs,
            inNumSearchedNodes
            );
      enumSelectionHeap->updateInterruptPriorityValue(*(pool[(solverThreadId)]->getSelectionHeapElement()), - static_cast<long double>(inNumSearchedNodes));
   }
   else
   {
      p = interruptingEnumSolvers.find(solverThreadId);
      if( p == interruptingEnumSolvers.end() )
      {
         THROW_LOGICAL_ERROR4("Invalid rank. Rank = ", rank, ", or invalid thread Id = ", threadId);
      }
   }
}


///
/// send message for specific solverType solvers
///
void
CMapLapParaSolverPool::sendMessageForSolvers(
      CMapLapParaSolverLocalComm *localComm,
      int tag,
      std::vector<bool> &hasSentManager,
      SolverType solverType
      )
{
   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm );

   for( unsigned int i = 0; i < nSolvers; i++ )
   {
      for( unsigned int j = 0; j < nThreadsPerRank; j++ )
      {
         unsigned int solverThreadId = getSolverThreadId(originRank+i, j);
         if( pool[( solverThreadId )]->isActive() &&
               pool[( solverThreadId )]->getSolverType() == solverType )
         {
            int rank = pool[( solverThreadId )]->getRank();
            int threadId = pool[( solverThreadId )]->getThreadId();
            if( !hasSentManager[solverThreadId] )
            {
               if( cmapLapParaComm )
               {
                  PARA_COMM_CALL(
                        cmapLapParaComm->send(&threadId, 1, UG::ParaINT, rank, tag);
                  );
               }
               else
               {
                  PARA_COMM_CALL(
                        localComm->send(&threadId, 1, UG::ParaINT, rank, tag);
                  );
               }
               hasSentManager[solverThreadId] = true;
            }
         }
      }
   }
}


///
/// interrupt all active solvers
///
void
CMapLapParaSolverPool::interruptAllSolvers(
      CMapLapParaSolverLocalComm *localComm,
      int tag
      )
{
   assert( tag == UG::TagInterruptRequest || tag == TagTimeLimitRequest );
   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm );

   ActiveSolvers           *activeSolvers;
   AscendingSelectionHeap  *selectionHeap;
   InterruptingSolvers     *interruptingSolvers;
   for( unsigned int i = 0; i < nSolvers; i++ )
   {
      for( unsigned int j = 0; j < nThreadsPerRank; j++ )
      {
         unsigned int solverThreadId = getSolverThreadId(originRank+i, j);
         if( pool[( solverThreadId )]->isActive() )
         {
            int rank = pool[( solverThreadId )]->getRank();
            int threadId = pool[( solverThreadId )]->getThreadId();
            if( cmapLapParaComm )
            {
               PARA_COMM_CALL(
                     cmapLapParaComm->send( &threadId, 1, UG::ParaINT, rank, tag);
               );
            }
            else
            {
               PARA_COMM_CALL(
                     localComm->send( &threadId, 1, UG::ParaINT, rank, tag);
               );
            }

            switch( pool[( solverThreadId )]->getSolverType() )
            {
            case DeepBkz:
            {
               activeSolvers = &activeDeepBkzSolvers;
               selectionHeap = deepBkzSelectionHeap;
               interruptingSolvers = &interruptingDeepBkzSolvers;
               break;
            }
            case Enum:
            {
               activeSolvers = &activeEnumSolvers;
               selectionHeap = enumSelectionHeap;
               interruptingSolvers = &interruptingEnumSolvers;
               break;
            }
            case Sieve:
            {
               activeSolvers = &activeSieveSolvers;
               selectionHeap = sieveSelectionHeap;
               interruptingSolvers = &interruptingSieveSolvers;
               break;
            }
            default:
               THROW_LOGICAL_ERROR2("CMapLapParaSolverPool::interruptAllSolvers: Invalid solver type = ",static_cast<int>(pool[( solverThreadId )]->getSolverType()));
            }

            auto p = activeSolvers->find(solverThreadId);
            if( p != activeSolvers->end() )
            {
               activeSolvers->erase(p);
               selectionHeap->deleteElement(pool[( solverThreadId )]);
               interruptingSolvers->insert(std::make_pair(solverThreadId,pool[( solverThreadId )]));
            }
            else
            {
               p = interruptingSolvers->find(solverThreadId);
               if( p != interruptingSolvers->end() )
               {
                  if( p->second->getRank() != rank ||
                        pool[( solverThreadId )]->getRank() != rank ||
                        pool[( solverThreadId )]->getThreadId() != threadId )
                  {
                     THROW_LOGICAL_ERROR4("Invalid rank. Rank = ", rank, ", or invalid thread Id = ", threadId);
                  }
               }
               else
               {
                  THROW_LOGICAL_ERROR2("Solver rank = ", rank);
               }
            }
         }
      }
   }
}


///
/// terminate all solvers
///
void
CMapLapParaSolverPool::terminateAllSolvers(
      CMapLapParaSolverLocalComm *localComm
      )
{
   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm );

   if( cmapLapParaComm )
   {
      for( unsigned int i = 0; i < nSolvers; i++ )
      {
         PARA_COMM_CALL(
               cmapLapParaComm->send(NULL, 0, UG::ParaBYTE, (originRank+i), UG::TagTerminateRequest);
         );
      }
   }
   else
   {
      // local solver termination
      for( unsigned int i = 0; i < nSolvers; i++ )
      {
         PARA_COMM_CALL(
               localComm->send(NULL, 0, UG::ParaBYTE, (originRank+i), UG::TagTerminateRequest);
         );
      }
   }
}


///
/// get basis of DeepBkz Solvers
/// @return deque of basis(Eigen::LatticeBasis<int>)
///
std::deque<std::shared_ptr<LatticeBasis<int>>>
CMapLapParaSolverPool::getBasisOfSolvers(
      )
{
   std::deque<std::shared_ptr<LatticeBasis<int>>> basisDeque;
   int numBasis = 0;
   for( auto x : activeDeepBkzSolvers )
   {
      // x->second is SolverPoolElement
      auto basis = x.second->getCurrentTask()->getBasis();
      basisDeque.push_back(basis);
      if( ++numBasis > 10000 ){ break; }
   }
   return basisDeque;
}


} // namespace ParaCMapLAP
