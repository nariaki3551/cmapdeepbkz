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

/**@file    cmapLapParaSolverPool.h
 * @brief   Solver pool.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_SOLVER_POOL_H__
#define __CMAP_LAP_PARA_SOLVER_POOL_H__


#include <unordered_map>
#include "ug/paraSolverPool.h"
#include "ug/paraComm.h"
#include "ug/paraDef.h"
#include "ug/paraSolverPool.h"
#include "ug/paraTask.h"
#include "cmapLapParaDef.h"
#include "cmapLapParaTaskPool.h"
#include "cmapLapParaSolverLocalComm.h"
#ifdef UG_WITH_ZLIB
#include "ug/gzstream.h"
#endif
namespace UG { class ParaParamSet; }
namespace UG { class ParaTimer; }
namespace ParaCMapLAP { class CMapLapParaSolverPoolElement; }


namespace ParaCMapLAP
{

struct DeepBkzSolverPoolElementData {
   double   enumCost;           ///< approximated time of full Enum [Pool] better a than b when a < b
   double   slopeGSA;           ///< slope of GSA line [Pool] better a than b when a > b
   double   topHalfSlopeGSA;    ///< slope of top-half GSA line [Pool] better a than b when a > b
   double   orthogonalFactor;   ///< orthogonal factor [Pool] better a than b when a < b
};                              ///< data that are used for CMapLap solver

struct EnumSolverPoolElementData {
   long int numSearchedNodes;   ///< number of searched nodes in the enumeration tree [Pool] better a than b when b > a
};                              ///< data that are used for Enum Solver

struct SievingSolverPoolElementData {

};                              ///< data that are used for Sieving Solver

union SolverPoolElementData {
   DeepBkzSolverPoolElementData deepBzkSolverData; ///< data for CMapLap solver
   EnumSolverPoolElementData    enumSolverData;    ///< data for Enum Solver
   SievingSolverPoolElementData sievingSolverData; ///< data for Sieving solver
};                                                 ///< union data type for data that are depending on solver type


typedef CMapLapParaSolverPoolElement * CMapLapParaSolverPoolElementPtr;

///
/// class CMapLapParaSolverPoolElement
/// (This class includes information about a Solver status)
///
class CMapLapParaSolverPoolElement
{

   int                   rank;                              ///< rank of the Solver
   int                   threadId;                          ///< thread ID
   SolverType            solverType;                        ///< type of Solver
   UG::SolverStatus      status;                            ///< status of the Solver
   long double           interruptPriorityValue;            ///< value related to interrupt priority. Smaller value has higher priority
   bool                  reseredForSieve;                   ///< indicate if Solver element is reserved for Sieve solver or not
   std::shared_ptr<CMapLapParaTask> currentTask;            ///< solving task
   CMapLapParaSolverPoolElementPtr *selectionHeapElement;   ///< pointer to selection heap element
   ///------------------------------------------------------------------
   /// the following values must be depending on Solver type
   ///------------------------------------------------------------------
   SolverPoolElementData solverPoolElementData;               ///< data that are depending on Solver type

public :

   ///
   /// constructor
   ///
   CMapLapParaSolverPoolElement(
         int inRank,
         int inThreadId
         ) : rank(inRank),
             threadId(inThreadId),
             solverType(Undefined),
             status(UG::Inactive),
             interruptPriorityValue(DBL_MAX),
             reseredForSieve(false),
             currentTask(nullptr),
             selectionHeapElement(nullptr)
   {
   }

   ///
   ///  destractor
   ///
   virtual ~CMapLapParaSolverPoolElement(
         )
   {
   }

   ///
   /// get rank of the Solver
   /// @return rank of Solver
   ///
   virtual int getRank(
         )
   {
      return rank;
   }

   ///
   /// get thread ID of the Solver
   /// @return thread ID of Solver
   ///
   virtual int getThreadId(
         )
   {
      return threadId;
   }

   ///
   /// get solver type of the Solver
   /// @return solver type of Solver
   ///
   virtual SolverType getSolverType(
         )
   {
      return solverType;
   }

   ///
   /// get current solving CMapLapParaNode
   /// @return pointer to CMapLapParaNode object
   ///
   virtual std::shared_ptr<CMapLapParaTask> getCurrentTask(
         )
   {
      return currentTask;
   }

   ///
   /// check if solver element is reserved for Sieve solver or not
   /// @return true if it is reserved for Sieve solver
   ///
   virtual bool isReservedForSieve(
         )
   {
      return reseredForSieve;
   }

   virtual void reserveForSieve(
         )
   {
      assert( reseredForSieve == false );
      reseredForSieve = true;
   }

   virtual void resetReservedForSieve(
         )
   {
      reseredForSieve = false;
   }

   ///
   /// extract current solving CMapLapParaNode
   /// @return pointer to CMapLapParaNode object
   ///
   virtual std::shared_ptr<CMapLapParaTask> extractCurrentTask(
         )
   {
      std::shared_ptr<CMapLapParaTask> task(std::move(currentTask));
      assert( currentTask = nullptr );
      return task;
   }

   ///
   /// activate this Solver
   ///
   virtual void activate(
         SolverType  inSolverType,                 ///< solver type
         std::shared_ptr<CMapLapParaTask> inTask   ///< pointer to CMapLapParaNode object to activate Solver
         )
   {
      status = UG::Active;
      solverType = inSolverType;
      currentTask = inTask;
      switch( solverType )
      {
      case DeepBkz:
         interruptPriorityValue = -(inTask->getEstimatedValue());
         solverPoolElementData.deepBzkSolverData.enumCost = inTask->getEstimatedValue();
         solverPoolElementData.deepBzkSolverData.slopeGSA = 0.0;
         solverPoolElementData.deepBzkSolverData.topHalfSlopeGSA = 0.0;
         solverPoolElementData.deepBzkSolverData.orthogonalFactor = 0.0;
         break;
      case Enum:
         solverPoolElementData.enumSolverData.numSearchedNodes = 0;
         break;
      case Sieve:
         // nothing to initialize
         break;
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaSolverPoolElement: Invalid solver type = ",static_cast<int>(solverType));
      }
   }

   ///
   ///  inactivate this Solver
   ///
   virtual void inactivate(
         )
   {
      status = UG::Inactive;
      solverType = Undefined;
      currentTask = nullptr;
   }

   ///
   /// kill this Solver
   /// @return pointer to CMapLapParaNode object assigned to this Solver
   ///
   virtual std::shared_ptr<CMapLapParaTask> died(
         )
   {
      status = UG::Dead;
      solverType = Undefined;
      std::shared_ptr<CMapLapParaTask> task = std::move(currentTask);
      assert( currentTask = nullptr );
      return task;
   }

   ///
   /// get Solver status
   /// @return Solver status
   ///
   virtual UG::SolverStatus getStatus(
         )
   {
      return status;
   }

   ///
   /// get interrupt priority value of Solver
   /// @return interrupt priority value
   ///
   virtual long double getInterruptPriorityValue(
         )
   {
      return interruptPriorityValue;
   }

   ///
   /// set enumeration cost on Solver
   ///
   virtual void setInterruptPriorityValue(
         long double inInterruptPriorityValue     ///< enumeration cost
         )
   {
      interruptPriorityValue = inInterruptPriorityValue;
   }

   ///
   /// get selection heap element
   /// @return pointer to CMapLapParaSolverPoolElementPtr
   ///
   virtual CMapLapParaSolverPoolElementPtr *getSelectionHeapElement(
         )
   {
      return selectionHeapElement;
   }


   ///
   /// set selection heap element
   ///
   virtual void setSelectionHeapElement(
         CMapLapParaSolverPoolElementPtr *inSelectionHeapElement   ///< selection heap element
         )
   {
      selectionHeapElement = inSelectionHeapElement;
   }

   ///
   /// check if this Solver is active or not
   /// @return true if this Solver is active, false otherwise
   ///
   virtual bool isActive(
         )
   {
      return ( status == UG::Active );
   }

   ///
   /// get approximated time of full Enum
   /// @return approximated time of full Enum
   ///
   virtual double getEnumCost(
         )
   {
      assert( solverType == DeepBkz );
      return solverPoolElementData.deepBzkSolverData.enumCost;
   }

   ///
   /// get slope of GSA line
   /// @return slope of GSA line
   ///
   virtual double getSlopeGSA(
         )
   {
      assert( solverType == DeepBkz );
      return solverPoolElementData.deepBzkSolverData.slopeGSA;
   }

   ///
   /// get slope of top-half GSA line
   /// @return slope of top-half GSA line
   ///
   virtual double getTopHalfSlopeGSA(
         )
   {
      assert( solverType == DeepBkz );
      return solverPoolElementData.deepBzkSolverData.topHalfSlopeGSA;
   }

   ///
   /// get orthogonal factor
   /// @return orthogonal factor
   ///
   virtual double getOrthogonalFactor(
         )
   {
      assert( solverType == DeepBkz );
      return solverPoolElementData.deepBzkSolverData.orthogonalFactor;
   }

   ///
   /// getter of basis of currentTask
   ///
   virtual std::shared_ptr<LatticeBasis<int>> getBasis(
         )
   {
      assert( solverType == DeepBkz );
      return currentTask->getBasis();
   }

   ///
   /// update DeepBkz solver data
   ///
   virtual void updateDeepBkzSolverData(
         int      inBlockSize,            ///< current blocksize
         int*     inBasis,                 ///<
         double   inEnumCost,             ///< approximated time of full Enum [Pool] better a than b when a < b
         double   inSlopeGSA,             ///< slope of GSA line [Pool] better a than b when a > b
         double   inTopHalfSlopeGSA,      ///< slope of top-half GSA line [Pool] better a than b when a > b
         double   inOrthogonalFactor      ///< orthogonal factor [Pool] better a than b when a < b
         )
   {
      assert( solverType == DeepBkz );
      interruptPriorityValue = -inEnumCost;
      currentTask->updateDeepBkz(inBlockSize, inBasis);
      solverPoolElementData.deepBzkSolverData.enumCost = inEnumCost;
      solverPoolElementData.deepBzkSolverData.slopeGSA = inSlopeGSA;
      solverPoolElementData.deepBzkSolverData.topHalfSlopeGSA = inTopHalfSlopeGSA;
      solverPoolElementData.deepBzkSolverData.orthogonalFactor = inOrthogonalFactor;
   }

   ///
   /// update Enum solver data
   ///
   virtual void updateEnumSolverData(
         int inDepth,
         int* inCoef,
         long int inNumSearchedNodes    ///< number of searched nodes in the enumeration tree [Pool] better a than b when b > a
         )
   {
      assert( solverType == Enum );
      interruptPriorityValue = static_cast<long double>(inNumSearchedNodes);
      currentTask->updateEnum(inDepth, inCoef);
      solverPoolElementData.enumSolverData.numSearchedNodes = inNumSearchedNodes;
   }


   ///
   /// get number of searched nodes
   /// @return number of searched nodes
   ///
   virtual long int getNumSearchedNodes(
         )
   {
      assert( solverType == Enum );
      return solverPoolElementData.enumSolverData.numSearchedNodes;
   }

};



///
/// class Selection Heap
///
class SelectionHeap
{
protected:
   std::size_t          maxHeapSize;               ///< maximum size of this heap
   std::size_t          heapSize;                  ///< current used heap size
   CMapLapParaSolverPoolElementPtr *heap;           ///< heap : contents are BbParaSolverPoolElementPtr

public:

   ///
   /// results of insert
   ///
   enum ResultOfInsert
   {
      SUCCEEDED,     ///< SUCCEEDED
      FAILED_BY_FULL ///< FAILED_BY_FULL
   };

   ///
   /// constructor
   ///
   SelectionHeap(
         std::size_t size                   ///< heap size
         )
   {
      maxHeapSize = size;
      heapSize = 0;
      heap = new CMapLapParaSolverPoolElementPtr[size+1]; /// index 0 is a sentinel
   }

   ///
   /// destructor of selection heap
   ///
   virtual ~SelectionHeap(
         )
   {
      delete[] heap;
   }

   ///
   /// insert BbParaSolverPoolElementPtr to Selection Heap
   /// @return succeeded or fail because of full
   ///
   virtual ResultOfInsert insert(
         CMapLapParaSolverPoolElementPtr solver    ///< Solver pool element to be inserted
         );

   ///
   /// obtain top priority BbParaSolverPoolElementPtr
   /// @return top priority BbParaSolverPoolElementPtr
   ///
   virtual CMapLapParaSolverPoolElementPtr top(
         ) const
   {
      return heap[1];
   }

   ///
   /// obtain top priority BbParaSolverPoolElementPtr
   /// @return top priority BbParaSolverPoolElementPtr
   ///
   virtual CMapLapParaSolverPoolElementPtr elem(
         int i
         ) const
   {
      assert( i > 0 );
      return heap[i];
   }

   ///
   /// remove top priority BbParaSolverPoolElementPtr from Selection Heap
   /// @return removed top priority BbParaSolverPoolElementPtr
   ///
   virtual CMapLapParaSolverPoolElementPtr remove(
         );

   ///
   /// resize Selection Heap
   ///
   virtual void resize(
         std::size_t size          ///< new size
         );

   ///
   /// get current used heap size
   /// @return used heap size
   ///
   virtual inline std::size_t getHeapSize(
         ) const
   {
      return heapSize;
   }

   ///
   /// get max heap size
   /// @return max heap size
   ///
   virtual inline std::size_t getMaxHeapSize(
         ) const
   {
      return maxHeapSize;
   }

   ///
   /// update selection heap by a new dual bound value of this Solver
   ///
   virtual void updateInterruptPriorityValue(
         CMapLapParaSolverPoolElementPtr solver,    ///< Solver pool element to be updated
         double value                        ///< dual bound value
         ) = 0;

   ///
   /// delete BbParaSolverPoolElementPtr from Selection Heap
   ///
   virtual void deleteElement(
         CMapLapParaSolverPoolElementPtr solver      ///< Solver pool element to be deleted
         ) = 0;

   ///
   /// up heap
   ///
   virtual void upHeap(
         std::size_t pos                      ///< start position to up heap
         ) = 0;

   ///
   /// down heap
   ///
   virtual void downHeap(
         std::size_t pos                      ///< start position to down heap
         ) = 0;

   //------------
   // for debug
   //------------
   ///
   /// stringfy of this object for debugging
   /// @return string to show inside of this object
   ///
   virtual const std::string toString(
         );

};

///
/// class AscendingSelectionHeap
///
class AscendingSelectionHeap : public SelectionHeap
{
public:

   ///
   /// constructor
   ///
   AscendingSelectionHeap(
         std::size_t size                        ///< heap size
         )
         : SelectionHeap(size)
   {
   }

   ///
   /// copy constructor
   ///
   AscendingSelectionHeap(
         const AscendingSelectionHeap &selectionHeap
         )
         : SelectionHeap(selectionHeap)
   {
   }

   ///
   /// destructor
   ///
   virtual ~AscendingSelectionHeap(
         )
   {
   }

   ///
   /// update selection heap by a new dual bound value of this Solver
   ///
   virtual void updateInterruptPriorityValue(
         CMapLapParaSolverPoolElementPtr solver,    ///< Solver pool element to be updated
         double value                        ///< dual bound value
         );

   ///
   /// delete BbParaSolverPoolElementPtr from Selection Heap
   ///
   virtual void deleteElement(
         CMapLapParaSolverPoolElementPtr solver      ///< Solver pool element to be deleted
         );

   ///
   /// up heap
   ///
   virtual void upHeap(
         std::size_t pos                      ///< start position to up heap
         );

   ///
   /// down heap
   ///
   virtual void downHeap(
         std::size_t pos                      ///< start position to down heap
         );

};



///
/// class CMapLapParaSolverPool
/// (Solver Pool basis class)
///
class CMapLapParaSolverPool : public UG::ParaSolverPool
{
   using InActiveSolvers      = std::unordered_map< int, CMapLapParaSolverPoolElementPtr >;
   using ActiveSolvers        = std::unordered_map< int, CMapLapParaSolverPoolElementPtr >;
   using InterruptingSolvers  = std::unordered_map< int, CMapLapParaSolverPoolElementPtr >;

protected:

   unsigned int                     nThreadsPerRank;              ///< maximum number of threads in a solver process
   long                             nDeepBkzTasks;                ///< total number of DeepBkz tasks assigned
   long                             nEnumTasks;                   ///< total number of Enum tasks assigned
   long                             nSieveTasks;                  ///< total number of Sieve tasks assinged
   InActiveSolvers                  inactiveSolvers;              ///< pointers to inactive Solvers
   ActiveSolvers                    activeDeepBkzSolvers;         ///< pointers to active DeepBkz Solvers
   InterruptingSolvers              interruptingDeepBkzSolvers;   ///< pointers to interruupting DeepBkz Solvers
   ActiveSolvers                    activeEnumSolvers;            ///< pointers to active Enum Solvers
   InterruptingSolvers              interruptingEnumSolvers;      ///< pointers to interruupting Enum Solvers
   ActiveSolvers                    activeSieveSolvers;           ///< pointers to active Sieve Solvers
   InterruptingSolvers              interruptingSieveSolvers;     ///< pointers to interruupting Sieve Solvers
   CMapLapParaSolverPoolElementPtr  *pool;                        ///< Solver pool indexed by Solver's rank
   AscendingSelectionHeap           *deepBkzSelectionHeap;        ///< pointers to active deepBkz Solvers in ascending order
   AscendingSelectionHeap           *enumSelectionHeap;           ///< pointers to active enum Solvers in ascending order
   AscendingSelectionHeap           *sieveSelectionHeap;          ///< pointers to active enum Solvers in ascending order

public:

   ///
   /// constructor
   ///
   CMapLapParaSolverPool(
         int inOriginRank,                   ///< origin rank of Solvers managed by this Solver pool
         UG::ParaComm *inParaComm,           ///< communicator used
         UG::ParaParamSet *inParaParams,     ///< paraParamSet used
         UG::ParaTimer *inParaTimer,         ///< timer used
         int inNThreadsPerRank               ///< maximum number of threads in a solver process
         ) : ParaSolverPool(inOriginRank, inParaComm, inParaParams, inParaTimer),
             nThreadsPerRank(inNThreadsPerRank),
             nDeepBkzTasks(0),
             nEnumTasks(0),
             nSieveTasks(0)
   {
      nSolvers = paraComm->getSize() - inOriginRank;

      pool = new CMapLapParaSolverPoolElementPtr[nSolvers*nThreadsPerRank];
      for( unsigned int i = 0; i < nSolvers; i++ )
      {
         for( unsigned int j = 0; j < nThreadsPerRank; j++ )
         {
            unsigned int solverThreadId = getSolverThreadId(originRank+i, j);
            pool[solverThreadId] = new CMapLapParaSolverPoolElement(originRank+i, j);
            inactiveSolvers.insert(std::make_pair(solverThreadId,pool[solverThreadId]));
         }
      }
      deepBkzSelectionHeap = new AscendingSelectionHeap(nSolvers*nThreadsPerRank);
      enumSelectionHeap    = new AscendingSelectionHeap(nSolvers*nThreadsPerRank);
      sieveSelectionHeap   = new AscendingSelectionHeap(nSolvers*nThreadsPerRank);
   }

   ///
   ///  destructor
   ///
   virtual ~CMapLapParaSolverPool(
         )
   {
      if( deepBkzSelectionHeap ) delete deepBkzSelectionHeap;
      if( enumSelectionHeap ) delete enumSelectionHeap;
      if( sieveSelectionHeap ) delete sieveSelectionHeap;
      for( unsigned int i = 0; i < (unsigned int)nSolvers*nThreadsPerRank; i++ )
      {
         delete pool[i];
      }
      if( pool ) delete[] pool;
   }

   ///
   /// clone object
   /// @return copied cmapLapParaSolverPool pointer
   ///
   virtual CMapLapParaSolverPool* clone(
         )
   {
      return new CMapLapParaSolverPool(*this);
   }

   ///
   /// get of solverThreadId
   /// @return threadd of (solverRank, threadId)
   ///
   virtual unsigned int getSolverThreadId(
         int solverRank,
         int threadId     ///< can be -1
         )
   {
      return ((solverRank - originRank )*nThreadsPerRank + std::max(threadId, 0));
   }

   ///
   /// get number of Solvers in this Solver pool
   /// @return number of Solvers
   ///
   virtual std::size_t getNSolvers(
         )
   {
      return nSolvers;
   }

   ///
   ///  get number of DeepBkz solver tasks
   ///  @return number of DeepBkz solver tasks
   ///
   virtual long getNDeepBkzTasks(
         )
   {
      return nDeepBkzTasks;
   }

   ///
   ///  add number of DeepBkz solver tasks
   ///
   virtual void addNDeepBkzTasks(
         int n       ///< number of new tasks
         )
   {
      nDeepBkzTasks += n;
   }

   ///
   /// get number of Enum solver tasks
   /// @return number of Enum solver tasks
   ///
   virtual long getNEnumTasks(
         )
   {
      return nEnumTasks;
   }

   ///
   ///  add number of Enum solver tasks
   ///
   virtual void addEnumTasks(
         int n       ///< number of new tasks
         )
   {
      nEnumTasks += n;
   }

   ///
   /// get number of Sieve solver tasks
   /// @return number of Sieve solver tasks
   ///
   virtual void addSieveTasks(
         int n       ///< number of new tasks
         )
   {
      nSieveTasks += n;
   }

   ///
   /// get number of active Solvers
   /// @return number of active Solvers
   ///
   virtual std::size_t getNumActiveSolvers(
         )
   {
      return (getNumActiveDeepBkzSolvers() + getNumActiveEnumSolvers() + getNumActiveSieveSolvers());
   }

   ///
   /// check if the specified Solver is active or not
   /// @return true if the specified Solver is active, false otherwise
   ///
   virtual bool isActive(
         int rank        ///< rank of the Solver
         )
   {
      bool active = false;
      for( unsigned int j = 0; j < nThreadsPerRank; j++ )
      {
         if( isSolverActive(rank, j) )
         {
            active = true;
            break;
         }
      }
      return active;
   }


   ///
   /// get number of inactive Solvers
   /// @return number of inactive Solvers
   ///
   virtual std::size_t getNumInactiveSolvers(
         )
   {
      return inactiveSolvers.size();
   }


   ///
   /// get number of active DeepBkz Solvers
   /// @return number of active DeepBkz Solvers
   ///
   virtual std::size_t getNumActiveDeepBkzSolvers(
         )
   {
      return activeDeepBkzSolvers.size() + interruptingDeepBkzSolvers.size();
   }

   ///
   /// get number of active Enum Solvers
   /// @return number of active Enum Solvers
   ///
   virtual std::size_t getNumActiveEnumSolvers(
         )
   {
      return activeEnumSolvers.size() + interruptingEnumSolvers.size();
   }

   ///
   /// get number of active Sieve Solvers
   /// @return number of active Sieve Solvers
   ///
   virtual std::size_t getNumActiveSieveSolvers(
         )
   {
      return activeSieveSolvers.size() + interruptingSieveSolvers.size();
   }

   ///
   /// check if the Solver specified by rank is active or not
   /// @return true if the Solver is active, false otherwise
   ///
   virtual bool isSolverActive(
         int rank      ///< rank of the Solver to be checked
         )
   {
      THROW_LOGICAL_ERROR2("CMapLapParaSolverPool::isSolverActive(rank) should not be called, rank = ", rank);
   }

   ///
   /// check if the Solver specified by rank is active or not
   /// @return true if the Solver is active, false otherwise
   ///
   virtual bool isSolverActive(
         int rank,      ///< rank of the Solver to be checked
         int threadId   ///< thread ID of the solver to be checked
         )
   {
      unsigned int solverThreadId = getSolverThreadId(rank, threadId);
      return pool[(solverThreadId)]->isActive();
   }

   ///
   /// get current solving CMapLapParaNode in the Solver specified by rank */
   /// @return pointer to CMapLapParaNode object
   ///
   virtual SolverType getSolverType(
         int rank,      ///< rank of the Solver to be checked
         int threadId   ///< thread ID of the solver to be checked
         )
   {
      unsigned int solverThreadId = getSolverThreadId(rank, threadId);
      return pool[(solverThreadId)]->getSolverType();
   }

   ///
   /// get current solving CMapLapParaNode in the Solver specified by rank */
   /// @return pointer to CMapLapParaNode object
   ///
   virtual UG::ParaTask *getCurrentTask(
         int rank      ///< rank of the Solver to be checked
         )
   {
      THROW_LOGICAL_ERROR2("CMapLapParaSolverPool::getCurrentNode(rank) should not be called, rank = ", rank);
   }

   ///
   /// get current solving CMapLapParaNode in the Solver specified by rank */
   /// @return pointer to CMapLapParaNode object
   ///
   virtual std::shared_ptr<CMapLapParaTask> getCurrentTask(
         int rank,      ///< rank of the Solver to be checked
         int threadId   ///< thread ID of the solver to be checked
         )
   {
      unsigned int solverThreadId = getSolverThreadId(rank, threadId);
      return  pool[(solverThreadId)]->getCurrentTask();
   }


   ///
   /// set the Solver specified by rank is interrupt requested
   ///
   virtual void interruptRequested(
         int rank     ///< rank of the Solver
         )
   {
      THROW_LOGICAL_ERROR1("CMapLapParaSolverPool::interruptRequested is not implemented");
   }

   ///
   /// check if the Solver specified by rank is interrupt requested or not
   /// @return return true if the Solver is interrupt requested, false otherwise
   ///
   virtual bool isInterruptRequested(
         int rank     ///< rank of the Solver
         )
   {
      THROW_LOGICAL_ERROR1("CMapLapParaSolverPool::isInterruptRequested is not implemented");
   }

   ///
   /// set the Solver specified by rank is terminate requested
   ///
   virtual void terminateRequested(
         int rank     ///< rank of the Solver
         )
   {
      THROW_LOGICAL_ERROR1("CMapLapParaSolverPool::terminateRequested is not implemented");
   }

   ///
   /// check if the Solver specified by rank is terminate requested or not
   /// @return return true if the Solver is terminate requested, false otherwise
   ///
   virtual bool isTerminateRequested(
         int rank     ///< rank of the Solver
         )
   {
      THROW_LOGICAL_ERROR1("CMapLapParaSolverPool::isTerminateRequested is not implemented");
   }

   ///
   /// set the Solver specified by rank is terminated
   ///
   virtual void terminated(
         int rank     ///< rank of the Solver
         )
   {
      THROW_LOGICAL_ERROR1("CMapLapParaSolverPool::terminated is not implemented");
   }

   ///
   /// check if the Solver specified by rank is terminated or not
   /// @return return true if the Solver is terminated, false otherwise
   ///
   virtual bool isTerminated(
         int rank     ///< rank of the Solver
         )
   {
      THROW_LOGICAL_ERROR1("CMapLapParaSolverPool::isTerminated is not implemented");
   }


#ifdef UG_WITH_ZLIB

   ///
   /// get copied task objects input queue
   /// @return the CMapLapParaTaskQueue pointer
   ///
   virtual CMapLapParaTaskQueue *getParaActiveTaskQueue(
         UG::ParaComm *comm      ///< communicator used
         )
   {
      CMapLapParaTaskQueue *cmapLapParaTaskQueue = new CMapLapParaTaskQueue();
      for( unsigned int i = 0; i < ((paraComm->getSize() - 1)*nThreadsPerRank); i++ )
      {
         if( pool[i]->getStatus() == UG::Active && pool[i]->getSolverType() != Sieve )
         {
            assert( getCurrentTask(pool[i]->getRank(), pool[i]->getThreadId()) );
            cmapLapParaTaskQueue->push(getCurrentTask(pool[i]->getRank(), pool[i]->getThreadId()));
         }
      }
      return cmapLapParaTaskQueue;
   }

   ///
   /// write cmapLapParaTask to checkpoint file
   /// @return the number of cmapLapParaTask saved
   ///
   virtual int writeActiveTasksToCheckpointFile(
         gzstream::ogzstream &out             ///< ogzstream to output
         )
   {
      int n = 0;
      for( unsigned int i = 0; i < ((paraComm->getSize() - 1)*nThreadsPerRank); i++ )
      {
         if( pool[i]->getStatus() == UG::Active && pool[i]->getSolverType() != Sieve )
         {
            assert( getCurrentTask(pool[i]->getRank(), pool[i]->getThreadId()) );
            getCurrentTask(pool[i]->getRank(), pool[i]->getThreadId())->write(out,1);   /// 1: means active
            n++;
         }
      }
      return n;
   }

#endif

   ///
   /// get an inactive Solver rank
   /// @return rank of inactivate Solver
   ///
   virtual int getInactiveSolverRankForThreadTask(
         int &threadId           ///< return thread ID
         );

   ///
   /// get an inactive Solver rank
   /// @return rank of inactivate Solver
   ///
   virtual int getInactiveSolverRankForRankTask(
         int &threadId           ///< return thread ID
         );

   ///
   /// this function is for restart, asssume that an inactive node exists
   /// @return true if the reservation of node for Sieve is succeeded else false
   ///
   virtual bool reserveForSieve(
         );

   ///
   /// activate a Solver with specified CMapLapParaTask which is sent within this Solver pool
   /// @return rank of Solver which is activated
   ///
   virtual int activateSolver(
         SolverType           solverType,             ///< solver type
         std::shared_ptr<CMapLapParaTask> paraTask,   ///< pointer to CMapLapParaTask object to be sent to a Solver
         int                  &threadId               ///< thread Id of solver activated (return value)
         );

   ///
   /// get interrupt DeepBkz solver
   /// @return interrupt DeepBkz solver rank, -1 is reserved for Sieve
   ///
   virtual int getInterruptDeepBkzSolver(
         int &threadId           ///< interrupt DeepBkz solver thread id (return value)
         );

   ///
   /// get interrupt Enum solver
   /// @return interrupt Enum solver rank
   ///
   virtual int getInterruptEnumSolver(
         int &threadId           ///< interrupt Enum solver thread id (return value)
         );

   ///
   /// interrupt DeepBkz Solvers
   ///
   virtual void interruptDeepBkzSolvers(
         int nInterrupt     ///< number of interrupt solvers
         );

   ///
   /// interrupt Enum Solvers
   ///
   virtual void interruptEnumSolvers(
         int nInterrupt     ///< number of interrupt solvers
         );

   ///
   /// interrupt DeepBkz and Enum Solvers to run a Sieve Solver
   ///
   virtual bool interruptDeepBkzAndEnumSolversToRunSieveSolver(
         );

   ///
   /// interrupt DeepBkz and Enum Solvers to run a Sieve Solver
   ///
   virtual bool interruptEnumAndDeepBkzSolversToRunSieveSolver(
         );

   ///
   /// check if interrupting solver exists
   ///
   virtual bool isThereInterruptingDeepBkzSolver(
         )
   {
      return ( interruptingDeepBkzSolvers.size() > 0 );
   }

   ///
   /// check if interrupting solver exists
   ///
   virtual bool isThereInterruptingEnumSolver(
         )
   {
      return ( interruptingEnumSolvers.size() > 0 );
   }

   ///
   /// inactivate the Solver specified by rank
   ///
   virtual void inactivateSolver(
         int rank,                                          ///< rank of the Solver to be inactivated
         int threadId,                                      ///< thread ID of the Solver to be inactivated
         std::shared_ptr<CMapLapParaTaskPool> paraNodePool  ///< pointer to CMapLapParaNodePool to change into collecting mode
         );

   ///
   /// getter of basis of current task
   ///
   virtual std::shared_ptr<LatticeBasis<int>> getBasis(
         int      rank,                ///< rank of the Solver
         int      threadId             ///< thread ID
         );

   ///
   /// update DeepBkz Solver status
   ///
   virtual void updateSolverStatus(
         int      rank,                ///< rank of the Solver
         int      threadId,            ///< thread ID
         int      inBlockSize,         ///< current blocksize
         int*     inBasis,              ///< current basis
         double   inEnumCost,          ///< approximated time of full Enum [Pool] better a than b when a < b
         double   inSlopeGSA,          ///< slope of GSA line [Pool] better a than b when a > b
         double   inTopHalfSlopeGSA,   ///< slope of top-half GSA line [Pool] better a than b when a > b
         double   inOrthogonalFactor   ///< orthogonal factor [Pool] better a than b when a < b
         );

   ///
   /// update Enum Solver status
   ///
   virtual void updateSolverStatus(
         int      rank,                ///< rank of the Solver
         int      threadId,            ///< thread ID
         int      inDepth,
         int*     inCoeffs,
         long int inNumSearchedNodes   ///< number of searched nodes in the enumeration tree
         );

   ///
   /// send message for specific solverType solvers
   ///
   virtual void sendMessageForSolvers(
         CMapLapParaSolverLocalComm *localComm,
         int tag,
         std::vector<bool> &hasSentManager,
         SolverType solverType
         );

   ///
   /// interrupt all active solvers
   ///
   virtual void interruptAllSolvers(
         CMapLapParaSolverLocalComm *localComm,
         int tag
         );

   ///
   /// terminate all solvers
   ///
   virtual void terminateAllSolvers(
         CMapLapParaSolverLocalComm *localComm
         );

   ///
   /// get basis of DeepBkz Solvers
   /// @return deque of basis(Eigen::LatticeBasis<int>)
   ///
   virtual std::deque<std::shared_ptr<LatticeBasis<int>>> getBasisOfSolvers(
         );
};

} // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_SOLVER_POOL_H__
