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

/**@file    cmapLapParaSolverLocalComm.h
 * @brief   CMapLapParaComm extension for CMAP-LAP
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_SOLVER_LOCAL_COMM_H__
#define __CMAP_LAP_PARA_SOLVER_LOCAL_COMM_H__


#include <vector>
#include <condition_variable>
#include <iostream>
#include <string>
#include <cassert>
#include "ug/paraSysTimer.h"
#ifdef _COMM_CPP11
#include "ug/paraTimerTh.h"
#else
#include "ug/paraTimerMpi.h"
#endif
#include "cmapLapParaTaskTh.h"
#include "cmapLapParaSolverStateTh.h"
#include "cmapLapParaCalculationStateTh.h"
#include "cmapLapParaSolverTerminationStateTh.h"
#include "cmapLapParaParamSet.h"
namespace ParaCMapLAP { class CMapLapParaBasis; }
namespace ParaCMapLAP { class CMapLapParaInstance; }
namespace ParaCMapLAP { class CMapLapParaPackedVector; }
namespace ParaCMapLAP { class CMapLapParaSolution; }
namespace ParaCMapLAP { class CMapLapParaSolverLocalComm; }
namespace ParaCMapLAP { class VectorElement; }


#define TAG_THREAD_TRACE( call, fromTo, sourceDest, tag ) \
   if( tagTraceFlag )  \
   {  \
      if( tag >= 0 ) \
      { \
         /* std::cout << " call = " << #call << ", Thread = " << getThread() << ", tag = " << tag << ", " << tagStringTable[tag] << std::endl; */ \
         *getOstream() << timer.getRTimeInterval() << " [Thread = " << getThreadId() << "] " << #call << " " << #fromTo  \
         << " " << sourceDest << " with Tag = " << getTagString(tag) << std::endl; \
      } \
      else \
      { \
         /* std::cout << " call = " << #call << ", Thread = " << getThread() << ", tag = " << tag << std::endl; */  \
         *getOstream() << timer.getRTimeInterval() << " [Thread = " << getThreadId() << "] " << #call << " " << #fromTo  \
         << " " << sourceDest << " as broadcast" << std::endl; \
      } \
}

namespace ParaCMapLAP
{

///
/// user defined transfer data types
///
static const int UG_USER_TYPE_FIRST                   = UG::TYPE_LAST + 1;
static const int ParaInstanceType                     = UG_USER_TYPE_FIRST +  0;
static const int ParaSolutionType                     = UG_USER_TYPE_FIRST +  1;
static const int ParaParamSetType                     = UG_USER_TYPE_FIRST +  2;
static const int ParaTaskType                         = UG_USER_TYPE_FIRST +  3;
static const int ParaSolverStateType                  = UG_USER_TYPE_FIRST +  4;
static const int ParaCalculationStateType             = UG_USER_TYPE_FIRST +  5;
static const int ParaSolverTerminationStateType       = UG_USER_TYPE_FIRST +  6;
static const int ParaRacingRampUpParamType            = UG_USER_TYPE_FIRST +  7;
static const int UG_USER_TYPE_LAST                    = UG_USER_TYPE_FIRST +  7;

#ifndef __CMAP_LAP_PARA_COMM_TH_H__
static const int CMAP_LAP_USER_TYPE_FIRST             = UG_USER_TYPE_LAST + 1;
static const int ParaBasisType                        = CMAP_LAP_USER_TYPE_FIRST + 0;
static const int CMapLapParaPackedVectorType          = CMAP_LAP_USER_TYPE_FIRST + 1;
static const int ParaCheckpointElementType            = CMAP_LAP_USER_TYPE_FIRST + 2;
static const int CMAP_LAP_USER_TYPE_LAST              = CMAP_LAP_USER_TYPE_FIRST + 2;
#endif

static const int ThreadTableSize = 751;     ///< size of thread table : this limits the number of threads

///
/// Class for message queue element
///
/// NOTE : For basic data types, this is copy of sender side memory.
///        When the memory is copied at receive function, the memory
///        have to be freed.
///        For user defined data type, this is the receiver side memory,
///        because it is better to allocate memory in the sender side for
///        mutex locking. Sender side functions have to allocate memory.
///        In this case, memory do not have to be freed. The memory is for
///        receiver side.
///
class LocalMessageQueueElement
{

   int source;                 ///< source thread thread of this message
   int count;                  ///< number of the data type elements
   int dataTypeId;             ///< data type id
   int tag;                    ///< tag of the message,  -1 : in case of broadcast message
   void *data;                 ///< data of the message
   LocalMessageQueueElement *next;  ///< point to next message queue element

public:

   ///
   /// default constructor of LocalMessageQueueElement
   ///
   LocalMessageQueueElement(
         )
         : source(-1),
           count(0),
           dataTypeId(-1),
           tag(-1),
           data(0),
           next(0)
   {
   }

   ///
   /// constructor of LocalMessageQueueElement
   ///
   LocalMessageQueueElement(
         int inSource,         ///< source thread thread of this message
         int inCount,          ///< number of the data type elements
         int inDataTypeId,     ///< data type id
         int inTag,            ///< tag of the message,  -1 : in case of broadcast message
         void *inData          ///< data of the message
         ) : source(inSource), count(inCount), dataTypeId(inDataTypeId), tag(inTag), data(inData), next(0)
   {
   }

   ///
   /// destructor of LocalMessageQueueElement
   ///
   virtual ~LocalMessageQueueElement(
         )
   {
   }

   ///
   /// getter of source thread
   /// @return thread of the source
   ///
   virtual int getSource(
         )
   {
      return source;
   }

   ///
   /// getter of the number of the data type elements
   /// @return the number of the data type elements
   ///
   virtual int getCount(
         )
   {
      return count;
   }

   ///
   /// getter of the data type id
   /// @return data type id
   ///
   virtual int getDataTypeId(
         )
   {
      return dataTypeId;
   }

   ///
   /// getter of the message tag
   /// @return tag of the message
   ///
   virtual int getTag(
         )
   {
      return tag;
   }

   ///
   /// getter of data
   /// @return pointer to the data
   ///
   virtual void *getData(
         )
   {
      return data;
   }

   ///
   /// getter of the pointer to the next LocalMessageQueueElement
   /// @return pointer to LocalMessageQueueElement
   ///
   virtual LocalMessageQueueElement *getNext(
         )
   {
      return next;
   }

   ///
   /// link to the next LocalMessageQueueElement
   ///
   virtual void link(
         LocalMessageQueueElement *nextElement  ///< pointer to LocalMessageQueueElement
         )
   {
      next = nextElement;
   }
};

///
/// Class of LocalMessageQueueTableElement
///
class LocalMessageQueueTableElement
{

   LocalMessageQueueElement *head;      ///< head of the message queue
   LocalMessageQueueElement *tail;      ///< tail of the message queue
   int                       size;      ///< number of the messages in queue

public:

   ///
   /// default constructor of LocalMessageQueueTableElement
   ///
   LocalMessageQueueTableElement(
           ) : head(0), tail(0), size(0)
   {
   }

   ///
   /// destructor of LocalMessageQueueTableElement
   ///
   virtual ~LocalMessageQueueTableElement(
         )
   {
      while( head )
      {
         LocalMessageQueueElement *next = head->getNext();
         delete head;
         head = next;
      }
   }

   ///
   /// check if the specified message exists or nor
   /// @return pointer to LocalMessageQueueElement, 0: no element
   ///
   virtual LocalMessageQueueElement *checkElement(
         int source,        ///< source thread
         int datatypeId,    ///< data type id
         int tag            ///< tag of the message
         )
   {
      LocalMessageQueueElement *ret = 0;
      for( LocalMessageQueueElement *current = head; current; current = current->getNext() )
      {
         if( current->getSource() == source
               && current->getDataTypeId() == datatypeId
               && current->getTag() == tag )
         {
            ret = current;
            break;
         }
      }
      return ret;
   }

   ///
   /// extracts a message
   /// @return pointer to the message
   ///
   virtual LocalMessageQueueElement *extarctElement(
         bool *sentMessage,  ///< for synchronization
         int source,         ///< source thread
         int datatypeId,     ///< data type id
         int tag             ///< tag of the message
         )
   {
      LocalMessageQueueElement *ret = 0;

      LocalMessageQueueElement *prev = head;
      LocalMessageQueueElement *current = head;
      while( current )
      {
         LocalMessageQueueElement *next = current->getNext();
         if( current->getSource() == source
               && current->getDataTypeId() == datatypeId
               && current->getTag() == tag )
         {
            ret = current;
            ret->link(0);
            if( current == head )
            {
               if( current == tail )
               {
                  head = 0;
                  tail = 0;
                  *sentMessage = false;
               }
               else
               {
                  head=next;
               }
            }
            else
            {
               if( current == tail )
               {
                  tail = prev;
               }
               prev->link(next);
            }
            size--;
            break;
         }
         else
         {
            prev = current;
            current = next;
         }
      }
      return ret;
   }

   ///
   /// extracts a message
   /// (This method is only for desctructor of CMapLapParaSolverLocalComm. No lock is necessary.)
   ///
   /// @return pointer to the message
   ///
   virtual LocalMessageQueueElement *extarctElement(
         bool              *sentMessage   ///< for synchronization
         )
   {
      if( head == 0 ) return 0;

      LocalMessageQueueElement *current = head;
      LocalMessageQueueElement *next = current->getNext();
      current->link(0);
      if( current == tail )
      {
         head = 0;
         tail = 0;
         *sentMessage = false;
      }
      else
      {
         head=next;
      }
      size--;
      return current;
   }

   ///
   /// enqueue a message
   ///
   virtual void enqueue(
         std::condition_variable& sentMsg,         ///< condition variable for synchronization
         std::mutex&              queueLockMutex,  ///< mutex for synchronization
         bool                     *sentMessage,    ///< flag for synchronization
         LocalMessageQueueElement      *newElement      ///< message queue element to enter
         )
   {
      if( tail == 0 )
      {
         head = tail = newElement;
      }
      else
      {
         tail->link(newElement);
         tail = newElement;
      }
      size++;
      *sentMessage = true;
      sentMsg.notify_one();
   }

   ///
   /// getter of head
   /// @return head element pointer
   ///
   virtual LocalMessageQueueElement *getHead(
         )
   {
      return head;
   }

   ///
   /// getter of size
   /// @return size of the message queue
   ///
   virtual int getSize(
         )
   {
      return size;
   }

   ///
   /// check if the queue is empty or not
   /// @return true, if it is empty, else false
   ///
   virtual bool empty(
         )
   {
      bool __empty = true;
      if( head ) __empty = false;
      assert( !__empty || size == 0 );
      return __empty;
   }

   ///
   /// wait for a message coming to a queue
   ///
   virtual void waitMessage(
         std::condition_variable& sentMsg,          ///< condition variable for synchronization
         std::mutex&              queueLockMutex,   ///< mutex for synchronization
         bool                     *sentMessage      ///< flag for synchronization
         )
   {
      std::unique_lock<std::mutex> lk(queueLockMutex);
      sentMsg.wait(lk, [&sentMessage]{return (*sentMessage == true);});
   }

   ///
   /// wait for a specified message coming to a queue
   ///
   virtual void waitMessage(
         std::condition_variable& sentMsg,         ///< condition variable for synchronization
         std::mutex&              queueLockMutex,  ///< mutex for synchronization
         bool                     *sentMessage,    ///< flag for synchronization
         int                      source,          ///< source thread of the message
         int                      tag              ///< tag of the message
         )
   {
      for(;;)
      {
         std::unique_lock<std::mutex> lk(queueLockMutex);
         sentMsg.wait(lk, [&sentMessage]{return (*sentMessage == true);});
         LocalMessageQueueElement *current = head;
         while( current )
         {
            LocalMessageQueueElement *next = current->getNext();
            if( current->getSource() == source
                  && current->getTag() == tag )
            {
               break;
            }
            current = next;
         }
         if( current ) break;
         *sentMessage = false;
      }
   }

   ///
   /// wait for a specified message coming to a queue
   ///
   virtual void waitMessage(
         std::condition_variable& sentMsg,         ///< condition variable for synchronization
         std::mutex&              queueLockMutex,  ///< mutex for synchronization
         bool                     *sentMessage,    ///< flag for synchronization
         int                      source,          ///< source thread of the message
         int                      datatypeId,      ///< data type id of the message
         int                      tag              ///< tag of the message
         )
   {
      for(;;)
      {
         std::unique_lock<std::mutex> lk(queueLockMutex);
         sentMsg.wait(lk, [&sentMessage]{return (*sentMessage == true);});
         LocalMessageQueueElement *current = head;
         while( current )
         {
            LocalMessageQueueElement *next = current->getNext();
            if( current->getSource() == source
                  && current->getDataTypeId() == datatypeId
                  && current->getTag() == tag )
            {
               break;
            }
            current = next;
         }
         if( current ) break;
         *sentMessage = false;
      }
   }

   ///
   /// wait for a specified message coming to a queue
   ///
   virtual void waitMessage(
         std::condition_variable& sentMsg,         ///< condition variable for synchronization
         std::mutex&              queueLockMutex,  ///< mutex for synchronization
         bool                     *sentMessage,    ///< flag for synchronization
         int                      source,          ///< source thread of the message
         int                      *tag             ///< tag of the message
         )
   {
      for(;;)
      {
         std::unique_lock<std::mutex> lk(queueLockMutex);
         sentMsg.wait(lk, [&sentMessage]{return (*sentMessage == true);});
         LocalMessageQueueElement *current = head;
         while( current )
         {
            LocalMessageQueueElement *next = current->getNext();
            if( current->getSource() == source )
            {
               *tag = current->getTag();
               break;
            }
            current = next;
         }
         if( current ) break;
         *sentMessage = false;
      }
   }

};

///
/// Class of LocalThreadsTableElement
///
class LocalThreadsTableElement
{

   int                thread;    ///< thread of this thread
   std::ostream       *tos;      ///< tag trace stream for this thread

public:

   ///
   /// default constructor of LocalThreadsTableElement
   ///
   LocalThreadsTableElement(
         )
         : thread(0),
           tos(0)
   {
   }

   ///
   /// constructor of LocalThreadsTableElement
   ///
   LocalThreadsTableElement(
         int       inThread,              ///< thread of this thread
         UG::ParaParamSet *paraParams     ///< UG parameter set
         )
         : thread(inThread),
           tos(0)
   {
      if( paraParams->getBoolParamValue(UG::TagTrace) )
      {
         if( paraParams->isStringParamDefaultValue(UG::TagTraceFileName) )
         {
            tos = &std::cout;
         }
         else
         {
            std::ostringstream s;
            s << paraParams->getStringParamValue(UG::TagTraceFileName) << inThread;
            std::ofstream  *ofs = new std::ofstream();
            ofs->open(s.str().c_str());
            tos = ofs;
         }
      }
   }

   ///
   /// destructor of LocalThreadsTableElement
   ///
   virtual ~LocalThreadsTableElement(
         )
   {
      if( tos )
      {
         if( tos != &std::cout )
         {
            std::ofstream *ofs = dynamic_cast<std::ofstream *>(tos);
            ofs->close();
            delete ofs;
         }
      }
   }

   ///
   /// getter of this thread thread
   /// @return the thread thread
   ///
   virtual int getThread(
         )
   {
      return thread;
   }

   ///
   /// setter of this thread thread
   ///
   virtual void setThread(
         int r     ///< thread to be set
         )
   {
      thread = r;
   }

   ///
   /// getter of tag trace stream of this thread
   /// @return ponter to the tag trace stream
   ///
   virtual std::ostream *getOstream(
         )
   {
      return tos;
   }
};


///
/// Communicator object for C++11 thread communications
///
class CMapLapParaSolverLocalComm : public UG::ParaComm
{
protected:
   int                        comSize;                         ///< communicator size : number of threads joined in this system
   bool                       tagTraceFlag;                    ///< indicate if tags are traced or not
   int                        **token;                         ///< index 0: token
                                                               ///< index 1: token color
                                                               ///<      -1: green
                                                               ///<     > 0: yellow ( termination origin solver number )
                                                               ///<      -2: red ( means the solver can terminate )
   UG::ParaSysTimer           timer;                           ///< timer
   static const char          *tagStringTable[];               ///< tag name string table
   LocalThreadsTableElement   *threadsTable[ThreadTableSize];  ///< threads table: index is thread thread
   static thread_local int    localThread;                       ///< local thread thread
   LocalMessageQueueTableElement   **messageQueueTable;             ///< message queue table
   bool                       *sentMessage;                    ///< sent message flag for synchronization
   std::mutex                 *queueLockMutex;                 ///< mutex for synchronization
   std::condition_variable    *sentMsg;                        ///< condition variable for synchronization

   std::mutex                 *tokenAccessLockMutex;           ///< mutex to access token
   std::mutex                 applicationLockMutex;            ///< mutex for applications
   std::mutex                 threadLockMutex;                 ///< mutex to access thread

   ///
   /// allocate memory and copy message
   ///
   virtual void *allocateMemAndCopy(
         const void* buffer,    ///< pointer to buffer of the message
         int count,             ///< the number of data element in the message
         const int datatypeId   ///< data type of each element in the message
         );

   ///
   /// copy message
   ///
   virtual void copy(
         void *dest,           ///< destination to copy the data
         const void *src,      ///< source of the data
         int count,            ///< the number of data element
         int datatypeId        ///< data type of each element
         );

   ///
   /// free memory
   ///
   virtual void freeMem(
         void* buffer,         ///< pointer to buffer of the message
         int count,            ///< the number of data element
         const int datatypeId  ///< data type of each element
         );

   ///
   /// free memory
   /// @return
   ///
   virtual bool freeStandardTypes(
         LocalMessageQueueElement *elem   ///< pointer to a message queue element
         );

   ///
   /// check if tag string table (for debugging) set up correctly
   /// @return true if tag string table is set up correctly, false otherwise
   ///
   virtual bool tagStringTableIsSetUpCorrectly(
         );

public:

   ///
   /// constructor of ParaComCPP11
   ///
   CMapLapParaSolverLocalComm(
         )
         : comSize(-1),
           tagTraceFlag(false),
           token(0),
           messageQueueTable(0),
           sentMessage(0),
           queueLockMutex(0),
           sentMsg(0),
           tokenAccessLockMutex(0)
   {
   }

   ///
   /// destructor of this communicator
   ///
   virtual ~CMapLapParaSolverLocalComm(
         );

   ///
   /// initializer of this communicator
   ///
   virtual void init(
         int argc,             ///< the number of arguments
         char **argv           ///< pointers to the arguments
         );

   ///
   /// get start time of this communicator
   /// (should not be used)
   /// @return start time
   ///
   virtual double getStartTime(
         )
   {
      return timer.getStartTime();
   }

   ///
   /// get rank of this process or this thread depending on run-time environment
   /// @return rank
   ///
   virtual int getRank(
         )
   {
      THROW_LOGICAL_ERROR1("*** getRank() is called in CMapLapParaSolverLocalComm class ***");
   }

   ///
   /// get thread ID of caller's thread
   /// @return thread ID of caller's thread
   ///
   virtual int getThreadId(
         );

   ///
   /// get size of this communicator, which indicates how many threads in a UG process
   /// @return the number of threads
   ///
   virtual int getSize(
         )
   {
      return comSize;
   }

   ///
   /// get size of the messageQueueTable
   /// @return if dest >= 0 then return the number of only messages waiting to send to dest,
   ///         othrewise return the number of all messages waiting to send.
   ///
   virtual int getNumOfMessagesWaitingToSend(
         int dest=-1
         )
   {
      if( dest >= 0 )
      {
         return messageQueueTable[dest]->getSize();
      }
      else
      {
         int nTotalMessages = 0;
         for( int i = 0; i < comSize; i++ )
         {
            nTotalMessages += messageQueueTable[i]->getSize();
         }
         return nTotalMessages;
      }
   }

   ///
   /// initializer for LoadCoordinator
   ///
   virtual void lcInit(
         UG::ParaParamSet *paraParams     ///< UG parameter set
         );

   ///
   /// initializer for main thread (threadId == 0)
   ///
   virtual void mainThreadInit(
         UG::ParaParamSet *paraParams     ///< UG parameter set
         )
   {
      lcInit(paraParams);
   }

   ///
   /// initializer for Solvers
   ///
   virtual void solverInit(
         UG::ParaParamSet *paraParams     ///< UG parameter set
         )
   {
   }

   ///
   /// initializer for a specific Solver
   ///
   virtual void solverInit(
         int thread,                      ///< thread of the Solver
         UG::ParaParamSet *paraParams     ///< UG parameter set
         );

   ///
   /// initializer for slave thread (threadId != 0)
   ///
   virtual void slaveThreadInit(
         int thread,                      ///< thread of the Solver
         UG::ParaParamSet *paraParams     ///< UG parameter set
         )
   {
      solverInit(thread, paraParams);
   }

   ///
   /// delete Solver from this communicator
   ///
   virtual void solverDel(
         int thread     ///< thread of the Solver
         )
   {
      THROW_LOGICAL_ERROR1("*** solverDel() is called in CMapLapParaSolverLocalComm class ***");
   }

   ///
   /// abort. How it works sometimes depends on communicator used
   ///
   virtual void abort(
         )
   {
      std::abort();
   }

   ///
   /// function to wait Terminated message
   /// (This function is not used currently)
   /// @return true when MPI communication is used, false when thread communication used
   ///
   virtual bool waitTerminatedMessage(
         )
   {
      return false;
   }

   ///
   /// wait token when UG runs with deterministic mode
   /// @return true, when token is arrived to the thread
   ///
   virtual bool waitToken(
         int thread        ///< thread to check if token is arrived
         )
   {
      THROW_LOGICAL_ERROR1("*** waitToken() is called in CMapLapParaSolverLocalComm class ***");
   }

   ///
   /// pass token to from the thread to the next
   ///
   virtual void passToken(
         int thread       ///< from this thread, the token is passed
         )
   {
      THROW_LOGICAL_ERROR1("*** passToken() is called in CMapLapParaSolverLocalComm class ***");
   }

   ///
   /// pass termination token from the thread to the next
   /// @return true, when the termination token is passed from this thread, false otherwise
   ///
   virtual bool passTermToken(
         int thread      ///< from this thread, the termination token is passed
         )
   {
      THROW_LOGICAL_ERROR1("*** passTermToken() is called in CMapLapParaSolverLocalComm class ***");
   }

   ///
   /// set received token to this communicator
   ///
   virtual void setToken(
         int thread,     ///< thread to set the token
         int *inToken  ///< token to be set
         )
   {
      THROW_LOGICAL_ERROR1("*** setToken() is called in CMapLapParaSolverLocalComm class ***");
   }

   ///
   /// get ostream pointer
   /// @return pointer to ostram
   ///
   virtual std::ostream *getOstream(
         );

   ///
   /// lock UG application to synchronize with other threads
   ///
   virtual void lockApp(
         )
   {
      applicationLockMutex.lock();
   }

   ///
   /// unlock UG application to synchronize with other threads
   ///
   virtual void unlockApp(
         )
   {
      applicationLockMutex.unlock();
   }

   ///
   /// lock thread
   ///
   virtual void lockThread()
   {
      threadLockMutex.lock();
   }

   ///
   /// unlock thread
   ///
   virtual void unlockThread()
   {
      threadLockMutex.unlock();
   }

   ///
   /// lock thread
   /// (for debugging)
   ///
   virtual void lockThread(
         char const *f,   ///< string to indicate what is locked
         int l            ///< a number to show something
         )
   {
      threadLockMutex.lock();
      // std::cout << "Thread Locked: " << f << ", " << l << std::endl;
   }

   ///
   /// unlock thread
   /// (for debugging)
   ///
   virtual void unlockThread(
         char const *f,    ///< string to indicate what is locked
         int l             ///< a number to show something
         )
   {
      threadLockMutex.unlock();
      // std::cout << "Thread Unlocked: " << f << ", " << l << std::endl;
   }

   ///
   /// create ParaTimer object
   /// @return pointer to ParaTimer object
   ///
   virtual UG::ParaTimer *createParaTimer(
         )
   {
#ifdef _COMM_CPP11
      return new UG::ParaTimerTh();
#else
      return new UG::ParaTimerMpi();
#endif
   }

   virtual UG::ParaInstance *createParaInstance(
         );

   virtual UG::ParaSolution *createParaSolution(
         );

   ///
   /// create ParaCalculationState object by default constructor
   /// @return pointer to ParaCalculationState object
   ///
   virtual UG::ParaCalculationState *createParaCalculationState(
         )
   {
      return new CMapLapParaCalculationStateTh();
   }

   ///
   /// create CMapLapParaCalculationState DeepBkz object
   /// @return pointer to CMapLapParaCalculationState object
   ///
   virtual UG::ParaCalculationState *createParaCalculationState(
         int          inTermState,          ///< termination status, 0: normal, -1: interrupted
         int          inThreadId,           ///< thread id
         int          inCurrentBlockSize,   ///< current DeepBkz block size
         int          inTour,               ///< number of DeepBkz loop
         double       inElapsedTime,        ///< elapsed time
         double       inShortestNorm,       ///< the shortest norm found
         double       inApproxFactor,       ///< approximated factor
         double       inHermiteFactor,      ///< hermite factor
         double       inRootHermiteFactor,  ///< (hermite factor)^(1/dim)
         double       inEnumCost,           ///< log of approximted nodes of enumeration tree with incumbent radius
         double       inEnumCostGH,         ///< log of approximted nodes of enumeration tree with GH radius
         double       inSlopeGSA,           ///< slope of GSA line
         double       inTopHalfSlopeGSA,    ///< slope of top-half GSA line
         double       inOrthogonalFactor    ///< orthogonal factor
         )
   {
      return new CMapLapParaCalculationStateTh(
            inTermState,
            inThreadId,
            inCurrentBlockSize,
            inTour,
            inElapsedTime,
            inShortestNorm,
            inApproxFactor,
            inHermiteFactor,
            inRootHermiteFactor,
            inEnumCost,
            inEnumCostGH,
            inSlopeGSA,
            inTopHalfSlopeGSA,
            inOrthogonalFactor
            );
   }


   ///
   /// create CMapLapParaCalculationState Enum object
   /// @return pointer to CMapLapParaCalculationState object
   ///
   virtual UG::ParaCalculationState *createParaCalculationState(
         int      inTermState,         ///< termination status, 0: normal, -1: interrupted
         int      inThreadId,          ///< thread id
         double   inElapsedTime,       ///< elapsed time
         double   inShortestNorm,      ///< the shortest norm found
         double   inApproxFactor,      ///< approximated factor
         double   inHermiteFactor,     ///< hermite factor
         double   inRootHermiteFactor, ///< (hermite factor)^(1/dim)
         long int inNSearch            ///< number of search nodes
         )
   {
      return new CMapLapParaCalculationStateTh(
            inTermState,
            inThreadId,
            inElapsedTime,
            inShortestNorm,
            inApproxFactor,
            inHermiteFactor,
            inRootHermiteFactor,
            inNSearch
            );
   }


   ///
   /// create CMapLapParaCalculationState Sieve object
   /// @return pointer to CMapLapParaCalculationState object
   ///
   virtual UG::ParaCalculationState *createParaCalculationState(
         int      inTermState,         ///< termination status, 0: normal, -1: interrupted
         int      inThreadId,          ///< thread id
         double   inElapsedTime,       ///< elapsed time
         int      inBlockSize,         ///< block size
         long int inNLoop,             ///< number of Sieve algorithm loop
         int      inListSize,          ///< size of List L
         int      inStackSize,         ///< size of Stack S
         int      inMaxListSize,       ///< maximum size of List L up to the point
         int      inNCollisions,       ///< number of collision
         double   inShortestNorm,      ///< the shortest norm found
         double   inApproxFactor,      ///< approximated factor
         double   inHermiteFactor,     ///< hermite factor
         double   inRootHermiteFactor  ///< (hermite factor)^(1/dim)
         )
   {
      return new CMapLapParaCalculationStateTh(
            inTermState,
            inThreadId,
            inElapsedTime,
            inBlockSize,
            inNLoop,
            inListSize,
            inStackSize,
            inMaxListSize,
            inNCollisions,
            inShortestNorm,
            inApproxFactor,
            inHermiteFactor,
            inRootHermiteFactor
            );
   }


   ///
   /// create ParaRacingRampUpParamSet object (Currently not used)
   /// @return pointer to ParaRacingRampUpParamSet object
   ///
   virtual UG::ParaRacingRampUpParamSet* createParaRacingRampUpParamSet(
         )
   {
      throw "** createParaRacingRampUpParamSet is called **";
   }

   ///
   /// create ParaNode object by default constructor
   /// @return pointer to ParaNode object
   ///
   virtual UG::ParaTask *createParaTask(
         )
   {
      return new CMapLapParaTaskTh();
   }

   ///
   /// create ParaNode DeepBkz object by constructor
   /// @return pointer to ParaNode object
   ///
   virtual UG::ParaTask *createParaNode(
         UG::TaskId inTaskId,                ///< task id
         UG::TaskId inGeneratorTaskId,       ///< generator task id
         double inEstimatedValue,            ///< estimated value
         int inThreadId,                     ///< thread ID
         int inBegin,                        ///< 1st index of block matrix
         int inEnd,                          ///< last index of block matrix
         int inBlockSize,                    ///< the current blocksize
         int inU,                            ///< the unimodular matrix size
         int inSeed,                         ///< the seed of randomize
         std::shared_ptr<LatticeBasis<int>> inBasis   ///< lattice basis
         )
   {
      return new CMapLapParaTaskTh(inTaskId, inGeneratorTaskId, inEstimatedValue, inThreadId,
            inBegin, inEnd, inBlockSize, inU, inSeed, inBasis);
   }

   ///
   /// create ParaNode Enum object by constructor
   /// @return pointer to ParaNode object
   ///
   virtual UG::ParaTask *createParaNode(
         UG::TaskId inTaskId,                ///< task id
         UG::TaskId inGeneratorTaskId,       ///< generator task id
         double   inEstimatedValue,          ///< estimated value
         int      inThreadId,                ///< thread ID
         int      inBegin,                   ///< 1st index of block matrix
         int      inEnd,                     ///< last index of block matrix
         int      inStart,                   ///< the index of current search node in the enumeration tree
         int      inLast,                    ///< search will be terminated when node reaches last-index depth node
         double   inProb,                    ///< the probability of the extreme pruning
         std::shared_ptr<LatticeVector<int>> inCoeffs,   ///< coefficients of basis vectors
         std::shared_ptr<LatticeBasis<int>> inBasis      ///< lattice basis
         )
   {
      return new CMapLapParaTaskTh(inTaskId, inGeneratorTaskId, inEstimatedValue, inThreadId,
            inBegin, inEnd, inStart, inLast, inProb, inCoeffs, inBasis);
   }

   ///
   /// create ParaSolverState object by default constructor
   /// @return pointer to ParaSolverState object
   ///
   virtual UG::ParaSolverState *createParaSolverState(
         )
   {
      return new CMapLapParaSolverStateTh();
   }

   ///
   /// create ParaSolverTerminationState object by default constructor
   /// @return pointer to ParaSolverTerminationState object
   ///
   virtual UG::ParaSolverTerminationState *createParaSolverTerminationState(
         )
   {
      return new CMapLapParaSolverTerminationStateTh();
   }

   virtual CMapLapParaInstance *createCMapLapParaInstance(
         char      *inProbFileName,
         char      *inProbName
         );

   virtual CMapLapParaInstance *createCMapLapParaInstance(
         char      *probName
         );

   virtual CMapLapParaSolution *createCMapLapParaSolution(
         int inThreadid,               ///< thread Id
         LatticeVector<int> inV,       ///< lattice vector
         double inObjValue             ///< objective function value
         );

   ///
   /// @brief create ParaBasis object
   ///
   virtual CMapLapParaBasis *createCMapLapParaBasis(
         int inThreadId,               ///< thread Id
         double inEnumCost,            ///< enumeration cost
         LatticeBasis<int> inBasis     ///< lattice basis
         );

   ///
   /// @brief create ParaBasis object
   ///
   virtual CMapLapParaBasis *createCMapLapParaBasis(
         );

   virtual CMapLapParaPackedVector *createCMapLapParaPackedVector(
         int      inThreadId,
         LatticeBasis<int> &inVectors,
         int      inRank,
         int      nThreadsPerRank
         );

   virtual CMapLapParaPackedVector *createCMapLapParaPackedVector(
         int inThreadId,
         std::vector<std::shared_ptr<VectorElement>> &inVectorElements
         );

   virtual CMapLapParaPackedVector *createCMapLapParaPackedVector(
         int inThreadId,
         int inNVectors
         );

   virtual CMapLapParaPackedVector *createCMapLapParaPackedVector(
         );

   ///
   /// broadcast function for standard ParaData types
   /// @return always 0 (for future extensions)
   ///
   virtual int bcast(
         void* buffer,            ///< point to the head of sending message
         int count,               ///< the number of data in the message
         const int datatypeId,    ///< data type in the message
         int root                 ///< root thread for broadcasting
         );

   ///
   /// send function for standard ParaData types
   /// @return always 0 (for future extensions)
   ///
   virtual int send(
         void* bufer,             ///< point to the head of sending message
         int count,               ///< the number of data in the message
         const int datatypeId,    ///< data type in the message
         int dest,                ///< destination to send the message
         const int tag            ///< tag of this message
         );

   ///
   /// receive function for standard ParaData types
   /// @return always 0 (for future extensions)
   ///
   virtual int receive(
         void* bufer,             ///< point to the head of receiving message
         int count,               ///< the number of data in the message
         const int datatypeId,    ///< data type in the message
         int source,              ///< source of the message coming from
         const int tag            ///< tag of the message
         );

   ///
   /// wait function for a specific tag from a specific source coming from
   /// @return always 0 (for future extensions)
   ///
   virtual void waitSpecTagFromSpecSource(
         const int source,         ///< source thread which the message should come from
         const int tag,            ///< tag which the message should wait
         int *receivedTag          ///< tag of the message which is arrived
         );

   ///
   /// probe function which waits a new message
   ///
   virtual bool probe(
         int *source,             ///< source thread of the message arrived
         int *tag                 ///< tag of the message arrived
         );

   ///
   /// iProbe function which checks if a new message is arrived or not
   /// @return true when a new message exists
   ///
   virtual bool iProbe(
         int *source,             ///< source thread of the message arrived
         int *tag                 ///< tag of the message arrived
         );

   ///
   /// User type send for created data type
   /// @return always 0 (for future extensions)
   ///
   virtual int uTypeSend(
         void* bufer,              ///< point to the head of sending message
         const int datatypeId,     ///< created data type
         int dest,                 ///< destination thread
         int tag                   ///< tag of the message
         );

   ///
   /// User type receive for created data type
   /// @return always 0 (for future extensions)
   ///
   virtual int uTypeReceive(
         void** bufer,             ///< point to the head of receiving message
         const int datatypeId,     ///< created data type
         int source,               ///< source thread
         int tag                   ///< tag of the message
         );

   ///
   /// get Tag string for debugging
   /// @return string which shows Tag
   ///
   virtual const char *getTagString(
         int tag                /// tag to be converted to string
         );

};

#define LOCK_THREAD( comm ) comm->lockThread(__FILE__, __LINE__)
#define UNLOCK_THREAD( comm ) comm->unlockThread(__FILE__, __LINE__)

} // namespace ParaCMapLAP

#endif  // __CMAP_LAP_PARA_SOLVER_LOCAL_COMM_H__
