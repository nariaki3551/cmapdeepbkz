/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*          This file is part of the program and software framework          */
/*                    UG --- Ubquity Generator Framework                     */
/*                                                                           */
/*  Copyright Written by Yuji Shinano <shinano@zib.de>,                      */
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

/**@file    paraCommPth.h
 * @brief   ParaComm extension for Pthreads communication
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_COMM_PTH_H__
#define __PARA_COMM_PTH_H__

#include <pthread.h>
#include <stdexcept>
#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cassert>
#include "paraDef.h"
#include "paraComm.h"
#include "paraSysTimer.h"
#include "paraParamSetTh.h"
#include "paraTimerTh.h"
#include "paraPthLock.h"
#include "paraPthCondVar.h"

#define HashEntry( tid ) ( hashCode(tid) % UG::HashTableSize )

#define TAG_TRACE( call, fromTo, sourceDest, tag ) \
   if( tagTraceFlag )  \
   {  \
      if( tag >= 0 ) \
      { \
         /* std::cout << " call = " << #call << ", Rank = " << getRank() << ", tag = " << tag << ", " << tagStringTable[tag] << std::endl; */ \
         *getOstream() << timer.getRTimeInterval() << " [Rank = " << getRank() << "] " << #call << " " << #fromTo  \
         << " " << sourceDest << " with Tag = " << getTagString(tag) << std::endl; \
      } \
      else \
      { \
         /* std::cout << " call = " << #call << ", Rank = " << getRank() << ", tag = " << tag << std::endl; */  \
         *getOstream() << timer.getRTimeInterval() << " [Rank = " << getRank() << "] " << #call << " " << #fromTo  \
         << " " << sourceDest << " as broadcast" << std::endl; \
      } \
}

namespace UG
{

///
/// user defined transfer data types
///
static const int UG_USER_TYPE_FIRST             = TYPE_LAST + 1;
static const int ParaInstanceType               = UG_USER_TYPE_FIRST +  0;
static const int ParaSolutionType               = UG_USER_TYPE_FIRST +  1;
static const int ParaParamSetType               = UG_USER_TYPE_FIRST +  2;
static const int ParaTaskType                   = UG_USER_TYPE_FIRST +  3;
static const int ParaSolverStateType            = UG_USER_TYPE_FIRST +  4;
static const int ParaCalculationStateType       = UG_USER_TYPE_FIRST +  5;
static const int ParaSolverTerminationStateType = UG_USER_TYPE_FIRST +  6;
static const int ParaRacingRampUpParamType      = UG_USER_TYPE_FIRST +  7;
// static const int ParaSolverDiffParamType        = USER_TYPE_FIRST +  8;
static const int UG_USER_TYPE_LAST              = UG_USER_TYPE_FIRST +  7;
//static const int ParaInitialStatType            = USER_TYPE_FIRST +  9;

static const int HashTableSize = 751;    ///< size of thread table : this limits the number of threads

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
class MessageQueueElement
{
   int source;                 ///< source thread rank of this message
   int count;                  ///< number of elements of the data type
   int dataTypeId;             ///< data type id
   int tag;                    ///< -1 : in case of broadcast message
   void *data;                 ///< NOTE : For basic data types, this is copy of sender side memory.
                               ///<        When the memory is copied at receive function, the memory
                               ///<        have to be freed.
                               ///<        For user defined data type, this is the receiver side memory,
                               ///<        because it is better to allocate memory in the sender side for
                               ///<        mutex locking. Sender side functions have to allocate memory.
                               ///<        In this case, memory do not hvae to freed. The memory is for
                               ///<        receiver side.
   MessageQueueElement *next;  ///< point to next message queue element

public:

   ///
   /// default constructor of MessageQueueElement
   ///
   MessageQueueElement(
         ) : source(-1), count(0), dataTypeId(-1), tag(-1), data(0), next(0)
   {
   }

   ///
   /// constructor of MessageQueueElement
   ///
   MessageQueueElement(
         int inSource,          ///< source thread rank of this message
         int inCount,           ///< number of the data type elements
         int inDataTypeId,      ///< data type id
         int inTag,             ///< tag of the message,  -1 : in case of broadcast message
         void *inData           ///< data of the message
         ) : source(inSource), count(inCount), dataTypeId(inDataTypeId), tag(inTag), data(inData), next(0)
   {
   }

   ///
   /// destructor of MessageQueueElement
   ///
   ~MessageQueueElement(
         )
   {
   }

   ///
   /// getter of source rank
   /// @return rank of the source
   ///
   int getSource(
         )
   {
      return source;
   }

   ///
   /// getter of the number of the data type elements
   /// @return the number of the data type elements
   ///
   int getCount(
         )
   {
      return count;
   }

   ///
   /// getter of the data type id
   /// @return data type id
   ///
   int getDataTypeId(
         )
   {
      return dataTypeId;
   }

   ///
   /// getter of the message tag
   /// @return tag of the message
   ///
   int getTag(
         )
   {
      return tag;
   }

   ///
   /// getter of data
   /// @return pointer to the data
   ///
   void *getData(
         )
   {
      return data;
   }

   ///
   /// getter of the pointer to the next MessageQueueElement
   /// @return pointer to MessageQueueElement
   ///
   MessageQueueElement *getNext(
         )
   {
      return next;
   }

   ///
   /// link to the next MessageQueueElement
   ///
   void link(
         MessageQueueElement *nextElement   ///< pointer to MessageQueueElement
         )
   {
      next = nextElement;
   }

};

///
/// Class of MessageQueueTableElement
///
class MessageQueueTableElement
{
   // bool                sentMessage;
   // Lock                   queueLock;
   // ConditionVariable   sentMsg;
   MessageQueueElement *head;       ///< head of the message queue
   MessageQueueElement *tail;       ///< tail of the message queue
   int                  size;      ///< number of the messages in queue

public:

   ///
   /// default constructor of MessageQueueTableElement
   ///
   MessageQueueTableElement(
           ) : head(0), tail(0), size(0)
   //      ) : sentMessage(false), head(0), tail(0)
   {
      // sentMsg.setLock(&queueLock);
   }

   ///
   /// destructor of MessageQueueTableElement
   ///
   ~MessageQueueTableElement(
         )
   {
      // LOCKED (&queueLock) {
         while( head )
         {
            MessageQueueElement *next = head->getNext();
            delete head;
            head = next;
         }
      //}
   }

   ///
   /// check if the specified message exists or nor
   /// @return pointer to MessageQueueElement, 0: no element
   ///
   MessageQueueElement *checkElement(
         int source,         ///< source rank
         int datatypeId,     ///< data type id
         int tag             ///< tag of the message
         )
   {
      MessageQueueElement *ret = 0;
      // LOCKED (&queueLock) {
         for( MessageQueueElement *current = head; current; current = current->getNext() )
         {
            if( current->getSource() == source
                  && current->getDataTypeId() == datatypeId
                  && current->getTag() == tag )
            {
               ret = current;
               break;
            }
         }
      // }
      return ret;
   }

   ///
   /// extracts a message
   /// @return pointer to the message
   ///
   MessageQueueElement *extarctElement(
         bool *sentMessage,   ///< for synchronization
         int source,          ///< source rank
         int datatypeId,      ///< data type id
         int tag              ///< tag of the message
         )
   {
      MessageQueueElement *ret = 0;

      // LOCKED (&queueLock) {
         MessageQueueElement *prev = head;
         MessageQueueElement *current = head;
         while( current )
         {
            MessageQueueElement *next = current->getNext();
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
      // }
      return ret;
   }

   ///
   /// extracts a message
   /// (This method is only for desctructor of ParaCommCPP11. No lock is necessary.)
   ///
   /// @return pointer to the message
   ///
   MessageQueueElement *extarctElement(
         bool              *sentMessage    ///< for synchronization
         )
   {
      if( head == 0 ) return 0;

      MessageQueueElement *current = head;
      MessageQueueElement *next = current->getNext();
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
   void enqueue(
         ConditionVariable *sentMsg,          ///< condition variable for synchronization
         bool              *sentMessage,      ///< mutex for synchronization
         MessageQueueElement *newElement      ///< message queue element to enter
         )
   {
      // LOCKED (&queueLock) {
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
         (*sentMsg).signal();
      // }
   }

   ///
   /// getter of head
   /// @return head element pointer
   ///
   MessageQueueElement *getHead(
         )
   {
      return head;
   }

   ///
   /// getter of size
   /// @return size of the message queue
   ///
   int getSize(
         )
   {
      return size;
   }

   ///
   /// check if the queue is empty or not
   /// @return true, if it is empty, else false
   ///
   bool isEmpty(
         )
   {
      bool empty = true;
      // LOCKED (&queueLock) {
         if( head ) empty = false;
      // }
      assert( !empty || size == 0 );
      return empty;
   }

   ///
   /// wait for a message coming to a queue
   ///
   void waitMessage(
         ConditionVariable *sentMsg,         ///< condition variable for synchronization
         bool              *sentMessage      ///< flag for synchronization
         )
   {
      CONDITIONVARIABLE_WAIT (sentMsg,  *sentMessage == true);
   }

   ///
   /// wait for a specified message coming to a queue
   ///
   void waitMessage(
         ConditionVariable *sentMsg,        ///< condition variable for synchronization
         bool              *sentMessage,    ///< flag for synchronization
         int               source,          ///< source rank of the message
         int               datatypeId,      ///< data type id of the message
         int               tag              ///< tag of the message
         )
   {
      for(;;)
      {
         // CONDITIONVARIABLE_WAIT (&sentMsg,  sentMessage == true);
         CONDITIONVARIABLE_WAIT (sentMsg,  *sentMessage == true);
         MessageQueueElement *current = head;
         while( current )
         {
            MessageQueueElement *next = current->getNext();
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
   void waitMessage(
         ConditionVariable *sentMsg,        ///< condition variable for synchronization
         bool *sentMessage,                 ///< flag for synchronization
         int source,                        ///< source rank of the message
         int *tag                           ///< tag of the message
         )
   {
      for(;;)
      {
         // CONDITIONVARIABLE_WAIT (&sentMsg,  sentMessage == true);
         CONDITIONVARIABLE_WAIT (sentMsg,  *sentMessage == true);
         MessageQueueElement *current = head;
         if( *tag == TagAny )
         {
            while( current )
            {
               MessageQueueElement *next = current->getNext();
               if( current->getSource() == source )
               {
                  *tag = current->getTag();
                  break;
               }
               current = next;
            }
         }
         else
         {
            while( current )
            {
               MessageQueueElement *next = current->getNext();
               if( current->getSource() == source
                     && current->getTag() == (*tag) )
               {
                  break;
               }
               current = next;
            }
         }
         if( current ) break;
         *sentMessage = false;
      }
   }

};

///
/// Class of ThreadsTableElement
///
class ThreadsTableElement
{

#ifdef _UG_NO_THREAD_LOCAL_STATIC
   pthread_t          tid;          ///< thread id of this thread
#endif
   int                 rank;        ///< rank of this thread
   std::ostream        *tos;        ///< tag trace stream for this thread
   ThreadsTableElement *next;       ///< next ThradsTableElement pointer

public:

   ///
   /// default constructor of ThreadsTableElement
   ///
   ThreadsTableElement(
         ) :
#ifdef _UG_NO_THREAD_LOCAL_STATIC
            tid(0),
#endif
            rank(0), tos(0), next(0)
   {
   }

   ///
   /// constructor of ThreadsTableElement
   ///
   ThreadsTableElement(
#ifdef _UG_NO_THREAD_LOCAL_STATIC
         pthread_t inTid,             ///< thread id
#endif
         int       inRank,            ///< thread rank
         ParaParamSet *paraParamSet   ///< UG parameter set
         ) :
#ifdef _UG_NO_THREAD_LOCAL_STATIC
            tid(inTid),
#endif
            rank(inRank), tos(0), next(0)
   {
      if( paraParamSet->getBoolParamValue(TagTrace) )
      {
         if( paraParamSet->isStringParamDefaultValue(TagTraceFileName) )
         {
            tos = &std::cout;
         }
         else
         {
            std::ostringstream s;
            s << paraParamSet->getStringParamValue(TagTraceFileName) << inRank;
            std::ofstream  *ofs = new std::ofstream();
            ofs->open(s.str().c_str());
            tos = ofs;
         }
      }
   }

   ///
   /// destructor of ThreadsTableElement
   ///
   ~ThreadsTableElement(
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

#ifdef _UG_NO_THREAD_LOCAL_STATIC
   ///
   /// get thread id
   /// @return thread id
   ///
   pthread_t getTid(
         )
   {
      return tid;
   }
#endif

   ///
   /// getter of this thread rank
   /// @return the thread rank
   ///
   int getRank(
         )
   {
      return rank;
   }

   ///
   /// setter of this thread rank
   ///
   void setRank(
         int r      ///< rank to be set
         )
   {
      rank = r;
   }

   ///
   /// getter of tag trace stream of this rank
   /// @return ponter to the tag trace stream
   ///
   std::ostream *getOstream(
         )
   {
      return tos;
   }

   ///
   /// get next element
   /// @return next element
   ///
   ThreadsTableElement *getNext(
         )
   {
      return next;
   }

#ifdef _UG_NO_THREAD_LOCAL_STATIC
   void link(
         ThreadsTableElement *nextElement
         )
   {
      next = nextElement;
   }
#endif

};

class ParaCalculationState;
class ParaTask;
class ParaSolverState;
class ParaSolverTerminationState;
class ParaDiffSubproblem;
class ParaInstance;
class ParaSolution;

///
/// Communicator object for pthreads thread communications
///
class ParaCommPth : public ParaComm
{
protected:
   int                      comSize;                          ///< communicator size : number of threads joined in this system
   bool                     tagTraceFlag;                     ///< indicate if tags are traced or not
   int                      **token;                          ///< index 0: token
                                                              ///< index 1: token color
                                                              ///<      -1: green
                                                              ///<     > 0: yellow ( termination origin solver number )
                                                              ///<      -2: red ( means the solver can terminate )
   ParaSysTimer             timer;                            ///< system timer
   static const char        *tagStringTable[];                ///< tag name string table
   static ThreadsTableElement  *threadsTable[HashTableSize];  ///< threads table: index is thread rank
#ifndef _UG_NO_THREAD_LOCAL_STATIC
   static __thread int localRank;                             ///< local thread rank
#endif
   MessageQueueTableElement **messageQueueTable;              ///< message queue table
   bool                     *sentMessage;                     ///< sent message flag for synchronization
   Lock                     *queueLock;                       ///< Lock for synchronization
   ConditionVariable        *sentMsg;                         ///< Condition variable for synchronization

   Lock                     *tokenAccessLock;                 ///< lock to access token
   Lock                     applicationLock;                  ///< lock for application
   Lock                     rankLock;                         ///< lock to access rank

   ///
   /// get hash code from thread id
   /// @return hash code
   ///
   unsigned int hashCode(
         pthread_t tid          ///< thread id
         );

   ///
   /// allocate memory and copy message
   ///
   void *allocateMemAndCopy(
         const void* buffer,    ///< pointer to buffer of the message
         int count,             ///< the number of data element in the message
         const int datatypeId   ///< data type of each element in the message
         );

   ///
   /// copy message
   ///
   void copy(
         void *dest,            ///< destination to copy the data
         const void *src,       ///< source of the data
         int count,             ///< the number of data element
         int datatypeId         ///< data type of each element
         );

   ///
   /// free memory
   ///
   void freeMem(
         void* buffer,          ///< pointer to buffer of the message
         int count,             ///< the number of data element
         const int datatypeId   ///< data type of each element
         );

   ///
   /// free memory
   /// @return
   ///
   bool freeStandardTypes(
         MessageQueueElement *elem   ///< pointer to a message queue element
         );

   ///
   /// check if tag string table (for debugging) set up correctly
   /// @return true if tag string table is set up correctly, false otherwise
   ///
   virtual bool tagStringTableIsSetUpCoorectly(
         );

public:

   ///
   /// constructor of ParaCommPth
   ///
   ParaCommPth(
         )
         : comSize(-1),
           tagTraceFlag(false),
           token(0),
           messageQueueTable(0),
           sentMessage(0),
           queueLock(0),
           sentMsg(0),
           tokenAccessLock(0)
   {
   }

   ///
   /// destructor of this communicator
   ///
   virtual ~ParaCommPth(
         );

   ///
   /// initializer of this communicator
   ///
   virtual void init(
         int argc,        ///< the number of arguments
         char **argv      ///< pointers to the arguments
         );

   ///
   /// get start time of this communicator
   /// (should not be used)
   /// @return start time
   ///
   double getStartTime(
         )
   {
      return timer.getStartTime();
   }

   ///
   /// get rank of caller's thread
   /// @return rank of caller's thread
   ///
   int getRank(
         );

   // int getRank( pthread_t tid);

   ///
   /// get size of this communicator, which indicates how many threads in a UG process
   /// @return the number of threads
   ///
   int getSize(
         )
   {
      return comSize;
   }

   ///
   /// get size of the messageQueueTable
   /// @return if dest >= 0 then return the number of only messages waiting to send to dest,
   ///         othrewise return the number of all messages waiting to send.
   ///
   int getNumOfMessagesWaitingToSend(
         int dest=-1
         )
   {
      if( dest >= 0 )
      {
         return messageQueueTable[dest]->getSize();
      }
      else
      {
         int nTotalMessages;
         for( int i = 0; i < ( comSize + 1 ); i++ )
         {
            nTotalMessages += messageQueueTable[i]->getSize();
         }
         return nTotalMessages;
      }
   }

   ///
   /// initializer for LoadCoordinator
   ///
   void lcInit(
         ParaParamSet *paraParamSet      ///< UG parameter set
         );

   ///
   /// initializer for Solvers
   ///
   void solverInit(
         ParaParamSet *paraParamSet      ///< UG parameter set
         )
   {
   }

   ///
   /// initializer for a specific Solver
   ///
   void solverInit(
         int rank,                       ///< rank of the Solver
         ParaParamSet *paraParamSet      ///< UG parameter set
         );

   ///
   /// reinitializer of a specific Solver
   ///
   void solverReInit(
         int rank,                       ///< rank of the Solver
         ParaParamSet *paraParamSet      ///< UG parameter set
         );

   ///
   /// delete Solver from this communicator
   ///
   void solverDel(
         int rank                        ///< rank of the Solver
         );

   ///
   /// abort. How it works sometimes depends on communicator used
   ///
   void abort(
         )
   {
      std::abort();
   }

   ///
   /// function to wait Terminated message
   /// (This function is not used currently)
   /// @return true when MPI communication is used, false when thread communication used
   ///
   bool waitTerminatedMessage(
         )
   {
      return false;
   }

   ///
   /// wait token when UG runs with deterministic mode
   /// @return true, when token is arrived to the rank
   ///
   bool waitToken(
         int rank                        ///< rank to check if token is arrived
         );

   ///
   /// pass token to from the rank to the next
   ///
   void passToken(
         int rank                        ///< from this rank, the token is passed
         );

   ///
   /// pass termination token from the rank to the next
   /// @return true, when the termination token is passed from this rank, false otherwise
   ///
   bool passTermToken(
         int rank                       ///< from this rank, the termination token is passed
         );

   ///
   /// set received token to this communicator
   ///
   void setToken(
         int rank,                      ///< rank to set the token
         int *inToken                   ///< token to be set
         );

   ///
   /// wait until thread id is registered to thread table
   ///
   void waitUntilRegistered(
         );

   ///
   /// notify that all solvers are registered
   ///
   void registedAllSolvers(
         );

   ///
   /// get ostream pointer
   /// @return pointer to ostram
   ///
   std::ostream *getOstream(
         );

   ///
   /// lock UG application to synchronize with other threads
   ///
   void lockApp(
         )
   {
      applicationLock.lock();
   }

   ///
   /// unlock UG application to synchronize with other threads
   ///
   void unlockApp(
         )
   {
      applicationLock.unlock();
   }

   ///
   /// lock rank
   ///
   void lockRank(
         )
   {
      rankLock.lock();
   }

   ///
   /// unlock rank
   ///
   void unlockRank(
         )
   {
      rankLock.unlock();
   }

   ///
   /// lock UG application to synchronize with other threads
   /// (for debugging)
   ///
   void lockApp(
         char const *f,        ///< string to indicate what is locked
         int l                 ///< a number to show something
         )
   {
      applicationLock.lock(f,l);
   }

   ///
   /// lock rank
   /// (for debugging)
   ///
   void lockRank(
         char const *f,        ///< string to indicate what is locked
         int l                 ///< a number to show something
         )
   {
      rankLock.lock(f,l);
      // std::cout << "Rank Locked: " << f << ", " << l << std::endl;
   }

   ///
   /// unlock rank
   /// (for debugging)
   ///
   void unlockRank(
         char const *f,        ///< string to indicate what is locked
         int l                 ///< a number to show something
         )
   {
      rankLock.unlock();
      // std::cout << "Rank Unlocked: " << f << ", " << l << std::endl;
   }

   ///
   /// create ParaTimer object
   /// @return pointer to ParaTimer object
   ///
   ParaTimer *createParaTimer(
         )
   {
      return new ParaTimerTh();
   }

   ///
   /// broadcast function for standard ParaData types
   /// @return always 0 (for future extensions)
   ///
   int bcast(
         void* buffer,            ///< point to the head of sending message
         int count,               ///< the number of data in the message
         const int datatypeId,    ///< data type in the message
         int root                 ///< root rank for broadcasting
         );

   ///
   /// send function for standard ParaData types
   /// @return always 0 (for future extensions)
   ///
   int send(
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
   int receive(
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
   void waitSpecTagFromSpecSource(
         const int source,       ///< source rank which the message should come from
         const int tag,          ///< tag which the message should wait
         int *receivedTag        ///< tag of the message which is arrived
         );

   ///
   /// probe function which waits a new message
   /// @return always true
   ///
   bool probe(
         int *source,            ///< source rank of the message arrived
         int *tag                ///< tag of the message arrived
         );

   ///
   /// iProbe function which checks if a new message is arrived or not
   /// @return true when a new message exists
   ///
   bool iProbe(
         int *source,            ///< source rank of the message arrived
         int *tag                ///< tag of the message arrived
         );

   ///
   /// User type send for created data type
   /// @return always 0 (for future extensions)
   ///
   int uTypeSend(
         void* bufer,             ///< point to the head of sending message
         const int datatypeId,    ///< created data type
         int dest,                ///< destination rank
         int tag                  ///< tag of the message
         );

   ///
   /// User type receive for created data type
   /// @return always 0 (for future extensions)
   ///
   int uTypeReceive(
         void** bufer,            ///< point to the head of receiving message
         const int datatypeId,    ///< created data type
         int source,              ///< source rank
         int tag                  ///< tag of the message
         );

   ///
   /// get Tag string for debugging
   /// @return string which shows Tag
   ///
   virtual const char *getTagString(
         int tag                /// tag to be converted to string
         );

};

#define DEF_PARA_COMM( para_comm, comm ) ParaCommPth *para_comm = dynamic_cast< ParaCommPth* >(comm)
#define LOCK_APP( comm ) comm->lockApp(__FILE__, __LINE__)
#define LOCK_RANK( comm ) comm->lockRank(__FILE__, __LINE__)
#define UNLOCK_RANK( comm ) comm->unlockRank(__FILE__, __LINE__)
}

#endif  // __PARA_COMM_PTH_H__
