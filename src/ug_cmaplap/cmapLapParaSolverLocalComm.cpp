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

/**@file    cmapLapParaSolverLocalComm.cpp
 * @brief   CMapLapParaSolverLocalComm extension for C++11 thread communication
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "cmapLapParaSolverLocalComm.h"
#include <cstring>
#include <cstdlib>
#ifndef _MSC_VER
#include <unistd.h>
#endif
#include "ug/paraTagDef.h"
#include "cmapLapParaTagDef.h"
#include "cmapLapParaInstanceTh.h"
#include "cmapLapParaSolutionTh.h"
#include "cmapLapParaPackedVectorTh.h"
#include "cmapLapParaBasisTh.h"
namespace ParaCMapLAP { class CMapLapParaBasis; }
namespace ParaCMapLAP { class CMapLapParaInstance; }
namespace ParaCMapLAP { class CMapLapParaPackedVector; }
namespace ParaCMapLAP { class CMapLapParaSolution; }
namespace ParaCMapLAP { class VectorElement; }


namespace ParaCMapLAP
{

std::mutex threadLockMutex;

//LocalThreadsTableElement *
//CMapLapParaSolverLocalComm::threadsTable[ThreadTableSize];

thread_local int
CMapLapParaSolverLocalComm::localThread = -1;     /*< local thread thread */


const char *
CMapLapParaSolverLocalComm::tagStringTable[] = {
   TAG_STR(UG::TagTask),
   TAG_STR(UG::TagTaskReceived),
   TAG_STR(TagDiffSubproblem),
   TAG_STR(TagRampUp),
   TAG_STR(UG::TagSolution),
   TAG_STR(TagIncumbentValue),
   TAG_STR(UG::TagSolverState),
   TAG_STR(UG::TagCompletionOfCalculation),
   TAG_STR(UG::TagNotificationId),
   TAG_STR(UG::TagTerminateRequest),
   TAG_STR(UG::TagInterruptRequest),
   TAG_STR(UG::TagTerminated),
   TAG_STR(TagRacingRampUpParamSets),
   TAG_STR(TagWinner),
   TAG_STR(UG::TagHardTimeLimit),
   TAG_STR(UG::TagAckCompletion),
   TAG_STR(TagToken),
#ifdef _COMM_CPP11
   TAG_STR(TagParaInstance),
   TAG_STR(TagCMapLapPackedVector),
   TAG_STR(TagBasisEnumCost),
   TAG_STR(TagTimeLimitRequest),
   TAG_STR(TagBasisRequest),
   TAG_STR(TagBasis),
   TAG_STR(TagVectorRequest),
   TAG_STR(TagUpdateNotificationInterval)
#endif
#ifdef _COMM_MPI_WORLD
   TAG_STR(TagCMapLapPackedVector),
   TAG_STR(TagBasisEnumCost),
   TAG_STR(TagTimeLimitRequest),
   TAG_STR(TagBasisRequest),
   TAG_STR(TagBasis),
   TAG_STR(TagVectorRequest),
   TAG_STR(TagUpdateNotificationInterval),
   TAG_STR(TagParaInstance),
   TAG_STR(TagTask1),
   TAG_STR(TagTask2),
   TAG_STR(TagTask3),
   TAG_STR(TagSolution1),
   TAG_STR(TagSolverState1),
   TAG_STR(TagCompletionOfCalculation1),
   TAG_STR(TagCMapLapPackedVector1),
   TAG_STR(TagBasis1)
#endif
};


bool
CMapLapParaSolverLocalComm::tagStringTableIsSetUpCorrectly(
      )
{
   // if( !UG::ParaComm::tagStringTableIsSetUpCorrectly() ) return false;
   // std::cout << "size = " << sizeof(tagStringTable)/sizeof(char*)
   //       << ", N_CMAP_LAP_TAGS = " << N_CMAP_LAP_TAGS << std::endl;
   return ( sizeof(tagStringTable)/sizeof(char*) == N_CMAP_LAP_TAGS );
}

const char *
CMapLapParaSolverLocalComm::getTagString(
      int tag                 /// tag to be converted to string
      )
{
   assert( tag >= 0 && tag < N_CMAP_LAP_TAGS );
   return tagStringTable[tag];
}

void
CMapLapParaSolverLocalComm::init( int argc, char *argv[] )
{
   // don't have to take any lock, because only LoadCoordinator call this function

   timer.start();

   /// Defulat comSize is the number of cores in a compute node
#ifdef _MSC_VER
   SYSTEM_INFO sysinfo;
   GetSystemInfo(&sysinfo);
   comSize = sysinfo.dwNumberOfProcessors; //includes logical cpu
#else
   comSize = sysconf(_SC_NPROCESSORS_CONF);
#endif

   for( int i = 1; i < argc; i++ )
   {
      if( strcmp(argv[i], "-ntpr") == 0 )
      {
         i++;
         if( i < argc )
            comSize = atoi(const_cast<const char*>(argv[i]));   // if -sth 0, then it is considered as use the number of cores system has
         else
         {
            std::cerr << "missing the number of solver threads after parameter '-ntpr'" << std::endl;
            exit(1);
         }
      }
   }

   ///
   /// if you add tag, you should add tagStringTale too
   ///
   // assert( sizeof(tagStringTable)/sizeof(char*) == N_TH_TAGS );
   assert( tagStringTableIsSetUpCorrectly() );

   ///
   /// initialize hashtable
   ///
   for(int i = 0; i < ThreadTableSize; i++ )
   {
      threadsTable[i] = 0;
   }

   messageQueueTable = new LocalMessageQueueTableElement *[comSize];
   sentMessage = new bool[comSize];
   queueLockMutex = new std::mutex[comSize];
   sentMsg = new std::condition_variable[comSize];
   for( int i = 0; i < comSize; i++ )
   {
      messageQueueTable[i] = new LocalMessageQueueTableElement;
      sentMessage[i] = false;
   }

}

void
CMapLapParaSolverLocalComm::lcInit(
      UG::ParaParamSet *paraParamSet
      )
{
   std::lock_guard<std::mutex> lock(threadLockMutex);
   assert( localThread == -1 );
   assert( threadsTable[0] == 0 );
   localThread = 0;
   threadsTable[0] = new LocalThreadsTableElement(0, paraParamSet);
   tagTraceFlag = paraParamSet->getBoolParamValue(UG::TagTrace);
}

void
CMapLapParaSolverLocalComm::solverInit(
      int thread,
      UG::ParaParamSet *paraParamSet
      )
{
   std::lock_guard<std::mutex> lock(threadLockMutex);
   assert( localThread == -1 );
   assert( threadsTable[thread] == 0 );
   localThread = thread;
   threadsTable[localThread] = new LocalThreadsTableElement(localThread, paraParamSet);
}

CMapLapParaSolverLocalComm::~CMapLapParaSolverLocalComm()
{
   std::lock_guard<std::mutex> lock(threadLockMutex);
   for(int i = 0; i < ThreadTableSize; i++ )
   {
      if( threadsTable[i] )
      {
         delete threadsTable[i];
      }
   }

   for(int i = 0; i < comSize; i++)
   {
      LocalMessageQueueElement *elem = messageQueueTable[i]->extarctElement(&sentMessage[i]);
      while( elem )
      {
         if( elem->getData() )
         {
            if( !freeStandardTypes(elem) )
            {
               ABORT_LOGICAL_ERROR2("Requested type is not implemented. Type = ", elem->getDataTypeId() );
            }
         }
         delete elem;
         elem = messageQueueTable[i]->extarctElement(&sentMessage[i]);
      }
      delete messageQueueTable[i];
   }
   delete [] messageQueueTable;

   if( sentMessage ) delete [] sentMessage;
   if( queueLockMutex ) delete [] queueLockMutex;
   if( sentMsg ) delete [] sentMsg;

}

int
CMapLapParaSolverLocalComm::getThreadId(
      )
{
   std::lock_guard<std::mutex> lock(threadLockMutex);
   if( localThread >= 0 ) return localThread;
   else return -1; // No ug threads
}

std::ostream *
CMapLapParaSolverLocalComm::getOstream(
      )
{
   std::lock_guard<std::mutex> lock(threadLockMutex);
   assert( threadsTable[localThread] );
   return threadsTable[localThread]->getOstream();
}

void *
CMapLapParaSolverLocalComm::allocateMemAndCopy(
      const void* buffer,
      int count,
      const int datatypeId
      )
{
   void *newBuf = 0;
   if( count == 0 ) return newBuf;

   switch(datatypeId)
   {
   case UG::ParaCHAR :
   {
      newBuf = new char[count];
      memcpy(newBuf, buffer, (unsigned long int)sizeof(char)*count);
      break;
   }
   case UG::ParaSHORT :
   {
      newBuf = new short[count];
      memcpy(newBuf, buffer, (unsigned long int)sizeof(short)*count);
      break;
   }
   case UG::ParaINT :
   {
      newBuf = new int[count];
      memcpy(newBuf, buffer, (unsigned long int)sizeof(int)*count);
      break;
   }
   case UG::ParaLONG :
   {
      newBuf = new long[count];
      memcpy(newBuf, buffer, (unsigned long int)sizeof(long)*count);
      break;
   }
   case UG::ParaUNSIGNED_CHAR :
   {
      newBuf = new unsigned char[count];
      memcpy(newBuf, buffer, (unsigned long int)sizeof(unsigned char)*count);
      break;
   }
   case UG::ParaUNSIGNED_SHORT :
   {
      newBuf = new unsigned short[count];
      memcpy(newBuf, buffer, (unsigned long int)sizeof(unsigned short)*count);
      break;
   }
   case UG::ParaUNSIGNED :
   {
      newBuf = new unsigned int[count];
      memcpy(newBuf, buffer, (unsigned long int)sizeof(unsigned int)*count);
      break;
   }
   case UG::ParaUNSIGNED_LONG :
   {
      newBuf = new unsigned long[count];
      memcpy(newBuf, buffer, (unsigned long int)sizeof(unsigned long)*count);
      break;
   }
   case UG::ParaFLOAT :
   {
      newBuf = new float[count];
      memcpy(newBuf, buffer, (unsigned long int)sizeof(float)*count);
      break;
   }
   case UG::ParaDOUBLE :
   {
      newBuf = new double[count];
      memcpy(newBuf, buffer, (unsigned long int)sizeof(double)*count);
      break;
   }
   case UG::ParaLONG_DOUBLE :
   {
      newBuf = new long double[count];
      memcpy(newBuf, buffer, (unsigned long int)sizeof(long double)*count);
      break;
   }
   case UG::ParaBYTE :
   {
      newBuf = new char[count];
      memcpy(newBuf, buffer, (unsigned long int)sizeof(char)*count);
      break;
   }
   case UG::ParaSIGNED_CHAR :
   {
      newBuf = new char[count];
      memcpy(newBuf, buffer, (unsigned long int)sizeof(char)*count);
      break;
   }
   case UG::ParaLONG_LONG :
   {
      newBuf = new long long[count];
      memcpy(newBuf, buffer, (unsigned long int)sizeof(long long)*count);
      break;
   }
   case UG::ParaUNSIGNED_LONG_LONG :
   {
      newBuf = new unsigned long long[count];
      memcpy(newBuf, buffer, (unsigned long int)sizeof(unsigned long long)*count);
      break;
   }
   case UG::ParaBOOL :
   {
      newBuf = new bool[count];
      memcpy(newBuf, buffer, (unsigned long int)sizeof(bool)*count);
      break;
   }
   default :
      THROW_LOGICAL_ERROR2("This type is not implemented. Type = ", datatypeId);
   }

   return newBuf;
}

void
CMapLapParaSolverLocalComm::copy(
      void *dest, const void *src, int count, int datatypeId
      )
{

   if( count == 0 ) return;

   switch(datatypeId)
   {
   case UG::ParaCHAR :
   {
      memcpy(dest, src, (unsigned long int)sizeof(char)*count);
      break;
   }
   case UG::ParaSHORT :
   {
      memcpy(dest, src, (unsigned long int)sizeof(short)*count);
      break;
   }
   case UG::ParaINT :
   {
      memcpy(dest, src, (unsigned long int)sizeof(int)*count);
      break;
   }
   case UG::ParaLONG :
   {
      memcpy(dest, src, (unsigned long int)sizeof(long)*count);
      break;
   }
   case UG::ParaUNSIGNED_CHAR :
   {
      memcpy(dest, src, (unsigned long int)sizeof(unsigned char)*count);
      break;
   }
   case UG::ParaUNSIGNED_SHORT :
   {
      memcpy(dest, src, (unsigned long int)sizeof(unsigned short)*count);
      break;
   }
   case UG::ParaUNSIGNED :
   {
      memcpy(dest, src, (unsigned long int)sizeof(unsigned int)*count);
      break;
   }
   case UG::ParaUNSIGNED_LONG :
   {
      memcpy(dest, src, (unsigned long int)sizeof(unsigned long)*count);
      break;
   }
   case UG::ParaFLOAT :
   {
      memcpy(dest, src, (unsigned long int)sizeof(float)*count);
      break;
   }
   case UG::ParaDOUBLE :
   {
      memcpy(dest, src, (unsigned long int)sizeof(double)*count);
      break;
   }
   case UG::ParaLONG_DOUBLE :
   {
      memcpy(dest, src, (unsigned long int)sizeof(long double)*count);
      break;
   }
   case UG::ParaBYTE :
   {
      memcpy(dest, src, (unsigned long int)sizeof(char)*count);
      break;
   }
   case UG::ParaSIGNED_CHAR :
   {
      memcpy(dest, src, (unsigned long int)sizeof(char)*count);
      break;
   }
   case UG::ParaLONG_LONG :
   {
      memcpy(dest, src, (unsigned long int)sizeof(long long)*count);
      break;
   }
   case UG::ParaUNSIGNED_LONG_LONG :
   {
      memcpy(dest, src, (unsigned long int)sizeof(unsigned long long)*count);
      break;
   }
   case UG::ParaBOOL :
   {
      memcpy(dest, src, (unsigned long int)sizeof(bool)*count);
      break;
   }
   default :
      THROW_LOGICAL_ERROR2("This type is not implemented. Type = ", datatypeId);
   }

}

void
CMapLapParaSolverLocalComm::freeMem(
      void* buffer,
      int count,
      const int datatypeId
      )
{

   if( count == 0 ) return;

   switch(datatypeId)
   {
   case UG::ParaCHAR :
   {
      delete [] static_cast<char *>(buffer);
      break;
   }
   case UG::ParaSHORT :
   {
      delete [] static_cast<short *>(buffer);
      break;
   }
   case UG::ParaINT :
   {
      delete [] static_cast<int *>(buffer);
      break;
   }
   case UG::ParaLONG :
   {
      delete [] static_cast<long *>(buffer);
      break;
   }
   case UG::ParaUNSIGNED_CHAR :
   {
      delete [] static_cast<unsigned char *>(buffer);
      break;
   }
   case UG::ParaUNSIGNED_SHORT :
   {
      delete [] static_cast<unsigned short *>(buffer);
      break;
   }
   case UG::ParaUNSIGNED :
   {
      delete [] static_cast<unsigned int *>(buffer);
      break;
   }
   case UG::ParaUNSIGNED_LONG :
   {
      delete [] static_cast<unsigned long *>(buffer);
      break;
   }
   case UG::ParaFLOAT :
   {
      delete [] static_cast<float *>(buffer);
      break;
   }
   case UG::ParaDOUBLE :
   {
      delete [] static_cast<double *>(buffer);
      break;
   }
   case UG::ParaLONG_DOUBLE :
   {
      delete [] static_cast<long double *>(buffer);
      break;
   }
   case UG::ParaBYTE :
   {
      delete [] static_cast<char *>(buffer);
      break;
   }
   case UG::ParaSIGNED_CHAR :
   {
      delete [] static_cast<char *>(buffer);
      break;
   }
   case UG::ParaLONG_LONG :
   {
      delete [] static_cast<long long *>(buffer);
      break;
   }
   case UG::ParaUNSIGNED_LONG_LONG :
   {
      delete [] static_cast<unsigned long long *>(buffer);
      break;
   }
   case UG::ParaBOOL :
   {
      delete [] static_cast<bool *>(buffer);;
      break;
   }
   default :
      THROW_LOGICAL_ERROR2("This type is not implemented. Type = ", datatypeId);
   }

}

bool
CMapLapParaSolverLocalComm::freeStandardTypes(
      LocalMessageQueueElement *elem   ///< pointer to a message queue element
      )
{
   if( elem->getDataTypeId() < UG_USER_TYPE_FIRST )
   {
      freeMem(elem->getData(), elem->getCount(), elem->getDataTypeId() );
   }
   else
   {
      switch( elem->getDataTypeId())
      {
      case ParaInstanceType:
      {
         delete reinterpret_cast<UG::ParaInstance *>(elem->getData());
         break;
      }
      case ParaSolutionType:
      {
         delete reinterpret_cast<UG::ParaSolution *>(elem->getData());
         break;
      }
      case ParaParamSetType:
      {
         delete reinterpret_cast<UG::ParaParamSet *>(elem->getData());
         break;
      }
      case ParaTaskType:
      {
         delete reinterpret_cast<UG::ParaTask *>(elem->getData());
         break;
      }
      case ParaSolverStateType:
      {
         delete reinterpret_cast<UG::ParaSolverState *>(elem->getData());
         break;
      }
      case ParaCalculationStateType:
      {
         delete reinterpret_cast<UG::ParaCalculationState *>(elem->getData());
         break;
      }
      case ParaSolverTerminationStateType:
      {
         delete reinterpret_cast<UG::ParaSolverTerminationState *>(elem->getData());
         break;
      }
      case ParaRacingRampUpParamType:
      {
         delete reinterpret_cast<UG::ParaRacingRampUpParamSet *>(elem->getData());
         break;
      }
      default:
      {
         return false;
      }
      }
   }
   return true;
}

int
CMapLapParaSolverLocalComm::bcast(
      void* buffer,
      int count,
      const int datatypeId,
      int root
      )
{
   if( getThreadId() == root )
   {
      for(int i=0; i < comSize; i++)
      {
         if( i != root )
         {
            send(buffer, count, datatypeId, i, -1);
         }
      }
   }
   else
   {
      receive(buffer, count, datatypeId, root, -1);
   }
   return 0;
}

int
CMapLapParaSolverLocalComm::send(
      void* buffer,
      int count,
      const int datatypeId,
      int dest,
      const int tag
      )
{
   {
      std::lock_guard<std::mutex> lock(queueLockMutex[dest]);
      messageQueueTable[dest]->enqueue(sentMsg[dest], queueLockMutex[dest], &sentMessage[dest],
            new LocalMessageQueueElement(getThreadId(), count, datatypeId, tag,
                  allocateMemAndCopy(buffer, count, datatypeId) ) );
   }
   TAG_THREAD_TRACE (Send, To, dest, tag);
   return 0;
}

int
CMapLapParaSolverLocalComm::receive(
      void* buffer,
      int count,
      const int datatypeId,
      int source,
      const int tag
      )
{
   int qThread = getThreadId();
   LocalMessageQueueElement *elem = 0;
   if( !messageQueueTable[qThread]->checkElement(source, datatypeId, tag) )
   {
      messageQueueTable[qThread]->waitMessage(sentMsg[qThread], queueLockMutex[qThread], &sentMessage[qThread], source, datatypeId, tag);
   }
   {
      std::lock_guard<std::mutex> lock(queueLockMutex[qThread]);
      elem = messageQueueTable[qThread]->extarctElement(&sentMessage[qThread],source, datatypeId, tag);
   }
   assert(elem);
   copy( buffer, elem->getData(), count, datatypeId );
   freeMem(elem->getData(), count, datatypeId );
   delete elem;
   TAG_THREAD_TRACE (Recv, From, source, tag);
   return 0;
}

void
CMapLapParaSolverLocalComm::waitSpecTagFromSpecSource(
      const int source,
      const int tag,
      int *receivedTag
      )
{
   /*
   // Just wait, iProbe and receive will be performed after this call
   messageQueueTable[getThread()]->waitMessage(source, datatypeId, tag);
   TAG_TRACE (Probe, From, source, tag);
   return 0;
   */
   int qThread = getThreadId();

   // LOCKED ( &queueLock[getThread()] )
   // {
   //   messageQueueTable[qThread]->waitMessage(sentMsg[qThread], queueLockMutex[qThread], &sentMessage[qThread], source, receivedTag);
   // }

   if( tag == UG::TagAny )
   {
      messageQueueTable[qThread]->waitMessage(sentMsg[qThread], queueLockMutex[qThread], &sentMessage[qThread], source, receivedTag);
   }
   else
   {
      messageQueueTable[qThread]->waitMessage(sentMsg[qThread], queueLockMutex[qThread], &sentMessage[qThread], source, tag);
      (*receivedTag) = tag;
   }
   TAG_THREAD_TRACE (Probe, From, source, *receivedTag);
   return;
}

bool
CMapLapParaSolverLocalComm::probe(
      int* source,
      int* tag
      )
{
   int qThread = getThreadId();
   messageQueueTable[qThread]->waitMessage(sentMsg[qThread], queueLockMutex[qThread],  &sentMessage[qThread]);
   LocalMessageQueueElement *elem = messageQueueTable[qThread]->getHead();
   *source = elem->getSource();
   *tag = elem->getTag();
   TAG_THREAD_TRACE (Probe, From, *source, *tag);
   return true;
}

bool
CMapLapParaSolverLocalComm::iProbe(
      int* source,
      int* tag
      )
{
   bool flag = false;
   int qThread = getThreadId();
   {
      std::lock_guard<std::mutex> lock(queueLockMutex[qThread]);
      flag = !(messageQueueTable[qThread]->empty());
      if( flag )
      {
         LocalMessageQueueElement *elem = messageQueueTable[qThread]->getHead();
         *source = elem->getSource();
         *tag = elem->getTag();
         TAG_THREAD_TRACE (Iprobe, From, *source, *tag);
      }
   }
   return flag;
}

int
CMapLapParaSolverLocalComm::uTypeSend(
      void* buffer,
      const int datatypeId,
      int dest,
      const int tag
      )
{
   {
      std::lock_guard<std::mutex> lock(queueLockMutex[dest]);
      messageQueueTable[dest]->enqueue(sentMsg[dest], queueLockMutex[dest], &sentMessage[dest],
            new LocalMessageQueueElement(getThreadId(), 1, datatypeId, tag, buffer ) );
   }
   TAG_THREAD_TRACE (Send, To, dest, tag);
   return 0;
}

int
CMapLapParaSolverLocalComm::uTypeReceive(
      void** buffer,
      const int datatypeId,
      int source,
      const int tag
      )
{
   int qThread = getThreadId();
   if( !messageQueueTable[qThread]->checkElement(source, datatypeId, tag) )
   {
      messageQueueTable[qThread]->waitMessage(sentMsg[qThread], queueLockMutex[qThread], &sentMessage[qThread], source, datatypeId, tag);
   }
   LocalMessageQueueElement *elem = 0;
   {
      std::lock_guard<std::mutex> lock(queueLockMutex[qThread]);
      elem = messageQueueTable[qThread]->extarctElement(&sentMessage[qThread], source, datatypeId, tag);
   }
   assert(elem);
   *buffer = elem->getData();
   delete elem;
   TAG_THREAD_TRACE (Recv, From, source, tag);
   return 0;
}

/*******************************************************************************
* transfer object factory
*******************************************************************************/

UG::ParaInstance*
CMapLapParaSolverLocalComm::createParaInstance(
      )
{
   return new CMapLapParaInstanceTh();
}

UG::ParaSolution*
CMapLapParaSolverLocalComm::createParaSolution(
      )
{
   return new CMapLapParaSolutionTh();
}

CMapLapParaInstance*
CMapLapParaSolverLocalComm::createCMapLapParaInstance(
      char      *inProbFileName,
      char      *inProbName
      )
{
    return new CMapLapParaInstanceTh(inProbFileName, inProbName);
}

CMapLapParaInstance*
CMapLapParaSolverLocalComm::createCMapLapParaInstance(
      char      *inProbName
      )
{
    return new CMapLapParaInstanceTh(inProbName);
}

CMapLapParaSolution*
CMapLapParaSolverLocalComm::createCMapLapParaSolution(
      int inThreadId,               ///< thread Id
      LatticeVector<int> inV,       ///< lattice vector
      double inObjValue             ///< objective function value
      )
{
   return new CMapLapParaSolutionTh(inThreadId, inV, inObjValue);
}


///
/// @brief create ParaBasis object
///
CMapLapParaBasis *
CMapLapParaSolverLocalComm::createCMapLapParaBasis(
      int inThreadId,            ///< thread Id
      double inEnumCost,         ///< enumeration cost
      LatticeBasis<int> inBasis  ///< lattice basis
      )
{
   return new CMapLapParaBasisTh(inThreadId, inEnumCost, inBasis);
}

///
/// @brief create ParaBasis object
///
CMapLapParaBasis *
CMapLapParaSolverLocalComm::createCMapLapParaBasis(
      )
{
   return new CMapLapParaBasisTh();
}


CMapLapParaPackedVector*
CMapLapParaSolverLocalComm::createCMapLapParaPackedVector(
      int      inThread,
      LatticeBasis<int> &inVectors,
      int      inRank,
      int      nThreadsPerRank
      )
{
   return new CMapLapParaPackedVectorTh(inThread, inVectors, inRank, nThreadsPerRank);
}

CMapLapParaPackedVector *
CMapLapParaSolverLocalComm::createCMapLapParaPackedVector(
      int inThreadId,
      std::vector<std::shared_ptr<VectorElement>> &inVectorElements
      )
{
   return new CMapLapParaPackedVectorTh(inThreadId, inVectorElements);
}

///
/// @brief create ParaPackedVector object
///
CMapLapParaPackedVector *
CMapLapParaSolverLocalComm::createCMapLapParaPackedVector(
      int inThreadId,
      int inNVectors
      )
{
   return new CMapLapParaPackedVectorTh(inThreadId, inNVectors);
}

CMapLapParaPackedVector *
CMapLapParaSolverLocalComm::createCMapLapParaPackedVector(
      )
{
  return new CMapLapParaPackedVectorTh();
}


} // namespace ParaCMapLAP
