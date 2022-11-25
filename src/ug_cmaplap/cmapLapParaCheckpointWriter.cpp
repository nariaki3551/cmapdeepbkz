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

/**@file    cmapLapParaCheckpointerWriter.cpp
 * @brief   Base class for CheckpointWriter.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "cmapLapParaCheckpointWriter.h"
#include "ug/paraComm.h"
#include "ug/paraDef.h"
#include "ug/paraParamSet.h"
#include "ug/paraTagDef.h"
#include "cmapLapParaSolution.h"
#include "cmapLapParaInstancePool.h"
#include "cmapLapParaLoadCoordinator.h"
#include "cmapLapParaShareDataPool.h"
#include "cmapLapParaSolverLocalComm.h"
#include "cmapLapParaTaskPool.h"


namespace ParaCMapLAP
{


///
/// constructor
///
CMapLapParaCheckpointWriter::CMapLapParaCheckpointWriter(
      CMapLapParaSolverLocalComm *inLocalComm,
      UG::ParaParamSet *inParaParamSet
      ) : localComm(inLocalComm),
          paraParams(inParaParamSet)
{
   /// creation of LC local solver
   paraTimer = std::unique_ptr<UG::ParaTimer>{localComm->createParaTimer()};
#ifdef _COMM_CPP11
   // In ParaTimerMpi, created time is start time
   paraTimer->init(localComm);
#endif
}


///
/// run this Checkpoint Writer
///
void
CMapLapParaCheckpointWriter::run(
      )
{
   std::cout
      << "*** "
      << "Time "
      << paraTimer->getElapsedTime()
      << " CMapLapParaCheckpointWriter"
      << " start"
      << std::endl;

   for(;;)
   {
      //  waitForAnyMessageFromLocal
      int source = -1;
      int tag = -1;
      localComm->probe(&source, &tag);
      if( tag == UG::TagTerminateRequest )
      {
         localComm->receive( NULL, 0, UG::ParaBYTE, 0, tag );
         std::cout
            << "*** "
            << "Time "
            << paraTimer->getElapsedTime()
            << " CMapLapParaCheckpointWriter::run"
            << " receive UG::TagTerminateRequest from LC "
            << std::endl;
         break;
      }
      else if( tag == UG::TagTask )
      {
         CheckpointElement *checkpointElement;
         localComm->uTypeReceive(
               (void **)&checkpointElement, ParaCheckpointElementType, 0, UG::TagTask);
         std::cout
            << "*** "
            << "Time "
            << paraTimer->getElapsedTime()
            << " CMapLapParaCheckpointWriter::run"
            << " receive currentCheckpointElement "
            << std::endl;
         updateCheckpointFiles( checkpointElement );
         localComm->send( NULL, 0, UG::ParaBYTE, 0, UG::TagCompletionOfCalculation );
      }
      else
      {
         THROW_LOGICAL_ERROR2( tag, " is invalid tag for checkpoint writer" );
      }
   }
   localComm->send( NULL, 0, UG::ParaBYTE, 0, UG::TagTerminated );  // destination is LC ( 0 )
   std::cout
      << "*** "
      << "Time "
      << paraTimer->getElapsedTime()
      << " CMapLapParaCheckpointWriter::run"
      << " send UG::TagTerminated to LC "
      << std::endl;
   return;
}


///
/// write checkpoint files
///
void
CMapLapParaCheckpointWriter::updateCheckpointFiles(
      CheckpointElement *checkpointElement
      )
{
   double startTime= paraTimer->getElapsedTime();

   time_t timer;
   char timeStr[30];


   // get checkpoint
   time(&timer);
   // make checkpoint time string
#ifdef _MSC_VER
   int bufsize = 256;
   ctime_s(timeStr, bufsize, &timer);
#else
   ctime_r(&timer, timeStr);
#endif
   for( int i = 0; timeStr[i] != '\0' && i < 26; i++ )
   {
      if( timeStr[i] == ' ') timeStr[i] = '_';
      if( timeStr[i] == '\n' ) timeStr[i] = '\0';
   }
   char *newCheckpointTimeStr = &timeStr[4];    // remove a day of the week

   if( strcmp(newCheckpointTimeStr, checkpointElement->lastCheckpointTimeStr) == 0 )
   {
      int l = strlen(newCheckpointTimeStr);
      newCheckpointTimeStr[l] = 'a';
      newCheckpointTimeStr[l+1] = '\0';
   }

   ///
   /// save tasks information
   ///
   char taskFileName[256];
   sprintf(taskFileName,"%s%s_%s_tasks.gz",
         paraParams->getStringParamValue(UG::CheckpointFilePath),
         checkpointElement->probName,newCheckpointTimeStr);
   gzstream::ogzstream checkpointTaskStream;
   checkpointTaskStream.open(taskFileName, std::ios::out | std::ios::binary);
   if( !checkpointTaskStream )
   {
      std::cout << "Checkpoint file for tasks cannot open. file name = " << taskFileName << std::endl;
      exit(1);
   }
   bool active = true;
   checkpointElement->paraActiveTaskQueue-> writeTasksToCheckPointFile(checkpointTaskStream,active);
   checkpointElement->paraDeepBkzTaskQueue->writeTasksToCheckPointFile(checkpointTaskStream,!active);
   checkpointElement->paraEnumTaskQueue->   writeTasksToCheckPointFile(checkpointTaskStream,!active);
   checkpointElement->paraSieveTaskQueue->  writeTasksToCheckPointFile(checkpointTaskStream,!active);
   checkpointTaskStream.close();


   ///
   /// save incumbent solution
   ///
   char solutionFileName[256];
   sprintf(solutionFileName,"%s%s_%s_solution.gz",
         paraParams->getStringParamValue(UG::CheckpointFilePath),
         checkpointElement->probName,newCheckpointTimeStr);
   gzstream::ogzstream checkpointSolutionStream;
   checkpointSolutionStream.open(solutionFileName, std::ios::out | std::ios::binary);
   if( !checkpointSolutionStream )
   {
      std::cout << "Checkpoint file for solution cannot open. file name = " << solutionFileName << std::endl;
      exit(1);
   }
   if( checkpointElement->globalBestCMapLapSolution )
   {
      checkpointElement->globalBestCMapLapSolution->writeCheckpointSolution(checkpointSolutionStream);
   }
   checkpointSolutionStream.close();     // empty solution file is necessary,
                                         // because it is removed next at the next checkpoint


   ///
   /// save InstancePool
   ///
   char instancePoolFileName[256];
   sprintf(instancePoolFileName,"%s%s_%s_instancePool.gz",
         paraParams->getStringParamValue(UG::CheckpointFilePath),
         checkpointElement->probName,newCheckpointTimeStr);
   gzstream::ogzstream checkpointInstancePoolStream;
   checkpointInstancePoolStream.open(instancePoolFileName, std::ios::out | std::ios::binary);
   if( !checkpointInstancePoolStream )
   {
      std::cout << "Checkpoint file for instancePool cannot open. file name = " << instancePoolFileName << std::endl;
      exit(1);
   }
   checkpointElement->basisElementQueue->writeBasisElementsToCheckpointFile(checkpointInstancePoolStream, -1);
   checkpointInstancePoolStream.close();


   ///
   /// save Incumbent Basis
   ///
   char basisTextFileName[256];
   sprintf(basisTextFileName,"%s%s_%s_basis.txt",
         paraParams->getStringParamValue(UG::CheckpointFilePath),
         checkpointElement->probName,newCheckpointTimeStr);
   checkpointElement->incumbentBasis->writeBasis(basisTextFileName);


   ///
   /// save ShareDataPool
   ///
   char shareDataPoolFileName[256];
   sprintf(shareDataPoolFileName,"%s%s_%s_shareDataPool.gz",
         paraParams->getStringParamValue(UG::CheckpointFilePath),
         checkpointElement->probName,newCheckpointTimeStr);
   gzstream::ogzstream checkpointShareDataPoolStream;
   checkpointShareDataPoolStream.open(shareDataPoolFileName, std::ios::out | std::ios::binary);
   if( !checkpointShareDataPoolStream )
   {
      std::cout << "Checkpoint file for ShareDataPool cannot open. file name = " << shareDataPoolFileName << std::endl;
      exit(1);
   }
   checkpointElement->vectorElementQueue->writeVectorElementsToCheckpointFile(
         checkpointShareDataPoolStream,
         -1    // nWrite number of written share data
         );
   checkpointShareDataPoolStream.close();


   ///
   /// save LoadCoordinator statistics
   ///
   char loadCoordinatorStatisticsFileName[256];
   sprintf(loadCoordinatorStatisticsFileName,"%s%s_%s_loadCoordinatorStatistics.txt",
         paraParams->getStringParamValue(UG::CheckpointFilePath),
         checkpointElement->probName,newCheckpointTimeStr);
   std::ofstream lcStaticFile(loadCoordinatorStatisticsFileName);
   lcStaticFile << checkpointElement->lcStatString;
   lcStaticFile.close();


   ///
   /// remove old checkfiles after second checkpointing
   ///
   if( checkpointElement->lastCheckpointTimeStr[0] != ' ' )  //  first checkpoint, should not remove
   {
      sprintf(taskFileName,"%s%s_%s_tasks.gz",
            paraParams->getStringParamValue(UG::CheckpointFilePath),
            checkpointElement->probName,checkpointElement->lastCheckpointTimeStr);
      sprintf(solutionFileName,"%s%s_%s_solution.gz",
            paraParams->getStringParamValue(UG::CheckpointFilePath),
            checkpointElement->probName,checkpointElement->lastCheckpointTimeStr);
      sprintf(basisTextFileName,"%s%s_%s_basis.txt",
            paraParams->getStringParamValue(UG::CheckpointFilePath),
            checkpointElement->probName,checkpointElement->lastCheckpointTimeStr);
      sprintf(instancePoolFileName,"%s%s_%s_instancePool.gz",
            paraParams->getStringParamValue(UG::CheckpointFilePath),
            checkpointElement->probName,checkpointElement->lastCheckpointTimeStr);
      sprintf(shareDataPoolFileName,"%s%s_%s_shareDataPool.gz",
            paraParams->getStringParamValue(UG::CheckpointFilePath),
            checkpointElement->probName,checkpointElement->lastCheckpointTimeStr);
      sprintf(loadCoordinatorStatisticsFileName,"%s%s_%s_loadCoordinatorStatistics.txt",
            paraParams->getStringParamValue(UG::CheckpointFilePath),
            checkpointElement->probName,checkpointElement->lastCheckpointTimeStr);
      if( remove(taskFileName) )
      {
         std::cout << "checkpoint task file cannot be removed: errno = " << strerror(errno) << std::endl;
         exit(1);
      }
      if ( remove(solutionFileName) )
      {
         std::cout << "checkpoint solution file cannot be removed: errno = " << strerror(errno) << std::endl;
         exit(1);
      }
      if ( remove(basisTextFileName) )
      {
         std::cout << "checkpoint Basis text file cannot be removed: errno = " << strerror(errno) << std::endl;
         exit(1);
      }
      if ( remove(instancePoolFileName) )
      {
         std::cout << "checkpoint InstancePool file cannot be removed: errno = " << strerror(errno) << std::endl;
         exit(1);
      }
      if ( remove(shareDataPoolFileName) )
      {
         std::cout << "checkpoint ShareDataPool file cannot be removed: errno = " << strerror(errno) << std::endl;
         exit(1);
      }
      if ( remove(loadCoordinatorStatisticsFileName) )
      {
         std::cout << "checkpoint LoadCoordinatorStatistics file cannot be removed: errno = " << strerror(errno) << std::endl;
         exit(1);
      }
   }


   ///
   /// remove after_checkpoint_solution if it exists
   ///
   char afterCheckpointingSolutionFileName[256];
   sprintf(afterCheckpointingSolutionFileName,"%s%s_after_checkpointing_solution.gz",
         paraParams->getStringParamValue(UG::CheckpointFilePath),
         checkpointElement->probName );
   gzstream::igzstream afterCheckpointingSolutionStream;
   afterCheckpointingSolutionStream.open(afterCheckpointingSolutionFileName, std::ios::in | std::ios::binary);
   if( afterCheckpointingSolutionStream  )
   {
      // afater checkpointing solution file exists
      afterCheckpointingSolutionStream.close();
      if ( remove(afterCheckpointingSolutionFileName) )
      {
         std::cout << "after checkpointing solution file cannot be removed: errno = " << strerror(errno) << std::endl;
         exit(1);
      }
   }

   // update last checkpoint time string
   checkpointElement->setLastCheckpointTimeStr(newCheckpointTimeStr);

   std::cout
      << "*** "
      << "Time "
      << paraTimer->getElapsedTime()
      << " CMapLapParaCheckpointWriter::updateCheckpointFiles"
      << " finish"
      << " lastCheckpointTimeStr <- "
      << checkpointElement->lastCheckpointTimeStr
      << std::endl;

   checkpointElement->writeCheckpointTime = paraTimer->getElapsedTime() - startTime;
}

} // namespace ParaCMapLAP
