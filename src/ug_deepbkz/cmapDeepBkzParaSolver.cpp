/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*          This file is part of the program and software framework          */
/*  CMAP-DeepBKZ --- Configurable Massively Parallel Solver for DeepBKZ      */
/*                                                                           */
/*  Copyright Written by Nariaki Tateiwa <n-tateiwa@kyudai.jp>,              */
/*                       Yuji Shinano <shinano@zib.de>,                      */
/*            Copyright (C) 2021 by Zuse Institute Berlin,                   */
/*            licensed under LGPL version 3 or later.                        */
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

/**@file    cmapDeepBkzParaSolver.cpp
 * @brief   ParaSolver extension for CMAP-DeepBKZ
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "cmapDeepBkzParaSolver.h"
#include "cmapDeepBkz.h"


using namespace ParaCMapLAP;


///
/// solve ParaTask
///
void
CMapDeepBkzParaSolver::solve(
      )
{
   terminationState = 0;
   double startElapsedTime = paraTimer->getElapsedTime();
   if( getSolverType() == DeepBkz )
   {
      nParaTasksDeepBkzReceived++;
      runDeepBkz();
      runningTimeDeepBkz += paraTimer->getElapsedTime() - startElapsedTime;
   }
   else
   {
      THROW_LOGICAL_ERROR2("CMapLapParaSolver::solve: Invalid solver type = ",static_cast<int>(getSolverType()));
   }
}


///
/// run CMAP-LAP Randomized DeepBkz
///
void
CMapDeepBkzParaSolver::runDeepBkz(
      )
{
   CMapLapParaTask *cmapLapParaTask = dynamic_cast<CMapLapParaTask *>(currentTask);

   // create Lattice object
   LatticeBasis<int> basis{*cmapLapParaTask->getBasis()};
   auto L = std::make_shared<LapTools::Lattice<int, double>>(basis);

   Config config(paraParams->getStringParamValue(CMapLapParamFilePath));
   int verbose = paraParams->getIntParamValue(DeepBkzVerbose);
   config.Quiet = ( verbose == 0 );
   L->setConfig(config);

   // randomize basis
   int begin = 0, end = L->m - 1;
   int randomizeSize = cmapLapParaTask->getU();
   if( randomizeSize >= 0 ){ begin = end - randomizeSize; }
   if( begin < 0 ){ begin = 0; }
   L->randomize(cmapLapParaTask->getSeed(), begin, end);

   // create DeepBkz object
   bool mergeBasisFromLC = ( paraParams->getIntParamValue(DimensionOfSharedLattice) > 0 );
   previousNotificationTime = paraTimer->getElapsedTime() - notificationInterval;
   CmapDeepBkz<int, double, double> reductionObj{
      L, this, getRank(), getThreadId(), verbose, mergeBasisFromLC};
   reductionObj.setNSendVectors(paraParams->getIntParamValue(DeepBkzNumOfSendVectorsToPool));

   // reduction (LLL, DeepLLL)
   reductionObj.lll();
   reductionObj.deeplll();

   // reduction (DeepBKZ)
   int startBlocksize      = paraParams->getIntParamValue(DeepBkzStartBlockSize);
   int endBlocksize        = paraParams->getIntParamValue(DeepBkzEndBlockSize);
   int intervalBlocksize   = paraParams->getIntParamValue(DeepBkzBlockSizeInterval);
   int blocksize;
   for( blocksize = startBlocksize; blocksize <= endBlocksize; blocksize += intervalBlocksize )
   {
      if( blocksize > endBlocksize ){ break; }
      reductionObj.deepbkz(blocksize);
   }


   // Notification message has to complete
   if( notificationProcessed ){ waitNotificationIdMessage(); }

   sendDeepBkzCalculationState(
         L->basis,
         blocksize,
         reductionObj.getNTour(),
         reductionObj.getRunningTime()
         );
}


///
/// constructor
///
CMapDeepBkzParaSolver::CMapDeepBkzParaSolver(
      int argc,
      char **argv,
      UG::ParaComm     *comm,
      CMapLapParaSolverLocalComm *inLocalComm,
      UG::ParaParamSet *inParaParamSet,
      UG::ParaInstance *inParaInstance,
      UG::ParaDeterministicTimer *inDetTimer,
      double timeOffset
      )
      :
      CMapLapParaSolver(argc, argv, comm, inLocalComm, inParaParamSet, inParaInstance, inDetTimer, timeOffset)
{
}


///
/// deconstructor
///
CMapDeepBkzParaSolver::~CMapDeepBkzParaSolver(
      )
{
}
