/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*          This file is part of the program and software framework          */
/*                 CMAP-TEST --- Test configure for CMAP-LAP                 */
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

/**@file    cmapTestParaSolver.cpp
 * @brief   CMapLapParaSolver extension for CMAP-TEST.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "cmapTestParaSolver.h"
#include "ug/paraDef.h"
#include "ug/paraParamSet.h"
#include "ug/paraTask.h"
#include "ug/paraTimer.h"
#include "cmapLapParaDef.h"
#include "cmapLapParaParamSet.h"
#include "cmapLapParaTask.h"
#include "Config.h"
#include "Def.h"
#include "Lattice.h"
#include "VectorElementPool.h"
#include "cmapTestDeepBkz.h"
#include "cmapTestEnumeration.h"
#include "cmapTestGaussSieve.h"
namespace ParaCMapLAP { class CMapLapParaSolverLocalComm; }
namespace UG { class ParaComm; }
namespace UG { class ParaDeterministicTimer; }
namespace UG { class ParaInstance; }


using namespace ParaCMapLAP;


///
/// solve ParaTask
///
void
CMapTestParaSolver::solve(
      )
{
   terminationState = 0;
   double startElapsedTime = paraTimer->getElapsedTime();
   switch ( getSolverType() )
   {
      case DeepBkz:
         nParaTasksDeepBkzReceived++;
         runDeepBkz();
         runningTimeDeepBkz += paraTimer->getElapsedTime() - startElapsedTime;
         break;
      case Enum:
         nParaTasksEnumReceived++;
         runEnum();
         runningTimeEnum += paraTimer->getElapsedTime() - startElapsedTime;
         break;
      case Sieve:
         nParaTasksSieveReceived++;
         runSieve();
         runningTimeSieve += paraTimer->getElapsedTime() - startElapsedTime;
         break;
      default:
         THROW_LOGICAL_ERROR2("CMapLapParaSolver::solve: Invalid solver type = ",static_cast<int>(getSolverType()));
   }
}


///
/// run CMAP-LAP Randomized DeepBkz
///
void
CMapTestParaSolver::runDeepBkz(
      )
{
   CMapLapParaTask *cmapLapParaTask = dynamic_cast<CMapLapParaTask *>(currentTask);
   auto L = std::make_shared<LapTools::Lattice<int, double>>(*cmapLapParaTask->getBasis());

   // load from param file
   LapTools::Config config(paraParams->getStringParamValue(CMapLapParamFilePath));
   int verbose = paraParams->getIntParamValue(DeepBkzVerbose);
   config.Quiet = ( verbose == 0 );
   L->setConfig(config);

   // randomize
   int begin = 0, end = L->m - 1;
   int randomizeSize = cmapLapParaTask->getU();
   if( randomizeSize >= 0 ){ begin = end - randomizeSize; }
   if( begin < 0 ){ begin = 0; }
   L->randomize(cmapLapParaTask->getSeed(), begin, end);

   // reduction
   bool mergeBasisFromLC = ( paraParams->getIntParamValue(DimensionOfSharedLattice) > 0 );
   previousNotificationTime = paraTimer->getElapsedTime() - notificationInterval;
   CmapDeepBkz<int, double, double> reductionObj{
      L, this, getRank(), getThreadId(), verbose, mergeBasisFromLC};
   reductionObj.setNSendVectors(paraParams->getIntParamValue(DeepBkzNumOfSendVectorsToPool));
   reductionObj.setNReceivedVectors(paraParams->getIntParamValue(DeepBkzNumOfReceiveVectorsFromPool));
   reductionObj.lll();
   reductionObj.deeplll();
   int startBlocksize      = paraParams->getIntParamValue(DeepBkzStartBlockSize);
   int endBlocksize        = paraParams->getIntParamValue(DeepBkzEndBlockSize);
   int blocksizeInterval   = paraParams->getIntParamValue(DeepBkzBlockSizeInterval);
   int blocksize;
   for( blocksize = startBlocksize; blocksize <= endBlocksize; blocksize += blocksizeInterval )
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
/// run CMAP-LAP Enum
///
void
CMapTestParaSolver::runEnum(
      )
{

   CMapLapParaTask *cmapLapParaTask = dynamic_cast<CMapLapParaTask *>(currentTask);
   auto L = std::make_shared<LapTools::Lattice<int, double>>(*cmapLapParaTask->getBasis());

   // load from param file
   LapTools::Config config(paraParams->getStringParamValue(CMapLapParamFilePath));
   int verbose = paraParams->getIntParamValue(EnumVerbose);
   config.Quiet = ( verbose == 0 );
   L->setConfig(config);

   // enumeration
   LatticeVector<int> v, coeffv;
   previousNotificationTime = paraTimer->getElapsedTime() - notificationInterval;
   CmapEnumeration<int, double, double> enumObj{
      L, this, getRank(), getThreadId(), verbose};
   enumObj.init();
   enumObj.setNSendVectors(paraParams->getIntParamValue(EnumNumOfSendVectorsToPool));
   enumObj.setSamplingDepth(2);
   enumObj.projectedEnum(v, coeffv);

   // Notification message has to complete
   if( notificationProcessed ){ waitNotificationIdMessage(); }

   sendEnumCalculationState(
         L->basis,
         enumObj.getRunningTime(),
         std::sqrt(enumObj.getBestVector().squaredNorm()),  // shortestNorm
         enumObj.getNSearch()
         );
}


///
/// run CMAP-LAP Sieve
///
void
CMapTestParaSolver::runSieve(
      )
{
   CMapLapParaTask *cmapLapParaTask = dynamic_cast<CMapLapParaTask *>(currentTask);
   auto L = std::make_shared<LapTools::Lattice<int, double>>(*cmapLapParaTask->getBasis());

   // load from param file
   LapTools::Config config(paraParams->getStringParamValue(CMapLapParamFilePath));
   int verbose = paraParams->getIntParamValue(SieveVerbose);
   config.Quiet = ( verbose == 0 );
   L->setConfig(config);

   // sieving
   previousNotificationTime = paraTimer->getElapsedTime() - notificationInterval;
   CmapGaussSieve<int, double> sieveObj{
      L, this, getRank(), getThreadId(), verbose};
   sieveObj.init(
         paraParams->getIntParamValue(SieveMaxListSize),    // listsize
         paraParams->getIntParamValue(SieveMaxStackSize),   // stacksize
         paraParams->getIntParamValue(SieveMaxCollision)    // maxCollision
         );
   sieveObj.setNSendVectors(paraParams->getIntParamValue(SieveNumOfSendVectorsToPool));
   sieveObj.setNReceiveVectors(paraParams->getIntParamValue(SieveNumOfReceiveVectorsFromPool));
   LatticeVector<int> shortestV;
   sieveObj.gaussSieve(shortestV);

   // Notification message has to complete
   if( notificationProcessed ){ waitNotificationIdMessage(); }

   // If the currentTask is a localTask, then finish before send calculation state
   if( cmapLapParaTask->isLocalTask() ) return;

   sendSieveCalculationState(
         L->basis,
         sieveObj.getRunningTime(),
         L->m,                         // blocksize
         sieveObj.getNLoop(),          // nLoop
         sieveObj.getList().size(),    // listsize
         sieveObj.getStack().size(),   // stacksize
         sieveObj.getList().size(),    // max listsize
         sieveObj.getNCollision(),     // number of collision
         std::sqrt(sieveObj.getBestVector().squaredNorm())  // shortestNorm
         );
}


///
/// constructor
///
CMapTestParaSolver::CMapTestParaSolver(
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
CMapTestParaSolver::~CMapTestParaSolver(
      )
{
}
