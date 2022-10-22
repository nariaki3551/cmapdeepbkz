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

#include <lattice/pbkz.hpp>
#include <DeepBKZ/lattice.h>
#include <DeepBKZ/DeepLLL.h>
#include <DeepBKZ/DeepBKZ.h>

#include "Lattice.h"
#include "cmapDeepBkz.h"


using namespace ParaCMapLAP;

///
/// @class CmapExDeepBkz
///
template<typename BasisFloat=int, typename GSFloat=double, typename EnumGSFloat=double>
class CmapExDeepBkz : public DeepBKZTool::lattice
{

using LatticePtr = std::shared_ptr<LapTools::Lattice<BasisFloat, GSFloat>>;

private:

   LatticePtr L;              ///< lattice
   LapTools::Config config;   ///< hyper parameters
   int rank;                  ///< rank
   int threadId;              ///< thread Id
   int verbose;               ///< verbose <= 0: none, 1: light, 2: medium, 3: heavy

   ParaCMapLAP::CMapLapParaSolver *cmapLapParaSolver; ///< communicator to LC
   bool mergeBasisFromLC;                             ///< merge basis from LoadCoordinator

public:

   ///
   /// Constructor
   ///
   CmapExDeepBkz(){};
   CmapExDeepBkz(
         LatticePtr inL,
         ParaCMapLAP::CMapLapParaSolver *inCmapLapParaSolver,
         int inRank=1,
         int inThreadId=0,
         int inVerbose=0,
         bool inMergeBasisFromLC=false
         )
      :
         rank(inRank),
         threadId(inThreadId),
         verbose(inVerbose),
         mergeBasisFromLC(inMergeBasisFromLC)
   {
      L = inL;
      config = inL->config;
      cmapLapParaSolver = inCmapLapParaSolver;
      SetDims(L->m, L->n);
      syncBasisToNTL();
   }


   ///
   /// @brief synchronize DeepBKZTool::lattice.basis to LapTools::Lattice.basis
   ///
   void syncBasisToEigen(
         )
   {
      for( int i = 0; i < NumRows; ++i )
         for( int j = 0; j < NumCols; ++j )
            L->basis.coeffRef(i, j) = basis[i][j];
      L->setGSO();
   }


   ///
   /// @brief synchronize LapTools::Lattice.basis to DeepBKZTool::lattice.basis
   ///
   void syncBasisToNTL(
         )
   {
      for( int i = 0; i < NumRows; ++i )
         for( int j = 0; j < NumCols; ++j )
            basis[i][j] = L->basis.coeff(i, j);
      SetGSO();
   }


   ///
   /// @brief communicate with LC
   /// @param[out] shouldAbort true if it should abort else false
   /// @return true if it communicated with LC else false
   ///
   bool communicate(
         bool &shouldAbort
         )
   {
      syncBasisToEigen();

      // double startTime = LapTools::Timer::getElapsedTime();

      if( cmapLapParaSolver->iReceiveIsNecessary() )
      {
         cmapLapParaSolver->iReceiveMessages();
      }

      // check interrupting
      if( cmapLapParaSolver->isInterrupted() )
      {
         shouldAbort = true;
         return true;
      }

      // check send-basis request
      if( cmapLapParaSolver->hasSendBasisRequest() )
      {
         ParaCMapLAP::LatticeBasis<int> basis = L->basis.template cast<int>();
         cmapLapParaSolver->sendBasis(basis, L->enumCost());
      }

      // check that it has received basis
      if( mergeBasisFromLC & cmapLapParaSolver->hasReceivedBasis() )
      {
         LapTools::LatticeBasis<BasisFloat> receivedSubBasis(cmapLapParaSolver->getReceivedBasis<BasisFloat>());
         LapTools::Lattice<BasisFloat, GSFloat> receivedSubLattice{receivedSubBasis};
         int index = 0;
         if( receivedSubLattice.isMoreReducedThan(*L, index) )
         {
            // double startMergeTime = LapTools::Timer::getElapsedTime();
            L->merge(receivedSubLattice);
            syncBasisToNTL();
            // this->mergeTime += LapTools::Timer::getElapsedTime() - startMergeTime;
            // updateBestObjectiveValue('U');
         }
      }

      // check to be necessary to send the status
      if( cmapLapParaSolver->notificationIsNecessary() )
      {
         sendStatus();
      }
      cmapLapParaSolver->iReceiveMessages();

      // this->runningTime += LapTools::Timer::getElapsedTime() - startTime;

      return true;
   }


   ///
   /// @brief communicate with LC
   /// @param[out] shouldAbort true if it should abort else false
   ///
   bool communicateInTour(
         bool &shouldAbort
         )
   {
      cmapLapParaSolver->iReceiveMessages();
      shouldAbort &= cmapLapParaSolver->isInterrupted();
      return true;
   }

   // ///
   // /// @brief update bestObjectiveValue
   // /// @param[in] sigh line-header character
   // ///
   // bool updateBestObjectiveValue(
   //       char sigh
   //       )
   // {
   //    if( LapTools::DeepBkz<BasisFloat, GSFloat, EnumGSFloat>::updateBestObjectiveValue(sigh) )
   //    {
   //       ParaCMapLAP::LatticeVector<int> u = L->basis.row(0).template cast<int>();
   //       cmapLapParaSolver->sendSolution(u, this->bestObjectiveValue);
   //       return true;
   //    }
   //    return false;
   // }

   ///
   /// @breif send SolverState of DeepBkz
   ///
   void sendStatus(
         )
   {
      ParaCMapLAP::LatticeBasis<int> _basis = L->basis.template cast<int>();
      cmapLapParaSolver->sendSolverState(
            _basis,
            -1,   // blocksize
            -1,   // number of tours
            -1    // running time
            );
   }

}; // class CmapExDeepBkz


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
      // runDeepBkz();
      runExDeepBkz();
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
   ParaCMapLAP::LatticeBasis<int> basis{*cmapLapParaTask->getBasis()};
   auto L = std::make_shared<LapTools::Lattice<int, double>>(basis);

   LapTools::Config config(paraParams->getStringParamValue(CMapLapParamFilePath));
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
/// run CMAP-LAP Randomized Extended DeepBkz
///
void
CMapDeepBkzParaSolver::runExDeepBkz(
      )
{
   CMapLapParaTask *cmapLapParaTask = dynamic_cast<CMapLapParaTask *>(currentTask);

   // create Lattice object
   ParaCMapLAP::LatticeBasis<int> basis{*cmapLapParaTask->getBasis()};
   auto L = std::make_shared<LapTools::Lattice<int, double>>(basis);

   LapTools::Config config(paraParams->getStringParamValue(CMapLapParamFilePath));
   int verbose = paraParams->getIntParamValue(DeepBkzVerbose);
   config.Quiet = ( verbose == 0 );
   L->setConfig(config);

   // randomize basis
   int begin = 0, end = L->m - 1;
   int randomizeSize = cmapLapParaTask->getU();
   if( randomizeSize >= 0 ){ begin = end - randomizeSize; }
   if( begin < 0 ){ begin = 0; }
   L->randomize(cmapLapParaTask->getSeed(), begin, end);

   // LLL reduction
   CmapDeepBkz<int, double, double> lllObj{
      L, this, getRank(), getThreadId(), 0, false};
   lllObj.deeplll();

   // create DeepBkz object
   bool mergeBasisFromLC = ( paraParams->getIntParamValue(DimensionOfSharedLattice) > 0 );
   previousNotificationTime = paraTimer->getElapsedTime() - notificationInterval;
   CmapExDeepBkz<int, double, double> reductionObj{
      L, this, getRank(), getThreadId(), verbose, mergeBasisFromLC};

   // reduction (DeepBKZ)
   int startBlocksize      = paraParams->getIntParamValue(DeepBkzStartBlockSize);
   int endBlocksize        = paraParams->getIntParamValue(DeepBkzEndBlockSize);
   int intervalBlocksize   = paraParams->getIntParamValue(DeepBkzBlockSizeInterval);
   int blocksize;
   int _start = 1, _end = L->n, gamma = L->n, abort = 4;
   double alpha = 0.99;
   for( blocksize = startBlocksize; blocksize <= endBlocksize; blocksize += intervalBlocksize )
   {
      if( blocksize > endBlocksize ){ break; }
      reductionObj.DeepBKZ(_start, _end, blocksize, alpha, gamma, abort);
   }

   // Notification message has to complete
   if( notificationProcessed ){ waitNotificationIdMessage(); }

   reductionObj.syncBasisToEigen();
   sendDeepBkzCalculationState(
         L->basis,
         blocksize,
         -1,
         -1
         // reductionObj.getNTour(),
         // reductionObj.getRunningTime()
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



///
/// instantiation
///
template class CmapExDeepBkz<int, double, double>;
template class CmapExDeepBkz<int, long double, double>;
template class CmapExDeepBkz<long int, long double, double>;
