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

/**@file    cmapDeepBkzParaLoadCoordinator.cpp
 * @brief   Load coordinator extension for CMAP-DeepBKZ.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "cmapDeepBkzParaLoadCoordinator.h"


using namespace ParaCMapLAP;


///
/// Constructor
///
CMapDeepBkzParaLoadCoordinator::CMapDeepBkzParaLoadCoordinator(
      UG::ParaComm *inComm,
      UG::ParaParamSet *inParaParamSet,
      UG::ParaInitiator *inParaInitiator,
      bool *inRacingSolversExist,
      UG::ParaTimer *inParaTimer,
      UG::ParaDeterministicTimer *detTimer,
      int inNThreadsInSolver
      )
   :
      CMapLapParaLoadCoordinator(
            inComm,
            inParaParamSet,
            inParaInitiator,
            inRacingSolversExist,
            inParaTimer,
            detTimer,
            inNThreadsInSolver
            ),
      numOfUnupdateSolverState(0)

{
   CMapLapParaParamSet *cmapLapParaParams = dynamic_cast<CMapLapParaParamSet *>(paraParams);
   int n = getCmapLapParaInitiator()->getDimension();
   int m = cmapLapParaParams->getIntParamValue(DimensionOfSharedLattice);
   dimensionOfSharedLattice = m;
   if( m > 0 )
   {
      auto instanceBasis = getCmapLapParaInitiator()->getBasis();
      LatticeBasis<int> subBasis(instanceBasis.block(0,0,m,n));
      sharedLattice.reset(new Lattice<int, double>(subBasis));
   }


   ///
   /// Setting of osLogSharedLattice
   ///
   {
      std::ostringstream s;
      s << paraParams->getStringParamValue(UG::LogSolvingStatusFilePath)
         << paraInitiator->getParaInstance()->getProbName() << "_LC" << paraComm->getRank() << ".shared_lattice.txt";
      ofsLogSharedLattice.open(s.str().c_str(), std::ios::app );
      if( !ofsLogSharedLattice )
      {
         std::cout << "Tabular vector pool log file cannot open : file name = " << s.str() << std::endl;
         exit(1);
      }
      osLogSharedLattice = &ofsLogSharedLattice;
   }

}


///
/// Destructors
///
CMapDeepBkzParaLoadCoordinator::~CMapDeepBkzParaLoadCoordinator(
      )
{
}


///
/// function to process TagSolverState message
/// @return always 0 (for extension)
///
int
CMapDeepBkzParaLoadCoordinator::processTagSolverState(
      int source,                                      ///< source solver rank
      int tag                                          ///< TagSolverState
      )
{
   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm);

   std::unique_ptr<CMapLapParaSolverState> solverState{dynamic_cast<CMapLapParaSolverState *>(cmapLapParaComm->createParaSolverState())};
   solverState->receive(cmapLapParaComm, source, tag);

   if( paraDetTimer
         && paraDetTimer->getElapsedTime() < solverState->getDeterministicTime() )
   {
      paraDetTimer->update( solverState->getDeterministicTime() - paraDetTimer->getElapsedTime() );
   }

   SolverType solverType = solverState->getSolverType();
   if( solverType != DeepBkz )
   {
      THROW_LOGICAL_ERROR2("CMapDeepBkzParaLoadCoordinator::activateSolver: Invalid solver type = ",static_cast<int>(solverType));
   }

   // for deepbkz
   int *basisArray;
   LatticeBasis<int> basis;

   // get basis
   if( solverType == DeepBkz )
   {
      std::shared_ptr<LatticeBasis<int>> prevBasis = getCmapLapParaSolverPool()->getBasis(
            source,
            solverState->getThreadId()
            );
      int n = getCmapLapParaInitiator()->getDimension();
      basisArray = solverState->getBasis();
      basis.setZero(n, n);
      int index = 0;
      for(int row = 0; row < n; ++row)
      {
         for(int col = 0; col < n; ++col)
         {
            basis(row, col) = basisArray[index++];
         }
      }
      int nNotChangeRows = 0;
      for( int i = 0; i < n; i++ )
      {
         if( basis.row(i) != prevBasis->row(i) ) break;
         nNotChangeRows++;
      }
      solverState->setNNotChangeRows(nNotChangeRows);
   }

   // update solver state
   getCmapLapParaSolverPool()->updateSolverStatus(
         source,
         solverState->getThreadId(),
         solverState->getCurrentBlockSize(),
         solverState->getBasis(),
         solverState->getEnumCost(),
         solverState->getSlopeGSA(),
         solverState->getTopHalfSlopeGSA(),
         solverState->getOrthogonalFactor()
         );
   getInstancePool()->insert(solverState->getEnumCost(), basis);

   // send-back global shortest vector
   if( getCmapLapParaInitiator()->getGlobalBestIncumbentSolution()
         && EPSLT(
            getCmapLapParaInitiator()->getGlobalBestIncumbentSolution()->getObjectiveFunctionValue(),
            solverState->getObjectiveFunctionValue(),
            DEFAULT_NUM_EPSILON)
            )
   {
      std::shared_ptr<CMapLapParaSolution> solution(
            dynamic_cast<CMapLapParaSolution *>(getCmapLapParaInitiator()->getGlobalBestIncumbentSolution()->clone(cmapLapParaComm))
            );
      solution->setThreadId(solverState->getThreadId());
      solution->send(cmapLapParaComm, source);
   }


   if( dimensionOfSharedLattice > 0 )
   {
      int n = getCmapLapParaInitiator()->getDimension();
      if( paraParams->getBoolParamValue(DynamicDimensionOfSharedLattice) &&
            numOfUnupdateSolverState >= getNSolvers() )
      {
         dimensionOfSharedLattice++;
         numOfUnupdateSolverState = 0;
      }
      int nShared = std::min(dimensionOfSharedLattice, n);
      std::unique_ptr<Lattice<int, double>> receivedLattice(new Lattice<int, double>{solverState->getBasis(), n, n});
      int index;
      if ( receivedLattice->isMoreReducedThan(*sharedLattice, index) )
      {
         // sharedLattice > receivedLattice
         *osLogSharedLattice << "Replace,Rank," << source << ",Thread," <<  solverState->getThreadId() << ",Time," << paraTimer->getElapsedTime()
            << ",unupdate," << numOfUnupdateSolverState << ",nShared," << nShared << ",index," << index << std::endl
            << receivedLattice->B.transpose()
            << std::endl;
         std::cout << "Replace,Rank," << source << ",Thread," <<  solverState->getThreadId() << ",Time," << paraTimer->getElapsedTime()
            << ",unupdate," << numOfUnupdateSolverState << ",nShared," << nShared << ",index," << index
            // << std::endl << receivedLattice->B.transpose()
            << std::endl;
         if( index == 0 )
         {
            double norm = std::sqrt(receivedLattice->basis.row(0).squaredNorm());
            std::cout
               << std::endl
               << "Norm = " << norm
               << ", AF = " << receivedLattice->approxFactor(norm)
               << ", RHF = " << receivedLattice->rootHermiteFactor(norm)
               << ", Vec = [" << receivedLattice->basis.row(0) << "]"
               << std::endl
               << std::endl;
         }

         // Eigen::VectorXd approxVector(n);
         // for( int i = 0; i < n; ++i )
         //    approxVector(i) = receivedLattice->approxFactor(-1, i);
         // *osLogSharedLattice << approxVector.transpose() << std::endl;

         sharedLattice = std::move(receivedLattice);
         numOfUnupdateSolverState = 0;
         dimensionOfSharedLattice = std::max(index+1, paraParams->getIntParamValue(DimensionOfSharedLattice));
         sharedLattice->resize(dimensionOfSharedLattice);
      }
      else if ( sharedLattice->isMoreReducedThan(*receivedLattice, index) )
      {
         // receivedLattice > sharedLattice
         // send-back sharedLattice
         std::shared_ptr<CMapLapParaBasis> paraBasis(
               cmapLapParaComm->createCMapLapParaBasis(solverState->getThreadId(), -1.0, sharedLattice->basis)
               );
         paraBasis->send(cmapLapParaComm, source);
         numOfUnupdateSolverState++;
         getLcts().nSendBackBasis++;
      }
      else if( sharedLattice->m < nShared )
      {
         /// sharedLattice == receivedLattice
         *osLogSharedLattice << "Replace,Rank," << source << ",Thread," <<  solverState->getThreadId() << ",Time," << paraTimer->getElapsedTime()
            << ",unupdate," << numOfUnupdateSolverState << ",nShared," << nShared << ",index,-1" << std::endl
            << receivedLattice->B.transpose()
            << std::endl;
         // Eigen::VectorXd approxVector(n);
         // for( int i = 0; i < n; ++i )
         //    approxVector(i) = receivedLattice->approxFactor(-1, i);
         // *osLogSharedLattice << approxVector.transpose() << std::endl;

         numOfUnupdateSolverState = 0;
         sharedLattice = std::move(receivedLattice);
         sharedLattice->resize(nShared);
      }
      else
      {
         numOfUnupdateSolverState++;
      }
   }

   if( logSolvingStatusFlag )
   {
      incrementNStatus();

      int threadId = solverState->getThreadId();
      *getOsLogSolvingStatus()
         << Logging::getSolverStateString(' ', paraTimer->getElapsedTime(), source, threadId, solverState->toStringLog())
         << std::endl;
      *getOsCsvLogSolvingStatus()
         << Logging::getSolverStateString(' ', paraTimer->getElapsedTime(), source, threadId, solverState->toStringLog(","), ",")
         << std::endl;
   }

   if( !paraParams->getBoolParamValue(NoWaitNotificationId) )
   {
      unsigned int notificationId[2];   /// notificationId[0]: notification ID
                                        /// notificationId[1]: thread ID
      notificationId[0] = solverState->getNotificaionId();
      notificationId[1] = solverState->getThreadId();
      PARA_COMM_CALL(
            cmapLapParaComm->send( notificationId, 2, UG::ParaUNSIGNED, source, UG::TagNotificationId)
            );
#ifdef _COMM_MPI_WORLD
      cmapLapParaComm->waitAllIsends();
#endif
   }
   return 0;

}
