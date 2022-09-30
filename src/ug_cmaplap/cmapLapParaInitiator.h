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

/**@file    cmapLapParaInitiator.h
 * @brief   ParaInitiator extension for CMapLap solver.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_INITIATOR_H__
#define __CMAP_LAP_PARA_INITIATOR_H__

#include <memory>
#include "ug/paraDef.h"
#include "ug/paraInitiator.h"
#include "cmapLapParaDef.h"
#include "cmapLapParaLattice.h"
#include "cmapLapParaComm.h"
#include "cmapLapParaInstance.h"
#include "cmapLapParaSolution.h"

#ifdef UG_WITH_UGS
#include "ugs/ugsParaCommMpi.h"
#endif


namespace ParaCMapLAP
{

enum FinalSolverState {
   InitialTasksGenerated,
   Aborted,
   HardTimeLimitIsReached,
   LowerBoundIsReached,
   ComputingWasInterrupted,
   ProblemWasSolved,
   RequestedSubProblemsWereSolved
};

///
/// Initiator class
///
class CMapLapParaInitiator : public UG::ParaInitiator
{
   Lattice<int, double> lattice;
   UG::ParaParamSet     *paraParams;
   std::shared_ptr<CMapLapParaInstance> instance;
   std::shared_ptr<CMapLapParaSolution> globalBestCMapLapSolution;   ///< this is for global one
   FILE                 *solutionFile;
   char                 *probFileName;
   char                 *probName;                    ///< probName for save
   char                 *basisprobname;               ///< probname for import basis
   char                 *logname;
   char                 *solutionFileName;
   int                  nSolutions;
   FinalSolverState     finalState;
   long long            nSolved;
   int                  nCMapLapThreads;
   int                  randomizeSeed;                ///< seed for randomization
#ifdef UG_WITH_UGS
   int                  seqNumber;
#endif

#ifdef UG_WITH_ZLIB
   gzstream::igzstream  checkpointTasksStream;        ///< gzstream for checkpoint tasks file
#endif

   ///
   /// print solution
   ///
   virtual void printSolution();

   ///
   /// print real
   ///
   virtual void printReal(
         FILE *file,
         double value,
         int width,
         int  precision
         );

public:
   ///
   /// constructor
   ///
   CMapLapParaInitiator(
         UG::ParaComm *inComm,
         UG::ParaTimer *inTimer
         )
         :  ParaInitiator(inComm, inTimer),
            paraParams(0),
            instance(nullptr),
            globalBestCMapLapSolution(nullptr),
            solutionFile(0),
            probFileName(0),
            probName(0),
            basisprobname(0),
            logname(0),
            solutionFileName(0),
            nSolutions(0),
            finalState(Aborted),
            nSolved(0),
            nCMapLapThreads(1),
            randomizeSeed(0)
#ifdef UG_WITH_UGS
            ,seqNumber(0)
#endif
   {
   }

   ///
   /// destructor
   ///
   virtual ~CMapLapParaInitiator(
         )
   {
      // if( probFileName ) delete probFileName; // Do not delete: argv[2]
      if( probName ) delete [] probName;

      // Close files
      fclose(solutionFile);
   }

   ///
   /// init function
   ///
   virtual int init(
         UG::ParaParamSet *paraParams,
         int          argc,
         char**       argv
         );



   virtual int reInit(
         int nRestartedRacing
         )
   {
      throw "** reInit was called **";
   }

   ///
   /// get instance
   ///
   virtual UG::ParaInstance *getParaInstance(
         )
   {
      return instance.get();
   }


   ///
   /// try to set incumbent solution
   ///
   virtual bool tryToSetIncumbentSolution(
         std::shared_ptr<CMapLapParaSolution> sol
         );

   ///
   /// send solver initialization message
   ///
   virtual void sendSolverInitializationMessage(
         );

   virtual CMapLapParaSolution *getGlobalBestIncumbentSolution(
         )
   {
      return globalBestCMapLapSolution.get();
   }

   ///
   /// get solution file name
   ///
   virtual char *getSolutionFileName(
         )
   {
      return solutionFileName;
   }

   ///
   /// get epsilon
   ///
   virtual double getEpsilon(
         )
   {
      return DEFAULT_NUM_EPSILON;
   }

   ///
   /// write solution
   ///
   virtual void writeSolution(
         const std::string& message
         );

   ///
   /// write ParaInstance
   ///
   virtual void writeParaInstance(
         const std::string& filename
         );

#ifdef UG_WITH_ZLIB
   ///
   /// write checkpoint solution
   ///
   virtual void writeCheckpointSolution(
         gzstream::ogzstream &out
         );

   ///
   /// write checkpoint solution
   ///
   virtual void writeCheckpointSolution(
         const std::string& filename
         );

   ///
   /// read solution from checkpoint file
   ///
   virtual double readSolutionFromCheckpointFile(
         char *afterCheckpointingSolutionFileName
         );
#endif

   ///
   /// read solution from checkpoint file
   ///
   virtual double readSolutionFile(
         char *solutionFileName
         );

   ///
   /// print solver version
   ///
   virtual void printSolverVersion(
         std::ostream *os  ///< output file (or NULL for standard output)
         );

   ///
   /// set final solver status
   ///
   virtual void setFinalSolverStatus(
         FinalSolverState status
         );

   ///
   /// output solution status
   ///
   virtual void outputFinalSolverStatistics(
         std::ostream *os, double time
         )
   {};
   virtual void outputFinalSolverStatistics(
         std::ostream *os,
         double time,
         double totalDeepBkzTime,
         double totalEnumTime,
         double totalSieveTime
         );

   ///
   /// Currently, not used virtual functions
   ///
   virtual void generateRacingRampUpParameterSets(
         int nParamSets, UG::ParaRacingRampUpParamSet **racingRampUpParamSets
         )
   {
      throw "** generateRacingRampUpParameterSets is called **";
   }
   virtual std::string getStatus(
         )
   {
      throw "** getStatus is called";
   }
   virtual void writeSolverParameters(
         std::ostream *os
         )
   {
      throw "** writeSolverParameters is called";
   }
   virtual void accumulateInitialStat(
         UG::ParaInitialStat *initialStat
         )
   {
      throw "** accumulateInitialStat is called";
   }


   ///
   /// getters
   //
   virtual Lattice<int, double> &getLattice()  { return lattice; }
   virtual int      getDimension()             { return lattice.n; }
   virtual double   getGH()                    { return lattice.GH; }
   virtual LatticeBasis<int> &getBasis()       { return lattice.basis; }
   virtual int      getNCMapLapThreads()       { return nCMapLapThreads; }
   virtual int      getRandomizeSeed()         { return randomizeSeed; }
};

} // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_INITIATOR_H__
