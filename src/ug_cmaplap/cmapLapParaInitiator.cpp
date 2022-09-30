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

/**@file    cmapLapParaInitiator.cpp
 * @brief   ParaInitiator extension for CMAP-LAP solver.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cmath>
#include <cassert>
#include <cfloat>
#include "cmapLapParaDef.h"
#include "cmapLapParaInstance.h"
#include "cmapLapParaInitiator.h"
#include "cmapLapParaLattice.h"
#include "cmapLapParaSolution.h"


namespace ParaCMapLAP
{

///
/// init function
///
int
CMapLapParaInitiator::init(
      UG::ParaParamSet *inParaParams,
      int     argc,
      char**  argv
      )
{
   int i;
#ifdef _COMM_PTH
   bool noUpgrade = false;
#endif
   paraParams = inParaParams;

   basisprobname = argv[2];

   ///
   /// Parse parameters
   ///
   for( i = 3; i < argc; ++i )   /// the first argument is runtime parameter file for ParaCMAP_LAP
                                 /// the second argument is problem file name
   {
      if( strcmp(argv[i], "-l") == 0 )
      {
         i++;
         if( i < argc )
            logname = argv[i];
         else
         {
            std::cerr << "missing log filename after parameter '-l'" << std::endl;
            exit(1);
         }
      }
      else if( strcmp(argv[i], "-w") == 0)
      {
#ifdef UG_WITH_ZLIB
         i++;
         if( i < argc )
         {
            prefixWarm = argv[i];
            char tasksFileName[256];
            sprintf(tasksFileName,"%s_tasks.gz",prefixWarm);
            checkpointTasksStream.open(tasksFileName, std::ios::in | std::ios::binary);
            if( !checkpointTasksStream.good() ){
                std::cerr << "ERROR: Opening file `" << tasksFileName << "' failed.\n";
                exit(1);
            }
         }
         else
         {
            std::cerr << "missing settings filename after parameter '-w'" << std::endl;
            exit(1);
         }
#else
         std::cerr << "Cannot work with parameter '-w' compiling without zlib" << std::endl;
         exit(1);
#endif
      }
      else if ( strcmp(argv[i], "-ntpr") == 0 )
      {
         i++;  // just omit this parameter and the following number.
         if( i < argc )
         {
            nCMapLapThreads = atoi(argv[i]);
         }
         else
         {
            std::cerr << "missing number after parameter '-ntpr'" << std::endl;
            exit(1);
         }
      }
      else if ( strcmp(argv[i], "-isol" ) == 0 )
      {
         i++;
         if( i < argc )
         {
            const char *extension = strrchr(argv[i], '.');
            if( extension && *(extension + 1) != '\0' )
            {
               extension += 1;
               if( strcmp(extension, "gz")  == 0 )
               {
                  readSolutionFromCheckpointFile(argv[i]);
                  continue;
               }
            }
            readSolutionFile(argv[i]);
         }
         else
         {
            std::cerr << "missing solution filename after parameter '-isol'" << std::endl;
            exit(1);
         }
      }
      else if ( strcmp(argv[i], "-fsol" ) == 0 )
      {
         i++;
         if( i < argc )
         {
            solutionFileName = argv[i];
         }
         else
         {
            std::cerr << "missing solution filename after parameter '-fsol'" << std::endl;
            exit(1);
         }
      }
      else if ( strcmp(argv[i], "-sth") == 0 )
      {
         i++;
      }
      // else if ( strcmp(argv[i], "-mfp") == 0 )
      // {
      //    i++;  // only solver use this
      // }
      else if ( strcmp(argv[i], "-o") == 0 )
      {
         i++;
         if( i < argc )
         {
            int nOverwriteFilenames = atoi(argv[i]);
            i += nOverwriteFilenames;
            if( i > argc )
            {
               std::cerr << "too many arguments after parameter '-o'" << std::endl;
               exit(1);
            }
         }
      }
      else if ( strcmp(argv[i], "-seed") == 0 )
      {
         i++;
         if( i < argc )
         {
            randomizeSeed = atoi(argv[i]);
         }
         else
         {
            std::cerr << "missing number after parameter '-seed'" << std::endl;
            exit(1);
         }
      }
#ifdef UG_WITH_UGS
      else if( strcmp(argv[i], "-ugsc") == 0 )
      {
         i++;
      }
#endif
      else
      {
         THROW_LOGICAL_ERROR3("invalid parameter <", argv[i], ">");
      }
   }

   ///
   /// Setup
   ///
   // initialize CMAP_LAP
   probFileName = argv[2];
   char tempName[1024];
   std::strcpy(tempName, probFileName);
   char *headOfName = tempName;
   for( int j = std::strlen(tempName); j > 0; j-- )
   {
      if( tempName[j] == '.' )
      {
         tempName[j] = '\0';
      }
      if( tempName[j] == '/' )
      {
         headOfName = &tempName[j+1];
         break;
      }
   }
   // probname <- saveFileName
   probName = new char[std::strlen(headOfName) + 1]; // メモリ確保
   std::strcpy(probName, headOfName);
   std::cout << probName << std::endl;

   ///
   /// Setup clock type
   ///
   if( paraParams->getRealParamValue(UG::TimeLimit) > 0.0 )
   {
      double timeRemains =  paraParams->getRealParamValue(UG::TimeLimit) - timer->getElapsedTime();
      std::cout << "The value of MAXTIME is now " << ceil(timeRemains) << std::endl;
   }

   if( basisprobname != NULL )
   {
      ///
      /// Load settings
      ///
      DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm );

      ///
      /// Start CMAP_LAP
      ///
      // Problem Creation
      instance = std::shared_ptr<CMapLapParaInstance>(cmapLapParaComm->createCMapLapParaInstance(probName, 0));
      std::cout << "*** CMapLap runs with " << nCMapLapThreads << " threads ***" << std::endl;

      // read input file
      lattice = Lattice<int, double>(std::string(basisprobname));

      std::ostringstream os;
      if( solutionFileName )
      {
         os << solutionFileName;
      }
      else
      {
         os << paraParams->getStringParamValue(UG::SolutionFilePath);
         os << instance->getProbName() << ".sol";
      }
      solutionFile = fopen(os.str().c_str(), "a");  // if solution file exists, append
      if( solutionFile == NULL )
      {
         THROW_LOGICAL_ERROR3("cannot open solution file <", os.str(), "> for writing");
      }

      //
      // Presolve the original instance
      //
   }
   else
   {
      std::cout << std::endl;
      std::cout << "syntax: " << argv[0] << "#solvers ppcmapLap_param_file problem_file_name "
                << "[-l <logfile>] [-q] [-sl <settings>] [-s <settings>] [-sr <root_settings>] [-w <prefix_warm>] [-sth <number>]" << std::endl;
      std::cout << "  -l <logfile>        : copy output into log file" << std::endl;
      std::cout << "  -q                  : suppress screen messages" << std::endl;
      std::cout << "  -sl <settings>      : load parameter settings (.set) file for LC presolving" << std::endl;
      std::cout << "  -s <settings>       : load parameter settings (.set) file for solvers" << std::endl;
      std::cout << "  -sr <root_settings> : load parameter settings (.set) file for root" << std::endl;
      std::cout << "  -w <prefix_warm>    : warm start file prefix ( prefix_warm_tasks.gz and prefix_warm_solution.txt are read )" << std::endl;
      std::cout << "  -sth <number>       : the number of solver threads used(FiberCMAP_LAP)" << std::endl;
      THROW_LOGICAL_ERROR1("invalid parameter");
   }

   return 0;
}

///
/// try to set globalBestCMapLapSolution to sol if sol is better than it
/// @return true if globalBestCMapLapSolution is replased else false
///
bool
CMapLapParaInitiator::tryToSetIncumbentSolution(
      std::shared_ptr<CMapLapParaSolution> sol
      )
{
   if( globalBestCMapLapSolution )
   {
      if( globalBestCMapLapSolution->getObjectiveFunctionValue() - 0.0001
            < sol->getObjectiveFunctionValue() )
      {
         return false;
      }
   }
   globalBestCMapLapSolution = sol;
   nSolutions++;
   if( !paraParams->getBoolParamValue(UG::Quiet) )
   {
      writeSolution("");
   }

   return true;
}

void
CMapLapParaInitiator::sendSolverInitializationMessage(
      )
{
   int warmStarted = 0;
   if( isWarmStarted() )
   {
      warmStarted = 1;
   }
   paraComm->bcast(&warmStarted,1, UG::ParaINT, 0);
}

void
CMapLapParaInitiator::writeSolution(
      const std::string& message
      )
{
   if( message == "Final Solution" )
   {
      fprintf(solutionFile, "[ Final Solution ]\n");
   }
   else
   {
      fprintf(solutionFile,"%s\n",message.c_str());
   }

   if( globalBestCMapLapSolution )   // solution keeps the best one
   {
      printSolution();
      fflush(solutionFile);
   }
   else
   {
      fprintf(solutionFile, "No Solution\n");
   }
}

void
CMapLapParaInitiator::writeParaInstance(
      const std::string& filename
      )
{
   // CMAP_LAP_CALL_ABORT( XPRSwriteprob(cmapLapProb, filename.c_str() , "lps") );
}


#ifdef UG_WITH_ZLIB
///
/// write checkpoint solution
///
void
CMapLapParaInitiator::writeCheckpointSolution(
      gzstream::ogzstream &out
      )
{
   if( globalBestCMapLapSolution )
   {
      globalBestCMapLapSolution->writeCheckpointSolution(out);
   }
}

///
/// write checkpoint solution
///
void
CMapLapParaInitiator::writeCheckpointSolution(
      const std::string& filename
      )
{
   gzstream::ogzstream checkpointSolutionStream;
   checkpointSolutionStream.open(filename.c_str(), std::ios::out | std::ios::binary);
   writeCheckpointSolution(checkpointSolutionStream);
   checkpointSolutionStream.close();
}

///
/// read solution from checkpoint file
///
double
CMapLapParaInitiator::readSolutionFromCheckpointFile(
      char *afterCheckpointingSolutionFileName
      )
{
   char tempSolutionFileName[256];
   sprintf(tempSolutionFileName,"%s_solution.gz", prefixWarm );
   gzstream::igzstream  checkpointSolutionStream;
   checkpointSolutionStream.open(tempSolutionFileName, std::ios::in | std::ios::binary);
   if( checkpointSolutionStream )
   {
      std::shared_ptr<CMapLapParaSolution> sol{dynamic_cast<CMapLapParaSolution*>(paraComm->createParaSolution())};
      if( sol->read(paraComm, checkpointSolutionStream) )
      {
         if( !globalBestCMapLapSolution )
         {
            globalBestCMapLapSolution = sol;
         }
         else
         {
            if( globalBestCMapLapSolution->getObjectiveFunctionValue() > sol->getObjectiveFunctionValue() )
            {
               globalBestCMapLapSolution = sol;
            }
         }
      }
   }
   checkpointSolutionStream.close();

   // check if after checkpoing solution file exists or not
   gzstream::igzstream afterCheckpointSolutionStream;
   afterCheckpointSolutionStream.open(afterCheckpointingSolutionFileName, std::ios::in | std::ios::binary);
   if( afterCheckpointSolutionStream )
   {
      // set up from after checkpointing solution file
      std::shared_ptr<CMapLapParaSolution> sol{dynamic_cast<CMapLapParaSolution*>(paraComm->createParaSolution())};
      if( sol->read(paraComm, afterCheckpointSolutionStream) )
      {
         if( !globalBestCMapLapSolution )
         {
            globalBestCMapLapSolution = sol;
            std::cout << "***** read after checkpoint solution is ****************" << std::endl;
            std::cout << globalBestCMapLapSolution->toString() << std::endl;
            std::cout << "from " << afterCheckpointingSolutionFileName << std::endl;
         }
         else
         {
            if( globalBestCMapLapSolution->getObjectiveFunctionValue() > sol->getObjectiveFunctionValue() )
            {
               globalBestCMapLapSolution = sol;
               std::cout << "***** read after checkpoint solution is ****************" << std::endl;
               std::cout << globalBestCMapLapSolution->toString() << std::endl;
               std::cout << "from " << afterCheckpointingSolutionFileName << std::endl;
            }
         }
      }
   }
   afterCheckpointSolutionStream.close();

   if( globalBestCMapLapSolution )
   {
      return globalBestCMapLapSolution->getObjectiveFunctionValue();
   }
   else
   {
      return DBL_MAX;
   }
}
#endif


///
/// read solution file
///
double
CMapLapParaInitiator::readSolutionFile(
      char *loadSolutionFileName
      )
{
   DEF_CMAP_LAP_PARA_COMM( cmapLapParaComm, paraComm);

   std::string loadSolutionFileNameString{loadSolutionFileName};
   Lattice<int, double> tmpLattice{loadSolutionFileName};
   int threadId = 0;
   double squaredNorm = tmpLattice.B[0];
   LatticeVector<int> v = tmpLattice.basis.row(0);
   std::shared_ptr<CMapLapParaSolution> sol{dynamic_cast<CMapLapParaSolution*>(
         cmapLapParaComm->createCMapLapParaSolution(threadId, v, squaredNorm)
         )};
   if( !globalBestCMapLapSolution )
   {
      globalBestCMapLapSolution = sol;
      std::cout << "***** read after checkpoint solution is ****************" << std::endl;
      std::cout << globalBestCMapLapSolution->toString() << std::endl;
      std::cout << "from " << loadSolutionFileName << std::endl;
   }
   else
   {
      if( globalBestCMapLapSolution->getObjectiveFunctionValue() > sol->getObjectiveFunctionValue() )
      {
         globalBestCMapLapSolution = sol;
         std::cout << "***** read initial solution is ****************" << std::endl;
         std::cout << "objective value " << sol->getObjectiveFunctionValue() << " from " << loadSolutionFileName << std::endl;
         std::cout << "from " << loadSolutionFileName << std::endl;
      }
   }

   if( globalBestCMapLapSolution )
   {
      return globalBestCMapLapSolution->getObjectiveFunctionValue();
   }
   else
   {
      return DBL_MAX;
   }
}


///
/// print solver version
///
void
CMapLapParaInitiator::printSolverVersion(
      std::ostream *os           ///< output file (or NULL for standard output)
      )
{
   char version[16] = "1.0";
   // CMAP_LAP_CALL_ABORT( XPRSgetversion(version) );
   if( os )
   {
      (*os) << "CMapLap-Optimizer version " << version << std::endl;
   }
   else
   {
      std::cout << "CMapLap-Optimizer version " << version << std::endl;
   }
}

///
/// set final solver status
///
void
CMapLapParaInitiator::setFinalSolverStatus(
      FinalSolverState state
      )
{
   finalState = state;
}

///
/// output solution status
///
void
CMapLapParaInitiator::outputFinalSolverStatistics(
      std::ostream *os,
      double time,
      double totalDeepBkzTime,
      double totalEnumTime,
      double totalSieveTime
      )
{
   if( os == 0 )
   {
      os = &std::cout;
   }

   if( finalState !=  Aborted )
   {
      *os << "CMapLap Status        : ";
   }

   switch ( finalState )
   {
   case InitialTasksGenerated:
      *os << "initial tasks were generated" << std::endl;
      break;
   case Aborted:
      *os << std::endl;
      break;
   case HardTimeLimitIsReached:
      *os << "solving was interrupted [ hard time limit reached ]" << std::endl;
      break;
   case LowerBoundIsReached:
      *os << "solving was interrupted [ lower bound reached ]" << std::endl;
      break;
   case ComputingWasInterrupted:
      *os << "solving was interrupted" << std::endl;
      break;
   case ProblemWasSolved:
      *os << "problem is solved" << std::endl;
      break;
   case RequestedSubProblemsWereSolved:
      *os << "requested subproblems are solved" << std::endl;
      break;
   default:
      THROW_LOGICAL_ERROR1("invalid final state");
   }

   double shortestNorm = 0.0;
   if( globalBestCMapLapSolution )
   {
      shortestNorm = std::sqrt(globalBestCMapLapSolution->getVector().squaredNorm());
   }
   *os << "Total Wall Time       : " << time << std::endl;
   *os << "  LC process          : " << 0 << std::endl;
   *os << "  solving             : " << time << std::endl;
   *os << "Total CPU Time        : " << totalDeepBkzTime + totalEnumTime + totalSieveTime << std::endl;
   *os << "  total DeepBkz time  : " << totalDeepBkzTime << std::endl;
   *os << "  total Enum    time  : " << totalEnumTime << std::endl;
   *os << "  total Sieve   time  : " << totalSieveTime << std::endl;
   *os << "Task                  : " << std::endl;
   *os << "  tasks (total)       : " << nSolved << std::endl;
   *os << "Best Objective        : " << shortestNorm << std::endl;
   *os << "  Approx. Factor      : " << lattice.approxFactor(shortestNorm) << std::endl;
   *os << "  Hermite Factor      : " << lattice.hermiteFactor(shortestNorm) << std::endl;
   *os << "  Root Hermite Factor : " << lattice.rootHermiteFactor(shortestNorm) << std::endl;
}

///
/// print solver version
///
void
CMapLapParaInitiator::printSolution(
      )
{
   assert( globalBestCMapLapSolution );

   fprintf(solutionFile, "objective value: ");
   printReal(solutionFile, globalBestCMapLapSolution->getObjectiveFunctionValue(), 20, 15);
   std::ostringstream s;
   s << globalBestCMapLapSolution->getVector().transpose();
   fprintf(solutionFile, "\n%s\n", s.str().c_str());
}


///
/// outputs a real number, or "+infinity", or "-infinity" to a file
///
void
CMapLapParaInitiator::printReal(
   FILE*   file,       ///< output file (or NULL for standard output)
   double  val,        ///< value to print
   int     width,      ///< width of the field
   int     precision   ///< number of significant digits printed
   )
{
   const static int MAXSTRLEN = 1024;
   char s[MAXSTRLEN];
   char strformat[MAXSTRLEN];

   if( ( REALABS(val) > CMAP_LAP_INFINITY) )
      (void) snprintf(s, 1024, "+infinity");
   else if( REALABS(-val) > CMAP_LAP_INFINITY )
      (void) snprintf(s, MAXSTRLEN, "-infinity");
   else
   {
      (void) snprintf(strformat, MAXSTRLEN, "%%.%dg", precision);
      (void) snprintf(s, MAXSTRLEN, (const char*)strformat, val);
   }
   (void) snprintf(strformat, MAXSTRLEN, "%%%ds", width);
   (void) fprintf(file, (const char*)strformat, s);
}

} // namespace ParaCMapLAP
