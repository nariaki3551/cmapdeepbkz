/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* MIT License                                                                     */
/*                                                                                 */
/* Copyright (c) 2022 Nariaki Tateiwa <n-tateiwa@kyudai.jp>                        */
/*                                                                                 */
/* Permission is hereby granted, free of charge, to any person obtaining a copy    */
/* of this software and associated documentation files (the "Software"), to deal   */
/* in the Software without restriction, including without limitation the rights    */
/* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell       */
/* copies of the Software, and to permit persons to whom the Software is           */
/* furnished to do so, subject to the following conditions:                        */
/*                                                                                 */
/* The above copyright notice and this permission notice shall be included in all  */
/* copies or substantial portions of the Software.                                 */
/*                                                                                 */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR      */
/* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        */
/* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     */
/* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          */
/* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,   */
/* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE   */
/* SOFTWARE.                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    synchronousCmapDeepBkzIdentity.cpp
 * @brief   main script for thread-parallel cmapdeepbkz.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <stdlib.h>
#include <algorithm>
#include <deque>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include <omp.h>

#include "Def.h"
#include "Reduction.h"
#include "DeepBkz.h"
#include "Lattice.h"
#include "Timer.h"
#include "synchronousCmapDeepBkzConfig.h"
#include "synchronousCmapDeepBkzSimilarity.h"


void usage(
      std::ostream *osLog = &std::cout
      )
{
   std::ostringstream s;
   s << "usage:" << std::endl;
   s << "    synchronous_cmapdeepbkz_identity                           " << std::endl;
   s << "       --dimension [dimension of identity basis]               " << std::endl;
   s << "       --param [param file]                                    " << std::endl;
   s << "       -b, --beta [beta (bkz;default 20)]                      " << std::endl;
   s << "       -c, --csv_output [csv log file path]                    " << std::endl;
   s << "       --num_solvers [number of solvers]                       " << std::endl;
   s << "       --num_shares [number of shared dimension]               " << std::endl;
   s << "       --num_threads [number of threads]                       " << std::endl;
   s << "       -o, --output [output file path]                         " << std::endl;
   s << "       --total_tour [total tour]                               " << std::endl;
   s << "       --random_seed [seed (randomize; default 0)]             " << std::endl;
   s << "       -q quiet output                                         " << std::endl;
   s << "       -h show help                                            " << std::endl;
   *osLog << s.str();
}


int main(int argc, char* argv[])
{

   if( argc < 3 )
   {
      usage();
      exit(1);
   }

   // initial values for options
   int verbose = 2;
   std::ostream *osLog = &std::cout;
   std::ofstream ofsLog;
   std::ofstream ofsCsvLog;
   std::ofstream ofsSimCsvLog;

   // set Config
   SynchronousCMapDeepBkzConfig config;
   config.setParams(argc, argv);
   int blocksize = config.beta;
   int numSolvers = config.NumSolvers;          // number of solvers
   int numShares = config.NumSharedDimension;   // number of shared dimension
   int numThreads = config.NumThreads;          // number of threads used
   int totalTours = config.TotalTours;          // number of total tours
   int seed = config.RandomizeSeed;

   // timer initialize
   LapTools::Timer::init();

   // set osLog, ofsLog, ofsCsvLog and ofsSimCsvLog
   if( config.Quiet )
   {
      ofsLog.open("/dev/null", std::ios::app);
      osLog = &ofsLog;
   }
   if( config.OutputLogFile )
   {
      ofsLog.open(config.LogFile);
      if( !ofsLog ){ std::cout << "can not open " << config.LogFile << std::endl; exit(1); }
      osLog = &ofsLog;
   }
   if( config.OutputCsvLogFile )
   {
      ofsCsvLog.open(config.CsvLogFile);
      if( !ofsCsvLog ){ std::cout << "can not open " << config.CsvLogFile << std::endl; exit(1); }
   }
   if( config.OutputCsvLogFile )
   {
      std::string SimCsvLogFile = "sim." + config.CsvLogFile;
      ofsSimCsvLog.open(SimCsvLogFile);
      if( !ofsSimCsvLog ){ std::cout << "can not open " << SimCsvLogFile << std::endl; exit(1); }
   }

   // read input file
   int dim = config.Dimension;
   LapTools::LatticeBasis<int> identity = 100 * MatrixXi::Identity(dim, dim);
   LapTools::Lattice<int, double> L(identity);
   L.setConfig(config);

   // display setting
   *osLog << "dimension   : " << dim                  << std::endl;
   *osLog << "outputfile  : " << config.OutputFile    << std::endl;
   *osLog << "paramfile   : " << config.ParamFile     << std::endl;
   *osLog << "csvlogfile  : " << config.CsvLogFile    << std::endl;
   *osLog << "beta        : " << blocksize            << std::endl;
   *osLog << "numSolvers  : " << numSolvers           << std::endl;
   *osLog << "numShares   : " << numShares            << std::endl;
   *osLog << "numThreads  : " << numThreads           << std::endl;
   *osLog << "totalTour   : " << totalTours           << std::endl;
   *osLog << "seed        : " << seed                 << std::endl;
   *osLog << std::endl;

   // initialize solver
   std::vector<LapTools::DeepBkz<int, double, double>> bkzObjs;
   bkzObjs.resize(numSolvers);
   for( int i = 0; i < numSolvers; ++i )
      bkzObjs[i].hasReduced = true;
   std::deque<Eigen::MatrixXi> basisDeque;

   // initialize restart flags and randomize seeds
   std::vector<bool> restart(numSolvers, false);
   std::vector<int>  seeds(numSolvers, seed);

   // initialize index for grassmann
   std::vector<int> grassmannIndexes(L.n);
   std::iota(grassmannIndexes.begin(), grassmannIndexes.end(), 0);
   BasisSimilarity::setGrassmannIndexes(grassmannIndexes);

   // header of logging
   ofsCsvLog << bkzObjs[0].getCsvLogHeader() << std::endl;
   ofsSimCsvLog << "tour,sigh,shortest,approx," << BasisSimilarity::outputSimilarityOfBasisHeader(L.n) << std::endl;

   // for logging
   double shortestNorm = L.shortestNorm();
   double approxFactor = L.approxFactor();

   // set global Lattice (global basis)
   LapTools::Lattice<int, double> globalLattice;
   if( numShares > 0 ){ globalLattice = L.copy(numShares); }


   // main loop
   for( int tour = 0; tour < totalTours; ++tour )
   {
      std::cout
         << "tour " << tour
         << " norm " << shortestNorm
         << " approx " << approxFactor
         << std::endl;

      for( int i = 0; i < numSolvers; ++i )
      {
         if( bkzObjs[i].idle() )
         {
            restart[i] = true;
            seeds[i] = seed++;
         }
      }

      // set instance
#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
      for( int i = 0; i < numSolvers; i++ )
      {
         if( restart[i] )
         {
            std::cout << "Solver " << i << " seed " << seeds[i] << " \r";
            std::cout.flush();
            auto _L = std::make_shared<LapTools::Lattice<int, double>>(L.basis);
            _L->setConfig(config);
            _L->randomize(seeds[i]);
            if( !bkzObjs[i].getLattice() )
            {
               bkzObjs[i] = LapTools::DeepBkz<int, double, double>{_L, i+1, 0, verbose};   // LapTools::Lattice, rank, threadId, verbose
               bkzObjs[i].setBlocksize(blocksize);
            }
            else
            {
               bkzObjs[i].resetLattice(_L);  // replace LapTools::Lattice
            }
         }
      }

      // run one tour
#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
      for( int i = 0; i < numSolvers; i++ )
      {
         if( !restart[i] )
         {
            bool shouldAbort;
            bkzObjs[i].tour(shouldAbort); // run one tour
         }
      }

      // for loggging
      for( int i = 0; i < numSolvers; ++i )
      {
         shortestNorm = std::min(shortestNorm, bkzObjs[i].getLattice()->shortestNorm());
	      approxFactor = std::min(approxFactor, bkzObjs[i].getLattice()->approxFactor());
      }

      // calculate similarity
      {
         basisDeque.clear();
         for( int i = 0; i < numSolvers; ++i )
            basisDeque.push_back(bkzObjs[i].getLattice()->basis);
         ofsSimCsvLog
            << tour
            << ",befor_merge,"
            << shortestNorm << ","
            << approxFactor << ","
            << BasisSimilarity::outputSimilarityOfBasis(basisDeque, tour, L.n, numThreads)
            << std::endl;
      }


      if( numShares > 0 )
      {
         // update global Lattice
         for( int i = 0; i < numSolvers; ++i )
         {
            int index;
            if( bkzObjs[i].getLattice()->isMoreReducedThan(globalLattice, index) )
            {
               std::cout
                  << "replace global LapTools::Lattice with " << i << " -th solver's LapTools::Lattice "
                  << std::endl
                  << globalLattice.B.head(numShares).transpose() << std::endl
                  << " <- " << std::endl
                  << bkzObjs[i].getLattice()->B.head(numShares).transpose()
                  << std::endl;
               globalLattice = bkzObjs[i].getLattice()->copy(numShares);
            }
         }

         // merge global LapTools::Lattice
#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
         for( int i = 0; i < numSolvers; ++i )
         {
            int index;
            if( !bkzObjs[i].idle() && globalLattice.isMoreReducedThan(*(bkzObjs[i].getLattice()), index) )
            {
               bkzObjs[i].getLattice()->merge(globalLattice);
            }
         }
      }

      // for loggging
      for( int i = 0; i < numSolvers; ++i )
      {
         shortestNorm = std::min(shortestNorm, bkzObjs[i].getLattice()->shortestNorm());
	      approxFactor = std::min(approxFactor, bkzObjs[i].getLattice()->approxFactor());
      }

      // calculate similarity
      {
         basisDeque.clear();
         for( int i = 0; i < numSolvers; ++i )
            basisDeque.push_back(bkzObjs[i].getLattice()->basis);
         ofsSimCsvLog
            << tour
            << ",after_merge,"
            << shortestNorm << ","
            << approxFactor << ","
            << BasisSimilarity::outputSimilarityOfBasis(basisDeque, tour, L.n, numThreads)
            << std::endl;
      }

      // write logging
      for( int i = 0; i < numSolvers; i++ )
      {
         char head = ' ';
         if( restart[i] ){ head = 'S'; restart[i] = false; }
         else if( bkzObjs[i].idle() ){ head= 'E'; }
         bkzObjs[i].outputLog(head, &ofsCsvLog);
      }
   }

   return 0;
}
