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

/**@file    seqcmaplap.cpp
 * @brief   main script for seqcmaplap.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <string>

#include "DeepBkz.h"
#include "ExDeepBkz.h"
#include "DeepLll.h"
#include "Enumeration.h"
#include "GaussSieve.h"
#include "Lattice.h"
// #include "RecursiveDeepBkz.h"
#include "SubEnumeration.h"
#include "Config.h"
#include "Def.h"
#include "Log.h"
#include "Timer.h"


using namespace LapTools;


void usage(
      std::ostream *osLog = &std::cout
      )
{
   std::ostringstream s;
   s << "usage:" << std::endl;
   s << "    seqcmaplap                                                 " << std::endl;
   s << "       -i, --input [basis file]                                " << std::endl;
   s << "       --param [param file]                                    " << std::endl;
   s << "       -a, --algorithm [algorithm]                             " << std::endl;
   s << "       -b, --beta [beta (bkz;default 20)]                      " << std::endl;
   s << "       -c, --csv_output [csv log file path]                    " << std::endl;
   s << "       -t, --time_limit [time limit[s]]                        " << std::endl;
   s << "       -l, --lower_bound [lower bound]                         " << std::endl;
   s << "       -o, --output [output file path]                         " << std::endl;
   s << "       --random_seed [seed (randomize; default 0)]             " << std::endl;
   s << "       -v, --verbose [int]                                     " << std::endl;
   s << "       -q quiet output                                         " << std::endl;
   s << "       -h show help                                            " << std::endl;
   s << "                                                               " << std::endl;
   s << "    floating setting                                           " << std::endl;
   s << "       --basis_float [<int, long_int>; default int]            " << std::endl;
   s << "       --gs_float [<double, long_double>; default double]      " << std::endl;
   s << "       --enum_gs_float [<double, long double>; default double] " << std::endl;
   s << "                                                               " << std::endl;
   s << "    algorithm:                                                 " << std::endl;
   s << "        deeplll, deepbkz(default)                              " << std::endl;
   s << "        exdeepbkz                                              " << std::endl;
   s << "        enum, subenum                                          " << std::endl;
   s << "        gausssieve                                             " << std::endl;
   s << "        data(show lattice data)                                " << std::endl;
   *osLog << s.str();
}


template<typename BasisFloat=int, typename GSFloat=double, typename EnumGSFloat=double>
int run(
      Config &config,
      int verbose
   )
{
   double startTime = Timer::getElapsedTime();

   // read input file
   auto L = std::make_shared<Lattice<BasisFloat, GSFloat>>(config.InputFile);
   L->setConfig(config);

   std::string algorithm = config.Algorithm;
   int rank = 1, thread = 0;
   std::shared_ptr<std::ofstream> ofsCsvLog;
   if( config.OutputCsvLogFile )
   {
      ofsCsvLog = std::shared_ptr<std::ofstream>(new std::ofstream(config.CsvLogFile, std::ios::app));
   }

   std::cout << "run " << algorithm << std::endl;
   if( algorithm == "enum" )
   {
      DeepLll<BasisFloat, GSFloat> deeplllObj{L, rank, thread, verbose};
      deeplllObj.deeplll();
      LatticeVector<BasisFloat> v = L->basis.row(0);
      LatticeVector<BasisFloat> coeffs(L->n); coeffs(0) = 1;
      Enumeration<BasisFloat, GSFloat, EnumGSFloat> enumObj{L, rank, thread, verbose};
      if( config.OutputCsvLogFile )
      {
         enumObj.setOsCsvLog(ofsCsvLog);
         *ofsCsvLog << enumObj.getCsvLogHeader() << std::endl;
      }
      enumObj.init();
      enumObj.projectedEnum(v, coeffs, -1, config.TimeLimit);
      std::cout
         << L->toSimpleString(v)
         << "Time                 : " << Timer::getElapsedTime()
         << std::endl;
      return 0;
   }
   else if( algorithm == "subenum" )
   {
      LatticeVector<BasisFloat> v = L->basis.row(0);
      LatticeVector<BasisFloat> coeffs(L->n); coeffs(0) = 1;
      int projDim = config.ProjectedDim;
      SubEnumeration<BasisFloat, GSFloat, EnumGSFloat> enumObj{L, rank, thread, verbose};
      enumObj.init(0, L->m-1, projDim);
      if( config.OutputCsvLogFile )
      {
         enumObj.setOsCsvLog(ofsCsvLog);
         *ofsCsvLog << enumObj.getCsvLogHeader() << std::endl;
      }
      enumObj.subEnum(v, coeffs, config.TimeLimit);
      std::cout
         << L->toSimpleString(v)
         << "Time                 : " << Timer::getElapsedTime() - startTime
         << std::endl;
      return 0;
   }
   else if( algorithm == "deeplll" )
   {
      if( verbose > 0 ){ std::cout << Log::getLogHeader() << std::endl; }
      DeepLll<BasisFloat, GSFloat> lllObj{L, rank, thread, verbose};
      if( config.OutputCsvLogFile )
      {
         lllObj.setOsCsvLog(ofsCsvLog);
         *ofsCsvLog << lllObj.getCsvLogHeader() << std::endl;
      }
      lllObj.deeplll();
   }
   else if( algorithm == "deepbkz" )
   {
      int blocksize = config.beta;
      DeepBkz<BasisFloat, GSFloat, EnumGSFloat> bkzObj{L, rank, thread, verbose};
      if( config.OutputCsvLogFile )
      {
         bkzObj.setOsCsvLog(ofsCsvLog);
         *ofsCsvLog << bkzObj.getCsvLogHeader() << std::endl;
      }
      bkzObj.deepbkz(blocksize, config.TimeLimit);
   }
   else if( algorithm == "exdeepbkz" )
   {
      int blocksize = config.beta;
      ExDeepBkz<BasisFloat, GSFloat, EnumGSFloat> bkzObj{L, rank, thread, verbose};
      if( config.OutputCsvLogFile )
      {
         bkzObj.setOsCsvLog(ofsCsvLog);
         *ofsCsvLog << bkzObj.getCsvLogHeader() << std::endl;
      }
      bkzObj.deepbkz(blocksize, config.TimeLimit);
   }
   else if( algorithm == "gausssieve" )
   {
      if( verbose > 0 ){ std::cout << Log::getLogHeader() << std::endl; }
      int listsize = -1, stacksize = -1, maxCollision = 10000;
      LatticeVector<BasisFloat> v;
      GaussSieve<BasisFloat, GSFloat> sieveObj{L, rank, thread, verbose};
      sieveObj.init(listsize, stacksize, maxCollision);
      if( config.OutputCsvLogFile )
      {
         sieveObj.setOsCsvLog(ofsCsvLog);
         *ofsCsvLog << sieveObj.getCsvLogHeader() << std::endl;
      }
      sieveObj.gaussSieve(v, config.TimeLimit);
      std::cout
         << L->toSimpleString(v)
         << "Time                 : " << Timer::getElapsedTime() - startTime
         << std::endl;
      return 0;
   }
   else if( algorithm != "data" )
   {
      std::cerr << " algorithm " << algorithm << " is incorrect." << std::endl;
      return 0;
   }

   std::cout
      << L->toSimpleString()
      << "Time                 : " << Timer::getElapsedTime() - startTime
      << std::endl;

   L->writeBasis(config.OutputFile);
   std::cout << "save basefile as " << config.OutputFile << std::endl;

   return 0;
}


int main(int argc, char* argv[])
{
   if( argc < 3 )
   {
      usage();
      exit(1);
   }

   Config config;
   if( !config.setParams(argc, argv) )
   {
      usage();
      exit(1);
   }

   int verbose = 1;
   for( int i = 0; i < argc; ++i )
   {
      if( std::strcmp(argv[i], "-v") == 0 || std::strcmp(argv[i], "--verbose") == 0 )
      {
         verbose = atoi(argv[++i]);
         break;
      }
   }

   // initialize config
   if( config.Quiet ){ verbose = 0; }

   // initialize timer
   Timer::init();

   // run main process
   if( config.BasisFloat == "int"
         && config.GSFloat == "double"
         && config.EnumGSFloat == "double" )
   {
      run<int, double, double>(config, verbose);
   }
   else if( config.BasisFloat == "int"
         && config.GSFloat == "long_double"
         && config.EnumGSFloat == "double" )
   {
      run<int, long double, double>(config, verbose);
   }
   else if( config.BasisFloat == "int"
         && config.GSFloat == "long_double"
         && config.EnumGSFloat == "long_double" )
   {
      run<int, long double, long double>(config, verbose);
   }
   else if( config.BasisFloat == "long_int"
         && config.GSFloat == "double"
         && config.EnumGSFloat == "double" )
   {
      run<long int, double, double>(config, verbose);
   }
   else if( config.BasisFloat == "long_int"
         && config.GSFloat == "long_double"
         && config.EnumGSFloat == "double" )
   {
      run<long int, long double, double>(config, verbose);
   }
   else if( config.BasisFloat == "long_int"
         && config.GSFloat == "long_double"
         && config.EnumGSFloat == "long_double" )
   {
      run<long int, long double, long double>(config, verbose);
   }

}



