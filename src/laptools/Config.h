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

/**@file    Config.h
 * @brief   Hyper parameters.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __LAPTOOLS_CONFIG_H__
#define __LAPTOOLS_CONFIG_H__

#include <iostream>
#include <fstream>
#include <cstring>
#include <unistd.h>
#include <cfloat>

namespace LapTools
{

class Config
{

public:
   ///
   /// Inputs
   ///
   std::string InputFile = "noset";    /// basis file path
   std::string Algorithm = "deepbkz";  /// algorithm name
   std::string OutputFile = "output";  /// saved basis file path
   std::string ParamFile;              /// parameter file name

   ///
   /// LogFile
   ///
   bool Quiet = false;              ///< if true, then do not output standard output
   std::string LogFile;             ///< logfile path
   std::string CsvLogFile;          ///< csv logfile path
   bool OutputLogFile    = false;   ///< whether output log in LogFile
   bool OutputCsvLogFile = false;   ///< whether output csv log in CsvLogFile

   ///
   /// Floating
   ///
   std::string BasisFloat  = "int";
   std::string GSFloat     = "double";
   std::string EnumGSFloat = "double";

   ///
   /// Global Parameter
   ///
   double TimeLimit = DBL_MAX;      ///< total timeLimit
   double LowerBound = 0;           ///< terminate when it finds the vector shoter than LowerBound.

   ///
   /// Randomize
   ///
   std::string RandomizeType = "Fplll";   ///< type of randome method
   int RandomizeSeed = 0;                 ///< Seed of the randomization
   double RandomizeScale = 1.2;           ///< parameter of randomized unimodular matrix
                                          ///< generate random coefficiences by this*random(0,1)

   ///
   /// GSO
   ///
   int GSOType = 0;  ///< type of caluculation method of GSO
                     ///< 0: CFA
                     ///< 1: Modified GSO

   ///
   /// DeepLll
   ///
   double eta = 0.501;  ///< parameter or sizeReduce
   double delta = 0.99; ///< parameter of swap bases
   int deep = -1;       ///< deep insertion limitter

   ///
   /// MLLL
   ///
   double MLLLdelta = 0.999;  ///< parameter of MLLL

   ///
   /// Enum
   ///
   double EnumTraverseNodesPerSecond = 3.3554432e7;   ///< number of nodes that can be searched per second (ITO)
   int    EnumPruningLowerBlockSize = 45;             ///< if the Enum's block size is smaller than it, then pruning is not done
   double EnumPruningGH = -1;                         ///< set sqR(upper bound of Enum) using Gaussian Heuristics
                                                      ///< if it is negative, then it is ignored
   bool   EnumPruningLowerBound = 1;                  ///< set sqR(upper bound of Enum) using upperbound when LowerBound is greather than 0
   bool   EnumPruningGaussEnum = 1;                   ///< update sqR(upper bound of Enum) using GaussEnum logic at first of Enum algorithm
   bool   EnumPruningUseCache = false;                ///< save and load curve info as cache files

   ///
   /// SubEnum
   ///
   int    ProjectedDim = 10;              ///< projected dimension, if this is negative, then this is ignored
   double lambda       = 1.05;            ///< projected dimension is determined by this paramter

   ///
   /// DeepBkz
   ///
   int    beta = 20;                ///< blocksize
   double EnumProbInBkz = 1.0;      ///< probability of Enum in Bkz
   double EnumTimeLimitInBkz = 1e9; ///< time limit of Enum in Bkz
   int    AutoAbortThreashold = 10; ///< auto abort if GSA Slope gets small in continuously AutoAbortThreashold times
                                    ///< if it is negative, then it is ignored
   ///
   /// RecursiveDeepBkz
   ///
   int RecursiveMaxTour = 200;
   double RecursiveLowerBeta = 1e100;

   ///
   /// Log Parameters
   ///
   bool   OutputUpdateVector       = 1;   ///< display update incumbent vector
   double LogIntervalDeepLll       = 1;   ///< main DeepLll
   double LogIntervalEnum          = 1;   ///< output log interval [s]
   double LogIntervalSubEnum       = 1;   ///< output log interval [s]
   double LogIntervalDeepBkz       = 1;   ///< DeepBkz
   double LogIntervalGaussSieve    = 1;   ///< Sieve

   ///
   /// constructor
   ///
   Config(){}
   Config( std::string paramFilePath )
   {
      load(paramFilePath);
   }

   ///
   /// @brief set parameters from argument
   /// @param[in] argc number of arguments
   /// @param[in] argv arguments
   /// @return true if success to read parameters else false
   ///
   virtual bool setParams(
         int argc,
         char **argv
         );


   ///
   /// @brief set parameter
   /// @param[in] argv arguments
   /// @param[in] i index of arguments
   /// @return true if success to read i-th argument else false
   ///
   virtual bool setParam(
         char **argv,
         int &i
         );


   ///
   /// @breif convert to string
   ///
   virtual std::string toString(
         );


   ///
   /// @brief laod parameters from file
   /// @param[in] paramFile parameter file path
   /// @return true if success to read paramters else false
   ///
   virtual bool load(
         std::string paramFile
         );


   ///
   /// @brief load parameter
   /// @param[in] name_str string of parameter
   /// @param[in] value string of parameter value
   /// @return true if success to read paramter else false
   ///
   virtual bool setParam(
         std::string name_str,
         char * value
         );

};  // class Config

}  // namespace LapTools

#endif  // __LAPTOOLS_CONFIG_H__
