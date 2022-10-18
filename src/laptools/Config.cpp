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

/**@file    Config.cpp
 * @brief   Hyper parameters.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "Config.h"
#include <sstream>


namespace LapTools
{


///
/// @brief set parameters from argument
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @return true if success to read parameters else false
///
bool
Config::setParams(
      int argc,
      char **argv
      )
{
   for( int i = 0; i < argc; ++i )
   {
      if( std::strcmp(argv[i], "-h") == 0 || std::strcmp(argv[i], "--help") == 0 )
      {
         return false;
      }
   }

   /// load parameter file
   for( int i = 0; i < argc; ++i )
   {
      if( std::strcmp(argv[i], "--param" ) == 0 )
      {
         load(std::string(argv[++i]));
         break;
      }
   }

   /// Parse parameters
   for( int i = 1; i < argc; ++i )   /// the first argument is runtime parameter file for ParaCMAP_LAP
                                     /// the second argument is problem file name
   {
      if( !setParam(argv, i) )
      {
         std::cout << "invalid parameter <" << argv[i] << ">" << std::endl;
         return false;
      }
   }
   return true;
}


///
/// @brief set parameter
/// @param[in] argv arguments
/// @param[in] i index of arguments
/// @return true if success to read i-th argument else false
///
bool
Config::setParam(
      char **argv,
      int &i
      )
{
   if( std::strcmp(argv[i], "-a") == 0 || std::strcmp(argv[i], "--algorithm") == 0 )
   {
      Algorithm = std::string(argv[++i]);
   }
   else if( std::strcmp(argv[i], "-b") == 0 || std::strcmp(argv[i], "--beta") == 0 )
   {
      beta = atoi(argv[++i]);
   }
   else if ( std::strcmp(argv[i], "-c") == 0 || std::strcmp(argv[i], "--csv_output") == 0 )
   {
      CsvLogFile = std::string(argv[++i]);
      OutputCsvLogFile = true;
   }
   else if ( std::strcmp(argv[i], "-i") == 0 || std::strcmp(argv[i], "--input") == 0 )
   {
      InputFile = std::string(argv[++i]);
   }
   else if ( std::strcmp(argv[i], "-l") == 0 || std::strcmp(argv[i], "--lower_bound") == 0 )
   {
      LowerBound = atof(argv[++i]);
   }
   else if ( std::strcmp(argv[i], "-o") == 0 || std::strcmp(argv[i], "--output") == 0 )
   {
      OutputFile = argv[++i];
   }
   else if( std::strcmp(argv[i], "--basis_float") == 0 )
   {
      BasisFloat = std::string(argv[++i]);
   }
   else if( std::strcmp(argv[i], "--gs_float") == 0 )
   {
      GSFloat = std::string(argv[++i]);
   }
   else if( std::strcmp(argv[i], "--enum_gs_float") == 0 )
   {
      EnumGSFloat = std::string(argv[++i]);
   }
   else if ( std::strcmp(argv[i], "-q") == 0 || std::strcmp(argv[i], "--quiet") == 0 )
   {
      Quiet = true;
   }
   else if ( std::strcmp(argv[i], "--random_seed") == 0 )
   {
      RandomizeSeed = atoi(argv[++i]);
   }
   else if ( std::strcmp(argv[i], "-t") == 0 || std::strcmp(argv[i], "--time_limit") == 0 )
   {
      TimeLimit = atof(argv[++i]);
   }
   else if ( std::strcmp(argv[i], "-v") == 0 || std::strcmp(argv[i], "--verbose") == 0 )
   {
      ++i;
   }
   else if ( std::strcmp(argv[i], "--param") == 0 )
   {
      ++i;
   }
   else
   {
      return false;
   }
   return true;
}


///
/// @breif convert to string
///
std::string
Config::toString(
      )
{
   std::ostringstream s;
   s << "# Config parameter" << std::endl;
   /// Inputs
   s << "InputFile                      = " << InputFile                      << std::endl;
   s << "Algorithm                      = " << Algorithm                      << std::endl;
   s << "OutputFile                     = " << OutputFile                     << std::endl;
   s << "ParamFile                      = " << ParamFile                      << std::endl;
   /// LogFile
   s << "Quiet                          = " << Quiet                          << std::endl;
   s << "LogFile                        = " << LogFile                        << std::endl;
   s << "CsvLogFile                     = " << CsvLogFile                     << std::endl;
   s << "OutputLogFile                  = " << OutputLogFile                  << std::endl;
   s << "OutputCsvLogFile               = " << OutputCsvLogFile               << std::endl;
   /// Floating
   s << "BasisFloat                     = " << BasisFloat                     << std::endl;
   s << "GSFloat                        = " << GSFloat                        << std::endl;
   s << "EnumGSFloat                    = " << EnumGSFloat                    << std::endl;
   /// Global Parameter
   s << "TimeLimit                      = " << TimeLimit                      << std::endl;
   s << "LowerBound                     = " << LowerBound                     << std::endl;
   /// Randomize
   s << "RandomizeType                  = " << RandomizeType                  << std::endl;
   s << "RandomizeSeed                  = " << RandomizeSeed                  << std::endl;
   s << "RandomizeScale                 = " << RandomizeScale                 << std::endl;
   /// GSO
   s << "GSOType                        = " << GSOType                        << std::endl;
   /// DeepLll
   s << "eta                            = " << eta                            << std::endl;
   s << "delta                          = " << delta                          << std::endl;
   s << "deep                           = " << deep                           << std::endl;
   /// MLLL
   s << "MLLLdelta                      = " << MLLLdelta                      << std::endl;
   /// Enum
   s << "EnumTraverseNodesPerSecond     = " << EnumTraverseNodesPerSecond     << std::endl;
   s << "EnumPruningLowerBlockSize      = " << EnumPruningLowerBlockSize      << std::endl;
   s << "EnumPruningGH                  = " << EnumPruningGH                  << std::endl;
   s << "EnumPruningLowerBound          = " << EnumPruningLowerBound          << std::endl;
   s << "EnumPruningGaussEnum           = " << EnumPruningGaussEnum           << std::endl;
   s << "EnumPruningUseCache            = " << EnumPruningUseCache            << std::endl;
   /// subEnum
   s << "ProjectedDim                   = " << ProjectedDim                   << std::endl;
   s << "lambda                         = " << lambda                         << std::endl;
   /// DeepBkz
   s << "beta                           = " << beta                           << std::endl;
   s << "EnumProbInBkz                  = " << EnumProbInBkz                  << std::endl;
   s << "EnumTimeLimitInBkz             = " << EnumTimeLimitInBkz             << std::endl;
   s << "AutoAbortThreashold            = " << AutoAbortThreashold            << std::endl;
   /// Recursive DeepBkz
   s << "RecursiveMaxTour               = " << RecursiveMaxTour               << std::endl;
   s << "RecursiveLowerBeta             = " << RecursiveLowerBeta             << std::endl;
   /// Log Parameter
   s << "OutputUpdateVector             = " << OutputUpdateVector             << std::endl;
   s << "LogIntervalDeepLll             = " << LogIntervalDeepLll             << std::endl;
   s << "LogIntervalEnum                = " << LogIntervalEnum                << std::endl;
   s << "LogIntervalSubEnum             = " << LogIntervalSubEnum             << std::endl;
   s << "LogIntervalDeepBkz             = " << LogIntervalDeepBkz             << std::endl;
   s << "LogIntervalGaussSieve          = " << LogIntervalGaussSieve          << std::endl;
   return s.str();
}


///
/// @brief laod parameters from file
/// @param[in] paramFile parameter file path
/// @return true if success to read paramters else false
///
bool
Config::load(
      std::string paramFile
      )
{
   if( paramFile.empty() )
   {
      return false;
   }
   std::ifstream fin(paramFile);
   if( !fin.is_open() )
   {
      std::cout << "can not open parameter file: " << paramFile << std::endl;
      return false;
   }

   FILE* fp = fopen(paramFile.c_str(), "r");
   char line[256], name[256], dummy[256], value[256];
   int lineNo = 0;
   std::string name_str;

   while( fgets(line, 256, fp) != NULL )
   {
      // skip head char is # or only \n
      lineNo++;
      if((line[0] == '#') || (!strcmp(line, "\n"))){
         continue;
      }
      sscanf(line, "%s %s %s\n", name, dummy, value);
      name_str = std::string(name);
      if( !setParam(name_str, value) )
      {
         std::cout << "param name " << name
            << " is invalid <line is " << lineNo
            << "; " << paramFile.c_str() << std::endl;
         fclose(fp);
         return false;
      }
   }
   fclose(fp);
   return true;
}


///
/// @brief load parameter
/// @param[in] name_str string of parameter
/// @param[in] value string of parameter value
/// @return true if success to read paramter else false
///
bool
Config::setParam(
      std::string name_str,
      char * value
      )
{
   if( name_str == "InputFile" ){
      InputFile = std::string(value);
   }
   else if( name_str == "Algorithm" ){
      Algorithm = std::string(value);
   }
   else if( name_str == "OutputFile" ){
      OutputFile = std::string(value);
   }
   else if( name_str == "ParamFile" ){
      ParamFile = std::string(value);
   }
   ///
   /// LogFile
   ///
   else if( name_str == "Quiet" ){
      Quiet = atoi(value);
   }
   else if( name_str == "LogFile" ){
      LogFile = std::string(value);
      OutputLogFile = true;
   }
   else if( name_str == "CsvLogFile" ){
      CsvLogFile = std::string(value);
      OutputCsvLogFile = true;
   }
   ///
   /// Floating
   ///
   else if( name_str == "BasisFloat" ){
      BasisFloat = std::string(value);
   }
   else if( name_str == "GSFloat" ){
      GSFloat = std::string(value);
   }
   else if( name_str == "EnumGSFloat" ){
      EnumGSFloat = std::string(value);
   }
   ///
   /// Global Parameter
   ///
   else if( name_str == "TimeLimit" ){
      TimeLimit = atof(value);
   }
   else if( name_str == "LowerBound" ){
      LowerBound = atof(value);
   }
   ///
   /// Randomize
   ///
   else if( name_str == "RandomizeType" ){
      RandomizeType = std::string(value);
   }
   else if( name_str == "RandomizeSeed" ){
      RandomizeSeed = atoi(value);
   }
   else if( name_str == "RandomizeScale" ){
      RandomizeScale = atof(value);
   }
   ///
   /// GSO
   ///
   else if( name_str == "GSOType" ){
      GSOType = atoi(value);
   }
   ///
   /// DeepLll
   ///
   else if( name_str == "eta" ){
      eta = atof(value);
   }
   else if( name_str == "delta" ){
      delta = atof(value);
   }
   else if( name_str == "deep" ){
      deep = atoi(value);
   }
   ///
   /// MLLL
   ///
   else if( name_str == "MLLLdelta" ){
      MLLLdelta = atof(value);
   }
   ///
   /// Enum
   ///
   else if( name_str == "EnumTraverseNodesPerSecond" ){
      EnumTraverseNodesPerSecond = atof(value);
   }
   else if( name_str == "EnumPruningLowerBlockSize" ){
      EnumPruningLowerBlockSize = atoi(value);
   }
   else if( name_str == "EnumPruningGH" ){
      EnumPruningGH = atof(value);
   }
   else if( name_str == "EnumPruningLowerBound" ){
      EnumPruningLowerBound = atoi(value);
   }
   else if( name_str == "EnumPruningGaussEnum" ){
      EnumPruningGaussEnum = atoi(value);
   }
   else if( name_str == "EnumPruningUseCache" ){
      EnumPruningUseCache = atoi(value);
   }
   ///
   /// SubEnum
   ///
   else if( name_str == "ProjectedDim" ){
      ProjectedDim = atoi(value);
   }
   else if( name_str == "lambda" ){
      lambda = atof(value);
   }
   ///
   /// DeepBkz
   ///
   else if( name_str == "beta" ){
      beta = atoi(value);
   }
   else if( name_str == "EnumProbInBkz" ){
      EnumProbInBkz = atof(value);
   }
   else if( name_str == "EnumTimeLimitInBkz" ){
      EnumTimeLimitInBkz = atof(value);
   }
   else if( name_str == "AutoAbortThreashold" ){
      AutoAbortThreashold = atof(value);
   }
   ///
   /// Recursive DeepBkz
   ///
   else if( name_str == "RecursiveMaxTour" ){
      RecursiveMaxTour = atof(value);
   }
   else if( name_str == "RecursiveLowerBeta" ){
      RecursiveLowerBeta = atof(value);
   }
   ///
   /// Log Parameters
   ///
   else if( name_str == "OutputUpdateVector" ){
      OutputUpdateVector = atoi(value);
   }
   else if( name_str == "LogIntervalDeepLll" ){
      LogIntervalDeepLll = atof(value);
   }
   else if( name_str == "LogIntervalEnum" ){
      LogIntervalEnum = atof(value);
   }
   else if( name_str == "LogIntervalSubEnum" ){
      LogIntervalSubEnum = atof(value);
   }
   else if( name_str == "LogIntervalDeepBkz" ){
      LogIntervalDeepBkz = atof(value);
   }
   else if( name_str == "LogIntervalGaussSieve" ){
      LogIntervalGaussSieve = atof(value);
   }
   else
   {
      return false;
   }
   return true;
}

}  // namespace LapTools
