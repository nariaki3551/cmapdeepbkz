/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*          This file is part of the program and software framework          */
/*      CMAP-ENUM --- Configurable Massively Parallel Solver for ENUM        */
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

/**@file    synchronousCmapDeepBkzConfig.h
 * @brief   Config extension for synchronousCmapdeepbkz.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SYNCHRONOUS_CMAP_DEEPBKZ_CONFIG_H__
#define __SYNCHRONOUS_CMAP_DEEPBKZ_CONFIG_H__

#include <iostream>
#include <fstream>
#include <cstring>
#include <unistd.h>
#include <cfloat>
#include "Config.h"


class SynchronousCMapDeepBkzConfig : public LapTools::Config
{

public:

   int NumSolvers = 1;           // number of solvers
   int NumSharedDimension = 0;   // number of shared dimension
   int NumThreads = 1;           // number of threads used
   int TotalTours = 100;         // number of total tours
   int Dimension = 80;           // dimension of lattice basis


   ///
   /// constructor
   ///
   SynchronousCMapDeepBkzConfig(): Config() {}
   SynchronousCMapDeepBkzConfig(
         std::string paramFilePath
         )
   {
      load(paramFilePath);
   }


   ///
   /// @brief set parameter
   /// @param[in] argv arguments
   /// @param[in] i index of arguments
   /// @return true if success to read i-th argument else false
   ///
   bool setParam(
         char **argv,
         int &i
         )
   {
      if( std::strcmp(argv[i], "-b") == 0 || std::strcmp(argv[i], "--beta") == 0 )
      {
         beta = atoi(argv[++i]);
      }
      else if ( std::strcmp(argv[i], "-c") == 0 || std::strcmp(argv[i], "--csv_output") == 0 )
      {
         CsvLogFile = std::string(argv[++i]);
         OutputCsvLogFile = true;
      }
      else if ( std::strcmp(argv[i], "--dimension") == 0 )
      {
         Dimension = std::atoi(argv[++i]);
      }
      else if ( std::strcmp(argv[i], "-i") == 0 || std::strcmp(argv[i], "--input") == 0 )
      {
         InputFile = std::string(argv[++i]);
      }
      else if ( std::strcmp(argv[i], "--num_solvers") == 0 )
      {
         NumSolvers = atoi(argv[++i]);
      }
      else if ( std::strcmp(argv[i], "--num_shares") == 0 )
      {
         NumSharedDimension = atoi(argv[++i]);
      }
      else if ( std::strcmp(argv[i], "--num_threads") == 0 )
      {
         NumThreads = atoi(argv[++i]);
      }
      else if ( std::strcmp(argv[i], "-o") == 0 || std::strcmp(argv[i], "--output") == 0 )
      {
         OutputFile = argv[++i];
      }
      else if ( std::strcmp(argv[i], "--total_tour") == 0 )
      {
         TotalTours = atoi(argv[++i]);
      }
      else if ( std::strcmp(argv[i], "--random_seed") == 0 )
      {
         RandomizeSeed = atoi(argv[++i]);
      }
      else if ( std::strcmp(argv[i], "-q") == 0 || std::strcmp(argv[i], "--quiet") == 0 )
      {
         Quiet = true;
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
   std::string toString(
         )
   {
      std::ostringstream s;
      s << Config::toString();
      s << "NumSolvers                     = " << NumSolvers                     << std::endl;
      s << "NumSharedDimension             = " << NumSharedDimension             << std::endl;
      s << "NumThreads                     = " << NumThreads                     << std::endl;
      s << "TotalTours                     = " << TotalTours                     << std::endl;
      s << "Dimension                      = " << Dimension                      << std::endl;
      return s.str();
   }


   ///
   /// @brief laod parameters from file
   /// @param[in] paramFile parameter file path
   /// @return true if success to read paramters else false
   ///
   bool load(
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
         if( !Config::setParam(name_str, value) )
         {
            if( name_str == "NumSolvers" )
            {
               NumSolvers = std::atoi(value);
            }
            else if( name_str == "NumSharedDimension" )
            {
               NumSharedDimension = std::atoi(value);
            }
            else if( name_str == "NumThreads" )
            {
               NumThreads = std::atoi(value);
            }
            else if( name_str == "TotalTours" )
            {
               TotalTours = std::atoi(value);
            }
            else if( name_str == "Dimension" )
            {
               Dimension = std::atoi(value);
            }
            else
            {
               std::cout << "param name " << name
                  << " is invalid <line is " << lineNo
                  << "; " << paramFile.c_str() << std::endl;
               fclose(fp);
               return false;
            }
         }
      }
      fclose(fp);
      return true;
   }

}; // class SynchronousCMapDeepBkzConfig

#endif  // __SYNCHRONOUS_CMAP_DEEPBKZ_CONFIG_H__
