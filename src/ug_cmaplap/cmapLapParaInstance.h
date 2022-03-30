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

/**@file    cmapLapParaInstance.h
 * @brief   ParaInstance extenstion for CMAP_LAP solver.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_INSTANCE_H__
#define __CMAP_LAP_PARA_INSTANCE_H__

#include <cstring>
#include <memory>
#include "ug/paraDef.h"
#include "ug/paraInstance.h"
#include "cmapLapParaDef.h"

namespace ParaCMapLAP
{


///
/// @class CMapLapParaInstance
///
class CMapLapParaInstance : public UG::ParaInstance, public std::enable_shared_from_this<CMapLapParaInstance>
{
protected:
   int           lProbName;               ///< length of problem name
   char          *probName;               ///< problem name


public:
   ///
   /// default constructor
   ///
   CMapLapParaInstance(
         ) : lProbName(0), probName(0)
   {
   }

   ///
   /// constructor
   ///
   CMapLapParaInstance(
         char      *inProbFileName,
         char      *inProbName
         ) : lProbName(0), probName(0)
   {
      int lProbFileName =  strlen(inProbFileName);
      char *tempProbName = new char[lProbFileName+1];
      strcpy(tempProbName, inProbFileName);
      int l = strlen(tempProbName) - 1 ;
      for( ; l > 0 && tempProbName[l] != '/' ; l--){}
      if( tempProbName[l] == '/' )
      {
         lProbName = strlen(&tempProbName[l+1]);
         probName = new char[lProbName+1];
         strcpy(probName, &tempProbName[l+1]);
      }
      else
      {
         lProbName = strlen(tempProbName);
         probName = new char[lProbName+1];
         strcpy(probName, tempProbName);
      }
      delete [] tempProbName;
   }

   ///
   /// constructor
   ///
   CMapLapParaInstance(
         char      *inProbName
         ) : lProbName(0), probName(0)
   {
      lProbName = strlen(inProbName);
      probName = new char[lProbName + 1];
      strcpy(probName, inProbName);
   }

   ///
   /// constructor
   ///
   CMapLapParaInstance(
         const CMapLapParaInstance& cmapLapParaInstance
         )
   {
      lProbName = cmapLapParaInstance.lProbName;
      probName = new char[lProbName + 1];
      strcpy(probName, cmapLapParaInstance.probName);
   }

   ///
   /// destractor
   ///
   virtual ~CMapLapParaInstance(
         )
   {
      if( probName ) delete [] probName;
   }

   ///
   /// stringfy ParaCalculationState
   ///
   virtual const std::string toString(
         )
   {
      std::ostringstream s;

      s << "lProbName = " << lProbName << std::endl;
      s << "probName = " << probName << std::endl;

      return s.str();
   }

   ///
   /// get problem name
   ///
   virtual const char *getProbName(
      )
   {
      return probName;
   }

};

} // namespace ParaCMapLAP

#if defined(_COMM_MPI_WORLD)
#include "cmapLapParaInstanceMpi.h"
#endif

#if defined(_COMM_PTH)
#include "cmapLapParaInstanceTh.h"
#endif

#endif  // __CMAP_LAP_PARA_INSTANCE_H__
