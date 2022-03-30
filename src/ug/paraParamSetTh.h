/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*          This file is part of the program and software framework          */
/*                    UG --- Ubquity Generator Framework                     */
/*                                                                           */
/*  Copyright Written by Yuji Shinano <shinano@zib.de>,                      */
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

/**@file    paraParamSetTh.h
 * @brief   
 * @author  Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PARA_PARAM_SET_TH_H__
#define __PARA_PARAM_SET_TH_H__
#include "paraDef.h"
#include "paraParamSet.h"

namespace UG
{

///
/// class ParaParamSetTh
///
class ParaParamSetTh : public ParaParamSet
{

public:

   ///
   /// constructor
   ///
   ParaParamSetTh(
         )
   {
   }

   ///
   /// constructor
   ///
   ParaParamSetTh(
         int inNParaParams
         )
         : ParaParamSet(inNParaParams)
   {
   }

   ///
   /// destructor
   ///
   ~ParaParamSetTh(
         )
   {
   }

   ///
   /// broadcast ParaParams
   /// @return always 0 (for future extensions)
   ///
   int bcast(
         ParaComm *comm,        ///< communicator used
         int root               ///< root rank for broadcast
         )
   {
      THROW_LOGICAL_ERROR1("bcast is called in ParaParamSetTh");
   }

};

}

#endif  // __PARA_PARAM_SET_TH_H__
