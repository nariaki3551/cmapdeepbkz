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

/**@file    cmapLapParaInstanceTh.cpp
 * @brief   CMapLapParaInstance extension for threads communication.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "cmapLapParaInstanceTh.h"
#include <cstring>
#include <memory>
#include "cmapLapParaSolverLocalComm.h"
#include "ug/paraTagDef.h"
#include "cmapLapParaTagDef.h"


namespace ParaCMapLAP
{

int
CMapLapParaInstanceTh::bcast(
      UG::ParaComm *comm,
      int root,
      int method
      )
{
#ifdef _COMM_CPP11
   CMapLapParaCommTh *commTh = dynamic_cast<CMapLapParaCommTh*>(comm);
   if( commTh )
   {
      // receiver side creation
      if( commTh->getRank() == root )
      {
         for( int i = 0; i < commTh->getSize(); i++ )
         {
            if( i != root )
            {
               PARA_COMM_CALL(
                  commTh->uTypeSend((void *)createDatatype(comm), ParaInstanceType, i, UG::TagParaInstance)
               );
            }
         }
      }
      else
      {
         std::unique_ptr<CMapLapParaInstanceTh> received;
         PARA_COMM_CALL(
            commTh->uTypeReceive((void **)&received, ParaInstanceType, root, UG::TagParaInstance)
         );
         lProbName    = received->lProbName;
         probName     = new char[lProbName + 1];
         strcpy(probName, received->probName);
      }
   }
   else
   {
      CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
      // receiver side creation
      if( localComm->getThreadId() == root )
      {
         for( int i = 0; i < localComm->getSize(); i++ )
         {
            if( i != root )
            {
               PARA_COMM_CALL(
                     localComm->uTypeSend((void *)createDatatype(comm), ParaInstanceType, i, UG::TagParaInstance)
               );
            }
         }
      }
      else
      {
         std::unique_ptr<CMapLapParaInstanceTh> received;
         PARA_COMM_CALL(
               localComm->uTypeReceive((void **)&received, ParaInstanceType, root, UG::TagParaInstance)
         );
         lProbName    = received->lProbName;
         probName     = new char[lProbName + 1];
         strcpy(probName, received->probName);
      }
   }
#else
   CMapLapParaSolverLocalComm *localComm = dynamic_cast<CMapLapParaSolverLocalComm *>(comm);
   // receiver side creation
   if( localComm->getThreadId() == root )
   {
      for( int i = 0; i < localComm->getSize(); i++ )
      {
         if( i != root )
         {
            PARA_COMM_CALL(
                  localComm->uTypeSend((void *)createDatatype(comm), ParaInstanceType, i, TagParaInstance)
            );
         }
      }
   }
   else
   {
      std::unique_ptr<CMapLapParaInstanceTh> received;
      PARA_COMM_CALL(
            localComm->uTypeReceive((void **)&received, ParaInstanceType, root, TagParaInstance)
      );
      lProbName    = received->lProbName;
      probName     = new char[lProbName + 1];
      strcpy(probName, received->probName);
   }
#endif
   return 0;

}

} // namespace ParaCMapLAP