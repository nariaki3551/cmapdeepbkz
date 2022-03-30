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

/**@file    cmapLapParaCheckpointWriter.h
 * @brief   Base class for CheckpointWriter.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_CHECKPOINT_WRITER_H__
#define __CMAP_LAP_PARA_CHECKPOINT_WRITER_H__

#include <memory>
#include "ug/paraTimer.h"
namespace ParaCMapLAP { class CMapLapParaSolverLocalComm; }
namespace ParaCMapLAP { class CheckpointElement; }
namespace UG { class ParaParamSet; }


namespace ParaCMapLAP
{

class CMapLapParaCheckpointWriter
{

   CMapLapParaSolverLocalComm *localComm;       ///< communicator
   std::unique_ptr<UG::ParaTimer> paraTimer;    ///< timer for this Checkpoint Writer
   UG::ParaParamSet *paraParams;                ///< ParaParamSet object

   ///
   /// write checkpoint files
   ///
   virtual void updateCheckpointFiles(
         CheckpointElement *checkpointElement
         );

public:
   ///
   /// constructor
   ///
   CMapLapParaCheckpointWriter(
         CMapLapParaSolverLocalComm *localComm,
         UG::ParaParamSet *paraParams
         );

   ///
   /// deconstructor
   ///
   virtual ~CMapLapParaCheckpointWriter(
         )
   {
   }

   ///
   /// run this Checkpoint Writer
   ///
   virtual void run(
         );

};

} // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_CHECKPOINT_WRITER_H__
