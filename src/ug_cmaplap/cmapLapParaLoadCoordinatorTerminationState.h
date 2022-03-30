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

/**@file    cmapLapParaLoadCoordinatorTerminationState.h
 * @brief   Load coordinator termination state.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_LOADCOORDINATOR_TERMINATION_STATE_H__
#define __CMAP_LAP_PARA_LOADCOORDINATOR_TERMINATION_STATE_H__

#include <cfloat>
#include <vector>
#include <numeric>
#include <cassert>
#include "ug/paraLoadCoordinatorTerminationState.h"
#ifdef UG_WITH_ZLIB
#include "ug/gzstream.h"
#endif
namespace UG { class ParaComm; }

namespace ParaCMapLAP
{

///
/// Class for LoadCoordinator termination state
/// which contains calculation state in a ParaLoadCoordinator
///
class CMapLapParaLoadCoordinatorTerminationState : public UG::ParaLoadCoordinatorTerminationState
{
public:

   // bool                  isCheckpointState;                   ///< indicate if this state is at checkpoint or not
   // int                   rank;                                ///< rank of this ParaLoadCoordinator
   ///
   /// Counters related to this ParaLoadCoordinator
   ///
   /// Define your status here

   ///
   ///  current dual bound value
   ///
   double                globalBestDualBoundValue;             ///< global best dual bound value (internal value)
   double                externalGlobalBestDualBoundValue;     ///< global best dual bound value (external value)
   ///
   ///  times of this LoadCoordinator
   ///
   double                idleTime;                             ///< idle time of this LoadCoordinator
   double                idleTimeToWaitIsend;                  ///< idle time to wait Isend messages
   std::vector<double>   processTimesOfMessageHandler;         ///< total time to proccess message handler functions
   std::vector<int>      processCallsOfMessageHandler;         ///< number of calls of message handler functions
   std::vector<double>   localProcessTimesOfMessageHandler;    ///< total time to proccess message handler functions in latest term
   std::vector<int>      localProcessCallsOfMessageHandler;    ///< number of calls of message handler functions in latest term
   int                   nHandlers;                            ///< number of handlers
   double                runningTime;                          ///< this ParaLoadCoordinator running time
   ///
   /// times of checkpoint
   ///
   int                   nWriteCheckpointing;                  ///< total number of writing checkpoint
   double                totalWriteCheckpointTime;             ///< total time used to write checkpoint
   int                   nCopyCheckpointing;                   ///< total number of coping checkpoint
   double                totalCopyCheckpointTime;              ///< total time used to copy checkpoint object
   ///
   /// other counter
   int                   nSendBackBasis;                       ///< number of send-back basis


   ///
   /// default constructor
   ///
   CMapLapParaLoadCoordinatorTerminationState(
         )
         : UG::ParaLoadCoordinatorTerminationState(),
           globalBestDualBoundValue(-DBL_MAX),
           externalGlobalBestDualBoundValue(-DBL_MAX),
           idleTime(0.0),
           idleTimeToWaitIsend(0.0),
           nHandlers(0),
           runningTime(0.0),
           nWriteCheckpointing(0),
           totalWriteCheckpointTime(0.0),
           nCopyCheckpointing(0),
           totalCopyCheckpointTime(0.0),
           nSendBackBasis(0)
   {
   }

   ///
   /// stringfy ParaCalculationState
   /// @return string to show inside of this object
   ///
   virtual std::string toString(
         );

   ///
   /// set processTimesOfMessageHandler and processCallsOfMessageHandler
   ///
   virtual void setProcessTimesOfMessageHandler(
         int inNHandlers
         )
   {
      nHandlers = inNHandlers;
      processTimesOfMessageHandler.resize(nHandlers, 0.0);
      processCallsOfMessageHandler.resize(nHandlers, 0);
      localProcessTimesOfMessageHandler.resize(nHandlers, 0.0);
      localProcessCallsOfMessageHandler.resize(nHandlers, 0);
   }


   ///
   /// reset localProcessTimesOfMessageHandler and localProcessCallsOfMessageHandler
   ///
   virtual void resetLocalProcessTimesOfMessageHandler(
         )
   {
      assert( nHandlers > 0 );
      std::fill(localProcessTimesOfMessageHandler.begin(), localProcessTimesOfMessageHandler.end(), 0.0);
      std::fill(localProcessCallsOfMessageHandler.begin(), localProcessCallsOfMessageHandler.end(), 0);
   }

   ///
   /// getter of total time of localProcessTimesOfMessageHandler
   ///
   virtual double totalLocalProcessTimesOfMessageHandler(
         )
   {
      return std::accumulate(localProcessTimesOfMessageHandler.begin(), localProcessTimesOfMessageHandler.end(), 0.0);
   }


   ///
   /// add process time into totalTimeToProcessMessageHandlers
   ///
   virtual void addProcessTime(
         int tag,
         double processTime
         )
   {
      processTimesOfMessageHandler[tag] += processTime;
      processCallsOfMessageHandler[tag]++;
      localProcessTimesOfMessageHandler[tag] += processTime;
      localProcessCallsOfMessageHandler[tag]++;
   }


#ifdef UG_WITH_ZLIB

   ///
   /// write to checkpoint file
   ///
   virtual void write(
         gzstream::ogzstream &out              ///< gzstream for output
         );

   ///
   /// read from checkpoint file
   ///
   virtual bool read(
         UG::ParaComm *comm,                       ///< communicator used
         gzstream::igzstream &in               ///< gzstream for input
         );

#endif

   ///
   /// @return string of each process time of message handler
   ///
   virtual std::string toStringProcessTimesOfMessageHandler(
         UG::ParaComm *comm,
         double elapsedTime,
         std::string delimiter=","
         );

   ///
   /// @return string of each process time of message handler
   ///
   virtual std::string toStringProcessTimesOfMessageHandlerHeader(
         UG::ParaComm *comm,
         std::string delimiter=","
         );

   ///
   /// @return string of each process time of message handler
   ///
   virtual std::string toStringLocalProcessTimesOfMessageHandlerHeader(
         UG::ParaComm *comm,
         std::string delimiter
         );

   ///
   /// @return string of each process time of message handler
   ///
   virtual std::string toStringLocalProcessTimesOfMessageHandler(
         UG::ParaComm *comm,
         std::string delimiter
         );

};

} // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_LOADCOORDINATOR_TERMINATION_STATE_H__
