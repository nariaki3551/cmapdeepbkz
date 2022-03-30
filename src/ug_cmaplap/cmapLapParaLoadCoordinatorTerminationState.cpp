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

/**@file    cmapLapParaLoadCoordinatorTerminationState.cpp
 * @brief   Load coordinator termination state.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "cmapLapParaLoadCoordinatorTerminationState.h"
#include <sstream>
#include <string>
#include "ug/paraComm.h"
#include "ug/gzstream.h"
#ifdef UG_WITH_ZLIB
#include "ug/gzstream.h"
#endif



namespace ParaCMapLAP
{

///
/// stringfy ParaCalculationState
///
std::string
CMapLapParaLoadCoordinatorTerminationState::toString(
      )
{
   std::ostringstream s;
   if( isCheckpointState )
      s << "######### LoadCoordinator Rank = " << rank << " is at checkpoint. #########" << std::endl;
   else
      s << "######### LoadCoordinator Rank = " << rank << " is terminated. #########" << std::endl;

   s << "#=== The number of ParaNodes received = " << nReceived << std::endl;
   s << "#=== The number of ParaNodes sent = " << nSent << std::endl;
   if( isCheckpointState )
   {
      s << "#=== Current global best dual bound value = " <<  externalGlobalBestDualBoundValue
        << "( internal value = " << globalBestDualBoundValue << " )" << std::endl;
   }

   s << "#=== Total Idle time to terminate this LoadCoordinator  = "
      << idleTime + idleTimeToWaitIsend + totalCopyCheckpointTime + totalWriteCheckpointTime << std::endl;
   s << "#=== Idle time to probe of this LoadCoordinator  = " << idleTime << std::endl;
   s << "#=== Idle time to wait Isend of this LoadCoordinator  = " << idleTimeToWaitIsend << std::endl;
   s << "#=== Total number of writing checkpoint files = " << nWriteCheckpointing << std::endl;
   s << "#=== Total elapsed time to write checkpoint files = " << totalWriteCheckpointTime << std::endl;
   s << "#=== Total number of coping checkpoint files = " << nCopyCheckpointing << std::endl;
   s << "#=== Total elapsed time to copy checkpoint files = " << totalCopyCheckpointTime << std::endl;
   return s.str();
}

#ifdef UG_WITH_ZLIB
void
CMapLapParaLoadCoordinatorTerminationState::write(
      gzstream::ogzstream &out
      )
{
   out.write(reinterpret_cast<char *>(&isCheckpointState), sizeof(bool));
   out.write(reinterpret_cast<char *>(&rank), sizeof(int));
   out.write(reinterpret_cast<char *>(&nWarmStart), sizeof(unsigned long long));
   out.write(reinterpret_cast<char *>(&nSent), sizeof(unsigned long long));
   out.write(reinterpret_cast<char *>(&nReceived), sizeof(unsigned long long));
   out.write(reinterpret_cast<char *>(&globalBestDualBoundValue), sizeof(double));
   out.write(reinterpret_cast<char *>(&externalGlobalBestDualBoundValue), sizeof(double));
   out.write(reinterpret_cast<char *>(&idleTime), sizeof(double));
   out.write(reinterpret_cast<char *>(&idleTimeToWaitIsend), sizeof(double));
   out.write(reinterpret_cast<char *>(&nHandlers), sizeof(int));
   out.write(reinterpret_cast<char *>(&runningTime), sizeof(double));
   out.write(reinterpret_cast<char *>(&nWriteCheckpointing), sizeof(int));
   out.write(reinterpret_cast<char *>(&totalWriteCheckpointTime), sizeof(double));
   out.write(reinterpret_cast<char *>(&nCopyCheckpointing), sizeof(int));
   out.write(reinterpret_cast<char *>(&totalCopyCheckpointTime), sizeof(double));
   // write process(Times/Calls)OfMessageHandler
   out.write(reinterpret_cast<char *>(processTimesOfMessageHandler.data()), nHandlers*sizeof(int));
   out.write(reinterpret_cast<char *>(processCallsOfMessageHandler.data()), nHandlers*sizeof(int));
}

bool
CMapLapParaLoadCoordinatorTerminationState::read(
      UG::ParaComm *comm,
      gzstream::igzstream &in
      )
{
   in.read(reinterpret_cast<char *>(&isCheckpointState), sizeof(bool));
   if( in.eof() ) return false;
   in.read(reinterpret_cast<char *>(&rank), sizeof(int));
   in.read(reinterpret_cast<char *>(&nWarmStart), sizeof(unsigned long long));
   in.read(reinterpret_cast<char *>(&nSent), sizeof(unsigned long long));
   in.read(reinterpret_cast<char *>(&nReceived), sizeof(unsigned long long));
   in.read(reinterpret_cast<char *>(&globalBestDualBoundValue), sizeof(double));
   in.read(reinterpret_cast<char *>(&externalGlobalBestDualBoundValue), sizeof(double));
   in.read(reinterpret_cast<char *>(&idleTime), sizeof(double));
   in.read(reinterpret_cast<char *>(&idleTimeToWaitIsend), sizeof(double));
   in.read(reinterpret_cast<char *>(&nHandlers), sizeof(int));
   in.read(reinterpret_cast<char *>(&runningTime), sizeof(double));
   in.read(reinterpret_cast<char *>(&nWriteCheckpointing), sizeof(int));
   in.read(reinterpret_cast<char *>(&totalWriteCheckpointTime), sizeof(double));
   in.read(reinterpret_cast<char *>(&nCopyCheckpointing), sizeof(int));
   in.read(reinterpret_cast<char *>(&totalCopyCheckpointTime), sizeof(double));
   // write process(Times/Calls)OfMessageHandler
   in.read(reinterpret_cast<char *>(processTimesOfMessageHandler.data()), nHandlers*sizeof(int));
   in.read(reinterpret_cast<char *>(processCallsOfMessageHandler.data()), nHandlers*sizeof(int));
   return true;
}


///
/// @return string of each process time of message handler
///
std::string
CMapLapParaLoadCoordinatorTerminationState::toStringProcessTimesOfMessageHandlerHeader(
         UG::ParaComm *comm,
         std::string delimiter
         )
{
   std::ostringstream sHead;

   sHead << "time" << delimiter <<  "idleTime" << delimiter << "idleTimeToWaitIsend";
   if( nHandlers > 0 )
   {
      for( int i = 0; i < nHandlers; i++ )
      {
         sHead << delimiter << comm->getTagString(i);
      }
      for( int i = 0; i < nHandlers; i++ )
      {
         sHead << delimiter << comm->getTagString(i) << "Call";
      }
   }
   sHead << delimiter << "numSendBackBasis";
   return sHead.str();
}


///
/// @return string of each process time of message handler
///
std::string
CMapLapParaLoadCoordinatorTerminationState::toStringProcessTimesOfMessageHandler(
         UG::ParaComm *comm,
         double elapsedTime,
         std::string delimiter
         )
{
   std::ostringstream sTime, sCall;

   sTime << elapsedTime << delimiter << idleTime << delimiter << idleTimeToWaitIsend;
   if( nHandlers > 0 )
   {
      for( int i = 0; i < nHandlers; i++ )
      {
         sTime << delimiter << processTimesOfMessageHandler[i];
         sCall << delimiter << processCallsOfMessageHandler[i];
      }
   }
   sCall << delimiter << nSendBackBasis;
   return sTime.str() + sCall.str();
}


///
/// @return string of each process time of message handler
///
std::string
CMapLapParaLoadCoordinatorTerminationState::toStringLocalProcessTimesOfMessageHandlerHeader(
         UG::ParaComm *comm,
         std::string delimiter
         )
{
   std::ostringstream sHead;

   sHead << "idleTime" << delimiter << "idleTimeToWaitIsend";
   if( nHandlers > 0 )
   {
      for( int i = 0; i < nHandlers; i++ )
      {
         sHead << delimiter << comm->getTagString(i);
      }
   }
   return sHead.str();
}


///
/// @return string of each process time of message handler
///
std::string
CMapLapParaLoadCoordinatorTerminationState::toStringLocalProcessTimesOfMessageHandler(
         UG::ParaComm *comm,
         std::string delimiter
         )
{
   std::ostringstream sTime;

   sTime << idleTime << delimiter << idleTimeToWaitIsend;
   if( nHandlers > 0 )
   {
      for( int i = 0; i < nHandlers; i++ )
      {
         sTime << delimiter << localProcessTimesOfMessageHandler[i];
      }
   }
   return sTime.str();
}
#endif

} // namespace ParaCMapLAP