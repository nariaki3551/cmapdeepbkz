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

/**@file    cmapLapParaIsendRequest.h
 * @brief   CMapLapISend cmapLapRequest data structure
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#ifndef __CMAP_LAP_PARA_ISEND_REQUEST_H__
#define __CMAP_LAP_PARA_ISEND_REQUEST_H__
#include <memory>
#include "ug/paraIsendRequest.h"
#include "cmapLapParaTaskMpi.h"
#include "cmapLapParaSolutionMpi.h"
#include "cmapLapParaBasisMpi.h"
#include "cmapLapParaPackedVectorMpi.h"
#include "cmapLapParaSolverStateMpi.h"
#include "cmapLapParaCalculationStateMpi.h"

namespace ParaCMapLAP
{

///
/// @class CMapLapParaIsendRequest
/// @brief Allocated memory for sending object
///
class CMapLapParaIsendRequest : public UG::ParaIsendRequest
{

private:

   std::shared_ptr<MPI_Request> cmapLapReq;
   std::shared_ptr<CMapLapParaCalculationState> cmapLapParaCalculationState;
   std::shared_ptr<CMapLapParaPackedVector> cmapLapParaPackedVector;
   std::shared_ptr<CMapLapParaSolution> cmapLapParaSolution;
   std::shared_ptr<CMapLapParaBasis> cmapLapParaBasis;
   std::shared_ptr<CMapLapParaSolverState> cmapLapParaSolverState;
   std::shared_ptr<CMapLapParaTask> cmapLapParaTask;

public:

   ///
   /// Constructor
   ///
   CMapLapParaIsendRequest(
         std::shared_ptr<MPI_Request>& inReq
         ) : UG::ParaIsendRequest(),
             cmapLapReq(inReq),
             cmapLapParaCalculationState(nullptr),
             cmapLapParaPackedVector(nullptr),
             cmapLapParaSolution(nullptr),
             cmapLapParaBasis(nullptr),
             cmapLapParaSolverState(nullptr),
             cmapLapParaTask(nullptr)
   {
   }
   CMapLapParaIsendRequest(
         std::shared_ptr<MPI_Request> inReq,
         std::shared_ptr<CMapLapParaCalculationState> inObjectPointer
         ) : CMapLapParaIsendRequest(inReq)
   {
      cmapLapParaCalculationState = inObjectPointer;
   }
   CMapLapParaIsendRequest(
         std::shared_ptr<MPI_Request> inReq,
         std::shared_ptr<CMapLapParaPackedVector> inObjectPointer
         ) : CMapLapParaIsendRequest(inReq)
   {
      cmapLapParaPackedVector = inObjectPointer;
   }
   CMapLapParaIsendRequest(
         std::shared_ptr<MPI_Request> inReq,
         std::shared_ptr<CMapLapParaSolution> inObjectPointer
         ) : CMapLapParaIsendRequest(inReq)
   {
      cmapLapParaSolution = inObjectPointer;
   }
   CMapLapParaIsendRequest(
         std::shared_ptr<MPI_Request> inReq,
         std::shared_ptr<CMapLapParaBasis> inObjectPointer
         ) : CMapLapParaIsendRequest(inReq)
   {
      cmapLapParaBasis = inObjectPointer;
   }
   CMapLapParaIsendRequest(
         std::shared_ptr<MPI_Request> inReq,
         std::shared_ptr<CMapLapParaSolverState> inObjectPointer
         ) : CMapLapParaIsendRequest(inReq)
   {
      cmapLapParaSolverState = inObjectPointer;
   }
   CMapLapParaIsendRequest(
         std::shared_ptr<MPI_Request> inReq,
         std::shared_ptr<CMapLapParaTask> inObjectPointer
         ) : CMapLapParaIsendRequest(inReq)
   {
      cmapLapParaTask = inObjectPointer;
   }
   
   ///
   /// deconstructor
   /// delete the object after sending, however, only paraTask is not deleted
   /// because the object is saved in Pool after sending.
   ///
   virtual ~CMapLapParaIsendRequest(
         )
   {
   }


   ///
   /// tester if cmapLapObjectPointer object has already been sent
   /// @return true if it has already been sent else false
   ///
   virtual bool test(
         )
   {
      assert( cmapLapReq );
      int flag = 0;
      MPI_CALL(
         MPI_Test(cmapLapReq.get(), &flag, MPI_STATUS_IGNORE)
      );
      return flag;
   }


   ///
   /// wait until all cmapLapObjectPointer objects compulete to send
   ///
   virtual void wait(
         )
   {
      assert( cmapLapReq );
      MPI_CALL(
         MPI_Wait(cmapLapReq.get(), MPI_STATUS_IGNORE)
      );
   }
};

} // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_ISEND_REQUEST_H__


