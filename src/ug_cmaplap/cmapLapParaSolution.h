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

/**@file    cmapLapParaSolution.h
 * @brief   ParaSolution extension for CMapLap solver.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_SOLUTION_H__
#define __CMAP_LAP_PARA_SOLUTION_H__

#include <cmath>
#include <sstream>
#include <memory>
#include "ug/paraSolution.h"
#include "cmapLapParaDef.h"
#ifdef UG_WITH_ZLIB
#include "ug/gzstream.h"
#endif
namespace UG { class ParaComm; }
namespace ParaCMapLAP { class CMapLapParaSolutionMpi; }
namespace ParaCMapLAP { class CMapLapParaSolutionTh; }


namespace ParaCMapLAP
{


///
/// CMapLapParaSolution class
///
class CMapLapParaSolution : public UG::ParaSolution, public std::enable_shared_from_this<CMapLapParaSolution>
{
protected:

   int threadId;                 ///< thread Id
   LatticeVector<int> v;         ///< lattice vector
   double projectedSquaredNorm;  ///< projected squared norm

public:

   friend CMapLapParaSolutionTh;
   friend CMapLapParaSolutionMpi;

   ///
   /// constructor
   ///
   CMapLapParaSolution(
         )
      :
         threadId(-1),
         projectedSquaredNorm(0.0)
   {
   }

   CMapLapParaSolution(
         int inThreadId,               ///< thread Id
         LatticeVector<int> &inV,      ///< lattice vector
         double inProjectedSquaredNorm ///< objective function value
         )
      :
         threadId(inThreadId),
         v(inV),
         projectedSquaredNorm(inProjectedSquaredNorm)
   {
   }

   ///
   /// destructor
   ///
   virtual ~CMapLapParaSolution(
         )
   {
   }

   /// get original objective function value
   /// double getOrgObjectiveFunctionValue(){ return norm; }

   ///
   /// get objective function value
   ///
   virtual double getObjectiveFunctionValue()
   {
      return projectedSquaredNorm;
   }

   virtual int getThreadId(
         )
   {
      return threadId;
   }

   virtual int setThreadId(
         int inThreadId
         )
   {
      threadId = inThreadId;
      return 0;
   }

   ///
   /// getter of v
   /// @return shortest vector
   ///
   virtual LatticeVector<int>& getVector(
         )
   {
      return v;
   }

   ///
   /// set objective function value
   ///
   virtual void setObjectiveFunctionValue(
         double val
         )
   {
      projectedSquaredNorm = val;
   }

#ifdef UG_WITH_ZLIB

   ///
   /// write checkpoint file
   ///
   virtual void writeCheckpointSolution(
         gzstream::ogzstream &out
         );

   ///
   /// user should implement write method
   ///
   virtual void write(
         gzstream::ogzstream &out
         );

   ///
   /// user should implement read method
   ///
   virtual bool read(
         UG::ParaComm *comm,
         gzstream::igzstream &in
         );

#endif

   ///
   /// stringfy solution
   ///
   virtual const std::string toString(
         )
   {
      std::ostringstream s;
      s << "threadId "  << threadId << std::endl;
      s << "obj = "     << projectedSquaredNorm << std::endl;
      s << "vector = "  << v.transpose();
      return s.str();
   }

};

} // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_SOLUTION_H__
