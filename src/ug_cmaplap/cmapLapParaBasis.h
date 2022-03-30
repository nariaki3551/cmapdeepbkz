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

/**@file    cmapLapParaBasis.h
 * @brief   Base class for communicating basis matrix.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_BASIS_H__
#define __CMAP_LAP_PARA_BASIS_H__

#include <sstream>
#include <memory>
#include "cmapLapParaDef.h"

#if UG_WITH_ZLIB
#include "ug/gzstream.h"
#endif
namespace UG { class ParaComm; }


namespace ParaCMapLAP
{

class CMapLapParaBasisTh;
class CMapLapParaBasisMpi;


///
/// CMapLapParaBasis class
///
class CMapLapParaBasis : public std::enable_shared_from_this<CMapLapParaBasis>
{

protected:

   int threadId;              ///< thread Id
   double enumCost;           ///< enumeration cost of basis
   LatticeBasis<int> basis;   ///< lattice basis

public:

   friend CMapLapParaBasisTh;
   friend CMapLapParaBasisMpi;

   ///
   /// constructor
   ///
   CMapLapParaBasis(
         )
         : threadId(-1)
   {
   }

   CMapLapParaBasis(
         int inThreadId,            ///< thread Id
         double inEnumCost,         ///< enumeration cost
         LatticeBasis<int> &inBasis ///< lattice basis
         ) : threadId(inThreadId),
             enumCost(inEnumCost),
             basis(inBasis)
   {
   }

   ///
   /// destructor
   ///
   virtual ~CMapLapParaBasis(
         )
   {
   }

   ///
   /// @brief getter of threadId
   /// @return thread id
   ///
   virtual int getThreadId(
         )
   {
      return threadId;
   }

   ///
   /// @brief setter of threadId
   /// @param[in] inThreadId thread id
   ///
   virtual void setThreadId(
         int inThreadId
         )
   {
      threadId = inThreadId;
   }

   ///
   /// @brief getter of basis
   ///
   virtual LatticeBasis<int> & getBasis(
         )
   {
      return basis;
   }

   ///
   /// @brief getter of basis
   ///
   virtual double getEnumCost(
         )
   {
      return enumCost;
   }


   virtual void send(
         UG::ParaComm *comm,
         int destination) = 0;

   virtual void receive(
         UG::ParaComm *comm,
         int source) = 0;


#ifdef UG_WITH_ZLIB
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
      s << "threadId "  << threadId       << std::endl;
      s << "nCols "     << basis.cols()  << std::endl;
      s << "nRows"      << basis.rows()  << std::endl;
      return s.str();
   }

};

} // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_BASIS_H__
