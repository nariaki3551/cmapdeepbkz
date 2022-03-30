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

/**@file    cmapLapParaSolution.cpp
 * @brief   ParaSolution extension for CMAP-LAP solver.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "cmapLapParaSolution.h"
#include <vector>
#include "cmapLapParaDef.h"
#ifdef UG_WITH_ZLIB
#include "ug/gzstream.h"
#endif
namespace UG { class ParaComm; }


namespace ParaCMapLAP
{

#ifdef UG_WITH_ZLIB
void
CMapLapParaSolution::writeCheckpointSolution(
      gzstream::ogzstream &out
      )
{
   write(out);
}

void
CMapLapParaSolution::write(
      gzstream::ogzstream &out
      )
{
   out.write(reinterpret_cast<char *>(&projectedSquaredNorm), sizeof(double));
   int n = v.size();
   out.write(reinterpret_cast<char *>(&n), sizeof(int));

   std::vector<int> vecArray(n);
   Eigen::Map<Eigen::VectorXi>(&vecArray[0], n) = v;
   out.write(reinterpret_cast<char *>(&vecArray[0]), sizeof(int)*n);
}

bool
CMapLapParaSolution::read(
      UG::ParaComm *comm,
      gzstream::igzstream &in
      )
{
   in.read(reinterpret_cast<char *>(&projectedSquaredNorm), sizeof(double));
   if( in.eof() ) return false;
   int n = 0;
   in.read(reinterpret_cast<char *>(&n), sizeof(int));

   std::vector<int> vecArray(n);
   in.read(reinterpret_cast<char *>(&vecArray[0]), sizeof(int)*n);
   v = Eigen::Map<Eigen::VectorXi>(&vecArray[0], n);
   return true;
}


} // namespace ParaCMapLAP

#endif
