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

/**@file    cmapLapParaBasis.cpp
 * @brief   Base class for communicating basis matrix.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "cmapLapParaBasis.h"
#include <vector>

namespace ParaCMapLAP
{


#ifdef UG_WITH_ZLIB
void
CMapLapParaBasis::write(
      gzstream::ogzstream &out
      )
{
   out.write(reinterpret_cast<char *>(&enumCost), sizeof(double));
   int nRows = basis.rows();
   int nCols = basis.cols();
   out.write(reinterpret_cast<char *>(&nRows), sizeof(int));
   out.write(reinterpret_cast<char *>(&nCols), sizeof(int));

   std::vector<int> basisArray(nRows*nCols);
   Eigen::Map<LatticeBasis<int>>(&basisArray[0], nRows, nCols) = basis;
   out.write(reinterpret_cast<char *>(&basisArray[0]), sizeof(int)*(nRows*nCols));
}

bool
CMapLapParaBasis::read(
      UG::ParaComm *comm,
      gzstream::igzstream &in
      )
{
   in.read(reinterpret_cast<char *>(&enumCost), sizeof(double));
   if( in.eof() ) return false;
   int nRows = 0;
   int nCols = 0;
   in.read(reinterpret_cast<char *>(&nRows), sizeof(int));
   in.read(reinterpret_cast<char *>(&nCols), sizeof(int));

   std::vector<int> basisArray(nRows*nCols);
   in.read(reinterpret_cast<char *>(&basisArray[0]), sizeof(int)*(nRows*nCols));
   basis = Eigen::Map<LatticeBasis<int>>(&basisArray[0], nRows, nCols);
   return true;
}

#endif


} // namespace ParaCMapLAP