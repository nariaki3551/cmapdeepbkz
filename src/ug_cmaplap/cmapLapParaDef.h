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

/**@file    cmapLapParaDef.h
 * @brief   Defines for CMAP-LAP
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_DEF_H__
#define __CMAP_LAP_PARA_DEF_H__

#include<cstdlib>
#include <eigen3/Eigen/Core>

namespace ParaCMapLAP
{

#define CMAP_LAP_CALL_ABORT(x) do                                         \
{                                                                        \
   int  _status_;                                                        \
   if( (_status_ = (x)) )                                                \
   {                                                                     \
      std::cout << "[CMAP_LAP CALL ERROR:" <<  __FILE__ << "] func = "    \
      << __func__ << ", line = " << __LINE__ << " - code = <"            \
      <<  _status_ << ">" << std::endl;                                  \
      abort();                                                           \
   }                                                                     \
}                                                                        \
while( 0 )

#define CMAP_LAP_INFINITY 1e+20  // values larger than this are considered infinity


template<typename BasisFloat=int>
using LatticeBasis = Eigen::Matrix<BasisFloat, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
template<typename BasisFloat=int>
using LatticeVector = Eigen::Matrix<BasisFloat, Eigen::Dynamic, 1>;

static const double infinity = 1e+75;
static const int ON_MEMORY_TRANSFER = 0;

enum SolverType {
   DeepBkz,
   Enum,
   Sieve,
   Undefined
};

} // namespace ParaCMapLAP

#endif
