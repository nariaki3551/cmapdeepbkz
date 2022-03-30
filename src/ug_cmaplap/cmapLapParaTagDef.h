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

/**@file    cmapLapParaTagDef.h
 * @brief   ug_cmaplap Tag definitions
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_TAG_DEF_H__
#define __CMAP_LAP_PARA_TAG_DEF_H__

#include "ug/paraTagDef.h"

namespace ParaCMapLAP
{

static const int TAG_CMAP_LAP_FIRST = UG::TAG_UG_LAST + 1;
//------------------------------------------------------------------------------------------------
// static const int TagCMapLapPackedVectorRequest          = TAG_CMAP_LAP_FIRST +  0;
static const int TagCMapLapPackedVector               = TAG_CMAP_LAP_FIRST +  0;
static const int TagBasisEnumCost                     = TAG_CMAP_LAP_FIRST +  1;
static const int TagTimeLimitRequest                  = TAG_CMAP_LAP_FIRST +  2;
static const int TagBasisRequest                      = TAG_CMAP_LAP_FIRST +  3;
static const int TagBasis                             = TAG_CMAP_LAP_FIRST +  4;
static const int TagVectorRequest                     = TAG_CMAP_LAP_FIRST +  5;
static const int TagUpdateNotificationInterval        = TAG_CMAP_LAP_FIRST +  6;
//-----------------------------------------------------------------------------------------------
static const int TAG_CMAP_LAP_BASE_LAST               = TAG_CMAP_LAP_FIRST +  6;
static const int N_CMAP_LAP_BASE_TAGS                 = TAG_CMAP_LAP_BASE_LAST - UG::TAG_UG_FIRST + 1;

#ifdef _COMM_MPI_WORLD

static const int TAG_CMAP_LAP_MPI_FIRST               = TAG_CMAP_LAP_BASE_LAST +  1;
//-----------------------------------------------------------------------------------------------
static const int TagParaInstance                      = TAG_CMAP_LAP_MPI_FIRST +  0;
static const int TagTask1                             = TAG_CMAP_LAP_MPI_FIRST +  1;
static const int TagTask2                             = TAG_CMAP_LAP_MPI_FIRST +  2;
static const int TagTask3                             = TAG_CMAP_LAP_MPI_FIRST +  3;
static const int TagSolution1                         = TAG_CMAP_LAP_MPI_FIRST +  4;
static const int TagSolverState1                      = TAG_CMAP_LAP_MPI_FIRST +  5;
static const int TagCompletionOfCalculation1          = TAG_CMAP_LAP_MPI_FIRST +  6;
static const int TagCMapLapPackedVector1              = TAG_CMAP_LAP_MPI_FIRST +  7;
static const int TagBasis1                            = TAG_CMAP_LAP_MPI_FIRST +  8;
//-----------------------------------------------------------------------------------------------
static const int TAG_CMAP_LAP_MPI_LAST                = TAG_CMAP_LAP_MPI_FIRST +  8;
static const int N_CMAP_LAP_MPI_TAGS                  = TAG_CMAP_LAP_MPI_LAST - UG::TAG_UG_FIRST + 1;
//-----------------------------------------------------------------------------------------------
static const int TAG_CMAP_LAP_LAST                    = TAG_CMAP_LAP_MPI_LAST;
static const int N_CMAP_LAP_TAGS                      = TAG_CMAP_LAP_LAST - UG::TAG_UG_FIRST + 1;

#endif

#if defined(_COMM_PTH) || defined (_COMM_CPP11)

static const int TAG_CMAP_LAP_TH_FIRST                  = TAG_CMAP_LAP_BASE_LAST + 1;
//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------
static const int TAG_CMAP_LAP_TH_LAST                   = TAG_CMAP_LAP_TH_FIRST - 1;     //  no tag
static const int N_CMAP_LAP_TH_TAGS                     = TAG_CMAP_LAP_TH_LAST - UG::TAG_UG_FIRST + 1;
//-----------------------------------------------------------------------------------------------
static const int TAG_CMAP_LAP_LAST                      = TAG_CMAP_LAP_TH_LAST;
static const int N_CMAP_LAP_TAGS                        = TAG_CMAP_LAP_LAST - UG::TAG_UG_FIRST + 1;

#endif

} // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_TAG_DEF_H__
