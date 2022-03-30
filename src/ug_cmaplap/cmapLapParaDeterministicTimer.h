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

/**@file    cmapLapParaDeterministicTimer.h
 * @brief   ParaDeterministicTimer extension for CMAP_LAP.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_DETERMINISITC_TIMER_H__
#define __CMAP_LAP_PARA_DETERMINISITC_TIMER_H__

#include "ug/paraDeterministicTimer.h"

namespace ParaCMapLAP
{

class CMapLapParaDeterministicTimer : public UG::ParaDeterministicTimer
{
   double current;
   int    normalizeFactor;
public:
   CMapLapParaDeterministicTimer() : current(0.0), normalizeFactor(1) {}
   virtual ~CMapLapParaDeterministicTimer() {}
   /**********************************************
    * if you want to set original initial time,  *
    * you can do it init()                       *
    **********************************************/
   virtual void normalize(UG::ParaComm *comm){ normalizeFactor = comm->getSize() - 1; }
   virtual void update(double value) { current += value; }
   virtual double getElapsedTime() { return current/normalizeFactor; }
};

} // namespace ParaCMapLAP

#endif  // __CMAP_LAP_PARA_DETERMINISTIC_TIMER_H__
