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

/**@file    cmapLapParaParamSet.h
 * @brief   Parameter set for CMAP-LAP.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_PARAM_SET_H__
#define __CMAP_LAP_PARA_PARAM_SET_H__
#if defined(_COMM_PTH) || defined(_COMM_CPP11)
#include "ug/paraParamSetTh.h"
#endif
#if defined(_COMM_MPI_WORLD)
#include "ug/paraParamSetMpi.h"
#endif
#include <stddef.h>
namespace UG { class ParaComm; }

#define OUTPUT_PARAM_VALUE_ERROR( msg1, msg2, msg3, msg4 ) \
   std::cout << "[PARAM VALUE ERROR] Param type = " << msg1 << ", Param name = " << msg2 \
     << ", Param value = " <<  msg3 <<  ": Param comment is as follows: " << std::endl \
     << msg4 << std::endl;  \
   return (-1)

namespace ParaCMapLAP
{

///
///  Bool parameters
///
static const int CMapLapParaParamsFirst               = UG::ParaParamsLast + 1;
static const int CMapLapParaParamsBoolFirst           = CMapLapParaParamsFirst;
//-------------------------------------------------------------------------
static const int WarmStartOnlyPool                    = CMapLapParaParamsBoolFirst + 0;
static const int NoWaitNotificationId                 = CMapLapParaParamsBoolFirst + 1;
static const int EnumPruningGH                        = CMapLapParaParamsBoolFirst + 2;
static const int EnumPruningLowerBound                = CMapLapParaParamsBoolFirst + 3;
static const int ShareIncumbentVector                 = CMapLapParaParamsBoolFirst + 4;
static const int LogShareDataPoolAll                  = CMapLapParaParamsBoolFirst + 5;
static const int LogShareDataPoolStat                 = CMapLapParaParamsBoolFirst + 6;
static const int CheckpointThreading                  = CMapLapParaParamsBoolFirst + 7;
static const int CheckpointReserving                  = CMapLapParaParamsBoolFirst + 8;
static const int OutputSimilarityOfBasis              = CMapLapParaParamsBoolFirst + 9;
static const int DynamicDimensionOfSharedLattice      = CMapLapParaParamsBoolFirst + 10;
static const int AutoAdjustmentNotificationInterval   = CMapLapParaParamsBoolFirst + 11;
//-------------------------------------------------------------------------
static const int CMapLapParaParamsBoolLast            = CMapLapParaParamsBoolFirst + 11;
static const int CMapLapParaParamsBoolN               = CMapLapParaParamsBoolLast - CMapLapParaParamsBoolFirst + 1;
///
/// Int parameters
///
static const int CMapLapParaParamsIntFirst            = CMapLapParaParamsBoolLast + 1;
//-------------------------------------------------------------------------
static const int MaxSizeOfMessageQueue                = CMapLapParaParamsIntFirst + 0;
static const int DeepBkzVerbose                       = CMapLapParaParamsIntFirst + 1;
static const int EnumVerbose                          = CMapLapParaParamsIntFirst + 2;
static const int SieveVerbose                         = CMapLapParaParamsIntFirst + 3;
static const int RandomizeRows                        = CMapLapParaParamsIntFirst + 4;
static const int DeepBkzReductionType                 = CMapLapParaParamsIntFirst + 5;
static const int DeepBkzStartBlockSize                = CMapLapParaParamsIntFirst + 6;
static const int DeepBkzEndBlockSize                  = CMapLapParaParamsIntFirst + 7;
static const int DeepBkzBlockSizeInterval             = CMapLapParaParamsIntFirst + 8;
static const int DeepBkzRecursiveMaxTour              = CMapLapParaParamsIntFirst + 9;
static const int DeepBkzRecursiveLowerBeta            = CMapLapParaParamsIntFirst + 10;
static const int DeepBkzNumOfSendVectorsToPool        = CMapLapParaParamsIntFirst + 11;
static const int DeepBkzNumOfReceiveVectorsFromPool   = CMapLapParaParamsIntFirst + 12;
static const int EnumPruningLowerBlockSize            = CMapLapParaParamsIntFirst + 13;
static const int EnumNumOfSendVectorsToPool           = CMapLapParaParamsIntFirst + 14;
static const int EnumSamplingDepth                    = CMapLapParaParamsIntFirst + 15;
static const int EnumLowerNumberOfDividedSearchTree   = CMapLapParaParamsIntFirst + 16;
static const int SubEnumProjectedDimension            = CMapLapParaParamsIntFirst + 17;
static const int SubEnumNSendVectorsToPool            = CMapLapParaParamsIntFirst + 18;
static const int SieveMaxCollision                    = CMapLapParaParamsIntFirst + 19;
static const int SieveMaxListSize                     = CMapLapParaParamsIntFirst + 20;
static const int SieveMaxStackSize                    = CMapLapParaParamsIntFirst + 21;
static const int SieveStartNVectors                   = CMapLapParaParamsIntFirst + 22;
static const int SieveNThreads                        = CMapLapParaParamsIntFirst + 23;
static const int SieveNumOfSendVectorsToPool          = CMapLapParaParamsIntFirst + 24;
static const int SieveNumOfReceiveVectorsFromPool     = CMapLapParaParamsIntFirst + 25;
static const int BlockSizeForLocalSolver              = CMapLapParaParamsIntFirst + 26;
static const int NumOfInitialDeepBkzSolvers           = CMapLapParaParamsIntFirst + 27;
static const int NumOfInitialEnumSolvers              = CMapLapParaParamsIntFirst + 28;
static const int NumOfInitialSieveSolvers             = CMapLapParaParamsIntFirst + 29;
static const int NumOfLocalSolversInLC                = CMapLapParaParamsIntFirst + 30;
static const int ShareDataPoolSize                    = CMapLapParaParamsIntFirst + 31;
static const int InstancePoolSize                     = CMapLapParaParamsIntFirst + 32;
static const int WriteSizeShareDataPool               = CMapLapParaParamsIntFirst + 33;
static const int LocalTaskNVectorLowerBound           = CMapLapParaParamsIntFirst + 34;
static const int LowerBoundOfInstancePoolShouldHold   = CMapLapParaParamsIntFirst + 35;
static const int ProjectedDimension                   = CMapLapParaParamsIntFirst + 36;
static const int DimensionOfSharedLattice             = CMapLapParaParamsIntFirst + 37;
static const int NumOfSamplingForSimilarityOfBasis    = CMapLapParaParamsIntFirst + 38;
//-------------------------------------------------------------------------
static const int CMapLapParaParamsIntLast             = CMapLapParaParamsIntFirst + 38;
static const int CMapLapParaParamsIntN                = CMapLapParaParamsIntLast - CMapLapParaParamsIntFirst + 1;
///
/// Longint parameters
///
static const int CMapLapParaParamsLongintFirst        = CMapLapParaParamsIntLast + 1;
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
static const int CMapLapParaParamsLongintLast         = CMapLapParaParamsLongintFirst - 1;  // No params -1
static const int CMapLapParaParamsLongintN            = CMapLapParaParamsLongintLast - CMapLapParaParamsLongintFirst + 1;
///
/// Real parameters
///
static const int CMapLapParaParamsRealFirst           = CMapLapParaParamsLongintLast + 1;
//-------------------------------------------------------------------------
static const int LowerBoundOfNorm                     = CMapLapParaParamsRealFirst   + 0;
static const int LowerBoundOfApproxFactor             = CMapLapParaParamsRealFirst   + 1;
static const int DeepBkzTaskMaxTimeLimit              = CMapLapParaParamsRealFirst   + 2;
static const int EnumTaskMaxTimeLimit                 = CMapLapParaParamsRealFirst   + 3;
static const int SieveTaskMaxTimeLimit                = CMapLapParaParamsRealFirst   + 4;
static const int EnumTotalProb                        = CMapLapParaParamsRealFirst   + 5;
static const int EnumTraversedNodesPerSeconds         = CMapLapParaParamsRealFirst   + 6;
static const int EnumPruningParameter                 = CMapLapParaParamsRealFirst   + 7;
static const int RandomizeScale                       = CMapLapParaParamsRealFirst   + 8;
static const int EnumStartTime                        = CMapLapParaParamsRealFirst   + 9;
static const int SieveStartTime                       = CMapLapParaParamsRealFirst   + 10;
static const int IntervalTimeOfAssignmentTableOutput  = CMapLapParaParamsRealFirst   + 11;
static const int IntervalTimeOfLogShareDataPool       = CMapLapParaParamsRealFirst   + 12;
static const int IReceiveInterval                     = CMapLapParaParamsRealFirst   + 13;
static const int ShareVectorsInterval                 = CMapLapParaParamsRealFirst   + 14;
static const int IntervalTimeOfOutputSimilarityOfBasis = CMapLapParaParamsRealFirst  + 15;
static const int LCTermTimeUpdateNotificationInterval = CMapLapParaParamsRealFirst   + 16;
static const int LCUpperIdleRatio                     = CMapLapParaParamsRealFirst   + 17;
static const int LCLowerIdleRatio                     = CMapLapParaParamsRealFirst   + 18;
//-------------------------------------------------------------------------
static const int CMapLapParaParamsRealLast            = CMapLapParaParamsRealFirst   + 18;
static const int CMapLapParaParamsRealN               = CMapLapParaParamsRealLast - CMapLapParaParamsRealFirst + 1;
///
/// Char parameters
///
static const int CMapLapParaParamsCharFirst           = CMapLapParaParamsRealLast + 1;
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
static const int CMapLapParaParamsCharLast            = CMapLapParaParamsCharFirst - 1;   // No params -1
static const int CMapLapParaParamsCharN                = CMapLapParaParamsCharLast - CMapLapParaParamsCharFirst + 1;
///
/// String parameters
///
static const int CMapLapParaParamsStringFirst         = CMapLapParaParamsCharLast    + 1;
//-------------------------------------------------------------------------
static const int CMapLapParamFilePath                 = CMapLapParaParamsStringFirst + 0;
//-------------------------------------------------------------------------
static const int CMapLapParaParamsStringLast          = CMapLapParaParamsStringFirst + 0;
static const int CMapLapParaParamsStringN             = CMapLapParaParamsStringLast - CMapLapParaParamsStringFirst + 1;
static const int CMapLapParaParamsLast                = CMapLapParaParamsStringLast;
static const int CMapLapParaParamsSize                = CMapLapParaParamsLast + 1;


///
/// class CMapLapParaParamSet
///
#if defined(_COMM_PTH) || defined(_COMM_CPP11)
class CMapLapParaParamSet : public UG::ParaParamSetTh
#endif
#if defined(_COMM_MPI_WORLD)
class CMapLapParaParamSet : public UG::ParaParamSetMpi
#endif
{

public:

   ///
   /// constructor
   ///
   CMapLapParaParamSet(
         );

   ///
   /// destructor
   ///
   virtual ~CMapLapParaParamSet(
         )
   {
   }

   ///
   /// read ParaParams from file
   ///
   void read(
         UG::ParaComm *comm,  ///< communicator used
         const char* filename ///< reading file name
         );

   ///
   /// get number of bool parameters
   /// @return size of parameter table
   ///
   virtual size_t getNumBoolParams(
         )
   {
      return (UG::ParaParamsBoolN + CMapLapParaParamsBoolN);
   }

   ///
   /// get number of int parameters
   /// @return size of parameter table
   ///
   virtual size_t getNumIntParams(
         )
   {
      return (UG::ParaParamsIntN + CMapLapParaParamsIntN);
   }

   ///
   /// get number of longint parameters
   /// @return size of parameter table
   ///
   virtual size_t getNumLongintParams(
         )
   {
      return (UG::ParaParamsLongintN + CMapLapParaParamsLongintN);
   }

   ///
   /// get number of real parameters
   /// @return size of parameter table
   ///
   virtual size_t getNumRealParams(
         )
   {
      return (UG::ParaParamsRealN + CMapLapParaParamsRealN);
   }

   ///
   /// get number of char parameters
   /// @return size of parameter table
   ///
   virtual size_t getNumCharParams(
         )
   {
      return (UG::ParaParamsCharN + CMapLapParaParamsCharN);
   }

   ///
   /// get number of string parameters
   /// @return size of parameter table
   ///
   virtual size_t getNumStringParams(
         )
   {
      return (UG::ParaParamsStringN + CMapLapParaParamsStringN);
   }

};

} // namespace ParaCMapLAP

#endif  // __CMAP_LAP_PARA_PARAM_SET_H__
