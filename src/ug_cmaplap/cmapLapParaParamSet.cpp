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

/**@file    cmapLapParaParamSet.cpp
 * @brief   Parameter set for CMAP-LAP.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <cfloat>
#include <climits>
#include "cmapLapParaParamSet.h"
namespace UG { class ParaComm; }


namespace ParaCMapLAP
{


CMapLapParaParamSet::CMapLapParaParamSet(
      )
#if defined(_COMM_PTH) || defined(_COMM_CPP11)
      : UG::ParaParamSetTh(CMapLapParaParamsSize)
#else
      : UG::ParaParamSetMpi(CMapLapParaParamsSize)
#endif
{
   ///
   /// bool params
   ///
   paraParams[WarmStartOnlyPool] = new UG::ParaParamBool(
         "WarmStartOnlyPool",
         "# WarmStartOnlyPool is TRUE: Load Coordnator reads only Basis and Vector pool when it uses warm start",
         false);
   paraParams[NoWaitNotificationId] = new UG::ParaParamBool(
         "NoWaitNotificationId",
         "# NoWaitNotificationId is TRUE: Solvers do not wait NotificationId from LC (LC doesnot send recieve tag of NoWaitNotificationId)",
         false);
   paraParams[EnumPruningGH] = new UG::ParaParamBool(
         "EnumPruningGH",
         "# EnumPruningGH is TRUE: pruning enumeration tree by GH, FALSE: no pruning enumeration tree by GH [Default value: FALSE]",
         false);
   paraParams[EnumPruningLowerBound] = new UG::ParaParamBool(
         "EnumPruningLowerBound",
         "# EnumPruningLowerBound is TRUE: pruning enumeration tree by Lower Bound, FALSE: no pruning enumeration tree by Lower Bound [Default value: FALSE]",
         false);
   paraParams[ShareIncumbentVector] = new UG::ParaParamBool(
         "ShareIncumbentVector",
         "# ShareIncumbentVector is TRUE: Each solver receives the incumbent solution, FALSE: Each solver does not receive the incumbent solution [Default value: TRUE]",
         true);
   paraParams[LogShareDataPoolAll] = new UG::ParaParamBool(
         "LogShareDataPoolAll",
         "# LogShareDataPoolAll is TRUE: write norm of vectors in vector pool [Default value: FALSE]",
         false);
   paraParams[LogShareDataPoolStat] = new UG::ParaParamBool(
         "LogShareDataPoolStat",
         "# LogShareDataPoolStat is TRUE: write statistics infomation of vectors in vector pool [Default value: FALSE]",
         false);
   paraParams[CheckpointThreading] = new UG::ParaParamBool(
         "CheckpointThreading",
         "# CheckpointThreading is TRUE: create therad for writing checkpoint [Default value: FALSE]",
         true);
   paraParams[CheckpointReserving] = new UG::ParaParamBool(
         "CheckpointReserving",
         "# CheckpointReserving is TRUE: copy reserved checkpoint objectes. Therefore when this is True, memory usage is increase [Default value: TRUE]",
         true);
   paraParams[OutputSimilarityOfBasis] = new UG::ParaParamBool(
         "OutputSimilarityOfBasis",
         "# OutputSimilarityOfBasis is TRUE: output similarity of basis whose each solver [Default value: FALSE]",
         false);
   paraParams[DynamicDimensionOfSharedLattice] = new UG::ParaParamBool(
         "DynamicDimensionOfSharedLattice",
         "# DynamicDimensionOfSharedLattice is TRUE: dimension of shared lattice will dynamically change in the session [Default value: FALSE]",
         false);
   paraParams[AutoAdjustmentNotificationInterval] = new UG::ParaParamBool(
         "AutoAdjustmentNotificationInterval",
         "# AutoAdjustmentNotificationInterval is TRUE: While adjusting the solver's notification interval so that the idle ratio of LC does not exceed LCUpperIdleRatio, allow it to fall below NotificationInterval. [Default value: FALSE]",
         false);


   ///
   /// int params
   ///
   paraParams[MaxSizeOfMessageQueue] = new UG::ParaParamInt(
         "MaxSizeOfMessageQueue",
         "# When the size of the message queue is greater than MaxSizeOfMessageQueue, Solver doesnot send any message other than paraSolution. If this parameter is -1, then this parameter is ignored. [Default value: 100][-1, INT_MAX]",
         100,
         -1,
         INT_MAX);
   paraParams[DeepBkzVerbose] = new UG::ParaParamInt(
         "DeepBkzVerbose",
         "# DeepBkzVerbose verbose level: [Default value: 0][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[EnumVerbose] = new UG::ParaParamInt(
         "EnumVerbose",
         "# EnumVerbose verbose level: [Default value: 0][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[SieveVerbose] = new UG::ParaParamInt(
         "SieveVerbose",
         "# SieveVerbose verbose level: [Default value: 0][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[RandomizeRows] = new UG::ParaParamInt(
         "RandomizeRows",
         "# Randomize the RandomizeRows rows from the bottom of basis: -1 - randomize overall of basis 0 no randomize[Default value -1][-1, INT_MAX]",
         -1,
         -1,
         INT_MAX);
   paraParams[DeepBkzReductionType] = new UG::ParaParamInt(
         "DeepBkzReductionType",
         "# Reduction type of DeepBkz; 1 - no-communication DeepBkz, 2 - parallel DeepBkz [Default value: 2][1, 2]",
         2,
         1,
         2);
   paraParams[DeepBkzStartBlockSize] = new UG::ParaParamInt(
         "DeepBkzStartBlockSize",
         "# DeepBkz reduction with blocksize from DeepBkzStartBlockSize to DeepBkzEndBlockSize at intervals of DeepBkzBlockSizeInterval after randomize [Default value: 30][2, INT_MAX]",
         30,
         2,
         INT_MAX);
   paraParams[DeepBkzEndBlockSize] = new UG::ParaParamInt(
         "DeepBkzEndBlockSize",
         "# DeepBkz reduction with blocksize from DeepBkzStartBlockSize to DeepBkzEndBlockSize at intervals of DeepBkzBlockSizeInterval after randomize [Default value: 30][2, INT_MAX]",
         30,
         2,
         INT_MAX);
   paraParams[DeepBkzBlockSizeInterval] = new UG::ParaParamInt(
         "DeepBkzBlockSizeInterval",
         "# DeepBkz reduction with blocksize from DeepBkzStartBlockSize to DeepBkzEndBlockSize at intervals of DeepBkzBlockSizeInterval after randomize [Default value: 30][1, INT_MAX]",
         5,
         1,
         INT_MAX);
   paraParams[DeepBkzRecursiveMaxTour] = new UG::ParaParamInt(
         "DeepBkzRecursiveMaxTour",
         "# Max number of tour of Recursive DeepBkz whose depth is greather than 0 [Default value: 200]",
         200,
         1,
         INT_MAX);
   paraParams[DeepBkzRecursiveLowerBeta] = new UG::ParaParamInt(
         "DeepBkzRecursiveLowerBeta",
         "# Performs a recursive DeepBkz subroutine when the Enum successes in DeepBkz tour and BlockSize is a greater than or equal to DeepBkzLowerBeta [Default value: 20][2, INT_MAX]",
         20,
         2,
         INT_MAX);
   paraParams[DeepBkzNumOfSendVectorsToPool] = new UG::ParaParamInt(
         "DeepBkzNumOfSendVectorsToPool",
         "# Number of the vectors sended in DeepBkz [Default value: 1][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[DeepBkzNumOfReceiveVectorsFromPool] = new UG::ParaParamInt(
         "DeepBkzNumOfReceiveVectorsFromPool",
         "# Number of the vectors received in DeepBkz [Default value: 10][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[EnumPruningLowerBlockSize] = new UG::ParaParamInt(
         "EnumPruningLowerBlockSize",
         "# Pruning the enumeration tree when a blocksize of Enum if greater than or equal to EnumPruningLowerBlockSize [Default value: 45][1, INT_MAX]",
         45,
         1,
         INT_MAX);
   paraParams[EnumNumOfSendVectorsToPool] = new UG::ParaParamInt(
         "EnumNumOfSendVectorsToPool",
         "# Number of sending vectors during Enumeration to pool of Load Coordnator  [Default value: 1][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[EnumSamplingDepth] = new UG::ParaParamInt(
         "EnumSamplingDepth",
         "# Sampling depth during Enumeration [Default value: 2][0, INT_MAX]",
         2,
         0,
         INT_MAX);
   paraParams[EnumLowerNumberOfDividedSearchTree] = new UG::ParaParamInt(
         "EnumLowerNumberOfDividedSearchTree",
         "# Lower number of the divided enumeration trees. [Default value: 1000][1, INT_MAX]",
         1000,
         1,
         INT_MAX);
   paraParams[SubEnumProjectedDimension] = new UG::ParaParamInt(
         "SubEnumProjectedDimension",
         "# Projected dimension of SubEnum, if it is 0 then run norml Enum [Default value: 0][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[SubEnumNSendVectorsToPool] = new UG::ParaParamInt(
         "SubEnumNSendVectorsToPool",
         "# Number of the vectors sended in SubEnum [Default value: 0][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[SieveMaxCollision] = new UG::ParaParamInt(
         "SieveMaxCollision",
         "# Sieve terminates the number of collision exceeds SieveMaxCollision, if SieveMaxCollision is -1, it is ignore [Default value: -1][-1, INT_MAX]",
         -1,
         -1,
         INT_MAX);
   paraParams[SieveMaxListSize] = new UG::ParaParamInt(
         "SieveMaxListSize",
         "# Maximum List L size, if SieveMaxListSize is -1, it is ignore [Default value: 10000000][-1, INT_MAX]",
         10000000,
         -1,
         INT_MAX);
   paraParams[SieveMaxStackSize] = new UG::ParaParamInt(
         "SieveMaxStackSize",
         "# Maximum Stack S size, if SieveMaxStackSize is -1, it is ignore [Default value: 10000000][-1, INT_MAX]",
         10000000,
         -1,
         INT_MAX);
   paraParams[SieveStartNVectors] = new UG::ParaParamInt(
         "SieveStartNVectors",
         "# Sieve algorithm start the number of vectors in ShareDataPool is greater than SieveStartNVectors [Default value: 100][-1, INT_MAX]",
         100,
         -1,
         INT_MAX);
   paraParams[SieveNThreads] = new UG::ParaParamInt(
         "SieveNThreads",
         "# Maximum number of threads of a Sieve solver, if value is 0, SieveNThreads is set to number of threads per rank in the machine [Default value: 32][0, INT_MAX]",
         32,
         0,
         INT_MAX);
   paraParams[SieveNumOfSendVectorsToPool] = new UG::ParaParamInt(
         "SieveNumOfSendVectorsToPool",
         "# Number of the vectors sended in Sieve [Default value: 10][0, INT_MAX]",
         10,
         0,
         INT_MAX);
   paraParams[SieveNumOfReceiveVectorsFromPool] = new UG::ParaParamInt(
         "SieveNumOfReceiveVectorsFromPool",
         "# Number of the vectors received in Sieve [Default value: 100][0, INT_MAX]",
         100,
         0,
         INT_MAX);
   paraParams[BlockSizeForLocalSolver] = new UG::ParaParamInt(
         "BlockSizeForLocalSolver",
         "# blocksize of DeepBkz for Local Solver [Default value: 10][1, 1000]",
         10,
         1,
         1000);
   paraParams[NumOfInitialDeepBkzSolvers] = new UG::ParaParamInt(
         "NumOfInitialDeepBkzSolvers",
         "# Number of DeepBkz solvers at the beginning of the run, if value is -1, NumOfInitialDeepBkzSolvers is set to number of solvers [Default value: -1][-1, INT_MAX]",
         -1,
         -1,
         INT_MAX);
   paraParams[NumOfInitialEnumSolvers] = new UG::ParaParamInt(
         "NumOfInitialEnumSolvers",
         "# Number of Enum solvers used  [Default value: 0][0, INT_MAX]",
         0,
         -1,
         INT_MAX);
   paraParams[NumOfInitialSieveSolvers] = new UG::ParaParamInt(
         "NumOfInitialSieveSolvers",
         "# Number of Sieve solvers used  [Default value: 0][0, INT_MAX]",
         0,
         -1,
         INT_MAX);
   paraParams[NumOfLocalSolversInLC] = new UG::ParaParamInt(
         "NumOfLocalSolversInLC",
         "# Number of local solvers in LoadCoordinator [Default value: 0][0, 999999999]",
         0,
         0,
         999999999);
   paraParams[ShareDataPoolSize] = new UG::ParaParamInt(
         "ShareDataPoolSize",
         "# Limit number of vectors stored in vector pool  [Default value: 100][0, INT_MAX]",
         100000,
         0,
         INT_MAX);
   paraParams[InstancePoolSize] = new UG::ParaParamInt(
         "InstancePoolSize",
         "# Limit number of basis stored in basis pool  [Default value: 50][1, INT_MAX]",
         50,
         1,
         INT_MAX);
   paraParams[WriteSizeShareDataPool] = new UG::ParaParamInt(
         "WriteSizeShareDataPool",
         "# maximum number of vector pool vectors that are written in every checkpoint  [Default value: 1000][1, INT_MAX]",
         1000,
         1,
         INT_MAX);
   paraParams[LocalTaskNVectorLowerBound] = new UG::ParaParamInt(
         "LocalTaskNVectorLowerBound",
         "# if number of vectors in vector pool is smaller than it, local task is not created [Default value: 10][2, INT_MAX]",
         10,
         2,
         INT_MAX);
   paraParams[LowerBoundOfInstancePoolShouldHold] = new UG::ParaParamInt(
         "LowerBoundOfInstancePoolShouldHold",
         "# Lower bound of number of basis that basis element pool tries to hold [Default value: 0][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[ProjectedDimension] = new UG::ParaParamInt(
         "ProjectedDimension",
         "# projected dimension of this solver [Default value: 0][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[DimensionOfSharedLattice] = new UG::ParaParamInt(
         "DimensionOfSharedLattice",
         "# dimension of always shared lattice [Default value: 0][0, INT_MAX]",
         0,
         0,
         INT_MAX);
   paraParams[NumOfSamplingForSimilarityOfBasis] = new UG::ParaParamInt(
         "NumOfSamplingForSimilarityOfBasis",
         "# number of sampling for calculation of basis similarity [Default value: 100][2, INT_MAX]",
         100,
         2,
         INT_MAX);

   ///
   /// longint params
   ///

   ///
   /// real params
   ///
   paraParams[LowerBoundOfNorm] = new UG::ParaParamReal(
         "LowerBoundOfNorm",
         "# Solver terminates when it finds the vector whose norm shorter than LowerBoundOfNorm [Default value: -1][0, DBL_MAX]",
         -1,
         0.0,
         DBL_MAX);
   paraParams[LowerBoundOfApproxFactor] = new UG::ParaParamReal(
         "LowerBoundOfApproxFactor",
         "# Solver terminates when it finds the vector whose approx factor shorter than LowerBoundOfApproxFactor [Default value: -1][0, DBL_MAX]",
         -1,
         0.0,
         DBL_MAX);
   paraParams[DeepBkzTaskMaxTimeLimit] = new UG::ParaParamReal(
         "DeepBkzTaskMaxTimeLimit",
         "# Max time limit of a DeepBkz para task. [Default value: DBL_MAX][0.0, DBL_MAX]",
         DBL_MAX,
         0.0,
         DBL_MAX);
   paraParams[EnumTaskMaxTimeLimit] = new UG::ParaParamReal(
         "EnumTaskMaxTimeLimit",
         "# Max time limit of a Enum para task. [Default value: DBL_MAX][0.0, DBL_MAX]",
         DBL_MAX,
         0.0,
         DBL_MAX);
   paraParams[SieveTaskMaxTimeLimit] = new UG::ParaParamReal(
         "SieveTaskMaxTimeLimit",
         "# Max time limit of a Sieve para task. [Default value: DBL_MAX][0.0, DBL_MAX]",
         DBL_MAX,
         0.0,
         DBL_MAX);
   paraParams[EnumTotalProb] = new UG::ParaParamReal(
         "EnumTotalProb",
         "# Total probability of finding the shortest vector in our solver. [Default value: 0.95][0, 1]",
         0.95,
         0.0,
         1.0);
   paraParams[EnumTraversedNodesPerSeconds] = new UG::ParaParamReal(
         "EnumTraversedNodesPerSeconds",
         "# of nodes in an enumeration tree that can be traversed per second. [Default value: 3.3554432e7][1.0, DBL_MAX]",
         3.3554432e7,
         1.0,
         DBL_MAX);
   paraParams[EnumPruningParameter] = new UG::ParaParamReal(
         "EnumPruningParameter",
         "# pruning parameter of Enum Task. [Default value: 1.0][-1.0, 1.0]",
         1.0,
         -1.0,
         1.0);
   paraParams[RandomizeScale] = new UG::ParaParamReal(
         "RandomizeScale",
         "# When we generate randomize matrix, we multiplied RandomizeScale for each elements of matrix. [Default value: 1.2][1.0, DBL_MAX]",
         1.2,
         1.0,
         DBL_MAX);
   paraParams[EnumStartTime] = new UG::ParaParamReal(
         "EnumStartTime",
         "# The start time of Enum solver. Negative value means that no use of Enum solvers [Default value: -1.0][-1.0, DBL_MAX]",
         -1.0,
         -1.0,
         DBL_MAX);
   paraParams[SieveStartTime] = new UG::ParaParamReal(
         "SieveStartTime",
         "# The start time of Enum solver. Negative value means that no use of Enum solvers [Default value: -1.0][-1.0, DBL_MAX]",
         -1.0,
         -1.0,
         DBL_MAX);
   paraParams[IntervalTimeOfAssignmentTableOutput] = new UG::ParaParamReal(
         "IntervalTimeOfAssignmentTableOutput",
         "# Interval time of assignment table output [Default value: 5.0][0.0, DBL_MAX]",
         5.0,
         0.0,
         DBL_MAX);
   paraParams[IntervalTimeOfLogShareDataPool] = new UG::ParaParamReal(
         "IntervalTimeOfLogShareDataPool",
         "# Interval time to write log of vector pool [Default: 600.0][0.0, DBL_MAX]",
         600.0,
         0.0,
         DBL_MAX);
   paraParams[IReceiveInterval] = new UG::ParaParamReal(
         "IReceiveInterval",
         "# An active Solver calls iReceive() when time elapsed from its previous notification. [Default: 1.0][0.0, DBL_MAX]",
         1.0,
         0.0,
         DBL_MAX);
   paraParams[ShareVectorsInterval] = new UG::ParaParamReal(
         "ShareVectorsInterval",
         "# An active Solver sends vectors when time elapsed from its previous notification. [Default: 1.0][0.0, DBL_MAX]",
         1.0,
         0.0,
         DBL_MAX);
   paraParams[IntervalTimeOfOutputSimilarityOfBasis] = new UG::ParaParamReal(
         "IntervalTimeOfOutputSimilarityOfBasis",
         "# Interval time to output the similarity of basis whose each solver. [Default: 600.0][0.0, DBL_MAX]",
         600.0,
         0.0,
         DBL_MAX);
   paraParams[LCTermTimeUpdateNotificationInterval] = new UG::ParaParamReal(
         "LCTermTimeUpdateNotificationInterval",
         "# Interval time to update the solver's NotificationInterval. [Default: 600.0][0.0, DBL_MAX]",
         600.0,
         0.0,
         DBL_MAX);
   paraParams[LCUpperIdleRatio] = new UG::ParaParamReal(
         "LCUpperIdleRatio",
         "# Upperbound of ratio of LC's idle. [Default: 1.0][0.0, 1.0]",
         1.0,
         0.0,
         1.0);
   paraParams[LCLowerIdleRatio] = new UG::ParaParamReal(
         "LCLowerIdleRatio",
         "# Lowerbound of ratio of LC's idle. [Default: 0.3][0.0, 1.0]",
         0.3,
         0.0,
         1.0);



   ///
   /// char params
   ///

   ///
   /// string params
   ///
   paraParams[CMapLapParamFilePath] = new UG::ParaParamString(
         "CMapLapParamFilePath",
         "# CMapLap parameter settings file name. Empty name use default settings. [Default: ]",
         "");

}

void
CMapLapParaParamSet::read(
      UG::ParaComm *comm,
      const char* filename
      )
{

   ParaParamSet::read(comm, filename);

   ///
   /// check parameter consistency
   ///
   if( (getIntParamValue(NumOfInitialDeepBkzSolvers) < 0)
         + (getIntParamValue(NumOfInitialEnumSolvers) < 0)
         + (getIntParamValue(NumOfInitialSieveSolvers) < 0) >= 2 )
   {
      THROW_LOGICAL_ERROR7(
            "Must not be more than two negative values among NumOfInitialDeepBkzSolvers, NumOfInitialEnumSolvers, and NumOfInitialSieveSolvers",
            "; NumOfInitialDeepBkzSolvers = ",
            getIntParamValue(NumOfInitialDeepBkzSolvers),
            "; NumOfInitialEnumSolvers = ",
            getIntParamValue(NumOfInitialEnumSolvers),
            "; NumOfInitialSieveSolvers = ",
            getIntParamValue(NumOfInitialSieveSolvers)
            );
   }

   if( getBoolParamValue(UG::Quiet) )
   {
      setBoolParamValue(UG::TagTrace, false);
      setBoolParamValue(UG::LogSolvingStatus, false);
      setBoolParamValue(UG::LogTasksTransfer, false);
   }

}

} // namespace ParaCMapLAP
