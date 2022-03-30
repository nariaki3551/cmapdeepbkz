/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* MIT License                                                                     */
/*                                                                                 */
/* Copyright (c) 2022 Nariaki Tateiwa <n-tateiwa@kyudai.jp>                        */
/*                                                                                 */
/* Permission is hereby granted, free of charge, to any person obtaining a copy    */
/* of this software and associated documentation files (the "Software"), to deal   */
/* in the Software without restriction, including without limitation the rights    */
/* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell       */
/* copies of the Software, and to permit persons to whom the Software is           */
/* furnished to do so, subject to the following conditions:                        */
/*                                                                                 */
/* The above copyright notice and this permission notice shall be included in all  */
/* copies or substantial portions of the Software.                                 */
/*                                                                                 */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR      */
/* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        */
/* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     */
/* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          */
/* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,   */
/* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE   */
/* SOFTWARE.                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    SubEnumeration.h
 * @brief   Base class for SubEnumeration.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __LAPTOOLS_SUB_ENUE_H__
#define __LAPTOOLS_SUB_ENUE_H__

#include <float.h>                              // for DBL_MAX
#include <cassert>                              // for assert
#include <iostream>                             // for ofstream, string
#include <memory>                               // for shared_ptr

#include "Def.h"
#include "Lattice.h"
#include "Log.h"
#include "Config.h"



namespace LapTools
{


///
/// @class SubEnumeration
///
/// @tparam BasisFloat floating point of basis matrix
/// @tparam GSFloat floating point of Gram-Schmidt matrix
/// @tparam EnumGSFloat floating point of Gram-Schmidt matrix for Enumeration
///
template<typename BasisFloat=int, typename GSFloat=double, typename EnumGSFloat=double>
class SubEnumeration
{


using LatticePtr = std::shared_ptr<Lattice<BasisFloat, GSFloat>>;


protected:

   LatticePtr L;                 ///< lattice
   Config   config;              ///< hyper parameters
   int      rank;                ///< solver-id
   int      threadId;            ///< solver-id
   int      verbose;             ///< verbose <= 0: none, 1: light, 2: medium, 3: heavy


private:

   double   timeLimit;           ///< time limit
   int      begin;               ///< begin index of basis for reduction
   int      end;                 ///< end index of basis for reduction
   int      projDim;             ///< the projected index ( begin + projDim <= end )

   int      nChildEnum;          ///< the number of finding vector in child enumeration
   double   bestObjectiveValue;  ///< squared norm of begin-th GSO vectors
   double   runningTime;         ///< time of just running algorthm

   std::shared_ptr<std::ofstream>   osCsvLog;   ///< ostream for csv log


public:

   ///
   /// @brief Constructor
   ///
   SubEnumeration(){}


   ///
   /// @brief Constructor
   /// @param[in] inL         Lattice
   /// @param[in] inRank      solver-id
   /// @param[in] inThreadId  solver-id
   /// @param[in] inVerbose   <= 0 : not, 1 : light, 2: heavy
   ///
   SubEnumeration(
         LatticePtr inL,
         int inRank=1,
         int inThreadId=0,
         int inVerbose=0
         )
      :
         rank(inRank),
         threadId(inThreadId),
         verbose(inVerbose),
         nChildEnum(0),
         bestObjectiveValue(DBL_MAX),
         runningTime(0.0),
         osCsvLog(nullptr)
   {
      L = inL;
      config = L->config;
      begin = 0;
      end = L->m-1;
   }


   ///
   /// @brief replace lattice object
   /// @param[in] inL   lattice
   ///
   void resetLattice(
         LatticePtr inL
         )
   {
      L = inL;
      config = L->config;
   }


   ///
   /// @brief getter of Lattice
   /// @return Lattice pointer
   ///
   LatticePtr getLattice(
         )
   {
      return L;
   }


   ///
   /// @brief set parameters
   /// @param[in] inBegin     begin index of enumeration
   /// @param[in] inEnd       end index of enumeration
   /// @param[in] inProjDim   projected dimension of subenum
   ///
   void init(
         int inBegin=0,
         int inEnd=-1,
         int inProjDim=0
       )
   {
      begin = inBegin;
      end = inEnd;
      if( end < 0 ) end = L->m - 1;
      projDim = inProjDim;
      assert( begin + projDim <= end );
      nChildEnum = 0;
   }


   ///
   /// @brief sub-Enumeration
   /// @param[out] resv          vector found by Enum
   /// @parem[out] coeffResv     coefficient of vector found by Enum
   /// @param[in] inTimeLimit    time limit
   ///
   bool subEnum(
         LatticeVector<BasisFloat> &resv,
         LatticeVector<BasisFloat> &coeffResv,
         double inTimeLimit=-1
         );


   ///
   /// @brief update bestObjectiveValue
   /// @param[in] sigh line-header character
   ///
   virtual void updateBestObjectiveValue(
         char sigh,
         LatticeVector<BasisFloat> &v
         );


   ///
   /// @brief createt log
   /// @param[in] sigh line-header character
   ///
   void outputLog(
      char sigh=' '
      );


   ///
   /// @brief setter of osCsvLog
   /// @param[in] inOsCsvLog
   ///
   void setOsCsvLog(
         std::shared_ptr<std::ofstream> inOsCsvLog
         )
   {
      osCsvLog = inOsCsvLog;
   }


   ///
   /// @brief get header of csv formed log
   ///
   std::string getCsvLogHeader(
      )
   {
      return Log::getCsvLogHeader();
   }


   ///
   /// @brief set verbose
   /// @param[in] inVerbose verbose level
   ///
   void setVerbose(
         int inVerbose
         )
   {
      verbose = inVerbose;
   }


   ///
   /// @brief get verbose
   /// @return inVerbose verbose level
   ///
   int getVerbose(
         )
   {
      return verbose;
   }


}; // class SubEnumeration


}  // namespace LapTools


#endif // __LAPTOOLS_SUB_ENUE_H__
