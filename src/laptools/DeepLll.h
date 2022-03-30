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

/**@file    DeepLll.h
 * @brief   Base class for DeepLLL.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __LAPTOOLS_DEEP_LLL_H__
#define __LAPTOOLS_DEEP_LLL_H__

#include <memory>
#include <float.h>
#include <iostream>
#include "Def.h"
#include "Config.h"
#include "Lattice.h"
#include "Log.h"
#include "Reduction.h"


namespace LapTools
{


///
/// @class DeepLll
///
/// @tparam BasisFloat floating point of basis matrix
/// @tparam GSFloat floating point of Gram-Schmidt matrix
///
template<typename BasisFloat=int, typename GSFloat=double>
class DeepLll : public Reduction<BasisFloat, GSFloat>
{

using LatticePtr = std::shared_ptr<Lattice<BasisFloat, GSFloat>>;


private:

   LatticePtr  L;                   ///< lattice
   Config      config;              ///< hyper parameters
   int         rank;                ///< rank
   int         threadId;            ///< thread Id
   int         verbose;             ///< verbose <= 0: none, 1: light, 2: medium, 3: heavy

   long int    nIter;               ///< number of iterations
   double      bestObjectiveValue;  ///< squared norm of 1-th GSO vectors
   double      runningTime;         ///< time of just running algorthm

   GSFloat           T__;        ///< allocated memory for GSO update
   VectorMu<GSFloat> P__;        ///< allocated memory for GSO update
   VectorMu<GSFloat> D__;        ///< allocated memory for GSO update
   VectorMu<GSFloat> S__;        ///< allocated memory for GSO update

   std::shared_ptr<std::ofstream>   osCsvLog;   ///< ostream for csv log


public:

   ///
   /// @brief Constructor
   ///
   DeepLll(){}


   ///
   /// @brief Constructor
   /// @param[in] inL         Lattice
   /// @param[in] inVerbose   <= 0 : not, 1 : light, 2: heavy
   /// @param[in] inRank      solver-id
   /// @param[in] inThreadId  solver-id
   ///
   DeepLll(
         LatticePtr inL,
         int inRank=1,
         int inThreadId=0,
         int inVerbose=0
         )
      :
         Reduction<BasisFloat, GSFloat>(inL, inRank, inThreadId, inVerbose),
         rank(inRank),
         threadId(inThreadId),
         verbose(inVerbose),
         nIter(0),
         bestObjectiveValue(DBL_MAX),
         runningTime(0.0),
         osCsvLog(nullptr)
   {
      L = inL;
      config = inL->config;
      P__.resize(inL->m);
      D__.resize(inL->m);
      S__.resize(inL->m);
   }


   ///
   /// @brief deconstructor
   ///
   virtual ~DeepLll(){}


   ///
   /// @brief DeepLll reduction algorithm
   ///        for Matrix[begin:end] and algorithm start at "start index";
   /// @param[in] start index of basis
   /// @param[in] begin index of basis
   /// @param[in] end   index of basis
   /// @return bool true: normal terminate, false: abnormal one
   ///
   virtual bool deeplll(
         int start=0,
         int begin=0,
         int end=-1
         );


   ///
   /// @brief communicate with LC
   /// @param[out] shouldAbort true if it should abort else false
   /// @note this is a dummy function for sequential execution
   ///
   virtual bool communicate(
         bool &shouldAbort
         )
   {
      return true;   // do nothing
   }


   ///
   /// @brief update GSO matrix with deep-insertion
   /// @param[in] i
   /// @param[in] k
   /// @details Update GSO information after executing DeepLll from row i to row k
   ///          {...,b_i,...,b_k,...} -> {...,b_i-1,b_k,b_i,...,b_k-1,b_k+1,...}
   /// @return true if B(i) >= B'(i) else false
   ///
   virtual bool GSOUpdateDeepLll(
         int i,
         int k
         );


   ///
   /// @brief update bestObjectiveValue
   /// @param[in] sigh line-header character
   ///
   virtual bool updateBestObjectiveValue(
         char sigh
         );


   ///
   /// @brief createt log
   /// @param[in] sigh line-header character
   ///
   virtual void outputLog(
      char sigh=' '
      );


   ///
   /// @brief setter of osCsvLog
   /// @param[in] inOsCsvLog
   ///
   virtual void setOsCsvLog(
         std::shared_ptr<std::ofstream> inOsCsvLog
         )
   {
      osCsvLog = inOsCsvLog;
   }


   ///
   /// @brief get header of csv formed log
   ///
   virtual std::string getCsvLogHeader(
      )
   {
      return Log::getCsvLogHeader();
   }


}; // class DeepLll


}  // namespace LapTools


#endif // __LAPTOOLS_DEEP_LLL_H__
