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

/**@file    DeepBkz.h
 * @brief   Base class for DeepBKZ.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __LAPTOOLS_DEEP_BKZ_H__
#define __LAPTOOLS_DEEP_BKZ_H__

#include <iostream>
#include <memory>
#include <float.h>
#include <cstddef>
#include <new>
#include "Def.h"
#include "Config.h"
#include "Lattice.h"
#include "Reduction.h"
#include "DeepLll.h"
#include "Enumeration.h"
#include "Utils.h"



namespace LapTools
{


///
/// @class DeepBkz
///
/// @tparam BasisFloat floating point of basis matrix
/// @tparam GSFloat floating point of Gram-Schmidt matrix
/// @tparam EnumGSFloat floating point of Gram-Schmidt matrix for Enumeration
///
template<typename BasisFloat=int, typename GSFloat=double, typename EnumGSFloat=double>
class DeepBkz : public Reduction<BasisFloat, GSFloat>
{

using LatticePtr = std::shared_ptr<Lattice<BasisFloat, GSFloat>>;


private:

   LatticePtr  L;                   ///< lattice
   Config      config;              ///< hyper parameters
   int         rank;                ///< rank
   int         threadId;            ///< thread Id
   int         verbose;    ///< verbose <= 0: none, 1: light, 2: medium, 3: heavy

   GSFloat                 T__;        ///< allocated memory for GSO update
   VectorMu<GSFloat>       nu__;       ///< allocated memory for GSO update
   VectorMu<GSFloat>       D__;        ///< allocated memory for GSO update
   VectorMu<GSFloat>       S__;        ///< allocated memory for GSO update
   LatticeBasis<BasisFloat> prevBasis; ///< basis of previous tour

   std::shared_ptr<std::ofstream>   osCsvLog;   ///< ostream for csv log


protected:

   int         blocksize;              ///< blocksize of BKZ
   double      timeLimit;              ///< time limit
   double      lowerBound;             ///< lowerbound of lattice vector norm

   DeepLll<BasisFloat, GSFloat>     lllObj;  ///< lll tool
   Enumeration<BasisFloat, GSFloat, EnumGSFloat> enumObj;   ///< enumeration tool

   int         nTour;                  ///< number of tours
   double      bestObjectiveValue;     ///< squared norm of the top of GSO vectors
   LatticeVector<BasisFloat> v;        ///< allocated memory for Enum
   LatticeVector<BasisFloat> coeffv;   ///< allocated memory for Enum

   int         nSlopeNoDecrease;       ///< number of times the GSA slope did not decrease in the recent past
   double      prevSlopeGSA;           ///< previous GSA Slope

   CountArray<size_t> hashCounter;     ///< counter of basis hash
   int         hashCounterSize;        ///< size of counter of basis hash

   double      runningTime;            ///< time of just running algorthm
   double      mergeTime;              ///< time of merging other lattice


public:

   bool        hasReduced;          ///< true if this lattice is already reduced with current reduction parameters
   bool        autoAborted;         ///< true if it is satisfied the comdition of auto abort;

   ///
   /// @brief Constructor
   ///
   DeepBkz(){}


   ///
   /// @brief Constructor
   /// @param[in] inL         Lattice
   /// @param[in] inVerbose   <= 0 : not, 1 : light, 2: heavy
   /// @param[in] inRank      solver-id
   /// @param[in] inThreadId  solver-id
   ///
   DeepBkz(
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
         T__(0.0),
         osCsvLog(nullptr),
         timeLimit(-1),
         lowerBound(-1),
         nTour(0),
         bestObjectiveValue(DBL_MAX),
         nSlopeNoDecrease(0),
         prevSlopeGSA(0),
         hashCounterSize(20),
         runningTime(0.0),
         mergeTime(0.0),
         hasReduced(false),
         autoAborted(false)
   {
      L = inL;
      config = inL->config;
      lllObj = DeepLll<BasisFloat, GSFloat>(inL, inRank, inThreadId);
      enumObj = Enumeration<BasisFloat, GSFloat, EnumGSFloat>(inL, inRank, inThreadId);
      v.resize(inL->n);
      coeffv.resize(inL->n);
      nu__.resize(inL->m);
      D__.resize(inL->m);
      S__.resize(inL->m);
      prevBasis = inL->basis;
      hashCounter.resize(hashCounterSize);
   }


   ///
   /// @brief deconstructor
   ///
   virtual ~DeepBkz(){}


   ///
   /// @brief replace lattice object
   /// @param[in] inL   lattice
   ///
   virtual void resetLattice(
         LatticePtr inL
         )
   {
      Reduction<BasisFloat, GSFloat>::resetLattice(inL);
      prevBasis = this->getLattice()->basis;
      lllObj.resetLattice(inL);
      enumObj.resetLattice(inL);
   }


   ///
   /// @brief return idle status
   /// @return true if it has reduced lattice basis or satisfies auto-abort condition else false
   ///
   virtual bool idle(
       )
   {
      return hasReduced || autoAborted;
   }


   ///
   /// @brief DeepBkz reduction algorithm
   ///        for Matrix[begin:end] and algorithm start at "start index";
   /// @param[in] inBlocksize block size;
   /// @param[in] inTimeLimit
   /// @return bool true: normal terminate, false: abnormal one
   /// @remark If lowerBound is positive value, algorithm terminates wheh it founds the lattice vector whose norm is less than lowerBound
   ///
   virtual bool deepbkz(
         int inBlocksize,
         double inTimeLimit=-1
         );


   ///
   /// @brief preprocess of DeepBkz main loop
   /// @param[out] shouldAbort true if it should abort else false
   ///
   virtual bool preprocess(
         bool &shouldAbort
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
      return true;
   }


   ///
   /// @brief communicate with LC
   /// @param[out] shouldAbort true if it should abort else false
   /// @note this is a dummy function for sequential execution
   ///
   virtual bool communicateInTour(
         bool &shouldAbort
         )
   {
      return true;
   }


   ///
   /// @brief run tour
   /// @param[out] shouldAbort true if it should abort else false
   /// @remark runningTime is updated in this function
   /// @remark autoAborted is updaetd in this function
   /// @remark If lowerBound is positive value, wheh it founds the lattice vector whose norm is less than lowerBound
   ///
   virtual bool tour(
         bool &shouldAbort
         );


   ///
   /// @brief one loop of algorithm
   /// @param[in] k
   /// @param[out] z
   /// @param[out] shouldAbort true if it should abort else false
   /// @return bool true: normal terminate, false: abnormal one
   ///
   virtual bool step(
         int k,
         int &z,
         bool &shouldAbort
         );


   ///
   /// @breif MLLL process after Enum found vector v during DeepBkz, and
   ///        v is inserted in the basis
   /// @param[in] k  begin index of block
   /// @param[in] l  end index of block
   /// @return
   ///
   virtual bool postProcess(
         int k,
         int l
         );


   ///
   /// @brief update mu and sqnorm for basis matrix
   ///        {b_0,...,b_k-1,v,b_k,...,b_l, b_l+1}
   /// @param[in] k index of basis
   /// @param[in] r tha last index whose coefficient is not zero
   ///
   virtual bool GSOupdateBkz(
         int k,
         int r
         );


   ///
   /// @brief DeepLll reduction algorithm
   /// @return bool true: normal terminate, false: abnormal one
   ///
   virtual bool deeplll(
         )
   {
      int lllVerbose = lllObj.getVerbose();
      lllObj.setVerbose(verbose);
      bool state = lllObj.deeplll();
      lllObj.setVerbose(lllVerbose);
      return state;
   }


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
      char sigh=' ',
      std::ofstream *os=nullptr
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
   /// @brief setter of lowerBound
   /// @param[in] inLowerBound
   ///
   virtual void setLowerBound(
         double inLowerBound
         )
   {
      lowerBound = inLowerBound;
   }


   ///
   /// @brief setter of lowerBound
   /// @param[in] inLowerBound
   ///
   virtual double getLowerBound(
         double inLowerBound
         )
   {
      return lowerBound;
   }


   ///
   /// @brief get header of csv formed log
   ///
   virtual std::string getCsvLogHeader(
      );


   ///
   /// @brief setter of blocksize
   ///
   virtual void setBlocksize(
         int inBlocksize
         )
   {
      blocksize = inBlocksize;
   }


   ///
   /// @brief getter of running time
   ///
   virtual double getRunningTime(
         )
   {
      return runningTime;
   };

   ///
   /// @brief getter of number of tour
   ///
   virtual double getNTour(
         )
   {
      return nTour;
   }


}; // class DeepBkz


}  // namespace LapTools


#endif // __LAPTOOLS_DEEP_BKZ_H__
