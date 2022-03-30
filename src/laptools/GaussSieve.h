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

/**@file    GaussSieve.h
 * @brief   Base class for GaussSieve.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __LAPTOOLS_GAUSS_SIEVE_H__
#define __LAPTOOLS_GAUSS_SIEVE_H__

#include <float.h>
#include <iostream>
#include <memory>
#include "Lattice.h"
#include "KleinSampler.h"
#include "VectorElement.h"
#include "VectorElementPool.h"
#include "Config.h"
#include "Def.h"


namespace LapTools
{


///
/// @class GaussSieve
///
/// @tparam BasisFloat floating point of basis matrix
/// @tparam GSFloat floating point of Gram-Schmidt matrix
///
template<typename BasisFloat=int, typename GSFloat=double>
class GaussSieve
{

using VectorElementType = VectorElementBase<BasisFloat>;
using VectorElementPtr = std::shared_ptr<VectorElementType>;
using VectorElementQueue = VectorElementBasePool<BasisFloat>;
using LatticePtr = std::shared_ptr<Lattice<BasisFloat, GSFloat>>;

private:

   LatticePtr  L;                      ///< lattice
   Config   config;                    ///< hyper parameters
   int      rank;                      ///< solver-id
   int      threadId;                  ///< solver-id
   int      verbose;                   ///< verbose level

   std::shared_ptr<std::ofstream>   osCsvLog;   ///< ostream for csv log


protected:
   int      maxCollision;       ///< terminate when the number of collisions becomes greater than it
   VectorElementQueue List;      ///< list
   VectorElementQueue Stack;     ///< stack
   KleinSampler<BasisFloat, GSFloat> kleinSampler; ///< sampler
   int      nCollision;          ///< number of collisions
   long int nLoop;               ///< number of loop
   LatticeVector<BasisFloat> bestVector;  ///< shortest vector found in gaussSieve

   double   timeLimit;           ///< time limit

   decltype(List.getIterBegin()) iterList;      ///< allocated memory for GaussReduce
   decltype(List.getIterBegin()) iterListBegin; ///< allocated memory for GaussReduce
   decltype(List.getIterBegin()) iterListEnd;   ///< allocated memory for GaussReduce

   bool     hasInitialized;      ///< flag for initialize

   double   bestObjectiveValue;  ///< squared norm of shortest vector

   double   startTime;           ///< time of start reduction
   double   runningTime;         ///< time of just running algorthm


public:

   ///
   /// @brief Constructor
   ///
   GaussSieve(){}


   ///
   /// @brief Constructor
   /// @param[in] inL         Lattice
   /// @param[in] inVerbose   <= 0 : not, 1 : light, 2: heavy
   /// @param[in] inRank      solver-id
   /// @param[in] inThreadId  solver-id
   ///
   GaussSieve(
         LatticePtr inL,
         int inRank=1,
         int inThreadId=0,
         int inVerbose=0
         )
      :
         rank(inRank),
         threadId(inThreadId),
         verbose(inVerbose),
         osCsvLog(nullptr),
         maxCollision(-1),
         nCollision(0),
         nLoop(0),
         timeLimit(-1),
         hasInitialized(false),
         bestObjectiveValue(DBL_MAX),
         startTime(0.0),
         runningTime(0.0)
   {
      L = inL;
      config = L->config;
      kleinSampler = KleinSampler<BasisFloat, GSFloat>(L);
   }


   ///
   /// @brief deconstructor
   ///
   virtual ~GaussSieve(){}


   ///
   /// @brief replace lattice object
   /// @param[in] inL   lattice
   ///
   virtual void resetLattice(
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
   virtual LatticePtr getLattice(
         )
   {
      return L;
   }


   ///
   /// @brief set size of list and statck
   /// @param[in] listsize
   /// @param[in] stacksize
   /// @param[in] inMaxCollision
   ///
   void init(
         int listsize,
         int stacksize,
         int inMaxCollision=-1
       )
   {
      List = VectorElementQueue(listsize);
      Stack = VectorElementQueue(stacksize);
      if( inMaxCollision > 0 ) maxCollision = inMaxCollision;
      hasInitialized = true;
   }


   ///
   /// @brief Gauss Sieve Algorithm
   /// @param[out] shortestV  vector
   /// @param[in] timeLimit   time limit
   ///
   virtual bool gaussSieve(
         LatticeVector<BasisFloat> &shortestV,
         double inTimeLimit=-1
         );


   ///
   /// @brief communicate with LC
   /// @param[out] shouldAbort true if it should abort else false
   ///
   virtual bool communicate(
         bool &shouldAbort
         )
   {
      return true;
   }


   ///
   /// @brief post process after reduce
   /// @param[in] v lattice vector has been reduced
   ///
   virtual void postProcessOfReduce(
         VectorElementPtr v
         )
   {}


   ///
   /// @brief Gauss Reduce Algorithm
   /// @details 1. vector p is reduced by v in List such that |v| <= |p|
   ///          2. vectors v in List are reduced by p such that |v| > |p|
   /// @param[in] p lattice vector which will be reduced
   /// @return false if p becomes 0-vector after reduce process else true
   ///
   virtual bool gaussReduce(
         VectorElementPtr p
         );


   ///
   /// @brief reduce v by w such that |w| < |v|
   /// @param[in] v lattice vector which will be reduced
   /// @param[in] w lattice vector
   /// @return true if v is reduced by w else false
   ///
   virtual bool reduce(
         VectorElementPtr v,
         VectorElementPtr w
         );


   ///
   /// @brief update bestObjectiveValue
   /// @param[in] sigh line-header character
   /// @param[in] shortest vector element
   ///
   virtual bool updateBestObjectiveValue(
         char sigh,
         VectorElementPtr v
         );


   ///
   /// @brief createt log
   /// @param[in] sigh line-header character
   ///
   virtual void outputLog(
      char sigh
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
   /// @brief getter of List
   ///
   virtual VectorElementQueue & getList(
         )
   {
      return List;
   }


   ///
   /// @brief getter of Stack
   ///
   virtual VectorElementQueue & getStack(
         )
   {
      return Stack;
   }


   ///
   /// @brief getter of number of loop
   ///
   virtual int getNLoop(
         )
   {
      return nLoop;
   }


   ///
   /// @brief getter of number of collisions
   ///
   virtual int getNCollision(
         )
   {
      return nCollision;
   }


   ///
   /// @brief getter of runnint time
   ///
   virtual double getRunningTime(
         )
   {
      return runningTime;
   }


   ///
   /// @brief getter of incumbent vector
   ///
   virtual LatticeVector<BasisFloat> & getBestVector(
         )
   {
      return bestVector;
   }


   ///
   /// @brief get header of csv formed log
   ///
   virtual std::string getCsvLogHeader(
      );


   ///
   /// @brief set verbose
   /// @param[in] inVerbose verbose level
   ///
   virtual void setVerbose(
         int inVerbose
         )
   {
      verbose = inVerbose;
   }


   ///
   /// @brief get verbose
   /// @return inVerbose verbose level
   ///
   virtual int getVerbose(
         )
   {
      return verbose;
   }


}; // class GaussSieve


}  // namespace LapTools


#endif // __LAPTOOLS_GAUSS_SIEVE_H__
