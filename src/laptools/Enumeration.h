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

/**@file    Enumeration.h
 * @brief   Base class for Enumeration.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __LAPTOOLS_ENUM_H__
#define __LAPTOOLS_ENUM_H__

#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include <deque>
#include <float.h>
#include <eigen3/Eigen/Core>
#include "Def.h"
#include "Config.h"
#include "Lattice.h"
#include "Log.h"


using namespace Eigen;


namespace LapTools
{


///
/// @class Enumeration
///
/// @tparam BasisFloat floating point of basis matrix
/// @tparam GSFloat floating point of Gram-Schmidt matrix
/// @tparam EnumGSFloat floating point of Gram-Schmidt matrix for Enumeration
///
template<typename BasisFloat=int, typename GSFloat=double, typename EnumGSFloat=double>
class Enumeration
{

using LatticePtr = std::shared_ptr<Lattice<BasisFloat, GSFloat>>;


private:

   LatticePtr  L;                      ///< lattice
   Config      config;                 ///< hyper parameters
   int         rank;                   ///< solver-id
   int         threadId;               ///< solver-id
   int         verbose;                ///< verbose <= 0: none, 1: light, 2: medium, 3: heavy

   bool        hasFinished;            ///< true if search has finished else false
   std::shared_ptr<std::ofstream>   osCsvLog;   ///< ostream for csv log


protected:

   MatrixMu<EnumGSFloat>   mu;         ///< memory for enumeration
   MatrixMu<EnumGSFloat>   sigma;      ///< memory for enumeration
   VectorMu<EnumGSFloat>   r;          ///< memory for enumeration
   VectorMu<EnumGSFloat>   rho;        ///< memory for enumeration
   LatticeVector<BasisFloat> coeffs;   ///< memory for enumeration
   VectorMu<EnumGSFloat>   c;          ///< memory for enumeration
   VectorXi    w;                      ///< memory for enumeration
   LatticeVector<BasisFloat> bestCoeffv;  ///< memory for enumeration
   LatticeVector<BasisFloat> bestVector;  ///< memory for enumeration
   int         lastNonZero;            ///< memory for enumeration
   int         k;                      ///< memory for enumeration
   int         numFixedCoeffs;         ///< number of fixed coefficients during search
   double      tmpScaler__;            ///< allocated memory for enumeration

   int         begin;                  ///< begin index of basis for enumeration
   int         end;                    ///< end index of basis for enumeration
   EnumGSFloat lowerBound;             ///< terminate when objective value becomes less than lowerBound
   double      timeLimit;              ///< time limit

   long int    nSearch;                ///< number of searched nodes
   EnumGSFloat bestObjectiveValue;     ///< squared norm of begin-th GSO vectors

   std::vector<EnumGSFloat> sqR;       ///< upperbound about coefficients of vector
   EnumGSFloat sqRadius;               ///< squared radius of enumeration search
   long double enumCost;               ///< number of enumeration tree nodes

   double      startTime;              ///< time of start reduction
   double      runningTime;            ///< time of just running algorthm
   double      nextLogTime;            ///< time for output next log


public:

   ///
   /// @brief Constructor
   ///
   Enumeration(){}


   ///
   /// @brief Constructor
   /// @param[in] inL         Lattice
   /// @param[in] inVerbose   <= 0 : not, 1 : light, 2: heavy
   /// @param[in] inRank      solver-id
   /// @param[in] inThreadId  solver-id
   ///
   Enumeration(
         LatticePtr inL,
         int inRank=1,
         int inThreadId=0,
         int inVerbose=0
         )
      :
         rank(inRank),
         threadId(inThreadId),
         verbose(inVerbose),
         hasFinished(false),
         osCsvLog(nullptr),
         numFixedCoeffs(0),
         timeLimit(-1),
         nSearch(0),
         bestObjectiveValue(DBL_MAX),
         sqRadius(0.0),
         enumCost(0.0),
         startTime(0.0),
         runningTime(0.0),
         nextLogTime(0.0)
   {
      L = inL;
      config = L->config;
      begin = 0;
      end = L->m-1;
      sqR.resize(L->m);
      config = L->config;
      sigma .resize(L->m+1, L->m);
      r     .resize(L->m+1);
      rho   .resize(L->m+1);
      coeffs.resize(L->m);
      c     .resize(L->m);
      w     .resize(L->m);
      bestCoeffv.resize(L->m);
   }


   ///
   /// @brief replace lattice object
   /// @param[in] inL   lattice
   ///
   virtual void resetLattice(
         LatticePtr inL
         )
   {
      L = inL; config = L->config;
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
   /// @brief set parameters of enumeration
   /// @param[in] inBegin  begin index of basis for reduction
   /// @param[in] inEnd    end   index of basis for reduction
   /// @param[in] inLowerBound
   /// @param[in] inTimeLimit
   ///
   virtual void init(
         int inBegin=0,
         int inEnd=-1,
         double inLowerBound=-1,
         double inTimeLimit=-1
         );


   ///
   /// @brief setter of lowerBound
   /// @param[in] inLowerBound lower bound of norm
   ///
   virtual void setLowerBound(
         double inLowerBound
         );


   ///
   /// set initial point of search
   /// @param[in] fixedCoeff coefficient array such that d-tail of coeffs equals fixedCoeff, where d = fixedCoeff.size()
   ///
   virtual void setFixedCoeffs(
         LatticeVector<BasisFloat> &fxedCoeffs
         );


   ///
   /// @brief Enum algorithm
   /// @details Vector resv satisifies norm(\pi_k(resv))^2 <= sqR_k for all k in be from begin to end
   ///          resv = sum( coeffs[i]*v_i for all begin <= i <= end )
   /// @param[out] resv          vector found by Enum
   /// @parem[out] coeffResv     coefficient of vector found by Enum
   /// @param[in] updateLimit    terminate when Enum update solution updateLimit times
   /// @param[in] inTimeLimit    time limit
   /// @return bool true if Enum found a shorter vector else false
   ///
   virtual bool projectedEnum(
         LatticeVector<BasisFloat> &resv,
         LatticeVector<BasisFloat> &coeffResv,
         int updateLimit=-1,
         double inTimeLimit=-1
         );


   ///
   /// @brief post process of search node
   ///
   virtual bool postProcess(
         )
   {
      return true;
   }


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
   /// @brief search just one node
   /// @return false if all nodes are searched else true
   ///
   virtual bool stepProjectedEnum(
      );


   ///
   /// @brief Calculate norm(pi_begin(v))^2
   /// @param[in] begin begin index of blockSize matrix
   /// @param[in] end last index of blockSize matrix
   /// @param[in] coefv coefficients of vector
   ///
   virtual double getProjectedSquaradNorm(
         int begin,
         int end,
         LatticeVector<BasisFloat> &coeffs
         );


   ///
   /// @brief update bestObjectiveValue
   /// @param[in] sigh line-header character
   /// @param[in] objectiveValue new objectiveValue
   ///
   virtual bool updateBestObjectiveValue(
         char sigh,
         EnumGSFloat objectiveValue
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


   ///
   /// @brief getter of number of searched nodes
   ///
   virtual long int getNSearch(
         )
   {
      return nSearch;
   }


   ///
   /// @brief getter of running time
   ///
   virtual double getRunningTime(
         )
   {
      return runningTime;
   }


   ///
   /// @brief update enumCost
   ///
   virtual void updateEnumCost(
         )
   {
      enumCost = L->enumCost(std::sqrt(sqRadius), begin, end);
   }

   ///
   /// @brief get progress of searching
   ///
   virtual double getProgress(
         )
   {
      return nSearch / enumCost;
   }


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
   /// @brief getter of verbose
   /// @return inVerbose verbose level
   ///
   virtual int getVerbose(
         )
   {
      return verbose;
   }


   ///
   /// @brief getter of best vector
   /// @return lattice vector
   ///
   virtual LatticeVector<BasisFloat> & getBestVector(
         );


   ///
   /// @brief getter of coeffs
   ///
   virtual LatticeVector<BasisFloat> & getCoeffs(
         )
   {
      return coeffs;
   }

   ///
   /// @brief getter of coeffs
   ///
   virtual int getDepth(
         )
   {
      return k;
   }


   ///
   /// @brief divide enumeration tree
   /// @param[in] N number of sub enumeration tree
   /// @note Depending on the input, it may only return a number of divisions smaller than N.
   ///
   virtual std::deque<LatticeVector<BasisFloat>> divide(
         int N
         );


   ///
   /// @brief calculation lower and upper bound of coefficient
   /// @param[in] w coefficient array such that d-tail of coeffs equals w, where d = w.size()
   /// @param[out] Lower lower bound of d+1-tail coefficients
   /// @param[out] Upper upper bound of d+1-tail coefficients
   ///
   virtual void range(
         LatticeVector<BasisFloat> &w,
         int &Lower,
         int &Upper
        );


   ///
   /// @brief check whether search has finished
   /// @return true if search has finished else false
   ///
   virtual bool finish(
         )
   {
      return hasFinished;
   }


}; // class Enumeration


}  // namespace LapTools


#endif // __LAPTOOLS_ENUM_H__
