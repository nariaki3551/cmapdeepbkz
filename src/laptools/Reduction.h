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

/**@file    Reduction.h
 * @brief   Base class for Lattice basis reduction of NTL library.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __LAPTOOLS_REDUCTION_H__
#define __LAPTOOLS_REDUCTION_H__

#include <memory>
#include "Lattice.h"
#include "Config.h"


namespace LapTools
{


///
/// @class Reduction
///
/// @tparam BasisFloat floating point of basis matrix
/// @tparam GSFloat floating point of Gram-Schmidt matrix
///
template<typename BasisFloat=int, typename GSFloat=double>
class Reduction
{

using LatticePtr = std::shared_ptr<Lattice<BasisFloat, GSFloat>>;


private:

   LatticePtr  L;          ///< lattice
   Config      config;     ///< hyper parameters
   int         rank;       ///< solver-id
   int         threadId;   ///< solver-id
   int         verbose;    ///< verbose <= 0: none, 1: light, 2: medium, 3: heavy

   double      delta;      ///< LLL parameter


public:

   ///
   /// @brief Constructor
   ///
   Reduction(){}


   ///
   /// @brief Constructor
   /// @param[in] inL         Lattice
   /// @param[in] inRank      solver-id
   /// @param[in] inThreadId  solver-id
   /// @param[in] inVerbose   <= 0 : not, 1 : light, 2: heavy
   ///
   Reduction(
         LatticePtr inL,
         int inRank=1,
         int inThreadId=0,
         int inVerbose=0
         )
      :
         rank(inRank),
         threadId(inThreadId),
         verbose(inVerbose),
         delta(0.99)
   {
      L = inL;
      config = L->config;

      // if( this->config.OutputCsvLogFile )
      // {
      //    ofsCsvLog.open(this->config.CsvLogFile, std::ios::app);
      //    ofsCsvLog << getCsvLogHeader() << std::endl;
      // }
   }


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
   /// @brief run LLL algorithm in NTL
   ///
   virtual bool lll(
       );


   ///
   /// @brief run bkz algorithm in NTL
   /// @param[in] beta blocksize
   ///
   virtual bool bkz(
         int beta
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


}; // class DeepLll


}  // namespace LapTools


#endif // __LAPTOOLS_REDUCTION_H__
