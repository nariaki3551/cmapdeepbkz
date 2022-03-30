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

/**@file    RecursiveDeepBkz.h
 * @brief   Base class for DeepBKZ with recursive.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __LAPTOOLS_RECURSIVE_DEEP_BKZ_H__
#define __LAPTOOLS_RECURSIVE_DEEP_BKZ_H__

#include <memory>                                     // for shared_ptr

#include "Lattice.h"
#include "DeepBkz.h"


namespace LapTools
{


///
/// @class RecursiveDeepBkz
///
/// @tparam BasisFloat floating point of basis matrix
/// @tparam GSFloat floating point of Gram-Schmidt matrix
/// @tparam EnumGSFloat floating point of Gram-Schmidt matrix for Enumeration
///
template<typename BasisFloat=int, typename GSFloat=double, typename EnumGSFloat=double>
class RecursiveDeepBkz : public DeepBkz<BasisFloat, GSFloat, EnumGSFloat>
{

using LatticePtr = std::shared_ptr<Lattice<BasisFloat, GSFloat>>;


private:

   int depth;     // recursive depth


public:



   ///
   /// @brief Constructor
   ///
   RecursiveDeepBkz(){}


   ///
   /// @brief Constructor
   /// @param[in] inL         Lattice
   /// @param[in] inVerbose   <= 0 : not, 1 : light, 2: heavy
   /// @param[in] inRank      solver-id
   /// @param[in] inThreadId  solver-id
   ///
   RecursiveDeepBkz(
         LatticePtr inL,
         int inRank=1,
         int inThreadId=0,
         int inVerbose=0
         )
      :
         DeepBkz<BasisFloat, GSFloat, EnumGSFloat>(inL, inRank, inThreadId, inVerbose),
         depth(0)
   {}


   ///
   /// @brief deconstructor
   ///
   virtual
   ~RecursiveDeepBkz(){}


   ///
   /// @brief one loop of algorithm
   /// @param[in] k
   /// @param[out] z
   /// @param[out] shouldAbort true if it should abort else false
   /// @return bool true: normal terminate, false: abnormal one
   ///
   bool
   step(
         int k,
         int &z,
         bool &shouldAbort
         );


   ///
   /// @brief setter of depth
   ///
   void
   setDepth(
         int inDepth
         )
   {
      depth = inDepth;
   }


}; // class RecursiveDeepBkz


}  // namespace LapTools


#endif // __LAPTOOLS_RECURSIVE_DEEP_BKZ_H__
