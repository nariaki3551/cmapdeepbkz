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

/**@file    RecursiveDeepBkz.cpp
 * @brief   Base class for DeepBKZ with recursive.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "RecursiveDeepBkz.h"

#include <math.h>         // for round
#include <algorithm>      // for min, max
#include "Config.h"       // for Config
#include "DeepLll.h"      // for DeepLll
#include "Enumeration.h"  // for Enumeration
#include "Lattice.h"      // for Lattice


namespace LapTools
{


///
/// @brief one loop of algorithm with recursive
/// @param[in] k
/// @param[out] z
/// @param[out] shouldAbort true if it should abort else false
/// @return bool return true if enumeration process is successed else false
///
template<typename BasisFloat, typename GSFloat, typename EnumGSFloat>
bool
RecursiveDeepBkz<BasisFloat, GSFloat, EnumGSFloat>::step(
      int k,
      int &z,
      bool &shouldAbort
      )
{
   int begin = 0, end = this->getLattice()->m - 1;
   int l = std::min(k+this->blocksize-1, end);
   int h = std::min(l+1, end);

   // preprocess
   if( this->timeLimit > 0 && this->runningTime >= this->timeLimit )
   {
      this->outputLog('T');
      shouldAbort = true;
   }
   if( shouldAbort ){ return true; }

   // sub-SVP
   this->enumObj.init(k, l, -1, (this->timeLimit > 0 ? this->timeLimit-this->runningTime : -1));
   bool enumSuccess = this->enumObj.projectedEnum(this->v, this->coeffv);

   // insert
   if( enumSuccess ){ this->postProcess(k, l); }

   // reduction for next sub-SVP
   if( enumSuccess )
   {
      if( this->blocksize > this->config.RecursiveLowerBeta )
      {
         int recursiveBlocksize = std::max(
               this->config.RecursiveLowerBeta,
               std::round(this->blocksize/2.0));
         RecursiveDeepBkz<BasisFloat, GSFloat, EnumGSFloat> bkzObj{
            this->L, this->rank, this->threadId, this->verbose};
         bkzObj.setDepth(depth+1);
         bkzObj.deepbkz(recursiveBlocksize);
      }
      z = 0;
   }
   else
   {
      this->lllObj.deeplll(h-1, begin, h);
      z++;
   }

   this->updateBestObjectiveValue('*');
   return enumSuccess;
}


///
/// instantiation
///
template class RecursiveDeepBkz<int, double, double>;
template class RecursiveDeepBkz<int, long double, double>;
template class RecursiveDeepBkz<int, long double, long double>;
template class RecursiveDeepBkz<long int, double, double>;
template class RecursiveDeepBkz<long int, long double, double>;
template class RecursiveDeepBkz<long int, long double, long double>;


}  // namespace LapTools
