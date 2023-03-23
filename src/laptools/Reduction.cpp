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

/**@file    Reduction.cpp
 * @brief   Base class for Lattice basis reduction of NTL library.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "Reduction.h"

#include <NTL/LLL.h>       // for BKZ_FP, BKZ_QP, LLL_FP, LLL_QP
#include <type_traits>     // for is_same
#include "NTL/mat_ZZ.h"    // for mat_ZZ
#include "Lattice.h"       // for Lattice


namespace LapTools
{


///
/// @brief run LLL algorithm in NTL
///
template<typename BasisFloat, typename GSFloat>
bool
Reduction<BasisFloat, GSFloat>::lll(
      )
{
   NTL::mat_ZZ Mbasis = L->toNLTMat();
   if( std::is_same<GSFloat, double>() )
      NTL::LLL_FP(Mbasis, delta, 0, 0, verbose);
   else
      NTL::LLL_QP(Mbasis, 0.99, 0, 0, verbose);
   L->fromNLTMat(Mbasis);
   return true;
}


///
/// @brief run bkz algorithm in NTL
/// @param[in] beta blocksize
///
template<typename BasisFloat, typename GSFloat>
bool
Reduction<BasisFloat, GSFloat>::bkz(
      int beta,
      long prune
      )
{
   NTL::mat_ZZ Mbasis = L->toNLTMat();
   if( std::is_same<GSFloat, double>() )
      NTL::BKZ_FP(Mbasis, delta, beta, prune, 0, verbose);
   else
      NTL::BKZ_QP(Mbasis, delta, beta, prune, 0, verbose);
   L->fromNLTMat(Mbasis);
   return true;
}


///
/// instantiation
///
template class Reduction<int, double>;
template class Reduction<int, long double>;
template class Reduction<long int, double>;
template class Reduction<long int, long double>;


}  // namespace LapTools
