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

/**@file    Def.h
 * @brief   Defines for Laptools.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __LAPTOOLS_DEF_H__
#define __LAPTOOLS_DEF_H__
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>
#include <cfloat>

#include <eigen3/Eigen/Core>

using namespace Eigen;

namespace LapTools
{


template<typename BasisFloat=int>
using LatticeBasis = Matrix<BasisFloat, Dynamic, Dynamic, RowMajor>;
template<typename BasisFloat=int>
using LatticeVector = Matrix<BasisFloat, Dynamic, 1>;

template<typename GSFloat=double>
using MatrixMu = Matrix<GSFloat, Dynamic, Dynamic, RowMajor>;
template<typename GSFloat=double>
using VectorMu = Matrix<GSFloat, Dynamic, 1>;

template<typename GSFloat=double>
using VectorB = Matrix<GSFloat, Dynamic, 1>;


#define SVP_MINEPSILON                   1e-5  ///< minimum value for any numerical epsilon


}  // namespace LapTools


#endif // __LAPTOOLS_DEF_H__
