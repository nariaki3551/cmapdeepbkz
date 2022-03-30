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

/**@file    Lattice.h
 * @brief   Base class for Lattice.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __LAPTOOLS_LATTICE_H__
#define __LAPTOOLS_LATTICE_H__

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stddef.h>
#include "NTL/mat_ZZ.h"
#include "Def.h"
#include "Config.h"


namespace LapTools
{


///
/// @class Lattice
///
/// @tparam BasisFloat floating point of basis matrix
/// @tparam GSFloat floating point of Gram-Schmidt matrix
/// @note Notation Variable Description
///       -----------------------------
///       bi  basis.row(i) the i-th basis vector
///       b*i -            the i-th GSO vector
///       Bi  B(i)         squared norm of b*i
///       mu  mu           GSO matrix
///
template<typename BasisFloat=int, typename GSFloat=double>
class Lattice
{

public:
   LatticeBasis<BasisFloat> basis;  ///< lattice basis
   MatrixMu<GSFloat>       mu;      ///< GSO coefficients
   VectorB<GSFloat>        B;       ///< squared norms of b*i

   int n;                           ///< dimension ( = the number of basis columns)
   int m;                           ///< the number of basis rows
   double GH;                       ///< the norm of the shortest vector by gaussian heuristic
   double logVolume;                ///< the logarithm of volume of the lattice

   Config config;                   ///< parameters

   VectorMu<GSFloat> r__;           ///< allcate for calculation GSO
   MatrixMu<GSFloat> p__;           ///< allcate for calculation GSO
   MatrixMu<GSFloat> q__;           ///< allcate for calculation GSO

   /// constructor
   Lattice(){}
   Lattice(std::string basefile){ readFile(basefile); }
   Lattice(LatticeBasis<BasisFloat> &latticeBasis){ constructFromBasis(latticeBasis); }
   Lattice(LatticeBasis<BasisFloat> *latticeBasis){ constructFromBasis(*latticeBasis); }
   Lattice(int *basisAarray, int row, int col){ constructFromArray(basisAarray, row, col); }
   void constructFromBasis(LatticeBasis<BasisFloat> &latticeBasis);
   void constructFromArray(int *basisAarray, int row, int col);
   void readFile(std::string basisfile);
   void setConfig(Config &inConfig){ config = inConfig; }

   // destructor
   virtual ~Lattice(){};

   // base function
   double projectedGH(int k=0, int l=-1);
   long double logProjectedVolume(int k=0, int l=-1);
   size_t hash();

   // projected function
   double projectedSqnorm(int k, LatticeVector<BasisFloat> &v);
   double projectedSqnormFromCoeff(int k, LatticeVector<BasisFloat> &coeffs);

   // GSO
   bool setGSO(int k=0, int l=-1);
   bool setGSO_CFA(int k, int l);
   bool setMGSO(int k, int l);

   // randomize
   bool randomize(int seed, int begin=0, int end=-1);

   // reduction
   bool sizeReduce(double eta);
   bool sizeReduce(int i, double eta);
   bool insertMlll(LatticeVector<BasisFloat> &v, int j, int k, int l);
   bool merge(Lattice &other);

   // metrics
   double shortestNorm();
   double approxFactor(double norm=-1, int k=0);
   double hermiteFactor(double norm=-1);
   double rootHermiteFactor(double norm=-1);
   double logOrthogonalityDefect();
   long double enumCost(double R=-1, int begin=0, int end=-1);
   double slopeGSA(int h=-1);

   // utils
   Lattice<BasisFloat, GSFloat> copy(int h=-1);
   NTL::mat_ZZ toNLTMat();
   void fromNLTMat(NTL::mat_ZZ &Mbasis);
   void resize(int rows);
   bool isMoreReducedThan(Lattice &other);
   bool isMoreReducedThan(Lattice &other, int &index);
   void outputBasis(std::ostream *os = &std::cout);
   void writeBasis(std::string filenam);
   std::string toSimpleString();
   std::string toSimpleString(LatticeVector<BasisFloat> &v);

};


}  // namespace LapTools

#endif // __LAPTOOLS_LATTICE_H__
