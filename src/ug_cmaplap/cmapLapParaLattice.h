/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*          This file is part of the program and software framework          */
/*  CMAP-LAP --- Configurable Massively Parallel Solver for Lattice Problems */
/*                                                                           */
/*  Copyright Written by Nariaki Tateiwa <n-tateiwa@kyudai.jp>,              */
/*                       Yuji Shinano <shinano@zib.de>,                      */
/*            Copyright (C) 2021 by Zuse Institute Berlin,                   */
/*            licensed under LGPL version 3 or later.                        */
/*            Commercial licenses are available through <licenses@zib.de>    */
/*                                                                           */
/* This code is free software; you can redistribute it and/or                */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>.     */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    cmapLapParaLattice.h
 * @brief   Base class for Lattice.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_LATTICE_H__
#define __CMAP_LAP_PARA_LATTICE_H__


#include <stddef.h>
#include <iostream>
#include <memory>
#include "NTL/mat_ZZ.h"

#include "cmapLapParaDef.h"


namespace ParaCMapLAP
{


template<typename GSFloat=double>
using MatrixMu = Eigen::Matrix<GSFloat, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template<typename GSFloat=double>
using VectorMu = Eigen::Matrix<GSFloat, Eigen::Dynamic, 1>;

template<typename GSFloat=double>
using VectorB = Eigen::Matrix<GSFloat, Eigen::Dynamic, 1>;


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

private:

   VectorMu<GSFloat> r__;           ///< allcate for calculation GSO
   MatrixMu<GSFloat> p__;           ///< allcate for calculation GSO
   MatrixMu<GSFloat> q__;           ///< allcate for calculation GSO


public:
   LatticeBasis<BasisFloat> basis;  ///< lattice basis
   MatrixMu<GSFloat>       mu;      ///< GSO coefficients
   VectorB<GSFloat>        B;       ///< squared norms of b*i

   int n;                           ///< dimension ( = the number of basis columns)
   int m;                           ///< the number of basis rows
   double GH;                       ///< the norm of the shortest vector by gaussian heuristic
   double logVolume;                ///< the logarithm of volume of the lattice

   int GSOType;                     ///< GSO method, 0: setGSO_CFA, 1: setMGSO
   double MLLLdelta;                ///< delta of MLLL in NTL

   /// constructor
   Lattice(){}
   Lattice(std::string basefile): GSOType(0), MLLLdelta(0.99){ readFile(basefile); }
   Lattice(LatticeBasis<BasisFloat> &latticeBasis): GSOType(0), MLLLdelta(0.99){ constructFromBasis(latticeBasis); }
   Lattice(LatticeBasis<BasisFloat> *latticeBasis): GSOType(0), MLLLdelta(0.99){ constructFromBasis(*latticeBasis); }
   Lattice(std::shared_ptr<LatticeBasis<BasisFloat>> latticeBasis): GSOType(0), MLLLdelta(0.99){ constructFromBasis(*latticeBasis); }
   Lattice(int *basisAarray, int row, int col): GSOType(0), MLLLdelta(0.99){ constructFromArray(basisAarray, row, col); }
   virtual void constructFromBasis(LatticeBasis<BasisFloat> &latticeBasis);
   virtual void constructFromArray(int *basisAarray, int row, int col);
   virtual void readFile(std::string basisfile);

   // destructor
   virtual ~Lattice(){}

   // base function
   virtual double projectedGH(int k=0, int l=-1);
   virtual long double logProjectedVolume(int k=0, int l=-1);
   virtual size_t hash();

   // GSO
   virtual bool setGSO(int k=0, int l=-1);
   virtual bool setGSO_CFA(int k, int l);
   virtual bool setMGSO(int k, int l);

   // randomize
   virtual bool randomize(int seed, int begin=0, int end=-1);

   // reduction
   virtual bool sizeReduce(double eta);
   virtual bool sizeReduce(int i, double eta);
   virtual bool insertMlll(LatticeVector<BasisFloat> &v, int j, int k, int l);
   virtual bool merge(Lattice &other);

   // metrics
   virtual double shortestNorm();
   virtual double approxFactor(double norm=-1, int k=0);
   virtual double hermiteFactor(double norm=-1);
   virtual double rootHermiteFactor(double norm=-1);
   virtual double logOrthogonalityDefect();
   virtual long double enumCost(double R=-1, int begin=0, int end=-1);
   virtual double slopeGSA(int h=-1);

   // utils
   virtual Lattice<BasisFloat, GSFloat> copy(int h=-1);
   virtual NTL::mat_ZZ toNLTMat();
   virtual void fromNLTMat(NTL::mat_ZZ &Mbasis);
   virtual void resize(int rows);
   virtual bool isMoreReducedThan(Lattice &other);
   virtual bool isMoreReducedThan(Lattice &other, int &index);
   virtual void outputBasis(std::ostream *os = &std::cout);
   virtual void writeBasis(std::string filenam);
   virtual std::string toSimpleString();

};


}  // namespace ParaCMapLAP

#endif // __CMAP_LAP_PARA_LATTICE_H__
