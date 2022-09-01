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

/**@file    Lattice.cpp
 * @brief   Base class for Lattice.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "Lattice.h"

#include <algorithm>
#include <utility>
#include <assert.h>
#include <float.h>
#include <memory>
#include "NTL/LLL.h"
#include "NTL/ZZ.h"
#include "NTL/matrix.h"
#include "NTL/tools.h"
#include <eigen3/Eigen/Core>
#include "Config.h"
#include "Def.h"
#include "GenUnimat.h"


namespace LapTools
{


///
/// @brief constructor
/// @param[in] latticeBasis basis
///
template<typename BasisFloat, typename GSFloat>
void
Lattice<BasisFloat, GSFloat>::constructFromBasis(
      LatticeBasis<BasisFloat> &latticeBasis
      )
{
   n = latticeBasis.cols();
   m = latticeBasis.rows();
   basis = latticeBasis;
   mu.resize(m, m);
   B.resize(m);
   r__.resize(m);
   setGSO();
   GH = projectedGH();
   logVolume = logProjectedVolume();
   config = Config();
}


///
/// @brief constructor
/// @param[in] basisArray element of basis of RowMajor
/// @param[in] row   number of rows of basis
/// @param[in] col   number of columns of basis
///
template<typename BasisFloat, typename GSFloat>
void
Lattice<BasisFloat, GSFloat>::constructFromArray(
      int *basisArray,
      int row,
      int col
      )
{
   LatticeBasis<BasisFloat> inBasis;
   inBasis.resize(row, col);
   int index = 0;
   for ( int i = 0; i < row; ++i )
   {
      for ( int j = 0; j < col; ++j )
      {
         inBasis(i, j) = basisArray[index++];
      }
   }
   constructFromBasis(inBasis);
}


///
/// @brief read basis from file
/// @param[in] basisfile basis file path
///
template<typename BasisFloat, typename GSFloat>
void
Lattice<BasisFloat, GSFloat>::readFile(
      std::string basisfile
      )
{
   std::ifstream fin(basisfile);
   if( !fin.is_open() )
   {
      std::cout << "can not open input file: " << basisfile << std::endl;
      exit(0);
   }
   NTL::Mat<int> inputBasis;
   fin >> inputBasis;

   n = inputBasis.NumCols();
   m = inputBasis.NumRows();

   LatticeBasis<BasisFloat> eigenInputBasis;
   eigenInputBasis.resize(m, n);
   for( int row = 0; row < m; ++row )
   {
      for( int col = 0; col < n; ++col )
      {
         eigenInputBasis(row, col) = inputBasis[row][col];
      }
   }
   constructFromBasis(eigenInputBasis);
}


///
/// @brief Calculate GH(L') where L' is projected lattice {pi_k(b_k), ..., pi_k(b_l)}
/// @details GH(L') = (vol(L')/nu(d))**(1/d)
///          = exp( (log(vol(L')) - log(nu(d))) / d )
///          where d is dimension of projected lattice (l-k+1).
///          vol(L') = \prod_{i=k}^{l} norm(b*i)
/// @param[in] k index of basis
/// @param[in] l index of basis
///
template<typename BasisFloat, typename GSFloat>
double
Lattice<BasisFloat, GSFloat>::projectedGH(
      int k,
      int l
      )
{
   if( l < 0 ) l = m - 1;
   int d = l-k+1;
   long double __logProjectedVolume = logProjectedVolume(k, l);
   long double logUnitSphereVolume = d*std::log(M_PI)*0.5 - std::log(std::tgamma(d*0.5+1.0));
   return std::exp( (__logProjectedVolume-logUnitSphereVolume)/d );
}


///
/// @brief log volume of the projected lattice L' = L({pi_k(b_k), ..., pi_k(b_l)})
/// @return log( vol(L') ) = log( prod( norm(b*i) for k <= i <= l ) )
///                        = sum( log( norm(b*i) for k <= i <= l ) )
/// @param[in] k index of basis
/// @param[in] l index of basis
///
template<typename BasisFloat, typename GSFloat>
long double
Lattice<BasisFloat, GSFloat>::logProjectedVolume(
      int k,
      int l
      )
{
   if( l < 0 ) l = m - 1;
   long double logVolume__ = 0.0;
   for( int i = k; i <= l; ++i ){ logVolume__ += std::log( B(i) ); };
   return logVolume__ * 0.5;
}


///
/// @brief hash of basis
/// @retrn hash value
///
template<typename BasisFloat, typename GSFloat>
size_t
Lattice<BasisFloat, GSFloat>::hash(
   )
{
   size_t seed = 0;
   for( int i = 0; i < basis.size(); ++i )
   {
     BasisFloat elem = *(basis.data() + i);
     seed ^= std::hash<BasisFloat>()(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
   }
   return seed;
}


///
/// @brief projected saured norm of vector v
/// @param[in] k index of basis
/// @param[in] vector lattice vector
/// @return pi_k(v) = sqnorm( sum( dot(v, b*i)/sqnorm(b*i) b*i for k <= i <= m) )
///
template<typename BasisFloat, typename GSFloat>
double
Lattice<BasisFloat, GSFloat>::projectedSqnorm(
      int k,
      LatticeVector<BasisFloat> &vector
      )
{
   // generate GSO matrix
   MatrixMu<GSFloat> GSO(m, m);
   for( int i = 0; i < m; i++ )
   {
      for( int j = 0; j < m; j++ ) { GSO(i, j) = basis(i,j); }
      for( int j = 0; j < i; j++ ) { GSO.row(i) -= mu(i,j) * GSO.row(j); }
   }

   // generate project vector
   LatticeVector<GSFloat> projectedV(n);
   projectedV.setZero();
   double inner = 0;
   for( int i = k; i < m; i++ )
   {
      inner = 0;
      for( int j = 0; j < m; j++ ){ inner += vector(j) * GSO(i, j); }
      projectedV += (inner / B(i)) * GSO.row(i);
   }
   return projectedV.squaredNorm();
}


///
/// @brief projected saured norm of vector v
/// @param[in] k index of basis
/// @param[in] coeffs coefficient array of lattice vector
/// @return pi_k(v), where v = sum( coeffs(i) * b(i) for all i )
///
template<typename BasisFloat, typename GSFloat>
double
Lattice<BasisFloat, GSFloat>::projectedSqnormFromCoeff(
      int k,
      LatticeVector<BasisFloat> &coeffs
      )
{
   double __projectedSqnorm = 0;
   for( int i = k; i < m; ++i )
   {
      double u = coeffs(i);
      for( int j = i+1; j < m; ++j )
      {
         u += mu(j,i) * coeffs(j);
      }
      __projectedSqnorm += u * u * B(i);
   }
   return __projectedSqnorm;
}


///
/// @brief compute and set orthogonalised basis
/// @return
/// @details set mu to GSO coefficient matrix, and
///          set B to squared norms of GSO vectors
///
template<typename BasisFloat, typename GSFloat>
bool
Lattice<BasisFloat, GSFloat>::setGSO(
      int k,
      int l
      )
{
   if( l == -1 ){ l = m-1; }
   if( config.GSOType == 0 )
      return setGSO_CFA(k, l);
   else
      return setMGSO(k, l);
}


///
/// @brief Calculate coefficients and square norm of
///        GSO b*k,...,b*l by Cholesky decomposition
/// @param[in] k index of basis
/// @param[in] l index of basis
/// @return
///
template<typename BasisFloat, typename GSFloat>
bool
Lattice<BasisFloat, GSFloat>::setGSO_CFA(
      int k,
      int l
      )
{
   int i, j, jj;
   mu.middleRows(k, l-k+1).setZero();
   for( i = k; i <= l; ++i )
      mu(i, i) = 1.0;

   B.segment(k, l-k+1).setZero();
   for( i = k; i <= l; ++i )
   {
      for( j = 0; j < i; ++j )
      {
         r__.coeffRef(j) = basis.row(i).dot(basis.row(j));
         for( jj = 0; jj < j; ++jj )
         {
            r__.coeffRef(j) -= mu.coeff(j,jj) * r__.coeff(jj);
         }
         mu.coeffRef(i,j) = r__.coeff(j)/B.coeff(j);
      }
      B.coeffRef(i) = basis.row(i).squaredNorm();
      for(j=1;j<=i;++j)
      {
         B.coeffRef(i) -= mu.coeff(i,j-1) * r__.coeff(j-1);
      }
   }
   return true;
}


///
/// @brief Calculate coefficients and square norm of
///        GSO b*k,...,b*l by Modified Gram-Schmidt Orthogonalization
/// @param[in] k index of basis
/// @param[in] l index of basis
/// @note r(i, j) = <bj, b*i> / norm(b*i) ( i < j )
///       r(i, i) = norm(b*i)
///       q.col(i) = b*_i / norm(b*i)
///                = (bi - sum(<bi, b*j>b*j for 0 <= j <= i-1) ) / norm(b*i)
///                = (bi - sum(r(j, i)b*j for 0 <= j <= i-1)   ) / norm(b*i)
/// @return false if it detects error for floating else true
///
template<typename BasisFloat, typename GSFloat>
bool
Lattice<BasisFloat, GSFloat>::setMGSO(
      int k,
      int l
      )
{
   int i, j;
   q__.setZero(n, m);
   p__.setZero(m, m);
   mu.middleRows(k, l-k+1).setZero();
   B.segment(k, l-k+1).setZero();

   for( i = k; i <= l; ++i )
   {
      // q__.col(i) = basis.row(i).cast<GSFloat>();
      for( j = 0; j < n; ++j )
         q__.coeffRef(j, i) = basis.coeff(i, j);
      for( j = 0; j < i; ++j )
      {
         p__.coeffRef(j, i) = q__.col(j).dot(q__.col(i));
         q__.col(i) -= p__.coeff(j, i) * q__.col(j);
      }
      p__.coeffRef(i, i) = q__.col(i).norm();
      q__.col(i) /= p__.coeff(i, i);
   }

   for( i = k; i <= l; ++i )
   {
      mu.coeffRef(i, i) = 1.0;
      B.coeffRef(i) = p__.coeff(i, i) * p__.coeff(i, i);
      if( B.coeffRef(i) < LDBL_EPSILON )
      {
         std::cout << "Length error in MGSO" << std::endl;
         return false;
      }
      for( j = 0; j < i; ++j )
      {
         mu.coeffRef(i, j) = p__.coeff(j, i) / p__.coeff(j, j);
      }
   }
   return true;
}


///
/// @brief shortest norm
///
template<typename BasisFloat, typename GSFloat>
double
Lattice<BasisFloat, GSFloat>::shortestNorm(
      )
{
   return std::sqrt(B(0));
}


///
/// @brief Transform b_i so that it is close to being
///        orthogonal to {b_0, ..., b_{m-1}} with keeping GSO vector.
/// @param[in] eta 0.501 -- 0.51
/// @return
///
template<typename BasisFloat, typename GSFloat>
bool
Lattice<BasisFloat, GSFloat>::sizeReduce(
      double eta
      )
{
   return sizeReduce(m-1, eta);
}


///
/// @brief Transform b_i so that it is close to being
///        orthogonal to {b_0, ..., b_{i-1}} with keeping GSO vector.
/// @param[in] i 0 -- m-1
/// @param[in] eta 0.501 -- 0.51
/// @return
///
template<typename BasisFloat, typename GSFloat>
bool
Lattice<BasisFloat, GSFloat>::sizeReduce(
      int i,
      double eta
      )
{
   int q;
   for( int j = i-1; j >= 0; j--)
   {
      if( std::abs(mu.coeff(i, j)) > eta )
      {
         q = std::round(mu.coeff(i,j));
         mu.row(i) -= q * mu.row(j);
         basis.row(i) -= q * basis.row(j);
      }
   }
   return true;
}


///
/// @brief get approximation factor
/// @param[in] norm
/// @param[in] k
/// @return norm / GH
///
template<typename BasisFloat, typename GSFloat>
double
Lattice<BasisFloat, GSFloat>::approxFactor(
      double norm,
      int k
      )
{
   assert( k >= 0 );
   if( k == 0 )
   {
      if( norm < 0 ) norm = shortestNorm();
      return norm / GH;
   }
   else
   {
      if( norm < 0 ) norm = std::sqrt(B[k]);
      return norm / projectedGH(k, m-1);
   }
}


///
/// @brief get hermite factor
/// @return norm / ( vol(L)^(1/m) ) = norm / exp( (1/m) log( vol(L) )
///                                 = exp( log(norm) - (1/m)log( vol(L) ) )
///
template<typename BasisFloat, typename GSFloat>
double
Lattice<BasisFloat, GSFloat>::hermiteFactor(
      double norm
      )
{
   if( norm < 0 ) norm = shortestNorm();
   return norm / std::exp( logVolume / m );
}


///
/// @brief get root hermite factor
/// @return ( hermiteFactor ) ^ (1/m)
///
template<typename BasisFloat, typename GSFloat>
double
Lattice<BasisFloat, GSFloat>::rootHermiteFactor(
      double norm
      )
{
   return std::pow( hermiteFactor(norm), (1.0/m) );
}


///
/// @brief get logarithm of orthogonality defect
/// @return prod( norm(bi) for 0 <= i <= m-1 ) / vol(L)
///         = exp( sum( log(norm(bi)) for 0 <= i <= m-1) - log(vol(L)) )
///
template<typename BasisFloat, typename GSFloat>
double
Lattice<BasisFloat, GSFloat>::logOrthogonalityDefect(
      )
{
   double logProdNorm = 0;
   for( int i = 0; i < m; ++i )
      logProdNorm += std::log( basis.row(i).squaredNorm() );
   return logProdNorm * 0.5 - logVolume;
}


///
/// @brief get enumeration cost of projceted lattice {pi_i(bi), ..., pi_i(bj)}, where i = begin, j = end
/// @details enumCost = 1/2 * sum( H_l  for 0 <= l <= m-1 ), where
///          H_l := nu(R, l+1) / vol(\pi_{m-1-l}(L)),
///          vol(\pi_{s}(L)) := prod( norm(b*i) for s <= i <= m-1 ) and
///          nu(R, s) is the volume of s-dimensional ball with R radius.
/// @note log(nu(R, s)) = s log(R) + log(nu(1, s)) ( because nu(R, s) = R^{s} * nu(1, s) ) and
///       log(vol(\pi_{s}(L))) = sum( log(norm(b*i)) for s <= i <= m-1 ), then we have
///       log(H_l) = (l+1) log(R) + log(nu(1, (l+1))) - sum( log(norm(b*i)) for m-1-l <= i <= m-1 )
/// @param[in] R radius; if it is negative, then R is set to shortest norm;
/// @param[in] begin
/// @param[in] end
///
template<typename BasisFloat, typename GSFloat>
long double
Lattice<BasisFloat, GSFloat>::enumCost(
      double R,
      int begin,
      int end
      )
{
   if( R < 0 ) R = std::sqrt(B(0));
   if( end < 0 ) end = m - 1;
   long double N = 0.0;
   long double __logProjectedVolume = 0.0;
   long double logSphereVolume = 0.0;
   long double logR = std::log(R);
   for( int i = m-1; i > m-1-begin ; --i ){ __logProjectedVolume += std::log(B(i)) * 0.5; }
   for( int l = begin; l <= end; ++l )
   {
      __logProjectedVolume += std::log(B(m-l-1)) * 0.5;
      logSphereVolume = (l+1) * logR + ( (l+1)*std::log(M_PI)*0.5 - std::log(std::tgamma((l+1)*0.5+1.0)) );
      N += std::exp( logSphereVolume - __logProjectedVolume );
   }
   return N * 0.5;
}


///
/// @brief slope of GSA({(i, log2(squared_norm(b*i))} for 0 <= i <= h })
/// @param[in] h index of basis; if it is negative, h is set to m
///
template<typename BasisFloat, typename GSFloat>
double
Lattice<BasisFloat, GSFloat>::slopeGSA(
      int h
      )
{
   if( h < 0 ) h = m;
   VectorB<GSFloat> x(h);
   double slope = 0.0;
   for( int i = 0; i < h; ++i ) x(i) = log2(B(i)) * 0.5;
   for( int i = 0; i < h; ++i ) slope += (2.0*i-h+1)*x(i);
   slope = (6.0*slope) / (h*(h-1)*(h+1));
   return slope;
}


///
/// @brief generate new lattice using basis[:h]
/// @param[in] h index of basis
///
template<typename BasisFloat, typename GSFloat>
Lattice<BasisFloat, GSFloat>
Lattice<BasisFloat, GSFloat>::copy(
      int h
      )
{
   if( h < 0 ) h = m - 1;
   LatticeBasis<BasisFloat> basis__ = basis.block(0, 0, h, n);
   Lattice <BasisFloat, GSFloat> copiedLattice{basis__};
   copiedLattice.setConfig(config);
   return copiedLattice;
}


///
/// @brief convert basis to NTL::mat_ZZ
/// @output NTL::mat_ZZ matrix
///
template<typename BasisFloat, typename GSFloat>
NTL::mat_ZZ
Lattice<BasisFloat, GSFloat>::toNLTMat(
      )
{
   NTL::mat_ZZ Mbasis;
   Mbasis.SetDims(m, n);
   for( int i = 0; i < m; i++ )
   {
      for( int j = 0; j < n; j++ )
      {
         Mbasis[i][j] = basis.coeff(i, j);
      }
   }
   return Mbasis;
}


///
/// @brief load basis from NTL::mat_ZZ
///
template<typename BasisFloat, typename GSFloat>
void
Lattice<BasisFloat, GSFloat>::fromNLTMat(
      NTL::mat_ZZ &Mbasis
      )
{
   // remove zero vectors
   // row from 0 to zeroIndex are zerovectors
   int zeroIndex = -1;
   bool zero = true;
   int ntlRows = Mbasis.NumRows();
   int ntlCols = Mbasis.NumCols();
   for( int i = 0; i < ntlRows; i++ )
   {
      for( int j = 0; j < ntlCols; j++ )
         if( NTL::to_int(Mbasis[i][j]) != 0 ){ zero = false; }
      if( zero == true )
         zeroIndex++;
      else
         break; // row i is not a zero vector
   }

   // resize
   ntlRows -= zeroIndex + 1;

   // check
   if( (ntlRows != m) || (ntlCols != n) )
   {
      std::cout
         << "invarid dimension, expected (" << m << ", " << n << "),"
         << "but got (" << ntlRows << ", " << ntlCols << ")"
         << std::endl;
      exit(0);
   }

   for( int i = 0; i < m; i++ )
   {
      for( int j = 0; j < n; j++ )
      {
         basis.coeffRef(i, j)
            = NTL::to_int(Mbasis[i+zeroIndex+1][j]);
      }
   }

   setGSO();
}


///
/// @brief resize of lattice basis
/// @param[in] rows number of resized basis
///
template<typename BasisFloat, typename GSFloat>
void
Lattice<BasisFloat, GSFloat>::resize(
      int rows
      )
{
   assert( rows <= m );
   if( rows == m ){ return; }
   m = rows;
   // error occures when directory resize by basis = basis.topRows(m)
   LatticeBasis<BasisFloat> __basis = basis;
   MatrixMu<GSFloat>       __mu     = mu;
   VectorB<GSFloat>        __B      = B;
   basis = __basis.topRows(m);
   mu    = __mu.topLeftCorner(m, m);
   B     = __B.head(m);
   // GH = projectedGH();
   // logVolume = logProjectedVolume();
}


///
/// @brief randomize basis[begin:end] and execute reduction
/// @param[in] seed        seed of unimodular randomization
/// @param[in] begin       index of basis
/// @param[in] end         index of basis
///
template<typename BasisFloat, typename GSFloat>
bool
Lattice<BasisFloat, GSFloat>::randomize(
      int seed,
      int begin,
      int end
      )
{
   assert( end == -1 || begin <= end );
   if( begin == end ){ return true; }
   if( end < 0 ) end = m - 1;
   int u = end - begin + 1;

   // generate unimodular matrix
   MatrixXi U(u, u);
   genUnimat(U, seed, config.RandomizeType, config.RandomizeScale);
   // unimodular conversion
   basis.block(end-u+1,0,u,n) = U.cast<BasisFloat>()*(basis.block(end-u+1,0,u,n));
   setGSO();
   return true;
}


///
/// @brief MLLL for part of basis matrix {b_j,...,b_k-1,v,b_k,...,b_l}
/// @param[in] v inserted vector
/// @param[in] j index of basis
/// @param[in] k index of basis
/// @param[in] l index of basis
///
template<typename BasisFloat, typename GSFloat>
bool
Lattice<BasisFloat, GSFloat>::insertMlll(
      LatticeVector<BasisFloat> &v,
      int j,
      int k,
      int l
      )
{
   int i, ii;
   NTL::mat_ZZ Mbasis;
   Mbasis.SetDims(l-j+2, n);
   for( i = j; i < k; i++ )
      for( ii = 0; ii < n; ++ii ){ Mbasis[i-j][ii] = basis(i,ii); }
   for( ii = 0; ii < n; ++ii )   { Mbasis[k-j][ii] = v(ii); }
   for( i = k; i <= l; i++ )
      for( ii = 0; ii < n; ++ii ){ Mbasis[i+1-j][ii] = basis(i,ii); }

   NTL::LLL_FP(Mbasis,config.MLLLdelta,0,0,0);

   // Mbasis[0] is zero vector
   for( i = j; i <= l; i++ )
      for( ii = 0; ii < n; ++ii ){ basis(i,ii) = NTL::to_int(Mbasis[i+1-j][ii]); }

   setGSO();
   return true;
}


///
/// @brief merge other latice
/// @param[in] other merged lattice
/// @remark runningTime and mergeTime are updated in this function
///
template<typename BasisFloat, typename GSFloat>
bool
Lattice<BasisFloat, GSFloat>::merge(
      Lattice<BasisFloat, GSFloat> &other
      )
{
   assert( n == other.n );
   VectorB<GSFloat> preSqnorm = B;
   LatticeVector<BasisFloat> v;
   for( int i = 0; i < std::min(m, other.m); i++ )
   {
      if( basis.row(i) != other.basis.row(i) )
      {
         v = other.basis.row(i);
         insertMlll(v, 0, 0, m-1);
      }
   }
   if( B != preSqnorm )
      return true;
   else
      return false;
}


///
/// @brief compair two lattice-basis
/// @details two basis are compared by lexicographic order of (B(0), ..., B(m-1))
/// @param[in] other lattice
/// @return bool of ( self <= other )
/// @note If the number of rows of basis is different, the comparison is done according to the smaller number of rows.
///
template<typename BasisFloat, typename GSFloat>
bool
Lattice<BasisFloat, GSFloat>::isMoreReducedThan(
      Lattice &other
      )
{
   int h = std::min(m, other.m);
   for( int i = 0; i < h; ++i )
   {
      if( B(i) < other.B(i) - 1.0 )
      {
         return true;
      }
      if( B(i) > other.B(i) )
         return false;
   }
   return false;
}


///
/// @brief compair two lattice-basis
/// @details two basis are compared by lexicographic order of (B(0), ..., B(m-1))
/// @param[in] other lattice
/// @param[out] index minimum index {i; B(i) < C(i)} where C is other lattice basis
/// @return bool of ( self <= other )
/// @note If the number of rows of basis is different, the comparison is done according to the smaller number of rows.
///
template<typename BasisFloat, typename GSFloat>
bool
Lattice<BasisFloat, GSFloat>::isMoreReducedThan(
      Lattice &other,
      int &index
      )
{
   int h = std::min(m, other.m);
   for( int i = 0; i < h; ++i )
   {
      if( B(i) < other.B(i) - 1.0 )
      {
         index = i;
         return true;
      }
      if( B(i) > other.B(i) )
         return false;
   }
   return false;
}


///
/// @brief output basis
/// @param[in] os
///
template<typename BasisFloat, typename GSFloat>
void
Lattice<BasisFloat, GSFloat>::outputBasis(
      std::ostream *os
      )
{
   *os << toNLTMat();
}


///
/// @brief write basis into textfile
/// @param[in] filename write file path
///
template<typename BasisFloat, typename GSFloat>
void
Lattice<BasisFloat, GSFloat>::writeBasis(
      std::string filename
      )
{
   std::ofstream basefile(filename);
   outputBasis(&basefile);
}


///
/// @brief stringfy lattice
///
template<typename BasisFloat, typename GSFloat>
std::string
Lattice<BasisFloat, GSFloat>::toSimpleString(
      )
{
   std::ostringstream s;
   s << std::endl
     << "Shortest vector      : "
     << basis.row(0)             << std::endl
     << "Shortest norm        : "
     << shortestNorm()           << std::endl
     << "Approximate factor   : "
     << approxFactor()           << std::endl
     << "Hermite factor       : "
     << hermiteFactor()          << std::endl
     << "Root Hermite factor  : "
     << rootHermiteFactor()      << std::endl
     << "Gaussian Heuristics  : "
     << GH                       << std::endl;
   return s.str();
}


///
/// @brief stringfy lattice
///
template<typename BasisFloat, typename GSFloat>
std::string
Lattice<BasisFloat, GSFloat>::toSimpleString(
      LatticeVector<BasisFloat> &v
      )
{
   std::ostringstream s;
   double norm = std::sqrt(v.squaredNorm());
   s << std::endl
     << "Shortest vector      : "
     << v.transpose()            << std::endl
     << "Shortest norm        : "
     << norm                     << std::endl
     << "Approximate factor   : "
     << approxFactor(norm)       << std::endl
     << "Hermite factor       : "
     << hermiteFactor(norm)      << std::endl
     << "Root Hermite factor  : "
     << rootHermiteFactor(norm)  << std::endl
     << "Gaussian Heuristics  : "
     << GH                       << std::endl;
   return s.str();
}


///
/// instantiation
///
template class Lattice<int, double>;
template class Lattice<int, long double>;
template class Lattice<long int, double>;
template class Lattice<long int, long double>;


}  // namespace LapTools
