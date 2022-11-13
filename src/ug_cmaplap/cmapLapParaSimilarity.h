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

/**@file    cmapLapParaSimilarity.h
 * @brief   Functions to calculate similarity of basis set.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_PARA_SIMILARITY_H__
#define __CMAP_LAP_PARA_SIMILARITY_H__

#include <vector>
#include <chrono>
#include <random>
#include <deque>
#include <set>
#include <cassert>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/SVD>


namespace ParaCMapLAP
{

namespace BasisSimilarity
{

/// index to calculate grassmann metric
std::vector<int> grassmannIndexes;


/// calculate grassmann metric
Eigen::MatrixXd scores_geodesic_metric;
Eigen::MatrixXd scores_chordal_metric;
Eigen::MatrixXd scores_fubini_study_metric;
Eigen::MatrixXd scores_chordal_2norm_metric;
Eigen::MatrixXd scores_chordal_fnorm_metric;
Eigen::MatrixXd scores_projection_2norm_metric;
Eigen::MatrixXd scores_max_metric;
Eigen::MatrixXd scores_mean_metric;
Eigen::VectorXd scores_num_upper_duplicate;
Eigen::VectorXd scores_num_all_duplicate;


///
/// @brief setter of grassmannIndexes
/// @param[in] inGrassmannIndexes
///
void setGrassmannIndexes(
      std::vector<int> inGrassmannIndexes
      )
{
   grassmannIndexes = inGrassmannIndexes;
}


///
/// @brief similarity log header
/// @param[in] n dimension
/// @return std::string
///
std::string outputSimilarityOfBasisHeader(
      int n
      )
{
   assert( grassmannIndexes.size() > 0 );

   std::ostringstream s;
   s  << "SimilarityStatus"
      << ",time"
      << ",numPairs";
   std::vector<std::string> metrics{"num_upper_duplicate", "num_all_duplicate"};
   std::vector<std::string> grassmannMetricNames{
      "grassmann_geodesic",
      "grassmann_chordal",
      "grassmann_fubini",
      "grassmann_chordal_2norm",
      "grassmann_chordal_fnorm",
      "grassmann_projection_2norm",
      "grassmann_max",
      "grassmann_mean"
   };
   for( auto name : grassmannMetricNames )
   {
      for( auto i : grassmannIndexes )
      {
         std::ostringstream metric;
         metric << name << "_" << i;
         metrics.push_back(metric.str());
      }
   }
   for( auto metric : metrics )
   {
      s  << "," << metric << "_min";
      s  << "," << metric << "_max";
      s  << "," << metric << "_average";
   }
   return s.str();
}


///
/// @brief sampling of vector
/// @tparam T population type
/// @param[in] population list
/// @param[out] sampled list
/// @param[in] sampled size
///
template<typename T> void sample(
      T &list,
      T &sampled_list,
      int size
      )
{
   assert( static_cast<int>(list.size()) >= size );

   std::vector<int> indexes(list.size(), 0);
   for ( size_t i = 0; i < list.size(); ++i )
      indexes[i] = i;

   std::random_device seed_gen;
   // std::mt19937 engine {seed_gen()};
   std::mt19937 engine {0};
   std::shuffle(indexes.begin(), indexes.end(), engine);

   sampled_list.resize(size);
   for ( int i = 0; i < size; ++i )
      sampled_list[i] = list[indexes[i]];
}


///
/// @brief sign function
/// @param[in] x
/// @return sign of x ( 1 or -1 )
///
int sign(
      int x
      )
{
   if( x >= 0 ) return 1;
   else         return -1;
}


///
/// @brief max(k; A(i) == B(i) || A(i) == -B(i) for all i <= k )
/// @param[in] basisA basis
/// @param[in] basisB basis
/// @return number of matches from the top of the basis vector
///
double num_upper_duplicate(
      LatticeBasis<int>& basisA,
      LatticeBasis<int>& basisB
      )
{
   double score = 0.0;
   for( int i = 0; i < basisA.rows(); i++ )
   {
      if( sign(basisA.coeff(i, 0))*basisA.row(i)
            == sign(basisB.coeff(i, 0))*basisB.row(i) )
         score += 1;
      else
         break;
   }
   return score;
}


///
/// @brief coun ( ( A(i) == B(i) || A(i) == -B(i) ) for all i )
/// @param[in] basisA basis
/// @param[in] basisB basis
/// @return number of basis vector overlaps
///
double num_all_duplicate(
      LatticeBasis<int>& basisA,
      LatticeBasis<int>& basisB
      )
{
   double score = 0.0;
   for( int i = 0; i < basisA.rows(); i++ )
   {
      if( sign(basisA.coeff(i, 0))*basisA.row(i)
            == sign(basisB.coeff(i, 0))*basisB.row(i) )
         score += 1;
   }
   return score;
}


/// @brief set Gram-Schmidt matrix
/// @parma[in] basis
/// @parma[out] GSO Gram-Schmidt of basis
///
void setGramSchmidt(
      LatticeBasis<int>& basis,
      Eigen::MatrixXd& GSO
      )
{
   GSO = LatticeBasis<int>(basis).cast<double>()
      .transpose()
      .householderQr()
      .householderQ()
      .transpose();
}

///
/// @brief Grassmann Geodesic Metric
/// @param[in] cc canonical correlations
/// @details sqrt(sum(theta_i^2 for 0 <= i < m))
///          = sqrt(sum(acos(cc(i))^2 for 0 <= i < m))
///
double grassmann_geodesic_metric(
      Eigen::VectorXd &cc
      )
{
   int m = cc.size();
   double d = 0.0;
   double theta = 0.0;
   for( int i = 0; i < m; ++i )
   {
      theta = std::acos(cc(i));
      d += theta * theta;
   }
   return std::sqrt(d);
}


///
/// @brief Grassmann Projected Metric
/// @param[in] cc canonicalCorrelations
/// @details 2^(-1/2) frobenius_norm(GSOA*GSOA^T - GSOB*GSOB^T)
///          = sqrt(m - sum(cos^2(theta_i) for 0 <= i < m))
///          = sqrt(m - sum(cc(i)^2) for 0 <= i <= m)
///
double grassmann_chordal_metric(
      Eigen::VectorXd &cc
      )
{
   int m = cc.size();
   double d = m;
   for( int i = 0; i < m; i++ )
      d -= cc(i) * cc(i);
   return std::sqrt(d);
}


///
/// @brief Grassmann Fubini-Study Metric
/// @param[in] cc canonical correlations
/// @details acos(prod(cos(theta_i) for 0 <= i < m))
///          = acos(prod(cc(i) for 0 <= i < m))
///
double grassmann_fubini_study_metric(
      Eigen::VectorXd &cc
      )
{
   return std::acos( cc.prod() );
}


///
/// @brief Grassmann Chordal 2-norm Metric
/// @param[in] cc canonical correlations
/// @details 2norm( U*GSOA^T - V*GSOB )
///          = infinity_norm(2*sin(theta_i/2) for 0 <= i < m)
///          = 2 * max(sin(theta_i/2) for 0 <= i < m)
///          = 2 * sqrt( max(sin^2(theta_i/2) for 0 <= i < m) ) (because 0 < theta_i < pi)
///          = 2 * sqrt( max((1-cc(i))/2 for 0 <= i < m) )
///
double grassmann_chordal_2norm_metric(
      Eigen::VectorXd &cc
      )
{
   int m = cc.size();
   double d = (1.0 - cc(0)) / 2.0;
   for( int i = 1; i < m; ++i )
   {
      d = std::max(d, (1.0 - cc(i))/2.0);
   }
   return 2.0 * std::sqrt(d);
}


///
/// @brief Grassmann Chordal frobenius-norm Metric
/// @param[in] cc canonical correlations
/// @details frobenius_norm( U*GSOA^T - V*GSOB )
///          = 2_norm(2*sin(theta_i/2) for 0 <= i < m)
///          = 2 * sqrt(sum(sin^2(theta_i/2) for 0 <= i < m)
///          = 2 * sqrt(sum((1-cc(i))/2 for 0 <= i < m)
///
double grassmann_chordal_fnorm_metric(
      Eigen::VectorXd &cc
      )
{
   int m = cc.size();
   double d = 0.0;
   for( int i = 0; i < m; ++i )
      d += (1.0 - cc(i)) / 2.0;
   return 2.0 * std::sqrt(d);
}


///
/// @brief Grassmann Projection 2-norm Metric
/// @param[in] cc canonicalCorrelations
/// @details 2norm( GSOA*GSOA^T - GSOB*GSOB^T )
///          infinity_norm( sin(theta_i) ) for 0 <= i < m) )
///          = sqrt( max( sin^2(theta_i) for 0 <= i < m) ) (because 0 <= theta_i <= pi)
///          = sqrt( max( 1 - cos^2(theta_i) for 0 <= i < m) )
///          = sqrt( max( 1 - cc(i)^2 for 0 <= i < m) )
///          = sqrt( 1 - cc(0)^2 )
///
double grassmann_projection_2norm_metric(
      Eigen::VectorXd &cc
      )
{
   // int m = cc.size();
   double d = 1.0 - cc(0) * cc(0);
   if( d < 1.0e-5 ) return 0;
   // for( int i = 1; i < m; ++i )
   // {
   //    d = std::max(d, 1.0 - cc(i) * cc(i));
   // }
   return std::sqrt(d);
}


///
/// @brief Grassmann max Metric
/// @param[in] cc canonicalCorrelations
/// @details min( sin(theta_i) ) for 0 <= i < m) )
///          = sqrt( min( sin^2(theta_i) for 0 <= i < m) ) (because 0 <= theta_i <= pi)
///          = sqrt( min( 1 - cos^2(theta_i) for 0 <= i < m) )
///          = sqrt( min( 1 - cc(i)^2 for 0 <= i < m) )
///          = sqrt( 1 - cc(m-1)^2 )
///
double grassmann_max_metric(
      Eigen::VectorXd &cc
      )
{
   int m = cc.size();
   double d = 1.0 - cc(m-1) * cc(m-1);
   if( d < 1.0e-5 ) return 0;
   return std::sqrt(d);
}


///
/// @brief Grassmann Mean Metric
/// @param[in] cc canonicalCorrelations
/// @details mean( sin(\theta_i)^2 for 0 <= i <= m )
///         = mean( 1 - cos(\theta_i)^2 for 0 <= i <= m )
///
double grassmann_mean_metric(
      Eigen::VectorXd &cc
      )
{
   int m = cc.size();
   double d = 0.0;
   for( int i = 0; i < m; i++ )
      d += 1.0 - cc(i)*cc(i);
   return d / static_cast<double>(m);
}


///
/// @brief output log
/// @param[in] basisDeque container of basis
/// @param[in] time
/// @param[in] n dimension of basis
/// @param[in] nSamples number of samples for calculation of basis similarity
/// @param[in] num_threads number threads
/// @param[in] verbose ( 0: not, 1: light, 2: medium, 3: heave )
///
std::string outputSimilarityOfBasis(
      std::deque<std::shared_ptr<LatticeBasis<int>>> &basisDeque,
      double time,
      int n,
      int nSamples=-1,
      int numThreads=1,
      int verbose=0
      )
{

   auto start = std::chrono::system_clock::now();
   auto elapsed = [&start](
         )
   {
      return std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now()-start
            ).count() / 1000.0;
   };

   int numBasis = basisDeque.size();
   int numPairs = numBasis * (numBasis-1) / 2;
   int k = -1;
   double tol = 1e-6;

   std::vector<std::pair<int, int>> combinations;
   combinations.resize(numPairs);
   k = -1;
   for ( int i = 0; i < numBasis; ++i )
   {
      for ( int j = i+1; j < numBasis; ++j )
      {
         combinations[++k] = std::make_pair(i, j);
      }
   }

   // sampling
   if( nSamples == -1 ){ nSamples = combinations.size(); }
   nSamples = std::min(nSamples, static_cast<int>(combinations.size()));
   decltype(combinations) sampledConbinations;
   sample(
         combinations,
         sampledConbinations,
         nSamples
         );

   // calculate GSO matrix
   std::vector<Eigen::MatrixXd> GSOList;
   GSOList.resize(basisDeque.size());
   std::set<int> calcGSOindexes;
   for( auto pair : sampledConbinations )
   {
      calcGSOindexes.insert(pair.first);
      calcGSOindexes.insert(pair.second);
   }
   for( auto _k : calcGSOindexes )
   {
      Eigen::MatrixXd GSO;
      setGramSchmidt(*basisDeque.at(_k), GSO);
      GSOList[_k] = GSO;
   }

   // calculate grassmann metric
   scores_geodesic_metric        .resize(nSamples, grassmannIndexes.size());
   scores_chordal_metric         .resize(nSamples, grassmannIndexes.size());
   scores_fubini_study_metric    .resize(nSamples, grassmannIndexes.size());
   scores_chordal_2norm_metric   .resize(nSamples, grassmannIndexes.size());
   scores_chordal_fnorm_metric   .resize(nSamples, grassmannIndexes.size());
   scores_projection_2norm_metric.resize(nSamples, grassmannIndexes.size());
   scores_max_metric             .resize(nSamples, grassmannIndexes.size());
   scores_mean_metric            .resize(nSamples, grassmannIndexes.size());
   scores_num_upper_duplicate    .resize(nSamples);
   scores_num_all_duplicate      .resize(nSamples);


#pragma omp parallel for num_threads(numThreads) schedule(static)
   for( k = 0; k < nSamples; ++k )
   {
      int i, j;
      std::tie(i, j) = sampledConbinations[k];

      // duplicate
      scores_num_upper_duplicate(k) = num_upper_duplicate(
               *(basisDeque.at(i)), *(basisDeque.at(j))
               );
      scores_num_all_duplicate(k) = num_all_duplicate(
               *(basisDeque.at(i)), *(basisDeque.at(j))
               );

      // grassmann
      // 1. generate sub GSO
      // 2. get canonical angles
      // 3. calculate metric
      for( int l = grassmannIndexes.size()-1; l > -1; --l )
      {
         int d = grassmannIndexes[l];
         Eigen::VectorXd canonicalCorrelations = Eigen::VectorXd::Ones(n-d);
         if( d >= n / 2.0 )
         {
            // GSO = [b*0; b*1; ...; b*n-1] -> [b*d; ...; b*(n-1)]
            Eigen::MatrixXd subGSOA{GSOList.at(i).block(d,0,n-d,n)};
            Eigen::MatrixXd subGSOB{GSOList.at(j).block(d,0,n-d,n)};
            // calculate canonical angles
            Eigen::JacobiSVD<Eigen::MatrixXd> SVD{subGSOA * subGSOB.transpose()};
            canonicalCorrelations.tail(n-d) = SVD.singularValues();
         }
         else if( d > 0 )
         {
            // GSO = [b*0; b*1; ...; b*n-1] -> [b*0; ...; b*d-1]
            Eigen::MatrixXd subGSOA{GSOList.at(i).block(0,0,d,n)};
            Eigen::MatrixXd subGSOB{GSOList.at(j).block(0,0,d,n)};
            // calculate canonical angles
            Eigen::JacobiSVD<Eigen::MatrixXd> SVD{subGSOA * subGSOB.transpose()};
            canonicalCorrelations.tail(d) = SVD.singularValues();
         }
         for( int ii = 0; ii < n-d; ++ii )
         {
            canonicalCorrelations(ii) = std::min(std::max(canonicalCorrelations(ii), -1.0), 1.0);
         }
         if( verbose > 0 )
         {
            std::cout << "k: " << k << " pair(" << i << ", " << j << ") d " << d << std::endl;
            std::cout << "canonicalCorrelations " << canonicalCorrelations.transpose() << std::endl;
         }
         if( (canonicalCorrelations.array() >= 1.0-tol).all() && (canonicalCorrelations.array() <= 1.0+tol).all() )
         {
            // GSOA == GSOB
            scores_geodesic_metric(k, l)           = 0.0;
            scores_chordal_metric(k, l)            = 0.0;
            scores_fubini_study_metric(k, l)       = 0.0;
            scores_chordal_2norm_metric(k, l)      = 0.0;
            scores_chordal_fnorm_metric(k, l)      = 0.0;
            scores_projection_2norm_metric(k, l)   = 0.0;
            scores_max_metric(k, l)                = 0.0;
            scores_mean_metric(k, l)               = 0.0;
         }
         else
         {
            // calculate metric
            scores_geodesic_metric(k, l)
               = grassmann_geodesic_metric(canonicalCorrelations);
            scores_chordal_metric(k, l)
               = grassmann_chordal_metric(canonicalCorrelations);
            scores_fubini_study_metric(k, l)
               = grassmann_fubini_study_metric(canonicalCorrelations);
            scores_chordal_2norm_metric(k, l)
               = grassmann_chordal_2norm_metric(canonicalCorrelations);
            scores_chordal_fnorm_metric(k, l)
               = grassmann_chordal_fnorm_metric(canonicalCorrelations);
            scores_projection_2norm_metric(k, l)
               = grassmann_projection_2norm_metric(canonicalCorrelations);
            scores_max_metric(k, l)
               = grassmann_max_metric(canonicalCorrelations);
            scores_mean_metric(k, l)
               = grassmann_mean_metric(canonicalCorrelations);
         }
      }
   }


   std::ostringstream s;
   s  << "SimilarityStatus"
      << "," << time
      << "," << scores_num_upper_duplicate.size();
      // upper_duplicate
   s  << "," << scores_num_upper_duplicate.minCoeff();
   s  << "," << scores_num_upper_duplicate.maxCoeff();
   s  << "," << scores_num_upper_duplicate.mean();
      // all_duplicate
   s  << "," << scores_num_all_duplicate.minCoeff();
   s  << "," << scores_num_all_duplicate.maxCoeff();
   s  << "," << scores_num_all_duplicate.mean();
   for( size_t i = 0; i < grassmannIndexes.size(); ++i )
   {
      // grassmann_geodesic_metric
      s << "," << scores_geodesic_metric.col(i).minCoeff();
      s << "," << scores_geodesic_metric.col(i).maxCoeff();
      s << "," << scores_geodesic_metric.col(i).mean();
   }
   for( size_t i = 0; i < grassmannIndexes.size(); ++i )
   {
      // grassmann_chordal_metric
      s << "," << scores_chordal_metric.col(i).minCoeff();
      s << "," << scores_chordal_metric.col(i).maxCoeff();
      s << "," << scores_chordal_metric.col(i).mean();
   }
   for( size_t i = 0; i < grassmannIndexes.size(); ++i )
   {
      // grassmann_fubini_study_metric
      s << "," << scores_fubini_study_metric.col(i).minCoeff();
      s << "," << scores_fubini_study_metric.col(i).maxCoeff();
      s << "," << scores_fubini_study_metric.col(i).mean();
   }
   for( size_t i = 0; i < grassmannIndexes.size(); ++i )
   {
      // grassmann_chordal_2norm
      s << "," << scores_chordal_2norm_metric.col(i).minCoeff();
      s << "," << scores_chordal_2norm_metric.col(i).maxCoeff();
      s << "," << scores_chordal_2norm_metric.col(i).mean();
   }
   for( size_t i = 0; i < grassmannIndexes.size(); ++i )
   {
      // grassmann_chordal_fnorm
      s << "," << scores_chordal_fnorm_metric.col(i).minCoeff();
      s << "," << scores_chordal_fnorm_metric.col(i).maxCoeff();
      s << "," << scores_chordal_fnorm_metric.col(i).mean();
   }
   for( size_t i = 0; i < grassmannIndexes.size(); ++i )
   {
      // grassmann_projection_2norm
      s << "," << scores_projection_2norm_metric.col(i).minCoeff();
      s << "," << scores_projection_2norm_metric.col(i).maxCoeff();
      s << "," << scores_projection_2norm_metric.col(i).mean();
   }
   for( size_t i = 0; i < grassmannIndexes.size(); ++i )
   {
      // grassmann_max_norm
      s << "," << scores_max_metric.col(i).minCoeff();
      s << "," << scores_max_metric.col(i).maxCoeff();
      s << "," << scores_max_metric.col(i).mean();
   }
   for( size_t i = 0; i < grassmannIndexes.size(); ++i )
   {
      // grassmann_mean_norm
      s << "," << scores_mean_metric.col(i).minCoeff();
      s << "," << scores_mean_metric.col(i).maxCoeff();
      s << "," << scores_mean_metric.col(i).mean();
   }

   // std::cout << "\r simlarity total " << elapsed() << " sec" << std::endl;
   return s.str();
}

}  // namespace BasisSimilarity


}  // namespace ParaCMapLAP


#endif // __CMAP_LAP_PARA_SIMILARITY_H__
