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

/**@file    cmapLapParaShareData.h
 * @brief   Base class for shared data amoung solvers.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP_SHARE_DATA_H__
#define __CMAP_LAP_SHARE_DATA_H__

#include <iostream>
#include <unordered_set>
#include "ug/gzstream.h"
#include "ug/paraParamSet.h"
#include "cmapLapParaDef.h"


namespace ParaCMapLAP
{

///
/// @class VectorElementBase
/// @brief Base class of ShareDataPool
///
template<typename BasisFloat=int>
class VectorElementBase
{


public:

   LatticeVector<BasisFloat> vector;  ///< vector (first coefficients is not negative)
   double __squaredNorm;              ///< squared norm of vector

   ///
   /// @brief constrctor
   ///
   VectorElementBase(){}
   VectorElementBase(
         LatticeVector<BasisFloat> &inVector
         )
   {
      init(inVector);
   }
   VectorElementBase(VectorElementBase<BasisFloat> *ve)
   {
      vector = LatticeVector<BasisFloat>(ve->vector);
      setSquaredNorm();
   }

   ///
   /// @brief initialize
   ///
   virtual void init(
      LatticeVector<BasisFloat> &inVector
      )
   {
      vector = LatticeVector<BasisFloat>(inVector);
      if( vector.coeff(0) < 0 ) vector = -vector;
      setSquaredNorm();
   }


   ///
   /// @brief deconstructor
   ///
   virtual ~VectorElementBase(){}


   ///
   /// @return size of vector
   ///
   virtual int size(
         )
   {
      return vector.size();
   }


   ///
   /// @return norm of vector
   ///
   virtual double norm(
         )
   {
      return std::sqrt(__squaredNorm);
   }


   ///
   /// @brief set squaredNorm from vector member variable
   ///
   virtual void setSquaredNorm(
         )
   {
      __squaredNorm = 0.0;
      for( int i = vector.size() - 1; i > -1; --i )
         __squaredNorm += vector.coeff(i) * vector.coeff(i);
   }

   ///
   /// @brief setter of squaredNorm
   /// @param[in] inSquaredNorm squaredNorm of vector member variable
   ///
   virtual void setSquaredNorm(
         double inSquaredNorm
         )
   {
      __squaredNorm = inSquaredNorm;
   }

   ///
   /// @brief setter of vector
   /// @param[in] vector lattice vector
   ///
   virtual void setVector(
         LatticeVector<BasisFloat> &inVector
         )
   {
      vector = inVector;
   }


   ///
   /// @breif getter of vector
   ///
   virtual LatticeVector<BasisFloat> & getVector(
         )
   {
      return vector;
   }

   ///
   /// @breif getter of squaredNorm
   ///
   virtual double & squaredNorm(
         )
   {
      return __squaredNorm;
   }


   ///
   /// @brief Stringfy Vector Element
   ///
   virtual std::string toString(
         )
   {
      std::ostringstream s;
      s << "vector: " << vector.transpose() << std::endl;
      s << "squaredNorm: " << __squaredNorm << std::endl;
      s << "norm: " << std::sqrt(__squaredNorm) << std::endl;
      s << std::endl;
      return s.str();
   }

   virtual std::string toSimpleString(
         )
   {
      std::ostringstream s;
      s << "vector: " << vector.head(5).transpose() << "..." << std::endl;
      s << "squaredNorm: " << __squaredNorm << std::endl;
      s << "norm: " << std::sqrt(__squaredNorm) << std::endl;
      s << std::endl;
      return s.str();
   }


}; // class VectorElementBase


///
/// @class VectorElement
/// @brief element of ShareDataPool
///
class VectorElement : public VectorElementBase<int>
{

public:

   SolverType solverType;        ///< solver type that finds this vector
   unsigned int solverThreadId;  ///< solver-ids that created this object
   std::unordered_set<unsigned int> solverThreadIds;  ///< the set of solver-ids already sent

   ///
   /// constrctor
   ///
   VectorElement(
         )
      :
         VectorElementBase<int>(),
         solverType(Undefined),
         solverThreadId(0)
   {}
   VectorElement(
         LatticeVector<int> &inVector,
         unsigned int inSolverThreadId,
         SolverType inSolverType=Undefined
         )
      :
         VectorElementBase<int>(inVector),
         solverType(inSolverType),
         solverThreadId(inSolverThreadId)
   {
      solverThreadIds.insert(solverThreadId);
   }
   VectorElement(
         VectorElement *ve
         )
      :
         VectorElementBase<int>(ve->vector),
         solverType(ve->solverType),
         solverThreadId(ve->solverThreadId),
         solverThreadIds(std::unordered_set<unsigned int>(ve->solverThreadIds))
   {
   }

   ///
   /// deconstructor
   ///
   virtual ~VectorElement(){}

   ///
   /// @brief getter of solverThreadId
   /// @return solverThreadId
   ///
   virtual unsigned int getSolverThreadId(
         )
   {
      return solverThreadId;
   }

   ///
   /// check this object already sent to a thread
   /// @return true: already sent, false: otherwise
   ///
   virtual bool isInSentSolverThreadId(
         unsigned int inSolverThreadId
         )
   {
      return( solverThreadIds.find(inSolverThreadId) != solverThreadIds.end() );
   }

   ///
   /// add thread id to id-set
   ///
   virtual void addSentSolverThreadId(
         unsigned int inSolverThreadId
         )
   {
      solverThreadIds.insert(inSolverThreadId);
   }

   ///
   /// erase thread id from id-set
   ///
   virtual void eraseSolverThreadId(
         unsigned int inSolverThreadId
         )
   {
      solverThreadIds.erase(inSolverThreadId);
   }

   ///
   /// clear thread id
   ///
   virtual void clearSolverThreadIds(
         )
   {
      solverThreadIds.clear();
   }

   ///
   /// setter of SolverType
   ///
   virtual void setSolverType(
         SolverType inSolverType
         )
   {
      solverType = inSolverType;
   }

   ///
   /// getter of SolverType
   /// @return solverType
   ///
   virtual SolverType getSolverType(
         )
   {
      return solverType;
   }


#ifdef UG_WITH_ZLIB
   ///
   /// write VectorElement in pool
   /// return number of solver has received this
   ///
   virtual int write(
         gzstream::ogzstream &out
         )
   {
      int dimension = this->vector.size();
      out.write(reinterpret_cast<char *>(&dimension), sizeof(int));

      std::vector<int> vec(dimension);
      Eigen::Map<LatticeVector<int>>(&vec[0], dimension) = this->vector;
      out.write(reinterpret_cast<char *>(&vec[0]), sizeof(int)*(dimension));

      out.write(reinterpret_cast<char *>(&this->squaredNorm()), sizeof(double));
      out.write(reinterpret_cast<char *>(&solverThreadId), sizeof(unsigned int));

      int sizeOfSolverThreadIds = solverThreadIds.size();
      out.write(reinterpret_cast<char *>(&sizeOfSolverThreadIds), sizeof(int));

      std::vector<int> ids(sizeOfSolverThreadIds);
      int k = 0;
      for( auto p = solverThreadIds.begin(); p != solverThreadIds.end() ; ++p )
      {
         ids[k++] = *p;
      }
      out.write(reinterpret_cast<char *>(&ids[0]), sizeof(int)*sizeOfSolverThreadIds);
      return sizeOfSolverThreadIds;
   }

   ///
   /// read
   ///
   virtual bool read(
         UG::ParaComm *comm,
         gzstream::igzstream &in
         )
   {
      int dimension = 0;
      in.read(reinterpret_cast<char *>(&dimension), sizeof(int));
      if( in.eof() ) return false;

      std::vector<int> vec(dimension);
      in.read(reinterpret_cast<char *>(&vec[0]), sizeof(int)*(dimension));
      this->vector = Eigen::Map<LatticeVector<int>>(&vec[0], dimension);

      in.read(reinterpret_cast<char *>(&this->squaredNorm()), sizeof(double));
      in.read(reinterpret_cast<char *>(&solverThreadId), sizeof(unsigned int));

      int sizeOfSolverThreadIds = 0;
      in.read(reinterpret_cast<char *>(&sizeOfSolverThreadIds), sizeof(int));

      std::vector<int> ids(sizeOfSolverThreadIds);
      in.read(reinterpret_cast<char *>(&ids[0]), sizeof(int)*sizeOfSolverThreadIds);
      for( int k = 0; k < sizeOfSolverThreadIds; k++ )
      {
         solverThreadIds.insert(ids[k]);
      }
      return true;
   }

#endif // WITH_UG_ZLIB


   ///
   /// Stringfy VectorElement
   ///
   virtual std::string toString(
         )
   {
      std::ostringstream s;

      s << "vector: "         << this->vector.transpose() << std::endl;
      s << "squaredNorm: "    << this->squaredNorm()      << std::endl;
      s << "solverThreadId: " << solverThreadId           << std::endl;
      s << "solverThreadIds: ";
      for( auto it = solverThreadIds.begin(); it != solverThreadIds.end(); it++ )
      {
          s << *it << " ";
      }
      s << std::endl;
      return s.str();
   }


   ///
   /// Stringfy simple VectorElement
   ///
   virtual std::string toSimpleString(
         )
   {
      std::ostringstream s;
      s << "vector: "         << this->vector.head(5).transpose() << "..." << std::endl;
      s << "squaredNorm: "    << this->squaredNorm() << std::endl;
      s << "solverThreadId: " << solverThreadId << std::endl;
      s << "solverThreadIds: ";
      for( auto it = solverThreadIds.begin(); it != solverThreadIds.end(); it++ )
      {
         s << *it << " ";
      }
      s << std::endl;
      return s.str();
   }
};


}  // namespace ParaCMapLAP


#endif // __CMAP_LAP_SHARE_DATA_H__
