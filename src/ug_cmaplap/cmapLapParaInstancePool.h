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

/**@file    cmapLapParaInstancePool.h
 * @brief   Basis Pool.
 * @author  Nariaki Tateiwa, Yuji Shinano
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __CMAP_LAP__POOL_H__
#define __CMAP_LAP__POOL_H__

#include <set>
#include <map>
#include <queue>
#include <cfloat>
#include <cmath>
#include <memory>

#include <NTL/LLL.h>

#include "ug/paraDef.h"
#include "cmapLapParaDef.h"


namespace ParaCMapLAP
{


///
/// @class BasisElement
/// @brief element for instance pool
///
class BasisElement
{

   double   enumCost;                        ///< enumeration cost of lattice basis
   std::shared_ptr<LatticeBasis<int>> basis; ///< lattice basis

public:

   ///
   /// constrctor
   ///
   BasisElement() : enumCost(0.0), basis(nullptr) {}
   BasisElement(
         double inEnumCost,
         LatticeBasis<int> &inBasis
         ) : enumCost(inEnumCost)
   {
      basis = std::make_shared<LatticeBasis<int>>(inBasis);
   }
   BasisElement(
         BasisElement* elem
         )
   {
      enumCost = elem->enumCost;
      basis = std::make_shared<LatticeBasis<int>>(*(elem->basis));
   }

   ///
   /// copy constructor
   ///
   BasisElement(
         const BasisElement &basisElement
         ) : enumCost(basisElement.enumCost)
   {
      if( basisElement.basis )
         basis = std::make_shared<LatticeBasis<int>>(*(basisElement.basis));
      else
         basis = nullptr;
   }

   ///
   /// deconstrctor
   ///
   virtual ~BasisElement(
         )
   {
   }

   virtual BasisElement *clone(
         )
   {
      return new BasisElement(*this);
   }

   virtual double getEnumCost(
         )
   {
      return enumCost;
   }

   virtual std::shared_ptr<LatticeBasis<int>> getBasis(
         )
   {
      return basis;
   }

   virtual double getShortestSquaredNorm(
         )
   {
      return basis->row(0).squaredNorm();
   }

#ifdef UG_WITH_ZLIB

   ///
   /// write BasisElement in NTL format
   ///
   virtual void writeBasis(
         char *basisTextFileName    ///< file name for saving an incumbent basis
         )
   {
      std::ofstream basefile(basisTextFileName);
      if( basis )
      {
         NTL::Mat<int> B;
         B.SetDims(basis->rows(), basis->cols());
         for(int row = 0; row < basis->rows(); row++)
         {
            for(int col = 0; col < basis->cols(); col++)
            {
               B[row][col] = basis->coeff(row, col);
            }
         }
         basefile << B;
      }
      basefile.close();
   }

   ///
   /// write BasisElement in pool
   /// return number of solver has received this
   ///
   virtual bool write(
         gzstream::ogzstream &out
         )
   {
      out.write(reinterpret_cast<char *>(&enumCost), sizeof(double));

      int nRows = basis->rows();
      int nCols = basis->cols();
      out.write(reinterpret_cast<char *>(&nRows), sizeof(int));
      out.write(reinterpret_cast<char *>(&nCols), sizeof(int));

      std::vector<int> basisArray(nRows*nCols);
      Eigen::Map<LatticeBasis<int>>(&basisArray[0], nRows, nCols) = *basis;
      out.write(reinterpret_cast<char *>(&basisArray[0]), sizeof(int)*(nRows*nCols));
      return true;
   }

   ///
   /// read
   ///
   virtual bool read(
         UG::ParaComm *comm,
         gzstream::igzstream &in
         )
   {
      in.read(reinterpret_cast<char *>(&enumCost), sizeof(double));
      if( in.eof() ) return false;

      int nRows, nCols;
      in.read(reinterpret_cast<char *>(&nRows), sizeof(int));
      in.read(reinterpret_cast<char *>(&nCols), sizeof(int));

      std::vector<int> basisArray(nRows*nCols);
      in.read(reinterpret_cast<char *>(&basisArray[0]), sizeof(int)*(nRows*nCols));
      basis = std::make_shared<LatticeBasis<int>>(Eigen::Map<LatticeBasis<int>>(&basisArray[0], nRows, nCols));
      return true;
   }
#endif

   ///
   /// Stringfy Basis Element
   ///
   virtual std::string toString(
         )
   {
     std::ostringstream s;
     s << "enumCost: " << enumCost << std::endl;
     s << "basis: " << basis << std::endl;
     return s.str();
   }
};


#ifdef UG_WITH_ZLIB

///
/// class BasisElementQueue
/// @note this container is used when checkpointing
///
class BasisElementQueue
{
private:

   std::queue< std::shared_ptr<BasisElement> > basisElementQueue;

public:
   BasisElementQueue(){};
   virtual ~BasisElementQueue()
   {
   }

   virtual int size()
   {
      return basisElementQueue.size();
   }

   ///
   /// add basisElement in the queue
   ///
   virtual void push(
         std::shared_ptr<BasisElement> basisElement
         )
   {
      basisElementQueue.push(basisElement);
   }

   ///
   /// write basisElements to checkpoint file
   /// @return number of basisElements written
   ///
   virtual int writeBasisElementsToCheckpointFile(
         gzstream::ogzstream &out,  ///< gzstream for output
         int nWrite                 ///< -1: write all elements
         )
   {
      int n = 0;
      while( !basisElementQueue.empty() )
      {
         if( nWrite > 0 && n >= nWrite ) break;
         basisElementQueue.front()->write(out);
         basisElementQueue.pop();
         n++;
      }
      return n;
   }
};
#endif



class BasisElementSortCriterion
{
public:
   virtual bool operator()(const std::shared_ptr<BasisElement>& n1, const std::shared_ptr<BasisElement>& n2) const
   {
      return (n1->getEnumCost() < n2->getEnumCost());
   }
};

class InstancePool
{

protected:

   size_t        limitUsageOfPool;  ///< limit usage of this pool
   std::set< std::shared_ptr<BasisElement>, BasisElementSortCriterion > ascendingPool;
   std::shared_ptr<BasisElement> incumbentBasis;      ///< basis has shortest vector and minimun Enum Cost

public:
   ///
   /// constructor
   ///
   InstancePool(
         double limit
         ) : limitUsageOfPool(limit), incumbentBasis(nullptr)
   {}

   ///
   /// deconstructor
   ///
   virtual ~InstancePool(
   )
   {
   }

   ///
   /// insert Basis Element pointer
   ///
   virtual bool insert(
         std::shared_ptr<BasisElement> basisElement
         )
   {
      // set incumbentBasis
      if( !incumbentBasis ||
            EPSLT( basisElement->getShortestSquaredNorm(),
                   incumbentBasis->getShortestSquaredNorm(),
                   1e-5) ||
            ( EPSEQ( basisElement->getShortestSquaredNorm(),
                     incumbentBasis->getShortestSquaredNorm(),
                     1e-5) &&
              basisElement->getEnumCost() < incumbentBasis->getEnumCost()
            )
        )
      {
         incumbentBasis = basisElement;  // share the ownership of basis pointer
      }

      bool inserted = true;
      auto findBasis = ascendingPool.find(basisElement);
      if( findBasis != ascendingPool.end() )
      {
         // there exists a vector in the pool whose norm is equal to the new vector
         return false;
      }
      else
      {
         // the new vector is not in the pool
         ascendingPool.insert(basisElement);
         if( limitUsageOfPool > 0 && ascendingPool.size() > limitUsageOfPool )
         {
            auto i = std::prev(ascendingPool.end());
            ascendingPool.erase(i);
         }
      }
      return inserted;
   }

   ///
   /// insert Basis Element data
   ///
   virtual bool insert(
         double inEnumCost,
         LatticeBasis<int> &inBasis
         )
   {
      return insert(std::shared_ptr<BasisElement>(new BasisElement(inEnumCost, inBasis)));
   }

   virtual double getMinEnumCost(
         )
   {
      if( ascendingPool.size() > 0 )
      {
         return (*ascendingPool.begin())->getEnumCost();
      }
      else
      {
         return DBL_MAX;
      }
   }


   ///
   /// @return true if instance pool has no entry else false
   ///
   virtual bool empty(
         )
   {
      return ( ascendingPool.size() == 0 );
   }


   ///
   /// extract minimum ecnum cost basis
   /// This function should be called after checking pool is not empty.
   ///
   virtual std::shared_ptr<LatticeBasis<int>> extractMinEnumCostBasis(
         )
   {
      assert( ascendingPool.size() > 0 );
      auto basis = (*ascendingPool.begin())->getBasis();
      ascendingPool.erase(ascendingPool.begin());
      return basis;
   }

   virtual LatticeBasis<int> getBasis(
         int index
         )
   {
      assert( static_cast<int>(ascendingPool.size()) > index );
      std::multimap<std::shared_ptr<BasisElement>, std::shared_ptr<BasisElement>, BasisElementSortCriterion >::iterator p;
      for( int i = 0; i < index; i ++ ) ++p;
      return LatticeBasis<int>(*p->second->getBasis());
   }

   virtual std::shared_ptr<BasisElement> getIncumbentBasis(
         )
   {
      return incumbentBasis;
   }

   virtual size_t size(
         )
   {
      return ascendingPool.size();
   }

#ifdef UG_WITH_ZLIB

   ///
   /// get copied basis elements
   /// @return the BasisElementQueue pointer
   ///
   virtual BasisElementQueue *getBasisElementQueue(
         UG::ParaComm *comm,     ///< communicator used
         int nWrite=-1           ///< -1: write all elements
         )
   {
      int n = 0;
      BasisElementQueue *basisElementQueue = new BasisElementQueue();
      for( const auto &p : ascendingPool )
      {
         if( nWrite > 0 && n >= nWrite ) break;
         basisElementQueue->push(p);
         n++;
      }
      return basisElementQueue;
   }

   ///
   /// write incumbent basis
   ///
   virtual void writeIncumbentBasis(
         char *basisTextFileName    ///< file name for saving an incumbent basis
         )
   {
      if( incumbentBasis )
      {
         incumbentBasis->writeBasis(basisTextFileName);
      }
   }

   ///
   /// write BasisElements in pool
   /// return the number of basis elements witten
   ///
   virtual int writeBasisElementsToCheckpointFile(
         gzstream::ogzstream &out,
         int nWrite                 ///< -1: write all elements
         )
   {
      int n = 0;
      for( const auto &p : ascendingPool )
      {
         if( nWrite > 0 && n >= nWrite) break;
         p->write(out);
         n++;
      }
      return n;
   }

   ///
   /// read BasisElements from Checkpoint file
   /// return the number of basis elements inserted
   ///
   virtual int readBasisElementsFromCheckpointFile(
         UG::ParaComm *comm,
         gzstream::igzstream &in
         )
   {
      int n = 1;  // It has an instance added in constructer of LoadCoordinator
      auto basisElement = std::make_shared<BasisElement>();
      while( basisElement->read(comm,in) )
      {
         if( insert(basisElement) ){ n++; }
         basisElement = std::make_shared<BasisElement>();
      }
      assert( n == static_cast<int>(size()) );
      return n;
   }
#endif

   ///
   /// Stringfy InstancePool
   ///
   virtual const std::string toString(
         )
   {
      std::ostringstream s;
      s << std::endl << std::endl;
      s << "===== Basis Pool (Size " << ascendingPool.size() << ") =====" << std::endl;
      for( const auto &p : ascendingPool )
      {
         s << p->toString();
         s << std::endl;
      }
      s << "===========================" << std::endl;
      return s.str();
   }
};

} // namespace ParaCMapLAP

#endif // __CMAP_LAP__POOL_H__
