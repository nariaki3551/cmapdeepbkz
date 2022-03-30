#include "cmapTestGaussSieve.h"
#include <iostream>
#include "Def.h"
#include "Timer.h"
#include "cmapLapParaDef.h"
#include "cmapLapParaSolver.h"


///
/// @brief reduce v by w such that |w| < |v|
/// @param[in] v lattice vector which will be reduced
///
template<typename BasisFloat, typename GSFloat>
void
CmapGaussSieve<BasisFloat, GSFloat>::postProcessOfReduce(
      VectorElementPtr v
      )
{
   if( SendStack.simpledInsertableTest(v) )
   {
      SendStack.insert(v);
      if( verbose > 2 )
      {
         std::cout
            << "sampling vector "
            << "[ #L " << this->List.size()
            << "; norm " << std::sqrt(v->squaredNorm())
            << "; approx " << L->approxFactor(v->norm()) << " ] "
            // << sendVector.transpose()
            << std::endl;
      }
      if( verbose > 2 && SendStack.size() % 1000 == 1 )
      {
         std::cout
            << "SENDSTACK," << LapTools::Timer::getElapsedTime()
            << "," << this->List.size()
            << SendStack.toStatString(L->GH)
            << std::endl;
      }
   }
}



///
/// @brief communicate with LC
/// @param[out] shouldAbort true if it should abort else false
///
template<typename BasisFloat, typename GSFloat>
bool
CmapGaussSieve<BasisFloat, GSFloat>::communicate(
      bool &shouldAbort
      )
{

   if( cmapLapParaSolver->iReceiveIsNecessary() )
   {
      cmapLapParaSolver->iReceiveMessages();
   }

   // check interrupting
   if( cmapLapParaSolver->isInterrupted() )
   {
      shouldAbort = true;
      return true;
   }

   // check to share vectors
   if( cmapLapParaSolver->shareVectorsIsNecessary() )
   {
      // send top-nSendVector vector
      int nRows = std::min(nSendVectors, static_cast<int>(SendStack.size()));
      if( nRows > 0 )
      {
         ParaCMapLAP::LatticeBasis<int> sendVectorsBasis(nRows, L->n);
         for( int i = 0; i < nRows; ++i )
         {
            sendVectorsBasis.row(i) = (SendStack.extract()->getVector()).template cast<int>();
         }
         if( verbose > 1 )
         {
            std::cout
               << "[Rank=" << rank << ", Thread=" << threadId << "] "
               << "send vectors";
            for( int i = 0; i < nRows; ++i )
            {
               std::cout << " " << std::sqrt(sendVectorsBasis.row(i).squaredNorm());
            }
            std::cout << std::endl;
         }
         cmapLapParaSolver->sendVectors(sendVectorsBasis);
      }

      cmapLapParaSolver->requestVectors(nReceiveVectors);
      if( cmapLapParaSolver->hasReceivedVectors() )
      {
         LapTools::LatticeBasis<BasisFloat> receivedVectors(nReceiveVectors, L->n);
         cmapLapParaSolver->getReceivedVectors<BasisFloat>(receivedVectors);
         int nReceivedVectors = receivedVectors.rows();
         for( int i = 0; i < nReceivedVectors; ++i )
         {
            LapTools::LatticeVector<BasisFloat> v = receivedVectors.row(i);
            auto vv = VectorElementPtr(new VectorElementType(v));
            this->Stack.insert(vv);
            if( vv->squaredNorm() < this->bestObjectiveValue )
            {
               updateBestObjectiveValue('U', vv);
            }
            if( verbose > 1 )
            {
               std::cout
                  << "[Rank=" << rank << ", Thread=" << threadId << "] "
                  << "receive (" << i+1 << "/" << nReceivedVectors << ")"
                  << ", norm: " << vv->norm()
                  << ", v: " << vv->vector.transpose()
                  << std::endl;
            }
         }
      }
   }

   // check to be necessary to send the status
   if( cmapLapParaSolver->notificationIsNecessary() )
   {
      sendStatus();
   }

   return true;
}


///
/// @brief update bestObjectiveValue
/// @param[in] sigh line-header character
/// @param[in] shortest vector element
///
template<typename BasisFloat, typename GSFloat>
bool
CmapGaussSieve<BasisFloat, GSFloat>::updateBestObjectiveValue(
      char sigh,
      VectorElementPtr v
      )
{
   if( LapTools::GaussSieve<BasisFloat, GSFloat>::updateBestObjectiveValue(sigh, v) )
   {
      ParaCMapLAP::LatticeVector<int> vv = v->vector.template cast<int>();
      cmapLapParaSolver->sendSolution(vv, vv.squaredNorm());
      return true;
   }
   return false;
}


///
/// @breif send SolverState of Sieve
///
template<typename BasisFloat, typename GSFloat>
void
CmapGaussSieve<BasisFloat, GSFloat>::sendStatus(
      )
{
   double shortestNorm = std::sqrt(this->bestObjectiveValue);
   cmapLapParaSolver->sendSolverState(
         this->runningTime,
         L->m, // blocksize
         this->nLoop,
         this->List.size(),
         this->Stack.size(),
         this->List.size(),
         this->nCollision,
         shortestNorm
         );
}


///
/// instantiation
///
template class CmapGaussSieve<int, double>;
template class CmapGaussSieve<int, long double>;
template class CmapGaussSieve<long int, double>;
template class CmapGaussSieve<long int, long double>;
