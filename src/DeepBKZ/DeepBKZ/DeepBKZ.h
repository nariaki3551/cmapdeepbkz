#ifndef DeepBKZ_H
#define DeepBKZ_H

#include "lattice.h"
#include "DeepLLL.h"
#include "Enum.h"
#include <time.h>

namespace DeepBKZTool
{

/**********
DeepBKZ
**********/
inline bool lattice::DeepBKZ(int start, int end, int b, FLOAT alpha, int gamma, int abort, int solver_id)
{
   bool res;
   int z, j, k, h, i, m, n = NumRows, l, tour = 0, N, ii, jj, bb, beta;
   double current;
   long double rho, rho1, slope;
   Vec<int> v; v.SetLength(n);
   Vec<double> R; R.SetLength(n);
   Vec<UINT> sv; sv.SetLength(n);
   mat_ZZ Basis; Basis.SetDims(n+1, n);
   mat_ZZ LL, VV; LL.SetDims(n, n);
   clock_t time_start = clock(), time_end;
   clock_t time_start1, time_end1;
   double total, time = 0.0;
   LatticeBasis<double> BB;
   double radius, prob;
   // int parallel = omp_get_num_procs();
   int parallel = 1;
   bool shouldAbort = false;
   long double tmp, alp, r;

   std::cout << "S rank " << std::setw(5) << right << solver_id
      << ": DeepBKZ-" << b
      << std::endl;
   SetGSO();
   res = DeepLLL(start, end, start+1, alpha, gamma); /* Initial DeepLLL */
   if (res != true) {
      std::cout << "E rank" << std::setw(5) << right << solver_id
         << ": DeepLLL error1 in DeepBKZ" << endl;
      return false;
   }
   communicate(shouldAbort);
   current = pow(B(start), 0.5);
   // std::cout << "  rank " << std::setw(5) << right << solver_id
   //    << ": Initial norm = " << current << std::endl;

	if (n <= 110) { bb = b+15; }
	else { bb = b + 10; }
	// bb = b + 15;
	tmp = bb/(2*M_PI*exp(1))*pow(M_PI*bb, 1.0/(bb));
	tmp = pow(tmp, 1.0/(2*(bb-1)));
	slope = -4*log(tmp);

   /* Main loop */
   GSA_slope(rho, start, end);
   if (rho > slope) {
      std::cout << "A rank " << std::setw(5) << right << solver_id
         << ": A large blocksize is required rho " << rho << " > slope " << slope
         << std::endl;
		return true;
   }

   tour = 0; N = 0;
   z = start-1; j = start-1;
   while( z < end-1 ){

      communicateInTour(shouldAbort);
      if (shouldAbort) {
         // std::cout << "A rank " << std::setw(5) << right << solver_id
         //    << ": terminate due to the request by Load Coordinator"
         //    << std::endl;
         return true;
      }

      if (j == end-1) {
         j = start-1;
         ++tour;
         GSA_slope(rho1, start, end);
         if (tour % 1 == 0) {
            std::cout << "  rank " << std::setw(5) << right << solver_id
               << ": block " << b << ", tour " << tour << ", norm " << current << ", GSA slope " << rho1 << std::endl;
         }

         /* Early termination */
         if (abort != 0) {
            if (0.999*fabs(rho) < fabs(rho1)) {
               ++N;
               if (N >= abort) {
                  time_end = clock();
                  total = static_cast<double>(time_end - time_start)/CLOCKS_PER_SEC;
                  std::cout << "  rank " << std::setw(5) << right << solver_id
                     << ": Auto-abort termination"
                     << ", Total time (seconds) = " << total
                     << ", DeepLLL time (seconds) = " << time << std::endl;
                  return true;
               }
            } else { N = 0; }
            rho = rho1;
         }

         communicate(shouldAbort);
         current = pow(B(start), 0.5);
         if (shouldAbort) {
            // std::cout << "A rank " << std::setw(5) << right << solver_id
            //    << ": terminate due to the request by Load Coordinator"
            //    << std::endl;
            return true;
         }

		   if (rho1 > slope) {
            std::cout << "A rank " << std::setw(5) << right << solver_id
               << ": A large blocksize is required rho " << rho1 << " > slope " << slope
               << std::endl;
		      return true;
		   }
      }
      ++j; k = min(j+b-1, end); h = min(k+1, end);
	  // ++j; k = min(j+min(b+j, 70)-1, end); h = min(k+1, end);

      /* Conversion to NTL */
      for (ii=1; ii<=n; ++ii) {
         for (jj=1; jj<=n; ++jj) {
            LL(ii, jj) = basis(ii, jj);
         }
      }
      BB = LL; BB.updateGSBasis();

      /* Enumeration */
      radius = pow(alpha*B(j), 0.50); prob = 1.0;
	   beta = k-j+1;
      if (beta > 45) {
	  	   if (beta >= 84) {
            //Heuristic parameter from experiments
            r = -40.60241266/beta/(beta-1) - 1.603898*log(beta)/(beta+113.3386);
            r = exp(r);
            double abar = 1.0 / (double)bkzconstants::ghconstant(beta) * exp(-0.25*(beta-1.0)*log(r));
            alp = abar * ((1.0+beta)/beta);
            prob = pow(alp, -beta)*2;
         } else {
            r = -1.50065628562105e-06 * beta*beta +  0.000444247990862759 * beta +  0.932324451217323;  //fitting
            double abar = 1.0 / (double)bkzconstants::ghconstant(beta) * exp(-0.25*(beta-1.0)*log(r));
            double minalpha = pow(1.38629436111989, 1.0/beta);
            minalpha = pow(2.77258872223978, 1.0/beta);
            abar = max(abar, minalpha);
            alp = abar * ((1.0+beta)/beta);
            prob = pow(alp, -beta)*2;
            // cout << "prob = " << prob << endl;
            if (beta<=55) prob = min(0.15+0.03*(55-beta),1.0);
         }
         if (beta <= 65) { prob *= 1.5; }

         // prob = 5.0/pow(1.05, k-j+1);
         tmp = 1.0;
         for (i=j; i<=k; ++i) { tmp *= B(i); }
         tmp = pow(tmp, 1.0/(k-j+1))*(k-j+1)/(2.0*M_PI*exp(1));
         tmp = sqrtl(tmp);
         // cout << "approx = " << radius/tmp << endl;
         if (radius > 1.05*tmp) { radius = 1.05*tmp; }
         // if (radius <= 0.90*tmp) { radius = 0; }
      }
      VV = ENUM(BB, radius, prob, enum_mode_find_shortest, 0, VL0, "istart="+to_stdstring(j) + " iend=" + to_stdstring(k) + "parallel=" + to_stdstring(parallel));

      /* Enumerate projected shortest vector */
      // for (i=1; i<=n; ++i) { R(i) = (double) alpha*B(j); }
      // res = Enum(v, R, j, k);

      time_start1 = clock();
      // if( res != true ) {
      if (VV.NumRows() == 0) {
         ++z;
         // res = DeepLLL(start, h, j, alpha, gamma);
      } else {
         z = start-1;
         // for (ii=1; ii<=n; ++ii) { sv(ii) = to_int(VV(1, ii)); }

         /* Insertion of a short lattice vector */
         for (i=1; i<=j-1; ++i) {
            for (l=1; l<=n; ++l) {Basis(i, l) = basis(i, l);}
         }
         for (l=1; l<=n; ++l) { Basis(j, l) = to_int(VV(1, l));}
         for (i=j+1; i<=n+1; ++i) {
            for (l=1; l<=n; ++l) {Basis(i, l) = basis(i-1, l);}
         }

         /* Modified LLL to remove the linear dependency  */
         NTL::LLL_QP(Basis, alpha, 0, 0, 0);

         /* Conversion from NTL::mat_ZZ to UINT */
         for (i=1; i<=n; ++i) {
            for (l=1; l<=n; ++l) {basis(i, l) = to_int(Basis(i+1, l));}
         }

         /* DeepLLL after insertion of short lattice vector */
         SetGSO();
         res = DeepLLL(start, end, j, alpha, gamma);
         res = SubDeepBKZ(start, end, min((int)round((b)/2.0), 40), 0.99, gamma, 300, solver_id);
#if 0
         // if (k != h) {
       if (end-j > 44) {
            if (b < 45) {
               // res = DeepLLL(start, end, j, alpha, gamma);
            } else {
               res = SubDeepBKZ(start, h, min((int)round(b-5)/2.0), 35), 0.99, gamma, 200, solver_id);
            }
            if (res != true) {
               cout << "SubDeepBKZ error in DeepBKZ" << endl;
               return false;
            }
         }
#endif

#if 1
         /* for debug */
         if (B(start) < 0.9999*current*current) {
            time_end = clock();
            total = static_cast<double>(time_end - time_start) / CLOCKS_PER_SEC;
            current = pow(B(start), 0.5);
            std::cout << "* rank " << std::setw(5) << right << solver_id
               << ": Time " << total << ", Norm = " << current
               << ", AF = " << approxFactor(current) << ", RHF = " << pow(hermiteFactor(current), 1.0/n)
               // << ", Vec = " << basis(start)
               << std::endl;
            sendSolution();
         }
#endif
      }
      time_end1 = clock();
      time += static_cast<double>(time_end1 - time_start1)/CLOCKS_PER_SEC;

      /* for debug */
      if (B(start) < 0.9999*current*current) {
         time_end = clock();
         total = static_cast<double>(time_end - time_start) / CLOCKS_PER_SEC;
         current = pow(B(start), 0.5);
         std::cout << "* rank " << std::setw(5) << right << solver_id
            << ": Time " << total << ", Norm = " << current
            << ", AF = " << approxFactor(current) << ", RHF = " << pow(hermiteFactor(current), 1.0/n)
            // << ", Vec = " << basis(start)
            << std::endl;
         sendSolution();
      }
   }

   time_end = clock();
   total = static_cast<double>(time_end - time_start)/CLOCKS_PER_SEC;
   std::cout << "F rank " << std::setw(5) << right << solver_id
      << ": Total time (seconds) = " << total
      << ": DeepLLL time (seconds) = " << time << std::endl;
   return true;
}


/**********
DeepBKZ
**********/
inline bool lattice::SubDeepBKZ(int start, int end, int b, FLOAT alpha, int gamma, int num, int solver_id)
{
   bool res;
   int z, j, k, h, i, m, n = NumRows, l, tour = 0, N, ii, jj;
   double current;
   long double rho, rho1;
   Vec<int> v; v.SetLength(n);
   Vec<double> R; R.SetLength(n);
   Vec<UINT> sv; sv.SetLength(n);
   mat_ZZ Basis;
   bool shouldAbort = false;

#if 0
   // cout << "\nDeepBKZ-" << b << endl;
   // SetGSO();
   res = DeepLLL(start, end, start+1, alpha, gamma); /* Initial DeepLLL */
   if (res != true) {
      std::cout << "rank " << std::setw(5) << right << solver_id
         << ": DeepLLL error1 in DeepBKZ" << std::endl;
      return false;
   }
#endif
   current = pow(B(start), 0.5);
   // cout << "Initial norm = " << current << endl;

   /* Main loop */
   GSA_slope(rho, start, end);
   tour = 0; N = 0;
   z = start-1; j = start-1;
   while( z < end-1 ){
      communicateInTour(shouldAbort);
      if (shouldAbort) return true;

      if (j == end-1) {
         j = start-1; ++tour;
         // if (tour % 50 == 0) {
         //    std::cout << "rank " << std::setw(5) << right <<
         //       solver_id << ": block " << b << ", sub-tour " << tour << ", norm " << current << ", GSA slope " << rho1 << std::endl;
         // }
         if (tour >= num) {
           return true;
         }

         GSA_slope(rho1, start, end);
         /* Early termination */
         if (fabs(rho) < fabs(rho1)) {
            ++N;
            if ( N >= 5 ) {
              // cout << "tours = " << tour << endl;
              // basis = tmp_basis; SetGSO();
              return true;
            }
         } else { N = 0; }
         rho = rho1;
      }
      ++j; k = min(j+b-1, end); h = min(k+1, end);

      /* Enumerate projected shortest vector */
      for (i=1; i<=n; ++i) { R(i) = (double) alpha*B(j); }
      res = Enum(v, R, j, k);

      if( res != true ) { ++z; }
      else {
         z = start-1;

         /* Construction of a short lattice vector */
         for (l=1; l<=n; ++l) { sv(l) = 0; }
         for (i=j; i<=k; ++i) {
            for (l=1; l<=n; ++l) { sv(l) += v(i)*basis(i, l);}
         }

         /* Insertion of a short lattice vector */
         Basis.SetDims(h+1, n);
         for (i=1; i<=j-1; ++i) {
            for (l=1; l<=n; ++l) {Basis(i, l) = basis(i, l);}
         }
         for (l=1; l<=n; ++l) { Basis(j, l) = sv(l);}
         for (i=j+1; i<=h+1; ++i) {
            for (l=1; l<=n; ++l) {Basis(i, l) = basis(i-1, l);}
         }

         /* Modified LLL to remove the linear dependency  */
         NTL::LLL_FP(Basis, alpha, 0, 0, 0);

         /* Conversion from NTL::mat_ZZ to UINT */
         for (i=1; i<=h; ++i) {
            for (l=1; l<=n; ++l) {basis(i, l) = to_int(Basis(i+1, l));}
         }

         /* DeepLLL after insertion of short lattice vector */
         SetGSO();
         // if (k != h) {
            res = DeepLLL(start, end, j, alpha, gamma);
            if (res != true) {
               std::cout << "rank " << std::setw(5) << right << solver_id
                  << ": DeepLLL error in DeepBKZ" << std::endl;
               return false;
            }
         // }
      }
#if 0
      /* for debug */
      if (B(start) < alpha*current*current) {
         current = pow(B(start), 0.5);
         std::cout << "rank " << std::setw(5) << right << solver_id
            << ": Norm = " << current << ", " << basis(start) << std::endl;
         sendSolution();
      }
#endif
   }
   return true;
}

}  // namespace DeepBKZTool

#endif // DeepBKZ_H
