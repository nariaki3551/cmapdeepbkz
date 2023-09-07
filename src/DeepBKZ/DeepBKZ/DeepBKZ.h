#ifndef DeepBKZ_H
#define DeepBKZ_H

#include "lattice.h"
#include "DeepLLL.h"
#include "Enum.h"
#include <time.h>
#include <cstdio>
#include "fplll/fplll.h"

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
   double sstart = 1; // same as DimensionOfSharedLattice
   long double av, av_old;

   // abort = 3;
   	// Screening
	av = 0.0;
	// for (jj=1; jj<=n-b+1; jj++) {
	for (jj=1; jj <= n-b+1; jj++) {
	// for (jj = 1; jj <= n-44; ++jj) {
		tmp = 1.0; k = min(jj+b-1, end);
		bb = k-jj+1;
		for (i=jj; i<=k; ++i) { tmp *= B(i); }
		tmp = pow(tmp, 0.50/bb);
		tmp *= pow(tgamma(0.50*bb + 1), 1.0/bb)/sqrtl(M_PI);
		// tmp = pow(tmp, 1.0/b)*b/(2.0*M_PI*exp(1));
		// tmp = sqrtl(tmp);
		av += pow(B(jj), 0.50)/tmp;
	}
	av /= (n-b+1);
	// av /= (n-44);
	// cout << "average = " << av << endl;
	av_old = av;
	// if (av <= 1.0) {return true;}

   std::cout << "rank " << std::setw(5) << right << solver_id
      << ": DeepBKZ-" << b << ", average:" << av << std::endl;
   SetGSO();
   res = DeepLLL(1, end, 2, alpha, gamma); /* Initial DeepLLL */
   if (res != true) {
      cout << "DeepLLL error1 in DeepBKZ" << endl;
      return false;
   }
   current = pow(B(1), 0.5);
   // std::cout << "rank " << std::setw(5) << right << solver_id
   //   << ": Initial norm = " << current << std::endl;

   tour = 0; N = 0;
   z = sstart-1; j = sstart-1;
   while( z < end-1 ){
      // communicateInTour(shouldAbort);
      // if (shouldAbort) return true;

	  communicate(shouldAbort);
      if (shouldAbort) return true;
	  current = pow(B(1), 0.5);
	  DeepLLL(1, end, 2, alpha, gamma);

      if (j == end-1) {

		// Screening
		av = 0.0;
		for (jj=1; jj<=n-b+1; ++jj) {
		// for (jj=1; jj <= n-44; ++jj) {
			tmp = 1.0; k = min(jj+b-1, end);
			bb = k-jj+1;
			for (i=jj; i<=k; ++i) { tmp *= B(i); }
			// tmp = pow(tmp, 1.0/b)*b/(2.0*M_PI*exp(1));
			// tmp = sqrtl(tmp);
			tmp = pow(tmp, 0.50/bb);
			tmp *= pow(tgamma(0.50*bb + 1), 1.0/bb)/sqrtl(M_PI);
			av += pow(B(jj), 0.50)/tmp;
		}
		av /= (n-b+1);
		// av /= (n-44);
		// cout << "average = " << av << endl;
		if (av_old <= av) {
			++N;
			// av_old = av;
			if (N >= 2) {
				time_end = clock();
				total = static_cast<double>(time_end - time_start)/CLOCKS_PER_SEC;
				std::cout << "rank " << std::setw(5) << right << solver_id
				 << ": Auto-abort termination by average"
				 << ", Total time (seconds) = " << total
				 << ", DeepLLL time (seconds) = " << time << std::endl;
				return true;
			}
		}
		else { N = 0; av_old = av; }

         j = sstart-1;
		 ++tour;
         GSA_slope(rho1, start, end);
         if (tour % 1 == 0) {
            std::cout << "rank " << std::setw(5) << right << solver_id
               << ": block " << b << ", tour " << tour << ", norm " << current << ", GSA slope " << rho1 << ", average " << av << std::endl;
         }
      }
      ++j; k = min(j+b-1, end); h = min(k+1, end);

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
      if (beta >= 45) {
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
		if (beta <= 60) { prob *= 1.3; }
		// if (60 < beta <= 70) { prob *= 1.1; }
		// prob *= 1.3;

		// prob = 5.0/pow(1.05, k-j+1);
		tmp = 1.0;
		for (i=j; i<=k; ++i) { tmp *= B(i); }
		// tmp = pow(tmp, 1.0/beta)*(beta)/(2.0*M_PI*exp(1));
		tmp = pow(tmp, 1.0/(2.0*beta));
		tmp *= pow(tgamma(0.50*beta + 1), 1.0/beta)/sqrtl(M_PI);
		// cout << "approx = " << radius/tmp << endl;
		if (radius > 1.08*tmp) { radius = 1.08*tmp; }
		if (radius <= av*tmp) { radius = 0.0; }
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
	  }
      else {
         z = sstart-1;
         // for (ii=1; ii<=n; ++ii) { sv(ii) = to_int(VV(1, ii)); }
		 // cout << VV(1) << endl;

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
		 // res = DeepLLL(start, end, j, alpha, gamma);
		 bb = min(40, ((int)round((b-10)/2.0)));
		 gamma = n;
		 // res = SubDeepBKZ(1, end, bb, 0.99, gamma, 300, solver_id);
		 // res = SubDeepBKZ(1, end, min((int)round(b-5)/2.0), 35)), 0.99, gamma, 300, solver_id);
#if 1
         // if (k != h) {
		 if (end-j >= 44) {
            if (b < 45) {
               res = DeepLLL(1, end, j, alpha, gamma);
            } else {
               res = SubDeepBKZ(1, end, bb, 0.99, gamma, 200, solver_id);
            }
            if (res != true) {
               cout << "SubDeepBKZ error in DeepBKZ" << endl;
               return false;
            }
         }
#endif

#if 1
		   /* for debug */
		   if (B(1) < alpha*current*current) {
			  time_end = clock();
			  total = static_cast<double>(time_end - time_start) / CLOCKS_PER_SEC;
			  current = pow(B(1), 0.5);
			  std::cout << "rank " << std::setw(5) << right << solver_id
				 << ": Time " << total << ", Norm = " << current << ", " << basis(start) << std::endl;
			  sendSolution();
		   }
#endif
      }
      time_end1 = clock();
      time += static_cast<double>(time_end1 - time_start1)/CLOCKS_PER_SEC;

#if 0
      /* for debug */
      if (B(1) < alpha*current*current) {
         time_end = clock();
         total = static_cast<double>(time_end - time_start) / CLOCKS_PER_SEC;
         current = pow(B(1), 0.5);
          std::cout << "rank " << std::setw(5) << right << solver_id
            << ": Time " << total << ", Norm = " << current << ", " << basis(1) << std::endl;
         sendSolution();
      }
#endif
   }

   time_end = clock();
   total = static_cast<double>(time_end - time_start)/CLOCKS_PER_SEC;
   std::cout << "rank " << std::setw(5) << right << solver_id
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
   double current = pow(B(1), 0.5);
   long double rho, rho1;
   Vec<int> v; v.SetLength(n);
   Vec<double> R; R.SetLength(n);
   Vec<UINT> sv; sv.SetLength(n);
   mat_ZZ Basis;
   bool shouldAbort = false;

   /* Main loop */
   GSA_slope(rho, start, end);
   tour = 0; N = 0;
   z = start-1; j = start-1;

   while( z < end-1 ){
	  if (b >= 30) {
		 communicate(shouldAbort);
		 if (shouldAbort) return true;
		 current = pow(B(1), 0.5);
		 DeepLLL(1, end, 2, alpha, gamma);
	  }

      if (j == end-1) {
		 // communicateInTour(shouldAbort);
         // if (shouldAbort) return true;

		 if (b < 30) {
			communicate(shouldAbort);
			if (shouldAbort) return true;
			current = pow(B(1), 0.5);
			DeepLLL(1, end, 2, alpha, gamma);
		 }
         j = start-1; ++tour;
		 if (b > 25) {
			 if (tour >= num) {
				// cout << "tour_limit:" << tour << endl;
				return true;
			}
			GSA_slope(rho1, start, end);
			/* Early termination */
			if (fabs(rho) < fabs(rho1)) {
			   ++N;
			   if ( N >= 5 ) {
				  // cout << "tour_termination:" << tour << endl;
				  return true;
			   }
			} else { N = 0; }
			rho = rho1;
		}
      }
      ++j; k = min(j+b-1, end); h = min(k+1, end);

      /* Enumerate projected shortest vector */
      for (i=1; i<=n; ++i) { R(i) = (double) 0.99*B(j); }
      res = Enum(v, R, j, k);

      if( res != true ) {
		++z;
		// res = DeepLLL(start, h, j, alpha, gamma);
	  } else {
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
         NTL::LLL_QP(Basis, alpha, 0, 0, 0);

         /* Conversion from NTL::mat_ZZ to UINT */
         for (i=1; i<=h; ++i) {
            for (l=1; l<=n; ++l) {basis(i, l) = to_int(Basis(i+1, l));}
         }

         /* DeepLLL after insertion of short lattice vector */
         SetGSO();
         // if (k != h) {
            res = DeepLLL(start, end, j, alpha, gamma);
			// res = DeepLLL(start, h, j, alpha, gamma);
            if (res != true) {
               std::cout << "rank " << std::setw(5) << right << solver_id
                  << ": DeepLLL error in DeepBKZ" << std::endl;
               return false;
          //  }
         }
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
