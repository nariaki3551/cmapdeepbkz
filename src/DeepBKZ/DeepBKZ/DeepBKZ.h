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
inline bool lattice::DeepBKZ(int start, int end, int b, FLOAT alpha, int gamma, int abort)
{
   bool res;
   int z, j, k, h, i, m, n = NumRows, l, tour = 0, N, ii, jj;
   double current;
   long double rho, rho1;
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

   cout << "\nDeepBKZ-" << b << endl;
   cout << "parallel =" << parallel << endl;
   SetGSO();
   res = DeepLLL(start, end, start+1, alpha, gamma); /* Initial DeepLLL */
   if (res != true) {
      cout << "DeepLLL error1 in DeepBKZ" << endl;
      return false;
   }
   current = pow(B(start), 0.5);
   cout << "Initial norm = " << current << endl;

   /* Main loop */
   GSA_slope(rho, start, end);
   tour = 0; N = 0;
   z = start-1; j = start-1;
   while( z < end-1 ){

      communicateInTour(shouldAbort);
      if (shouldAbort) return true;

      if (j == end-1) {
         j = start-1; ++tour;
         GSA_slope(rho1, start, end);
         if (tour % 1 == 0) {
            cout << "tour = " << tour << ": GSA slope = " << rho1 << endl;
         }
         /* Early termination */
         if (abort != 0) {
            if (fabs(rho) < fabs(rho1)) {
               ++N;
               if (N >= abort) {
                  cout << "Auto-abort termination" << endl;
                  time_end = clock();
                  total = static_cast<double>(time_end - time_start)/CLOCKS_PER_SEC;
                  cout << "Total time (seconds) = " << total << endl;
                  cout << "DeepLLL time (seconds) = " << time << endl;
                  return true;
               }
            } else { N = 0; }
            rho = rho1;
         }

         communicate(shouldAbort);
         if (shouldAbort) return true;
      }
      ++j; k = min(j+b-1, end); h = min(k+1, end);

      /* Conversion to NTL */
      for (ii=1; ii<=n; ++ii) {
         for (jj=1; jj<=n; ++jj) {
            LL(ii, jj) = basis(ii, jj);
         }
      }
      BB = LL;
      BB.updateGSBasis();

      /* Enumeration */
      radius = pow(alpha*B(j), 0.50); prob = 1.0;
      if (k-j+1 > 45) { prob = 4.0/pow(1.05, k-j+1); }
      VV = ENUM(BB, radius, prob, enum_mode_find_shortest, 0, VL0, "istart="+to_stdstring(j) + " iend=" + to_stdstring(k) + "parallel=" + to_stdstring(parallel));

      /* Enumerate projected shortest vector */
      // for (i=1; i<=n; ++i) { R(i) = (double) alpha*B(j); }
      // res = Enum(v, R, j, k);

      time_start1 = clock();
      // if( res != true ) {
      if (VV.NumRows() == 0) { ++z;}
      else {
         z = start-1;
         for (ii=1; ii<=n; ++ii) {
            sv(ii) = to_int(VV(1, ii));
         }
         /* Insertion of a short lattice vector */
         for (i=1; i<=j-1; ++i) {
            for (l=1; l<=n; ++l) {Basis(i, l) = basis(i, l);}
         }
         for (l=1; l<=n; ++l) { Basis(j, l) = sv(l);}
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
         if (k != h) {
            if (b < 40) {
               res = DeepLLL(start, end, j, alpha, gamma);
            } else {
               res = SubDeepBKZ(start, end, roundf(b/2.0), 0.99, gamma, 5*b);
            }
            if (res != true) {
               cout << "SubDeepBKZ error in DeepBKZ" << endl;
               return false;
            }
         }
#if 0
         if (h < end) {
         // if (end-j > 70) {
            // res = DualDeepBKZ(25, 0.99, j+1, end, n, 5);
            res = DualDeepLLL(0.99, j+1, end, end-1, 20);
            if (res != true) {
               cout << "DualDeepBKZ error in DeepBKZ" << endl;
               return false;
            }
         }
#endif
      }
      time_end1 = clock();
      time += static_cast<double>(time_end1 - time_start1)/CLOCKS_PER_SEC;

      /* for debug */
      if (B(start) < alpha*current*current) {
         time_end = clock();
         total = static_cast<double>(time_end - time_start) / CLOCKS_PER_SEC;
         current = pow(B(start), 0.5);
         cout << "Time " << total << " Norm = " << current << endl;
         cout << basis(start) << endl;
      }
   }

   time_end = clock();
   total = static_cast<double>(time_end - time_start)/CLOCKS_PER_SEC;
   cout << "Total time (seconds) = " << total << endl;
   cout << "DeepLLL time (seconds) = " << time << endl;
   return true;
}


/**********
DeepBKZ
**********/
inline bool lattice::SubDeepBKZ(int start, int end, int b, FLOAT alpha, int gamma, int num)
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

   // cout << "\nDeepBKZ-" << b << endl;
   // SetGSO();
   res = DeepLLL(start, end, start+1, alpha, gamma); /* Initial DeepLLL */
   if (res != true) {
      cout << "DeepLLL error1 in DeepBKZ" << endl;
      return false;
   }
   current = pow(B(start), 0.5);
   // cout << "Initial norm = " << current << endl;

   /* Main loop */
   // GSA_slope(rho, start, end);
   tour = 0; N = 0;
   z = start-1; j = start-1;
   while( z < end-1 ){

      communicateInTour(shouldAbort);
      if (shouldAbort) return true;

      if (j == end-1) {
         j = start-1; ++tour;
         if (tour >= num) {
            return true;
         }
      }
      ++j; k = min(j+b-1, end); h = min(k+1, end);

      /* Enumerate projected shortest vector */
      for (i=1; i<=n; ++i) { R(i) = (double) alpha*B(j); }
      res = Enum(v, R, j, k);

      if( res != true ) { ++z;}
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
         NTL::LLL_QP(Basis, alpha, 0, 0, 0);

         /* Conversion from NTL::mat_ZZ to UINT */
         for (i=1; i<=h; ++i) {
            for (l=1; l<=n; ++l) {basis(i, l) = to_int(Basis(i+1, l));}
         }

         /* DeepLLL after insertion of short lattice vector */
         SetGSO();
         if (k != h) {
            res = DeepLLL(start, end, j, alpha, gamma);
            if (res != true) {
               cout << "DeepLLL error in DeepBKZ" << endl;
               return false;
            }
         }
      }
      /* for debug */
      if (B(start) < alpha*current*current) {
         current = pow(B(start), 0.5);
         cout << "Norm = " << current << endl;
         cout << basis(start) << endl;
      }
   }
   return true;
}

}  // namespace DeepBKZTool

#endif // DeepBKZ_H
