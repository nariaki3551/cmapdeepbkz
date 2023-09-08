#ifndef ENUM_H
#define ENUM_H

#include "lattice.h"

namespace DeepBKZTool
{

/*****************************************************************************
Pruned Enumeration (Latest Date: January 29, 2018 by Yasuda)
Input
	Vec<long long>	(output) & result:
	Lattice			(input)
	Vec<double>		(input) SR: Squared bounds
    int				(input) g
    int	 			(input) h
Output
	Shorter vector of Projected lattice L(pi_g(b_g),...,pi_g(b_h))
	true(1)
	false(0)
Reference:
- "Lattice Enumeration using Extreme Pruning" by Gama-Nguyen-Regev (EUROCRYPT 2011)
- Appendix B, Algorithm 2(Pruned Enumeration)
***********************************************************************************/
inline bool lattice::Enum(Vec<int> &result, Vec<double> SR, int g, int h)
{
	/* Assume that GSO information is already input */
	int i, j, k, last_nonzero, n = NumRows, flag = 0;
	Mat<double> sigma; sigma.SetDims(n+1, n);
	Vec<int> r; r.SetLength(n+1);
	Vec<double> rho; rho.SetLength(n+1);
	Vec<int> v; v.SetLength(n);
	Vec<double> c; c.SetLength(n);
	Vec<int> w; w.SetLength(n);
	double dtmp;

	/* Step 1 */
	for (i = 1; i <= n+1; ++i) {
		for (j = 1; j <= n; ++j) { sigma(i, j) = 0.0; }
	}
	for (i = 0; i <= n; ++i) { r[i] = i; }

	/* Step 2: partial norm rho */
	for (i = 1; i <= n+1; ++i) { rho(i) = 0.0; }

	/* Step 3: current combination of integers for a lattice vector */
	for (i = 1; i <= n; ++i) { v(i) = 0; }
	v(g) = 1;

	/* Step 4: centers */
	for (i = 1; i <= n; ++i) { c(i) = 0.0; }

	/* Step 5: jumps */
	for (i = 1; i <= n; ++i) { w(i) = 0; }

	/* Step 6: largest index i for which v_i != 0 */
	last_nonzero = g;

#if 1
	int n_sub = h-g+1, p, t;
	long double templ, temp;
	Vec<double> eta; eta.SetLength(n_sub);
	if (n_sub >= 35) {
		p = 15;
		templ = sqrtl(B(g));
		for (t = 2; t<= n_sub; ++t) {
			temp = tgammal(0.50*(t+1.0));
			eta[t-1] =  powl(temp*pow(2.0, -1.0*p)*templ, 2.0/(t - 1.0))/M_PI;
			templ *= sqrtl(B(g+t-1));
		}
		for (i = 1; i<n_sub; ++i) {
			SR[i] -= eta[i];
			if (SR[i] <= DBL_EPSILON) {
				return false;
			}
		}
	}
#endif

	k = g;
	while (true) {
		/* Step 9: compute squared norm of current node */
		dtmp = v(k) - c(k);
		dtmp *= dtmp;
		rho(k) = rho(k+1) + dtmp * B(k);
		// if (rho(k) <= SR(n+1-k)) {
		if (rho(k) <= SR[k-g]) {
			if (k == g) {
				/* A solution is found */
				++flag;
				for (i=1; i<=n; ++i) { result(i) = v(i); }
				for (i=1; i<=n; ++i) { SR(i) = min(0.99*rho(g), SR(i)); }
				// SR = 0.99*rho(g); /* SR update */
			} else {
				--k;
				r[k - 1] = max(r[k - 1], r[k]);
				for (i = r[k]; i >= k+1; --i) {
					sigma(i, k) = sigma(i + 1, k) + v(i) * mu(i, k);
				}
				/* Steps 17 and 18 */
				c(k) = -sigma(k + 1, k);
				v(k) = round(c(k));
				w(k) = 1;
			}
		} else {
			++k;
			if (k == h+1) {
				if (flag > 0) {
					// cout << result << endl;
					return true;
				} else {
					/* There is no solution */
					for (i = 1; i <= n; ++i) { result(i) = 0; }
					return false;
				}
			}
			r[k - 1] = k;
			if (k >= last_nonzero) {
				last_nonzero = k;
				++v(k);
			} else {
				if (v(k) > c(k)) { v(k) -= w(k); }
				else { v(k) += w(k); }
				++w(k);
			}
		}
	}
}

}  // namespace DeepBKZTool

#endif //ENUM_H
