#ifndef DEEPLLL_H
#define DEEPLLL_H

#include "lattice.h"

/******************
size-reduce at k
******************/
inline void lattice::size_reduce(int k) 
{	
	int i, t, n=NumRows;
	UINT r;
	
	if (k >= 2) {
		for (t = k-1; t >= 1; --t) {
			if (fabsl(mu(k, t)) > ETA) {
				r = (UINT) round(mu(k, t));
				/* basis update */
				for (i=1; i<=n; ++i) { basis(k, i) -= r*basis(t, i); }
				/* GS-mu update */
				for (i=1; i<=t; ++i) { mu(k, i) -= r*mu(t, i); }
			}
		}
	}
}

/*****************************************************
DeepLLL (with deep insertion restriction) 
- Input
	int start		: start index
	int end			: end index
	int h			: deep start index
	FLOAT alpha		: reduction parameter
	int gamma		: insertion restriction parameter
- Output 
DeepLLL-reduced basis
******************************************************/
inline bool lattice::DeepLLL(int start, int end, int h, FLOAT alpha, int gamma) 
{
	int k = h, i, j, l, s, t, flag, n = NumRows;
	int ii, t_1, s_1;
	bool res; 
	
	while (k <= end) {
		/* Size-reduce at k */
		size_reduce(k);
		
		res = InnerProduct(tmp, basis(k), basis(k));
		if (res == false) {
			cout << "InnerProduct error in DeepLLL" << endl; 
			return false; 
		}
		i = 1; flag = 0;
		while (i < k) {
			/* Check of deep exchange condition */
			if ( i>=start && tmp < alpha*B(i) ) {
				/* Check of deep insertion restriction */
				if ( (i-start+1) <= gamma || (k-i+1) <= gamma) { 
					flag = 1;
				}
			}
			if (flag == 0) {
				tmp -= mu(k, i)*mu(k, i)*B(i);
				++i;
			} else {	
				/* Deep insertion at (i, k) */
				tmpp = basis(k);
				for (t=k; t>i; --t) { basis(t) = basis(t-1); }
				basis(i) = tmpp; 
				
				/* (Efficient) GSO update after deep insertion */
				if ( flag == 1 ) {
					D(k) = B(k);
					BB(k) = B(k);
					DD(k) = 1.0 / D(k);
					ii = i-1;

					for (l = k-1; l != ii; --l) {
						BB(l) = mu(k, l) * B(l);
						D(l) = D(l + 1) + mu(k, l) * BB(l);
						DD(l) = 1.0 / D(l);
						B(l + 1) = D(l + 1) * B(l) * DD(l);
					}
					B(i) = D(i);
					
					t_1 = k - 1;
					tmp1 = mu(k, t_1) * DD(k);
					for (s = k + 1; s <= n; ++s) {
						tmpld(s) = BB(k) * mu(s, k);
						mu(s, k) = mu(s, t_1) - tmp1 * tmpld(s);
					}

					for (t = k - 1; t > i; --t) {
						t_1 = t - 1;
						tmp1 = mu(k, t_1) * DD(t);
						for (s = k + 1; s <= n; ++s) {
							tmpld(s) += BB(t) * mu(s, t);
							mu(s, t) = mu(s, t_1) - tmp1 * tmpld(s);
						}
						for (s = k; s > t + 1; --s) {
							s_1 = s - 1;
							tmpld(s) += BB(t) * mu(s_1, t);
							mu(s, t) = mu(s_1, t_1) - tmp1 * tmpld(s);
						}
						s_1 = t + 1;
						tmpld(s_1) = BB(t);
						mu(s_1, t) = mu(t, t - 1) - tmp1 * tmpld(s_1);
					}

					for (s = k + 1; s <= n; ++s) {mu(s, i) = (tmpld(s) + BB(i) * mu(s, i)) * DD(i);}
					ii = i + 1;
					for (s = k; s > ii; --s) {mu(s, i) = (tmpld(s) + BB(i) * mu(s - 1, i)) * DD(i);}

					mu(ii, i) = BB(i) * DD(i);
					for(t = 1; t < i; ++t) {mut(t) = mu(k, t);}
					for(s = k; s != i; --s) {
						for(t = 1; t < i; ++t) {
							mu(s, t) = mu(s - 1, t);
						}
					}
					for(t = 1; t < i; ++t) {mu(i, t) = mut(t);}
				} 	
				/* Update index k */
				k = max(i, start+1);
				size_reduce(k);
				break;
			}
		}			
		++k;
	}

#if 0
	/* for debug */
	for (i=2; i<end; ++i) {
		for (j = 1; j<i; ++j) {
			if (abs(mu(i, j)) > ETA) {
				cout << "Size-reduce error in DeepLLL" << endl;
				cout << "i=" << i << ", j=" << j << ": " << "mu(i, j) = " << mu(i, j) << endl;
				return false;
			}
		}
	}
#endif

	return true; 
}


#endif
