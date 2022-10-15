#include "bkzlibwrapper.hpp"

//Size reduction

#ifdef define_normal
template <typename T> inline  int SizeReduce(LatticeBasis<T>& B,mat_ZZ*U,int j,int bound=1,int opt=gsupdate,int redtype=normalgsreduction) {
#endif
    
#ifdef define_small
template <typename GFLOAT,typename T> inline int SizeReduce(SmallLatticeBasis<GFLOAT,T>& B,pbkzmatrix<int>* U,int j,int bound=1,int opt=gsupdate,int redtype=normalgsreduction) {
#endif
//start of the subroutine

#ifdef define_normal
            B.updateGSBasis(1,j);
#endif
#ifdef define_small
            B.updateGSBasis_noapprox(1,j);
#endif

            
        int m = B.L.NumCols();
        T prev_maxbias = -1;
        T maxbias;
        int bigcount = 0;
        bool updateflag=false;
        do {
            maxbias = 0;
            char flag_bigmu = false;
            int gloop = 0;
            do {
                T bias,absbias;
                maxbias = 0;
                int maxflag = 1;    
                for (int i=j-1;i>=bound;i--) {
                    bias = round(B.gs.mu[j][i]);
                    if (abs(B.gs.mu[j][i]) < 0.501) bias = 0;    //for precision-loss
                    absbias = abs(bias);   //need to blush up heuristic constants
        #ifdef define_normal
                    if (absbias > 1e+6) {
                        flag_bigmu = true;
                        bigcount++;
                        gloop += maxflag;
                        maxflag = 0;
                        if (gloop>100) {
                            B.GScomputed = j-1;
                            B.updateGSBasis(j,j);
                            B.gs.displaymu(j,j);
                            gloop = 0;
                            return -1;  //infinite loop?
                        }
                    } else {
                        flag_bigmu = false;
                        if (redtype == lazyreduction) {
                            if (i < j-5) {
                                bias = 0;
                                absbias = 0;
                            }
                        }
                    }
            #endif
                    maxbias = max(maxbias,absbias);
                    if (bias!=0) {
                        T* mult_i = B.gs.mu[i].data()+1;
                        T* mult_j = B.gs.mu[j].data()+1;
                        if (1) {
                            if (bias==1) {
                                for (int k=1;k<i;k++) (*mult_j++) -= (*mult_i++);
                            } else 
                            if (bias==-1) {
                                for (int k=1;k<i;k++) (*mult_j++) += (*mult_i++);
                            } else {
                                for (int k=1;k<i;k++) (*mult_j++) -= (*mult_i++) * bias;
                            }                   
                        }
        #ifdef define_normal
                        if (bias!=0) {
                            ZZ zzbias;
                            int ibias; 
                            //Here, bias is in [-10^7,10^7]
                            if (flag_bigmu==true) {
                                zzbias = conv<ZZ>(to_intstring(B.gs.mu[j][i]).c_str());
                            } else {
                                double dbias = conv_to_double(bias);
                                //double dbias = bias;
                                ibias = round(dbias);
                                zzbias = ibias;
                            }
                            RowTransformwrap(B.L[j-1],B.L[i-1],zzbias);    //L[j-1] <- L[j-1] - bias*L[i-1]
                            if (U!=0) RowTransformwrap((*U)[j-1],(*U)[i-1],zzbias); 
                            updateflag = true;
                        }
#endif
        #ifdef define_small
                        if (bias!=0) {
                            addvector<T>(B.L[j-1],-bias,B.L[i-1]);
                            if (U!=0) {
                                int ibias;
                                ibias = static_cast<int>(bias);
                                //conv(ibias,bias);
                                addvector((*U)[j-1], (-ibias),(*U)[i-1]);
                            }
                            updateflag = true;
                        }
        #endif
#ifdef define_normal
                        if (flag_bigmu==true) {
                             if (bias!=0) *mult_j -= bias;
                        } else {                        
                             if (bias!=0)  *mult_j -= bias;
                        }
#endif
#ifdef define_small
                       if (bias!=0)   *mult_j -= bias;
#endif
                    }
                }   //end of for
#ifdef define_normal
                if (flag_bigmu==true) {
                    B.GScomputed = j-1;
                    B.updateGSBasis(j,j);
                    updateflag = false;
                }
#endif
                if (prev_maxbias > 0) {
                    if (prev_maxbias <= maxbias) {
                        //some precision error?
                        maxbias = -1;   //error flag;
                        break;
                    }
                }
                prev_maxbias = maxbias;
            } while (flag_bigmu==true);
        }   while (maxbias > 100000);   //the constant need to polish

        if ((opt==gsupdate) && (updateflag == true)) {
            B.GScomputed = j-1;
#ifdef define_normal
            B.updateGSBasis(j,j);
#endif
#ifdef define_small
            B.updateGSBasis_noapprox(j,j);
#endif
        } else {
            B.updateapprox(j,j);
            B.GScomputed = j;
        }
        
        return 0;
}

#define stnormal 0x00
#define stgreedy 0x01


#ifdef define_normal
template <typename T> int local_LLL(LatticeBasis<T>& B, mat_ZZ* U, double delta,int istart,int iend,int strategy,int vl,int& swapcount,T zerolimit=1) {
#endif

#ifdef define_small
template <typename T2,typename T> int local_LLL_small(SmallLatticeBasis<T2,T>& B,pbkzmatrix<int>* U,T delta,int istart,int iend,int strategy,int vl,int& swapcount,T2 zerolimit=1e-4) {
//Note: for a dependent basis, a vector shorter than zerolimit is regarded as the zero vector
#endif

    
#ifdef define_normal
    B.updateGSBasis(1,istart);
#endif

#ifdef define_small
    B.updateGSBasis_noapprox(1,istart);
#endif
    
    int j = max(1,istart);
    swapcount=0;
    double ss = gettimeofday_sec();
    
    int prev_swap_pos=0;
    int same_swap_count=0;
    int achieved = istart;
    while (1) {
        do {

#ifdef define_normal
            B.updateGSBasis(j,j);
#endif
#ifdef define_small
            B.updateGSBasis_noapprox(j,j);
#endif
            //if (SizeReduce(B,U,j,gsupdate,lazyreduction)<0) {
            if (SizeReduce(B,U,j,1,gsupdate)<0) {
                return -1;
            }
            
            if (vl>=3) {
                ccout << "LLL: j=" << j << "; ";
                B.gs.displayc(1,j);
            }

#ifdef define_normal
            T len = LengthOf<T>(B.L[j-1]) ;
            if (len<zerolimit) {
#endif
#ifdef define_small
            if (LengthOf(B.L[j-1])<zerolimit) {
#endif
                
                for (int i=j;i<=iend-1;i++) {
                    swap(B.L[i-1],B.L[i]);
                    if (U!=0) swap((*U)[i-1],(*U)[i]);
                }
                B.GScomputed = j-1;
 #ifdef define_normal
                B.updateGSBasis(j,j);
#endif
#ifdef define_small
                B.updateGSBasis_noapprox(j,j);
#endif
                iend--;
                if (j > iend) break;
            } else {
                break;
            }
        } while (1);        
        if (j > iend) break;
        
        //check Lovasz condition


        bool approxlovaszcond = true;
        if (j>=istart + 1) {
            if (B.gs.cd[j-1] > B.gs.cd[j]*2 ) {
                approxlovaszcond = false;
            } else 
            if (B.gs.cd[j-1] * delta > B.gs.cd[j] + B.gs.mu[j][j-1] * B.gs.mu[j][j-1] * B.gs.cd[j-1]) {
                approxlovaszcond = false;
            }
        }
                
        if ( approxlovaszcond == false ) {
            swap(B.L[j-1],B.L[j-2]);    //indexes are shifted
            if (U!=0) swap((*U)[j-1],(*U)[j-2]);
            B.GScomputed = j-2;
            if (j==prev_swap_pos) {
                if (same_swap_count++ > 10000) {
                    //infinite loop?
                    ccout << "infinite_loop: " << j << endl;
                    B.GScomputed = 0;
#ifdef define_normal
            B.updateGSBasis(1,j-2);
#endif
#ifdef define_small
            B.updateGSBasis_noapprox(1,j-2);
#endif
                    same_swap_count=0;
                }
            } else {
                    same_swap_count=0;
            }
            prev_swap_pos=j;

            if (strategy==stnormal) j--;
            swapcount++;
        } else {
            if (strategy==stnormal) j++;
        }

        if ((j > achieved) || (swapcount%1000==0) || (vl>=3)) {
            achieved = max(j,achieved);
            if (vl>=2) {
                ccout << "small_LLL: time=" << gettimeofday_sec() - ss << " computed " << achieved << "/" << iend << " #swap=" <<  swapcount <<"       \r";
                ccout.flush();
            }
            for (int i=1;i<j;i++) {
                if (B.gs.c[i] <= 0) {
                    ccout << endl;
                    ccout << "Norm error? at i=" << i << endl;
                    ccout << endl;
                    B.gs.displayc(1,j);
                    exit(0);
                    return -1;
                } 
            }
        }
        if (j > iend) break;
    }
        

    return iend;    //if iend is changed B[iend+1...] are zero vectors
}

