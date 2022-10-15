#ifndef _inc_sampling_boost_cpp
#define _inc_sampling_boost_cpp

//Sampling from uniform and normal distributions

namespace sampling_tools{
    //basic samplier
   boost::normal_distribution<> n01(0.0, 1.0);
   boost::uniform_01<> u01;
   boost::random::independent_bits_engine<boost::mt19937, 180, boost::multiprecision::cpp_int> rand_gen;
   std::vector<boost::mt19937> mt;

   
    using namespace boost::multiprecision;
    using namespace boost::random;

    inline void resizemt(int mtindex) {
        if (mt.size() <= mtindex+1) {
            #pragma omp critical
            {
                if (mt.size() <= mtindex+1) {
                    int oldsize = mt.size();
                    mt.resize(mtindex+1);
                    for (int i=oldsize;i<mtindex+1;i++) mt[i].seed(i);
                }
            }
        }
    }

    template <typename DFLOAT> DFLOAT random(int mtindex=0) {
        resizemt(mtindex);
        DFLOAT a = (DFLOAT)u01(mt[mtindex]); 
        return a;
    }

    template <typename DFLOAT> DFLOAT unique(DFLOAT a,DFLOAT b,int mtindex=0) {
        resizemt(mtindex);
        return a + (b-a) * u01(mt[mtindex]);
    }

    template <typename DFLOAT> DFLOAT normal(int mtindex=0) {
        resizemt(mtindex);
        DFLOAT a = (DFLOAT)n01(mt[mtindex]); 
        return a;
    }

    template <typename DFLOAT> DFLOAT normal(DFLOAT mean,DFLOAT sigma,int mtindex=0) {
        resizemt(mtindex);
        return mean + sigma * normal<DFLOAT>(mtindex);
    }
    
    template <typename DFLOAT> bool checkR(DFLOAT* a,DFLOAT* R,int n,int istart=1) {
        //Check sum of a[istart]^2 + ... + a[k]^2 < R[k] for all k=istart..n-1
        DFLOAT s;
        s = a[istart]*a[istart];
        if (s>R[istart]) return false;
        for (int i=istart+1;i<n;i++) {
            s += a[i]*a[i];
            if (s>R[i]) return false;
        }
        return true;
    }
    
    template <typename DFLOAT> bool isContainedDKLR(DFLOAT* y,DFLOAT* L,DFLOAT* R,int n,int l) {
        //Check: y[k] in [0,1] for k=1,...,l and
        //          L[k] < \sum y[i] < R[k] for k=1,...,l
        int i;
        DFLOAT sum=0;
        for (i=1;i<=l;i++) {
            if (y[i]<0) return false;
            if (y[i]>1) return false;
            sum += y[i];
            if (sum>R[i]) return false;
            if (sum<L[i]) return false;
        }
        for (i=l+1;i<=n;i++) {
            if (y[i]<0) return false;
            if (y[i]>1) return false;
        }
        return true;
    }

    template <typename DFLOAT> void SamplefromDKLR(DFLOAT* y,DFLOAT*L,DFLOAT* R,int n,int l,int mythread) {
        //Uniform sampling from the object (y[1],...,y[n]) s.t. 
        //  y[k] in [0,1] for k=1,...,n and L[k] < \sum_{i=1,....,k} y[i] < R[k] for k=1,...,n
        // by the hit-and-run algorithm
        DFLOAT* d = y + n + 1;
        DFLOAT tl,tu;
        DFLOAT s,ys;
        DFLOAT rr=0;
        int j=0;
        
        s = 0;
        //setting initial point in the object
        while (isContainedDKLR(y,L,R,n,l) == false) {
            y[1] = unique(L[1],R[1]);
            s = y[1];
            for (int i=2;i<=l;i++) {
                y[i] = unique<DFLOAT>(max(L[i]-s,0.0),max(0.95*R[i]-s,0.0),mythread);
                s += y[i];
            }
            for (int i=l+1;i<=n;i++) y[i] = unique<DFLOAT>(0,1,mythread);
            j++;
            if (j>10) {
                j=100;
            }
        }

        do {
            rr += 1.0;
            if (rr>3) {
                rr = 0.01;
                do {
                    y[1] = unique<DFLOAT>(L[1],R[1],mythread);
                    s = y[1];
                    for (int i=2;i<=l;i++) {
                        y[i] = unique<DFLOAT>(max(L[i]-s,0.0),0.9*R[i]-s,mythread);
                        s += y[i];
                    }
                    for (int i=l+1;i<=n;i++) y[i] = unique<DFLOAT>(0,1,mythread);
                    j++;
                    if (j>10) {
                        ccout << "SamplefromDKLR in samplingboost.cpp: index j error?" << endl; 
                        for (int i=1;i<=n;i++) {
                            ccout << "RL[" << i << "]=" << R[i] << " " << L[i] << endl;
                        }
                    }
                } while (isContainedDKLR(y,L,R,n,l) == false);
                rr = 0;
            }

            //random direction
            s = 0;
            for (int i=1;i<=n;i++) {
                d[i] = normal<DFLOAT>(mythread);
                s += d[i]*d[i];
            }
            s = sqrt(s);
            for (int i=1;i<=n;i++) {
                d[i] /= s;
            }

            tl = -999999;
            tu = 999999;

            s = 0;
            ys = 0;
            for (int i=1;i<=l;i++) {
                if (d[i]<0) {
                    tu = min(tu,-y[i]/d[i]);
                }
                if (d[i]>0) {
                    tl = max(tl,-y[i]/d[i]);
                }
                s += d[i];
                ys += y[i];
                if (s<0) {
                    tl = max(tl,(R[i]-ys)/s);
                    tu = min(tu,(L[i]-ys)/s);
                }
                if (s>0) {
                    tu = min(tu,(R[i]-ys)/s);
                    tl =  max(tl,(L[i]-ys)/s);
                }
            }

            for (int i=l+1;i<=n;i++) {
                if (d[i]<0) {
                    tu = min(tu,-y[i]/d[i]);
                }
                if (d[i]>0) {
                    tl = max(tl,-y[i]/d[i]);
                }
                if (d[i]<0) {
                    tl = max(tl,(1-y[i])/d[i]);
                }
                if (d[i]>0) {
                    tu = min(tu,(1-y[i])/d[i]);
                }
            }
            tl = unique<DFLOAT>(tl,tu,mythread);
            for (int i=1;i<=n;i++) y[i] += tl * d[i]; 
        } while (isContainedDKLR(y,L,R,n,l) == false);

    }
    
    template <typename DFLOAT> void SamplefromCpKLR(DFLOAT* yD,DFLOAT* yC,DFLOAT* L,DFLOAT* R,int n,int mythread) {
        int i;
        DFLOAT s,t = 0,u;
        SamplefromDKLR(yD,L,R,n,n,mythread);
        for (i=1;i<=n;i++) {
            t = unique<DFLOAT>(0.0,2*3.1415926535897932,mythread);
            s = sqrt(yD[i]);
            yC[2*i] = s * cos(t);
            yC[2*i-1] = s * sin(t);
        }
    }

    template <typename DFLOAT> void SamplefromCpKLR_plane(DFLOAT* yD,DFLOAT* yC,DFLOAT* L,DFLOAT* R,int n,int mythread) {
        int i;
        DFLOAT s,t = 0,u;
        SamplefromDKLR(yD,L,R,n-1,n-1,mythread);
        u = R[n];
        for (i=1;i<=n-1;i++) {
            t = unique<DFLOAT>(0.0,2*3.1415926535897932,mythread);
            s = sqrt(yD[i]);
            u -= yD[i];
            yC[2*i] = s * cos(t);
            yC[2*i-1] = s * sin(t);
        }
        s = sqrt(u);
        yC[2*i] = s * cos(t);
        yC[2*i-1] = s * sin(t);
    }

    template <typename DFLOAT> inline void Samplefromsphere(DFLOAT* y,int n,DFLOAT radius=1.0,int mtindex=0) {
        resizemt(mtindex);
        DFLOAT sum=0;
        for (int i=1;i<=n;i++) {
            y[i] = normal<DFLOAT>(mtindex);
            sum += y[i]*y[i];
        }
        sum=sqrt(sum);
        for (int i=1;i<=n;i++) {
            y[i] *= radius / sum ;
        }
    }
}

#define define_even
    #include "samplingboostcore.cpp"
#undef define_even

#define define_odd
    #include "samplingboostcore.cpp"
#undef define_odd


//Sampling from objects
namespace sampling_tools{

    int initialized = 0;
    
    template <typename DFLOAT> DFLOAT ApproxRatio(std::vector<DFLOAT>& rd,int istart,int iend,int parallel=1,int sampling_mult=10) {
        
        int dim = iend-istart+1;
        int numsamples = sampling_mult * dim * dim;
        if (dim%2==0) {
            return ApproxRatio_even<double,DFLOAT>(rd,istart,iend,numsamples,parallel);
        } else {
            return ApproxRatio_odd<double,DFLOAT>(rd,istart,iend,numsamples,parallel);
        }
    }
    
};

#endif
