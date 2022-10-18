#ifndef _inc_volume_and_probability_wrapper
#define _inc_volume_and_probability_wrapper


#include "memorytools.cpp"
#include "bkzconstants.hpp"
using namespace boost::multiprecision;

#include <lattice/samplingboost.cpp>


#define opt_upper 0x01
#define opt_lower 0x02

#define opt_volume 0x01 
#define opt_volume_prob 0x02
#define opt_surface 0x03
#define opt_surface_prob 0x04
#define opt_gaussian_prob 0x05

typedef std::vector<std::vector<std::vector<mpfr_float> > > EStemp;

#define lexical_convert

namespace EvenSimplex {

    std::vector<EStemp*> Flarge;
    
    template <typename T> int GetESprecision(int dim,T uprob) {
        //uprob is guaranteed probability bound
        int prec = 1 + dim / 2.1;   //for full-enum
        if (uprob <= 0.1) {
            prec = min(prec,(int)(2 + dim/3.6));
        }
        if (prec < 4) return 4;
        return prec;
    }
    
    EStemp* GetEStemp(int dim,int prec) {
        int precisionback = mpfr_float::default_precision() ;
        
        mpfr_float::default_precision(prec);
        int slot = prec;
        if (slot <0) {
            ccout << "GetEStemp: precision error" << endl;
            return NULL;
        }
        int cs = Flarge.size();
        if (slot >= cs) {
            Flarge.resize(slot+1);
            for (int i=cs;i<=slot;i++) Flarge[i] = NULL;
        }
        if (Flarge[slot]==NULL) Flarge[slot] = new EStemp();
        EStemp* ret = Flarge[slot];

        ret->resize(dim+2);
        for (int i=0;i<=dim;i++) (*ret)[i].resize(dim+2);
        mpfr_float::default_precision(precisionback); 
        return ret;
    }
}

template <typename T,typename T2> T EvalAt(T* F,T2 a,int n) {
    //Assume a be a small number
    T ret,t;
    ret = 0;
    if (a==0) return F[0];

    int i;
    t = 1;
    for (i=0;i<n;i++) {
        ret += t * F[i];
        t *= T(a);
    }
    return ret;
}

template <typename T,typename T2> void Integ(T* ret, T* F,T2 low,int n) {    
    if (n<0) return;
    int i;
    for (i=n-1;i>=0;i--) {
        ret[i+1] = F[i] / (i+1);
    }
    ret[0] = 0;
    ret[0] = -EvalAt(ret,low,n+1);
    return;
}

template <typename T> void diffpoly(T* ret,T* F,int n) {    
 
    for (int i=0;i<n;i++) {
        ret[i] = F[i+1] * (i+1);
    }
    
}
    
//#include <boost/math/special_functions/binomial.hpp>
template <typename CFLOAT> void SetBinomial(CFLOAT *bin,int deg,CFLOAT a) {
    //Set bin = (x+a)^deg
    for (int i=0;i<=deg;i++) {
        bin[i] = pow(a,deg-i) * (CFLOAT)bkzconstants::binary_coeff(deg,i);
    }
}

template <typename DFLOAT,typename CFLOAT> DFLOAT EvenSimplexVolumeWrapperMain(std::vector<DFLOAT>& rd,int dim,int option,int parallel=1,int usedprecision=-1,DFLOAT uprob=1,DFLOAT* diag=0) {

    //Computing the volume of n-dimensional trancated simplex
    //{ (x1,...,xn) : \sum[i=1,...,l] x_i < R_l for l=1,...,n}
    //F[n,n] is the volume
    //pb is the probability that Pr{x <- surface of full-simplex}[x is in trancated simplex]

    if (dim==1) return rd[1];

    //SetPrecision
    int prec;
    int precisionback;

    bool mpfr_flag = false; 
    if ((typeid(CFLOAT) != typeid(double)) && (typeid(CFLOAT) != typeid(long double)))  {
        if (usedprecision==-1) {
            prec = EvenSimplex::GetESprecision(dim,uprob);
        } else {
            prec = usedprecision;
        }
        precisionback = mpfr_float::default_precision() ;
        mpfr_float::default_precision(prec);
        mpfr_flag = true;
    }
    
    
    CFLOAT** F = (CFLOAT**)shared_memory::allocate2<CFLOAT>(14002,dim+1,dim+1);
    CFLOAT* bin = (CFLOAT*)shared_memory::allocate1<CFLOAT>(3332,dim+1);
    CFLOAT* CR = (CFLOAT*)shared_memory::allocate1<CFLOAT>(3333,dim+1);
    for (int i=1;i<=dim;i++) CR[i] = boost::lexical_cast<CFLOAT>(rd[i]);
    for (int i=0;i<dim;i++) F[0][i] = 0;
    

    if (option==opt_volume_prob) {
        F[0][dim] = 1.0;   //computes F[d,d] * d!, probability
    } else 
    if (option==opt_surface_prob) {
        F[0][dim] = 1.0/dim;   //computes F[d,d] * d!, probability
    }
    if (option==opt_volume) {
        F[0][dim] =  (CFLOAT)1.0/(CFLOAT)bkzconstants::factorial(dim); 
    } else 
    if (option==opt_surface) {
        F[0][dim] =  (CFLOAT)1.0/(CFLOAT)bkzconstants::factorial(dim); 
    }

    for (int j=1;j<=dim;j++) {
        CFLOAT ht = F[j-1][dim-(j-1)];
        CFLOAT pw = 1;
        for (int k=dim-j+1;k>=0;k--) {
#ifdef lexical_convert
            F[j][k] = F[j-1][k] - ht * pw * boost::lexical_cast<CFLOAT>(bkzconstants::binary_coeff(dim-j+1,k));
#else
            F[j][k] = F[j-1][k] - ht * pw * (CFLOAT)bkzconstants::binary_coeff(dim-j+1,k);
#endif
            if (k>0) pw *= (-CR[j]);
        }
    }

    for (int i=1;i<dim;i++) {
        //When CFLOAT=mpfr_float, the above code is extremely slow
#ifdef lexical_convert
        F[i][0] = F[i][dim-i] * boost::lexical_cast<CFLOAT>(bkzconstants::factorial(dim-i));
#else
        F[i][0] = F[i][dim-i] * (CFLOAT)bkzconstants::factorial(dim-i);
#endif
    }
    if (diag!=0) {
        for (int i=1;i<=dim;i++) diag[i] = boost::lexical_cast<DFLOAT>(F[i][0]); 
    }
    DFLOAT ret;

    if ((option==opt_volume) || (option==opt_volume_prob)) {
        ret = boost::lexical_cast<DFLOAT>(F[dim][0]);
    }
    if ((option==opt_surface) || (option==opt_surface_prob)) {
        ret = boost::lexical_cast<DFLOAT>(F[dim-1][0]);
    }

    //if CFLOAT==mpfr_float
    if (mpfr_flag == true) {
        mpfr_float::default_precision(precisionback); 
    }
    return ret;
}
      
template <typename DFLOAT> DFLOAT EvenSimplexVolumeWrapper(std::vector<DFLOAT> rd,int dim,int option,int parallel,DFLOAT uprob=1,DFLOAT* diag=0) {
    //need to recompute thresholds 
    if (uprob > 0.1) { 
        if (dim<=35) {
            return EvenSimplexVolumeWrapperMain<DFLOAT,double>(rd,dim,option,parallel,-1,uprob,diag);
        } else 
        if (dim<=40) {
            return EvenSimplexVolumeWrapperMain<DFLOAT,long double>(rd,dim,option,parallel,-1,uprob,diag);
        } else
        if (dim<=62) {
            return EvenSimplexVolumeWrapperMain<DFLOAT,float10>(rd,dim,option,parallel,-1,uprob,diag);
        } else
        if (dim<=95) {
            return EvenSimplexVolumeWrapperMain<DFLOAT,float30>(rd,dim,option,parallel,-1,uprob,diag);
        } else
        if (dim<=165) {
            return EvenSimplexVolumeWrapperMain<DFLOAT,float50>(rd,dim,option,parallel,-1,uprob,diag);
        } else
        if (dim<=215) {
            return EvenSimplexVolumeWrapperMain<DFLOAT,float80>(rd,dim,option,parallel,-1,uprob,diag);
        } else
        if (dim<=265) {
            return EvenSimplexVolumeWrapperMain<DFLOAT,float100>(rd,dim,option,parallel,-1,uprob,diag);
        } else
        if (dim<=360) {
            return EvenSimplexVolumeWrapperMain<DFLOAT,float150>(rd,dim,option,parallel,-1,uprob,diag);
        } else 
        if (dim<=460) {
            return EvenSimplexVolumeWrapperMain<DFLOAT,float200>(rd,dim,option,parallel,-1,uprob,diag);
        } else 
        if (dim<=700) {
            return EvenSimplexVolumeWrapperMain<DFLOAT,float300>(rd,dim,option,parallel,-1,uprob,diag);
        } else {
  
            int usedprecision = EvenSimplex::GetESprecision(dim,uprob);
            bkzconstants::factorial(dim+1);
            DFLOAT ret;
            #pragma omp critical
            {
                ret = EvenSimplexVolumeWrapperMain<DFLOAT,mpfr_float>(rd,dim,option,parallel,usedprecision,uprob,diag);
            }
            return ret;
            }
    } else {
        if (dim<=54) {
            return EvenSimplexVolumeWrapperMain<DFLOAT,double>(rd,dim,option,parallel,-1,uprob,diag);
        } else 
        if (dim<=70) {
            return EvenSimplexVolumeWrapperMain<DFLOAT,long double>(rd,dim,option,parallel,-1,uprob,diag);
        } else 
        if (dim<=100) {
            return EvenSimplexVolumeWrapperMain<DFLOAT,float10>(rd,dim,option,parallel,-1,uprob,diag);
        } else 
        if (dim<=160) {
            return EvenSimplexVolumeWrapperMain<DFLOAT,float30>(rd,dim,option,parallel,-1,uprob,diag);
        } else 
        if (dim<=260) { 
            return EvenSimplexVolumeWrapperMain<DFLOAT,float50>(rd,dim,option,parallel,-1,uprob,diag);
        } else 
        if (dim<=440) {
            return EvenSimplexVolumeWrapperMain<DFLOAT,float100>(rd,dim,option,parallel,-1,uprob,diag);
        } else 
        if (dim<=620) {
            return EvenSimplexVolumeWrapperMain<DFLOAT,float150>(rd,dim,option,parallel,-1,uprob,diag);
        } else 
        if (dim<=800) {
            return EvenSimplexVolumeWrapperMain<DFLOAT,float200>(rd,dim,option,parallel,-1,uprob,diag);
        } else 
        if (dim<=1160) {
            return EvenSimplexVolumeWrapperMain<DFLOAT,float300>(rd,dim,option,parallel,-1,uprob,diag);
        } else {
            int usedprecision = EvenSimplex::GetESprecision(dim,uprob);
            DFLOAT ret;
            bkzconstants::factorial(dim+1);
            #pragma omp critical
            {
                ret =  EvenSimplexVolumeWrapperMain<DFLOAT,mpfr_float>(rd,dim,option,parallel,usedprecision,uprob,diag);
            }
            return ret;
        }
    }
 
}

namespace pruning_func {

    template <typename DFLOAT> DFLOAT RigidBoundProb(std::vector<DFLOAT>& rd,int istart,int dim,int option,int parallel,DFLOAT uprob=1) {
        //istart=1 is default parameter
        int n = dim/2+1;
    
        std::vector<DFLOAT> tt;
        tt.resize(n+1);

        if (option==opt_lower) {
            n = dim/2;
            if (dim%2==1) n = (dim-1)/2;
            for (int i=1;i<=n;i++) tt[i] = rd[2*i-1+(istart-1)];
        } else {
            n = dim/2;
            if (dim%2==1) n = (dim-1)/2;
            for (int i=1;i<=n;i++) tt[i] = rd[2*i+(istart-1)];
            if (dim%2==1) {
                n++;
                tt[n] = 1.0;
            }
        }
        return EvenSimplexVolumeWrapper<bkzfloat>(tt,n,opt_surface_prob,parallel,uprob);
    }

    template <typename DFLOAT> DFLOAT RigidBoundVolume(std::vector<DFLOAT>& rd,int istart,int dim,int option,int parallel,DFLOAT uprob=1) {
        //istart=1 is default parameter
        int n = dim/2+1;
        std::vector<DFLOAT> tt;
        tt.resize(n+1);

        if (option==opt_lower) {
            n = dim/2;
            if (dim%2==1) n = (dim-1)/2;
            for (int i=1;i<=n;i++) tt[i] = rd[2*i-1+(istart-1)];
        } else {
            n = dim/2;
            if (dim%2==1) n = (dim-1)/2;
            for (int i=1;i<=n;i++) tt[i] = rd[2*i+(istart-1)];
            if (dim%2==1) {
                n++;
                tt[n] = 1.0;
            }
        }
        return EvenSimplexVolumeWrapper<bkzfloat>(tt,n,opt_volume_prob,parallel,uprob);
    }

    template <typename DFLOAT,typename IFLOAT> DFLOAT RigidBoundCost(std::vector<DFLOAT>& lcost,std::vector<DFLOAT>& rd,std::vector<IFLOAT>& c,DFLOAT clim,int istart,int dim,int option,int parallel=1,DFLOAT uprob=1) {
        //pruning_func[1,...,n]
        //c[1,...,n]
        DFLOAT sum = 0;
        DFLOAT t = 1;
        DFLOAT t2 = 1;
        
        lcost.resize(dim/2+1);
        
        int prec = EvenSimplex::GetESprecision(dim/2,uprob);
        
        std::vector<DFLOAT> R2;
        R2.resize(1+dim/2);

        int shift = 0;
        if (option==opt_lower) shift = 1;
        
        int oddflag = 0;
        if (dim%2==1) oddflag = 1;
        
        for (int i=1;i<=dim/2;i++) {
            R2[i] = rd[(istart-1)+2*i-shift]/rd[istart+dim-1];
        }
        DFLOAT* diag = (DFLOAT*)shared_memory::allocate1<DFLOAT>(3330,dim/2+1);
        EvenSimplexVolumeWrapper<bkzfloat>(R2,dim/2,opt_volume,parallel,uprob,diag) ;
        t2 = 1;
        sum = 0;
        for (int i=1;i<=dim/2;i++) {
            t = diag[i]*(DFLOAT)bkzconstants::factorial(i);  //probability of the area/simplex
            lcost[i] = t;
            t *= bkzconstants::vol_unit_ball(i*2); //volume of ball
            t2 *= clim * clim / (DFLOAT)c[istart+dim-1-(i-1)*2] / (DFLOAT)c[istart+dim-2-(i-1)*2];
            DFLOAT tt = t*t2;
            sum += tt*2 ;    //Summing up depth[2*i] and depth[2*i+1]
        }
        sum *= 2;   //Actual numbe of nodes is doubled since ENUM have to check a node is inside of cylinder-intersection
        return sum*0.5; 
    }


    template <typename DFLOAT> void RigidBoundCostFactor(std::vector<DFLOAT>& lcost,std::vector<DFLOAT>& rd,int istart,int dim,int option,int parallel,DFLOAT uprob=1) {
        //pruning_func[1,...,n]
        //c[1,...,n]
        DFLOAT sum = 0;
        DFLOAT t = 1;
        DFLOAT t2 = 1;
        
        lcost.resize(dim/2+1);
        
        int prec = EvenSimplex::GetESprecision(dim/2,uprob);
        
        std::vector<DFLOAT> R2;
        R2.resize(1+dim/2);
        for (int k=2;k<=dim;k+=2) {
            bool exec_condition = false;
            if ((option==opt_upper) && (rd[k]!=0)) exec_condition = true;
            if ((option==opt_lower) && (rd[k-1]!=0)) exec_condition = true;
            if (exec_condition) {
                int shift = 0;
                if (option==opt_lower) shift = 1;
                for (int i=1;i<=k/2;i++) {
                    R2[i] = rd[(istart-1)+2*i-shift]/rd[(istart-1)+k-shift];
                }
                t = EvenSimplexVolumeWrapper<bkzfloat>(R2,k/2,opt_volume_prob,parallel,1.0) ;
                t *= pow(rd[(istart-1)+k-shift],0.5*k);
                lcost[k/2] = t;;
            }
        }
    }

    //Below rd[1..n] is squared values, 0 <= rd[i] <= 1 for all i
    //Monotonically increasing
    template <typename DFLOAT> DFLOAT Rigid_upper_prob(std::vector<DFLOAT>& rd,int istart,int dim,int parallel=1,DFLOAT uprob=1) {
        if (dim==1) return rd[1];
        return RigidBoundProb(rd,istart,dim,opt_upper,parallel,uprob);
    }

    template <typename DFLOAT> DFLOAT Rigid_lower_prob(std::vector<DFLOAT>& rd,int istart,int dim,int parallel=1,DFLOAT uprob=1) {
        if (dim==1) return rd[1];
        return RigidBoundProb(rd,istart,dim,opt_lower,parallel,uprob);
    }

    template <typename DFLOAT> DFLOAT Rigid_upper_volume(std::vector<DFLOAT>& rd,int istart,int dim,int parallel=1,DFLOAT uprob=1) {
        if (dim==1) return rd[1];
        return RigidBoundVolume(rd,istart,dim,opt_upper,parallel,uprob);
    }

    template <typename DFLOAT> DFLOAT Rigid_lower_volume(std::vector<DFLOAT>& rd,int istart,int dim,int parallel=1,DFLOAT uprob=1) {
        if (dim==1) return rd[1];
        return RigidBoundVolume(rd,istart,dim,opt_lower,parallel,uprob);
    }

    template <typename DFLOAT,typename IFLOAT> DFLOAT Rigid_upper_cost(std::vector<DFLOAT>& lcost,std::vector<DFLOAT>& rd,std::vector<IFLOAT>& c,DFLOAT clim,int istart,int dim,int parallel=1,DFLOAT uprob=1) {
        return RigidBoundCost(lcost,rd,c,clim,istart,dim,opt_upper,parallel,uprob);
    }

    template <typename DFLOAT,typename IFLOAT> DFLOAT Rigid_upper_cost(std::vector<DFLOAT>& rd,std::vector<IFLOAT>& c,DFLOAT clim,int istart,int dim,int parallel=1,DFLOAT uprob=1) {
        std::vector<DFLOAT> lcost;
        return RigidBoundCost(lcost,rd,c,clim,istart,dim,opt_upper,parallel,uprob);
    }

    template <typename DFLOAT> void Rigid_upper_costfactor(std::vector<DFLOAT>& lcost,std::vector<DFLOAT>& rd,int istart,int dim,int parallel=1,DFLOAT uprob=1) {
        RigidBoundCostFactor(lcost,rd,istart,dim,opt_upper,parallel,uprob);
    }

    template <typename DFLOAT> void Rigid_lower_costfactor(std::vector<DFLOAT>& lcost,std::vector<DFLOAT>& rd,int istart,int dim,int parallel=1,DFLOAT uprob=1) {
        RigidBoundCostFactor(lcost,rd,istart,dim,opt_lower,parallel,uprob);
    }
    
    template <typename DFLOAT,typename IFLOAT> DFLOAT Rigid_upper_cost_general(std::vector<DFLOAT>& lcost,std::vector<DFLOAT>& rd,std::vector<IFLOAT>& c,int istart,int dim,int parallel=1,DFLOAT uprob=1) {
        DFLOAT clim=0;
        for (int i=istart;i<istart+dim;i++) clim = max(clim,rd[i]); 
        std::vector<DFLOAT> rr;
        rr.resize(istart+dim+1);
        for (int i=istart;i<istart+dim;i++) rr[i] = rd[i]/clim; 
        return RigidBoundCost(lcost,rr,c,(DFLOAT)sqrt(clim),istart,dim,opt_upper,parallel,uprob);
    }

    template <typename DFLOAT,typename IFLOAT> DFLOAT Rigid_upper_cost_general(std::vector<DFLOAT>& rd,std::vector<IFLOAT>& c,int istart,int dim,int parallel=1,DFLOAT uprob=1) {
        std::vector<DFLOAT> lcost;
        DFLOAT clim=0;
        for (int i=istart;i<istart+dim;i++) clim = max(clim,rd[i]); 
        std::vector<DFLOAT> rr;
        rr.resize(istart+dim+1);
        for (int i=istart;i<istart+dim;i++) rr[i] = rd[i]/clim; 
        return RigidBoundCost(lcost,rr,c,(DFLOAT)sqrt(clim),istart,dim,opt_upper,parallel,uprob);
    }

    template <typename DFLOAT,typename IFLOAT> DFLOAT Rigid_lower_cost(std::vector<DFLOAT>& lcost,std::vector<DFLOAT>& rd,std::vector<IFLOAT>& c,DFLOAT clim,int istart,int dim,int parallel=1,DFLOAT uprob=1) {
        return RigidBoundCost(lcost,rd,c,clim,istart,dim,opt_lower,parallel,uprob);
    }

    template <typename DFLOAT,typename IFLOAT> DFLOAT Rigid_lower_cost(std::vector<DFLOAT>& rd,std::vector<IFLOAT>& c,DFLOAT clim,int istart,int dim,int parallel=1,DFLOAT uprob=1) {
       std::vector<DFLOAT> lcost;
       return RigidBoundCost(lcost,rd,c,clim,istart,dim,opt_lower,parallel,uprob);
    }

    template <typename DFLOAT> DFLOAT Approx_prob(std::vector<DFLOAT>& rd,int istart,int dim,DFLOAT uprob=1,int parallel=1,int sampling_mult=10) {
        if (dim==1) return rd[1];
        DFLOAT upper = RigidBoundProb(rd,istart,dim,opt_upper,parallel,uprob);
        DFLOAT ratio = sampling_tools::ApproxRatio(rd,istart,istart+dim-1,parallel,sampling_mult);
        return upper * ratio;
    }
}

#endif