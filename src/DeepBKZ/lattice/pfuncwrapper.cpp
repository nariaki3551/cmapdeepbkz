#ifndef _inc_pfuncwrapper_cpp
#define _inc_pfuncwrapper_cpp

#include "bkzlibwrapper.hpp"
#include "volandprobwrapper.cpp"

#define pf_fromtable 0x00
#define pf_lower 0x01
#define pf_systematic 0x02

#define pf_crossentropy 0x05
#define pf_crossentropy_exact 0x06


#define target_prob 0x01
#define target_volume 0x02

struct PruningFunction {
    //structure to represent pruning function
    //rd[i] is (c.R[i])^2
    //rd[1..n] is data , in increasing form
    //cost[i..n] saves cost factors in some of subroutines; these are from simplex volumes
    std::vector<bkzfloat> rd;
    std::vector<bkzfloat> costfactor;   //Cost of C_k

    PruningFunction& operator=(const PruningFunction& A) {
        rd.resize(A.rd.size());
        for (int i=0;i<A.rd.size();i++) rd[i] = A.rd[i];
        costfactor.resize(A.costfactor.size());
        for (int i=0;i<A.costfactor.size();i++) costfactor[i] = A.costfactor[i];
        return *this;
    }
    
    void display() {
        ccout << "[";
        for (int i=1;i<rd.size();i++) {
            ccout << rd[i] << " ";
        }
        ccout << "]" << endl;
    }
    void savetofile(std::string fname) {
        std::ofstream of;
        of.open(fname.c_str(),ios::trunc);
        for (int i=1;i<rd.size();i++) {
            of << i << "\t" << rd[i] << endl;
        }
        of.close();
    }
    void loadfromfile(std::string fname) {
        std::ifstream ifs;
        ifs.open(fname.c_str());
        do {
            if (ifs.eof()==true) break;
            int i;
            double rr;
            ifs >> i;
            ifs >> rr;
            rd.resize(i+1);
            rd[i] = rr;
        } while (1);
        ifs.close();
    }
    
};

template <typename T> T r_to_delta(T r,int dim) {
    return exp(-log(r)*(dim-1.0)/(4.0*dim)); 
}

template <typename T> T delta_to_r(T delta,int dim) {
    return pow(delta,-4.0*dim/(dim-1)); 
}

template <typename T> double Detect_GSAr(std::vector<T>& cd,int istart,int iend,int option=INPUT_SQUARED) {
    //Estimating delta by least square
    T xs,ys,xxs;
    double r;
    int i;
    int dim = iend-istart+1;
    xs = dim*(dim-1.0)/2.0;
    xxs = dim*(dim-1.0)*(2.0*dim-1.0)/6.0;
    ys = 0;
    for (int i=istart;i<=iend;i++) {
        ys += log(cd[i]);
    }
    r = 0;
    for (int i=istart;i<=iend;i++) r += (double)log(cd[i]) * (i-istart);
    if (option==INPUT_SQUARED) {
        ys *= 0.5;
        r *= 0.5;
    }
    r = (1.0*dim*r - (double)(xs*ys)) / (1.0*dim*(double)xxs-(double)(xs*xs));
    r = exp(2.0*r);
    return r;
}


template <typename T> bkzfloat Detect_delta(std::vector<T>& cd,int istart,int iend,int option=INPUT_SQUARED) {
    
    int dim = iend - istart + 1;
    
    bkzfloat r = Detect_GSAr(cd,istart,iend,option);
    bkzfloat ret = r_to_delta(r,dim);  //convert to root Hermite Factor
    
    //range check
    if ((ret<=0.5) || (ret>=1.5) || (boost::math::isnan<bkzfloat>(ret)==true)) {
        if (dim>10) {
            ccout << "delta-error: " << ret << " " << r <<  endl;
            ccout << istart << " : " << iend << endl;
            for (int i=istart;i<=iend;i++) {
                ccout << "cd[" << i << "]=" << cd[i] << endl;
            }
        }
    }
    return  ret;
}

template <typename T,typename T2,typename T3> T FullENUMCostbase(std::vector<T2>& c,std::vector<T>& partcost,int istart,int iend,T3 radius,int opt=INPUT_NONSQUARED) {
    //c[1..n] is |b*i|
    T ret = 0;
    T t = 1;
    partcost.resize(iend+1);
    for (int i=iend;i>=istart;i--) {
        if (opt==INPUT_SQUARED) {
            t *= (T)(radius / sqrt(c[i]));
        } 
        if (opt==INPUT_NONSQUARED) {
            t *= (T)(radius / c[i]);
        } 
        partcost[i] = t * (T)bkzconstants::vol_unit_ball(iend-i+1);
        ret += partcost[i];
    }    
    return ret*0.5;
}

template <typename T,typename T2,typename T3> T FullENUMCost(std::vector<T2>& c,int istart,int iend,T3 radius,int opt=INPUT_NONSQUARED) {
    std::vector<T> partcost;   //dummy
    return FullENUMCostbase(c,partcost,istart,iend,radius,opt);
}

template <typename T,typename T2,typename T3> T FullENUMCost(LatticeBasis<T2>& B,T3 radius) {
    B.updateGSBasis();
    return FullENUMCost<T,T2,T3>(B.gs.c,1,B.gs.c.size()-1,radius);
}

template <typename T,typename T2> T FullENUMCost(std::vector<T2>& c,int istart,int iend,int opt=INPUT_NONSQUARED) {
    T radius = lattice_tools::LatticeGH(c,istart,iend,opt);
    //ccout << "radius=" << radius << endl;
    return FullENUMCost<T>(c,istart,iend,radius,opt);

}

#define _enable_pfcache
//int chit = 0;
//int cmiss=0;

namespace pruning_func {

    int initpf = 0;
    std::string pfsavedir;
    std::vector<std::vector<double> > dtable;     //table of pruning function 
    
    void init_pruning_func() {
        
        if (initpf>0) return;
        
        if (FileExists("pfdata2table.dat")==true) {
            loaddoublevector(dtable,"pfdata2table.dat");
        } else {
            ccout << "pfdata2table.dat not found!" << endl;
            exit(0);
        }
        pfsavedir =  makefullpath(ReadConf("bkz.conf","pfcache"));
        mkdirRecursive(pfsavedir.c_str(), 0777);
        initpf = 1;
    }
    
    
    template <typename T> void getprobUL(T p0,T& pu,T& pl) {
        if (p0>1e-1)  { pu = 1e-0; pl=1e-1; return; }
        if (p0>1e-2)  { pu = 1e-1; pl=1e-2; return; }
        if (p0>1e-3)  { pu = 1e-2; pl=1e-3; return; }
        if (p0>1e-4)  { pu = 1e-3; pl=1e-4; return; }
        if (p0>1e-6)  { pu = 1e-4; pl=1e-6; return; }
        if (p0>1e-9)  { pu = 1e-6; pl=1e-9; return; }
        if (p0>1e-12) { pu = 1e-9; pl=1e-12; return; }
        if (p0>1e-24) { pu = 1e-12; pl=1e-24; return; }
        if (p0>1e-48) { pu = 1e-24; pl=1e-48; return; }
        pu = 1e-48; pl=0; return;
    }

    template <typename T> void getdeltaUL(T d0,T& du,T& dl) {

        if (d0 < 1.005) d0 = 1.005;
        if (d0 > 1.02) d0 = 1.02;
        if (d0<1.01) {du=1.01;dl=1.005; return;}
        if (d0<1.015) {du=1.015;dl=1.01; return;}
        if (d0<=1.02) {du=1.02;dl=1.015; return;}
        return;
    }

    void getdimUL(int d0,int& du,int& dl) {
    
        d0 = max(20,d0);
        d0 = min(300,d0);
        if (d0%20==0) {du=dl=d0; return;};
        
        dl = int(d0/20)*20;
        du = dl + 20;
    }
    
    
    template <typename T> int getpf2probindex(T prob) {
        int pow = round( boost::lexical_cast<double>(-log(prob)/log(10.0)));
        int ret = 0;
        if ((1<=pow) && (pow<=4)) ret = pow - 1;
        else if (pow==6) ret = 4;
        else if (pow==9) ret = 5;
        else if (pow==12) ret = 6;
        else if (pow==24) ret = 7;
        else if (pow==48) ret = 8;
        return ret;
    }
    int getpf2dimindex(int dim) {
        return (dim-20)/20;
    }
    template <typename T> int getpf2deltaindex(T delta) {
        int ret = (int)( boost::lexical_cast<double>((delta-1.000+0.0025)/0.005)); 
        return ret;
    }

    template <typename T> int getpf2index(T prob,int dim,T delta) {

        int i,j,k;
        i = getpf2probindex(prob);
        j = getpf2dimindex(dim);
        k = getpf2deltaindex(delta);
        return i + j * 9 + k * 15 * 9;
    }

    template <typename T> void getfuncfromtable(std::vector<std::vector<double> >& pf,std::vector<double>& R,T prob,int dim,T delta) {
        int j;
        R.resize(dim+3);
        if (prob==1) {
            for (j=0;j<dim+3;j++) R[j] = 1.0;
            return;
        }
        if (prob==0) {
            for (j=0;j<dim+3;j++) R[j] = 0.0;
            return;
        }
        int i = getpf2index(prob,dim,delta);
        R.resize(pf[i].size());
        for (j=0;j<(int)pf[i].size();j++) R[j] = pf[i][j];
    }

    template <typename DFLOAT> DFLOAT SetApproximatePruningFunction(std::vector<DFLOAT>& rd,int dim,DFLOAT delta,DFLOAT tprob) {

        pruning_func::init_pruning_func();
        using namespace pruning_func;

        double startcputime = clock();
        DFLOAT pp;      //final probability
        //Input: cd[1..n] = |b*i|^2
        //       probability
        //Output: rd[1..n]: R[i]^2, increasing form
        delta = round(delta * 20000) / 20000.0;
        
#ifdef _enable_pfcache
        //Cache memory
        std::size_t ch = std::hash<double>{}(dim);
        ch += std::hash<double>{}((double)delta);
        ch += std::hash<double>{}((double)boost::math::log1p(tprob));
        int cindex = ch%20000;
        
        DFLOAT* pfcache = (DFLOAT*)shared_memory::allocate1<DFLOAT>(50000+cindex,3+dim);
        if (pfcache[0] == ch) {
            for (int i=1;i<=dim;i++) rd[i] = pfcache[i+1];
            return pfcache[1];
        }
#endif        
        
        //For small dimension, it reuses the function for dim=20
        if (dim<20) {
            if (dim<=10) {
                int i;
                for (i=1;i<=dim;i++) rd[i] = 1.0;
                return 1.0;
            } else {
                DFLOAT pp;
                rd.resize(21);
                pp = SetApproximatePruningFunction(rd,20,delta,tprob);
                rd.resize(dim+1);
                return pp;
            }            
        }
        DFLOAT pu,pl;
        DFLOAT du = 0,dl = 0;
        double a;   //adjusting coefficient in odd dimension
        int dimu,diml;

        getprobUL(tprob,pu,pl);
        getdeltaUL(delta,du,dl);
        getdimUL(dim,dimu,diml);

        std::vector<double> Ruuu,Ruul,Rulu,Rull;        //{u,l}^3 is corresponding to prob,dim,delta resp.
        std::vector<double> Rluu,Rlul,Rllu,Rlll;        //{u,l}^3 is corresponding to prob,dim,delta resp.
        DFLOAT d0 = delta;      //for readability
        //Read basic function data from the table
        //Note: Ruuu and other values are stored in [0..n-1]

        //for upper prob
        getfuncfromtable(dtable,Ruuu,pu,dimu,du);
        getfuncfromtable(dtable,Ruul,pu,dimu,dl);
        a = boost::lexical_cast<double>((d0-dl) / (du-dl));

        for (int i=0;i<(int)Ruuu.size();i++) Ruuu[i] = Ruul[i] * (1.0-a) + Ruuu[i] * a;   
            if (diml==dimu) {
            Rulu.resize(Ruuu.size());    
            for (int i=0;i<(int)Ruuu.size();i++) Rulu[i] = Ruuu[i];   
        } else {
            getfuncfromtable(dtable,Rulu,pu,diml,du);
            getfuncfromtable(dtable,Rull,pu,diml,dl);
            a = boost::lexical_cast<double>((d0-dl) / (du-dl));
            for (int i=0;i<(int)Rulu.size();i++) Rulu[i] = Rull[i] * (1.0-a) + Rulu[i] * a;   
        }

        if (pl>0) {
            //the same operation to the lower prob
            getfuncfromtable(dtable,Rluu,pl,dimu,du);
            getfuncfromtable(dtable,Rlul,pl,dimu,dl);
            a = boost::lexical_cast<double>((d0-dl) / (du-dl));
            for (int i=0;i<(int)Rluu.size();i++) Rluu[i] = Rlul[i] * (1.0-a) + Rluu[i] * a;   

            if (diml==dimu) {
                Rllu.resize(Rluu.size());    
                for (int i=0;i<(int)Rluu.size();i++) Rllu[i] = Rluu[i];   
            } else {
                getfuncfromtable(dtable,Rllu,pl,diml,du);
                getfuncfromtable(dtable,Rlll,pl,diml,dl);
                a = boost::lexical_cast<double>((d0-dl) / (du-dl));
                for (int i=0;i<(int)Rllu.size();i++) Rllu[i] = Rlll[i] * (1.0-a) + Rllu[i] * a;   
            }
        }
        //Here Ruuu and Rulu are interpolated function for upper prob
        //     Rluu and Rllu are interpolated function for lower prob
        //Next, interpolate the dimension factor

        //for upper prob
        std::vector<DFLOAT> RU,RL;
        RU.resize(dim+1);
        RL.resize(dim+1);

        double au,al = 0;
        int shift=0;    //=1 if dim is odd
        if (dim%2==1) shift = 1;
        if (diml==dimu) {
            if (diml==dim) {
                for (int i=0;i<dim;i++) {
                    if (pl>0) RL[i+1] = Rllu[i+2];
                    RU[i+1] = Ruuu[i+2];
                }
            } else {
                for (int i=0;i<dim;i++) {
                    if (pl>0) RL[i+1] = Rllu[ i*diml/dim+2];
                    RU[i+1] = Ruuu[i*dimu/dim+2];
                }
            }
            au = Ruuu[shift];
            if (pl>0) al = Rllu[shift];
        } else {
            a = boost::lexical_cast<double>((1.0 * dim-diml) / (dimu-diml));
            double il;
            double id;
            for (int i=0;i<dim;i++) {
                RL[i+1] = 0;
                RU[i+1] = 0;

                //Decide RL[i+1]
                il = i * (dimu-1.0-0.001) / (dim-1.0);
                id = il - (int)il;
                id = 1.0 - id;
                if (pl>0) RL[i+1] += a * (id * Rluu[2+(int)il]  +  (1.0-id) * Rluu[2+(int)il+1]); 
                RU[i+1] += a * (id * Ruuu[2+(int)il]  +  (1.0-id) * Ruuu[2+(int)il+1]); 

                il = i * (diml-1.0-0.001) / (dim-1.0);   //-0.001 is an implementing technique to avoid an error at i=dim-1
                id = il - (int)il;
                id = 1.0 - id;
                if (pl>0) RL[i+1] += (1.0-a) * (id * Rllu[2+(int)il]  +  (1.0-id) * Rllu[2+(int)il+1]); 
                RU[i+1] += (1.0-a) * (id * Rulu[2+(int)il]  +  (1.0-id) * Rulu[2+(int)il+1]); 
            }
            au = a * Ruuu[shift] + (1.0-a) * Rulu[shift];
            if (pl>0) al = a * Rluu[shift] + (1.0-a) * Rllu[shift];
        }

        //generate final function by binary search
        if (pl>0) {
            for (int i=2;i<=dim;i++) RL[i] = max(RL[i],RL[i-1]);
            for (int i=1;i<=dim;i++) RL[i] = min(RL[i],boost::lexical_cast<DFLOAT>(1.0));
            for (int i=dim/2;i>=1;i--) {
                if (RL[i]<0) RL[i] = RL[i+1];
            }
        }

        for (int i=2;i<=dim;i++) RU[i] = max(RU[i],RU[i-1]);
        for (int i=1;i<=dim;i++) RU[i] = min(RU[i],boost::lexical_cast<DFLOAT>(1.0));
        for (int i=dim/2;i>=1;i--) {
            if (RU[i]<0) RU[i] = RU[i+1];
        }

        //Assume 0<RU,RL<1
        if (pl>0) for (int i=1;i<=dim;i++) RL[i] = log(RL[i]);
        for (int i=1;i<=dim;i++) RU[i] = log(RU[i]);

        double f = 0.5;
        double df = 0.5;
        while (1) {
            if (pl>0) {
                for (int i=1;i<=dim;i++) rd[i] = exp(f*RU[i] + (1.0-f) * RL[i]);
            } else {
                for (int i=1;i<=dim;i++) rd[i] = exp((1.5-f)*RU[i]);
            }
            //modify
            for (int i=2;i<=dim;i++) rd[i] = max(rd[i],rd[i-1]);
            for (int i=1;i<=dim;i++) rd[i] = min(rd[i],boost::lexical_cast<DFLOAT>(1.0));
            for (int i=dim/2;i>=1;i--) {
                if (rd[i]<0) rd[i] = rd[i+1];
            }
            a = f * au + (1.0-f) * al;
            if (pl==0) a = au;
            pp = a * Rigid_upper_prob(rd,1,dim) + (1.0-a) * Rigid_lower_prob(rd,1,dim);

            if (pp>tprob) f-=df;
            if (pp<tprob) f+=df;
            df *= 0.5;
            if (fabs(pp-tprob)/tprob < 0.01) break;
            if (df<1e-6) break;
        }

#ifdef _enable_pfcache
        pfcache[0] = ch;
        pfcache[1] = pp;
        for (int i=1;i<=dim;i++) pfcache[i+1] = rd[i];
#endif
        return pp;
    }

    bkzfloat Rigid_upper_prob(PruningFunction& PF ) {
        std::vector<bkzfloat> r;
        int n = PF.rd.size()-1;
        r.resize(n+1);
        for (int i=1;i<=n;i++) r[i] = PF.rd[i];
        return Rigid_upper_prob(r,1,n);
    }

    bkzfloat Rigid_lower_prob(PruningFunction& PF ) {
        std::vector<bkzfloat> r;
        int n = PF.rd.size()-1;
        r.resize(n+1);
        for (int i=1;i<=n;i++) r[i] = PF.rd[i];
        return Rigid_lower_prob(r,1,n);
    }

    bkzfloat Rigid_upper_volume(PruningFunction& PF ) {
        std::vector<bkzfloat> r;
        int n = PF.rd.size()-1;
        r.resize(n+1);
        for (int i=1;i<=n;i++) r[i] = PF.rd[i];
        return Rigid_upper_volume(r,1,n);
    }

    bkzfloat Rigid_lower_volume(PruningFunction& PF ) {
        std::vector<bkzfloat> r;
        int n = PF.rd.size()-1;
        r.resize(n+1);
        for (int i=1;i<=n;i++) r[i] = PF.rd[i];
        return Rigid_lower_volume(r,1,n);
    }


    bkzfloat Approx_prob(PruningFunction& PF,int parallel=1,int sampling_mult=10 ) {
        std::vector<bkzfloat> r;
        int n = PF.rd.size()-1;
        r.resize(n+1);
        for (int i=1;i<=n;i++) r[i] = PF.rd[i];
        return pruning_func::Approx_prob(r,1,n,(bkzfloat)1.0,parallel,sampling_mult);
    }

    void Rigid_upper_costfactor(PruningFunction& PF) {
        Rigid_upper_costfactor(PF.costfactor,PF.rd,1,PF.rd.size()-1);
    }

    void Rigid_lower_costfactor(PruningFunction& PF) {
        Rigid_lower_costfactor(PF.costfactor,PF.rd,1,PF.rd.size()-1);
    }

    template <typename T> bkzfloat Rigid_upper_cost(PruningFunction& PF,LatticeBasis<T>& B,bkzfloat clim) {
        int n = PF.rd.size()-1;
        B.updateGSBasis();
        std::vector<bkzfloat> lcost;
        bkzfloat cost = Rigid_upper_cost(lcost,PF.rd,B.gs.c,clim,1,n);
        copyvec(PF.costfactor,lcost);
        return cost;
    }

    template <typename T> bkzfloat Rigid_lower_cost(PruningFunction& PF,LatticeBasis<T>& B,bkzfloat clim) {
        int n = PF.rd.size()-1;
        B.updateGSBasis();
        std::vector<bkzfloat> lcost;
        bkzfloat cost = Rigid_lower_cost(lcost,PF.rd,B.gs.c,clim,1,n);
        copyvec(PF.costfactor,lcost);
        return cost;
    }

    template <typename T> bkzfloat Rigid_upper_cost(PruningFunction& PF,std::vector<T>& c,bkzfloat clim) {
        int n = PF.rd.size()-1;
        std::vector<bkzfloat> lcost;
        bkzfloat cost = Rigid_upper_cost(lcost,PF.rd,c,clim,1,n);
        copyvec(PF.costfactor,lcost);
        return cost;
    }

    template <typename T> bkzfloat Rigid_lower_cost(PruningFunction& PF,std::vector<T>& c,bkzfloat clim) {
        int n = PF.rd.size()-1;
        std::vector<bkzfloat> lcost;
        bkzfloat cost =  Rigid_lower_cost(lcost,PF.rd,c,clim,1,n);
        copyvec(PF.costfactor,lcost);
        return cost;
    }

    void fromtable(PruningFunction& PF,bkzfloat prob,int dim,bkzfloat delta) {
        pruning_func::init_pruning_func();
        std::vector<double> R;
        getfuncfromtable(dtable,R,prob,dim,delta);
        PF.rd.resize(1+dim);
        for (int i=1;i<=dim;i++) PF.rd[i] = R[i+1];
    }
    
}

#include "optimizePFCEboost.cpp"

#define pf_nooption 0x00
#define pf_lightoptimize 0x01 
#define pf_CEoptimize 0x02

#define PPFbase 0x01
#define PPFfast 0x02

#include "monotone_opt.cpp"


namespace pruning_func {

    double SetLBPruningFunction(PruningFunction& PF,int istart,int iend,bkzfloat prob,int vl=0) {

        int dim = iend-istart+1;
        PF.rd.resize(iend+1);

        ValueTable<bkzfloat> tt;
        tt.direction = -1;
        bkzfloat range;
        bkzfloat logprob = log(prob);
        bkzfloat alpha;
        bkzfloat pp = prob;
        while (1) {
            if (tt.data.size()==0) {
                alpha = 1.0;
            } else 
            if (tt.data.size()==1) {
                alpha = 0.1;
            } else {
                alpha = tt.nextcandidate(logprob,range);
                if (range < 1e-6) break;
            }
            for (int k=istart;k<iend;k++) {
                int kk = k-istart+1;
                bkzfloat rr = pow(ibeta_inv<bkzfloat>(kk,dim-kk,prob),alpha);
                PF.rd[k] = min(1.0,rr.convert_to<double>());
            }
            PF.rd[iend] = 1.0;
            pp = pruning_func::Rigid_lower_prob(PF.rd,istart,dim);
            tt.add(alpha,log(pp));
        }
        return pp.convert_to<double>();
    }

    double SetPowerPruningFunction(PruningFunction& PF,int istart,int iend,bkzfloat& prob,int iceil,int parallel,int opttarget,int vl=0,char opt=PPFfast) {
        //Setting PF[istart,iend]=i/(dim-ceil)^alpha

        extern double pfupperalpha(int dim,int iceil,int pi,int parallel,int opttarget,int vl);
        ValueTable<bkzfloat> tt;
        tt.direction = -1;
        
        bkzfloat alpha = 0.0;
        bkzfloat astep = 1.0; 
        PF.rd.resize(iend+1);
        if (prob >=1) {
            for (int i=istart;i<=iend;i++) PF.rd[i] = 1.0;
        } 
        if (prob < 0) {
            ccout << "SetPowerPruningFunction: probability error: " << prob << endl;
            exit(0);
        }
        
        int dim = iend-istart+1;
        int j = dim - iceil;
        int phase = 0;

        if (opt==PPFfast) {
            //Calling alpha range 
            bkzfloat logprob = -log(prob) / 2.30258509299;
            logprob = floor(logprob);
            int pi = logprob.convert_to<int>();
            alpha = pfupperalpha(dim,iceil,pi,parallel,opttarget,vl);
            astep = pfupperalpha(dim,iceil,pi+1,parallel,opttarget,vl);
            tt.add(alpha,- 2.30258509299 * pi);
            tt.add(astep,- 2.30258509299 *(pi+1));
            astep = 0.5 * (astep - alpha);
            phase = 1;
            if (abs(prob-exp(- 2.30258509299 * pi))/prob <  1e-6 ) {
                for (int i=istart;i<=iend;i++) {
                    PF.rd[i] = pow((bkzfloat)(1.0*(i-istart+1)/j),alpha);
                    if (PF.rd[i]>=1.0) PF.rd[i] = 1.0;
                }        
                return alpha.convert_to<double>();
            }
            if (abs(prob-exp(- 2.30258509299 * (pi+1)))/prob <  1e-6 ) {
                for (int i=istart;i<=iend;i++) {
                    PF.rd[i] = pow((bkzfloat)(1.0*(i-istart+1)/j),astep);
                    if (PF.rd[i]>=1.0) PF.rd[i] = 1.0;
                }        
                return astep.convert_to<double>();
            }
        } else {
            //Base case
            tt.add(0,1);
        }
        bkzfloat range;
        bkzfloat logprob = log(prob);
        
        while (1) {
            if (tt.data.size()>=2) {
                alpha = tt.nextcandidate(logprob,range);
                if (range < 1e-6) break;
            } else {
                alpha = 1.0;
            }
            for (int i=istart;i<=iend;i++) {
                PF.rd[i] = pow((bkzfloat)(1.0*(i-istart+1)/j),alpha);
                if (PF.rd[i]>=1.0) PF.rd[i] = 1.0;
            }        
            bkzfloat pp = pruning_func::Rigid_lower_prob(PF.rd,istart,dim,parallel);
            if (vl>=VL2) ccout << "alpha=" << alpha << " prob=" << pp << endl;
            if (pp<0) {
                ccout << "Rigid_lower_prob returns error value: pp=" << pp << " alpha=" << alpha << " prob=" << pp << endl;
                ccout << "dim=" << dim << endl;
                PF.display();
                PF.savetofile("bugpf.txt");
                exit(0);  //bug??
            }
            if (abs(pp-prob)/prob < 1e-6) break;
            tt.add(alpha,log(pp));
        }
        return alpha.convert_to<double>();
    }

    template <typename T> bkzfloat SetPruningFunction(std::vector<T>& c, PruningFunction& PF,int istart,int iend,bkzfloat radius,bkzfloat prob, int sharpness=0,int parallel=1,int vl=0,std::string stringoptions="") {
        //input radius is non-squared form
        //returned value is probability

        bkzfloat cost;
        bkzfloat retprob;

        int n = c.size()-1;
        PF.rd.resize(n+1);

        char opttarget = target_prob;
        std::map<std::string,std::string> options;
        if (stringoptions!="") {
            ExtractOptions(options,stringoptions);
            if (options["optimize"] == "volume") {
                opttarget = target_volume;
            }
        }
        bkzfloat fullnumber;
        if (opttarget == target_volume) {
            T gh = lattice_tools::LatticeGH<T>(c,istart,iend);
            fullnumber = boost::lexical_cast<bkzfloat>(pow(radius/gh,iend-istart+1)*0.5);
            prob = prob / fullnumber;   //volume ratio
        }

        if (prob>=1) {
            for (int i=0;i<=n;i++) PF.rd[i] = 1.0;
            return 1.0;
        }
        if (prob <= 0) {
            //probability error?
            ccout << "negative probability: p=" << prob << endl;
            exit(0);
        }

        //Todo: set precision
        int prec = 2 * n;   //temporal

        if (sharpness==0) {
            bkzfloat delta = Detect_delta(c,istart,iend,INPUT_NONSQUARED);
            std::vector<bkzfloat> rdtemp;
            rdtemp.resize(1+n);
            retprob  = SetApproximatePruningFunction(rdtemp,iend-istart+1,delta,prob);
            for (int i=1;i<=iend-istart+1;i++) PF.rd[ istart + i-1] = rdtemp[i];
        }
        
        if (sharpness==pf_lower) {
            pruning_func::SetLBPruningFunction(PF,istart,iend,prob,0);
        }

        if (sharpness>=pf_systematic) {
            bkzfloat prevcost;
            bkzfloat mincost=-1;
            int dim = iend-istart+1;
            COptTable<int,double> CT;
            CT.direction = -1;
            int range;
            int index;
            int loop=0;
            
            while (1) {
                if (CT.data.size()==0) index = 0;
                if (CT.data.size()==1) {
                    index = dim/3;
                    if (dim>=150) index = pow(dim,0.75);
                }
                if (CT.data.size()==2) {
                    index = dim/6;
                    if (dim>=150) index =  pow(dim,0.75)*0.5;
                }
                if (CT.data.size()>=3) {
                    index = CT.nextvalue(range);
                    if (range<=2) break;
                    while (CT.valexist(index)==true) {
                        index = CT.nextvaluerandom(range);
                    }
                }
                bkzconstants::savecachetable();

                double alpha = pruning_func::SetPowerPruningFunction(PF,istart,iend,prob,index,parallel,opttarget,VL0);
                cost = pruning_func::Rigid_upper_cost(PF.rd,c,radius,istart,iend-istart+1,parallel);
                if (mincost < 0) {
                    mincost = cost;
                }
                mincost = min(mincost,cost);
                bkzfloat logcost = log(cost);
                CT.add(index,logcost.convert_to<double>());
                if (loop++>1000) break;
            }
            index = CT.minkey();
            double alpha = pruning_func::SetPowerPruningFunction(PF,istart,iend,prob,index,parallel,opttarget,VL0);
            cost = pruning_func::Rigid_upper_cost(PF.rd,c,radius,istart,iend-istart+1,parallel);
            //ccout << "cost=" << cost << endl;
        }
        
        if (sharpness>=pf_crossentropy) {

            int timemax = iend-istart+1;
            std::vector<T> cshift;
            int n = iend-istart+1;
            cshift.resize(n+1);
            for (int i=istart;i<=iend;i++) {
                cshift[i-istart+1] = c[i];
            }
            std::vector<bkzfloat> rdtemp;
            rdtemp.resize(1+n);
            for (int i=1;i<=iend-istart+1;i++) rdtemp[i] = PF.rd[ istart + i-1];
            bkzfloat gslower = (bkzfloat)c[iend]*(bkzfloat)c[iend]/radius/radius;
            for (int i=1;i<=iend-istart+1;i++) rdtemp[i] = max(PF.rd[ istart + i-1],gslower);
            //for (int i=1;i<=iend-istart+1;i++) ccout << i << " " << rdtemp[i] << endl; 
            int pfinit=0;   //use input pf as the initial
            cost = pruning_func::optimize_pruning_functionCE(rdtemp,radius,cshift,prob,sharpness,pfinit,vl,-1,opttarget);
            for (int i=1;i<=iend-istart+1;i++) PF.rd[ istart + i-1] = rdtemp[i];
    
        }
        return cost;
    }


    template <typename T> bkzfloat SetPruningFunction(LatticeBasis<T>& B, PruningFunction& PF,int istart,int iend,bkzfloat radius,bkzfloat prob, int sharpness=0,int parallel=1,int vl=0,std::string options="") {
        return SetPruningFunction(B.gs.c,PF,istart,iend,radius, prob,sharpness,parallel,vl,options);
    }   
    
    template <typename T> bkzfloat SetPruningFunction(LatticeBasis<T>& B, PruningFunction& PF,bkzfloat radius,bkzfloat prob, int sharpness=0,int parallel=1,int vl=0,std::string options="") {
        return SetPruningFunction(B,PF,1,B.dim,radius,prob,sharpness,parallel,vl,options);
    }

    double pfupperalpha(int dim,int iceil,int pi,int parallel=1,int opttarget=target_prob,int vl=0) {
        //range check
        bkzconstants::loadcachetable();
        
        if ((dim < 0) || (dim>4096)) {
            ccout << "pfupperalpha: dimension error: " << dim << endl;
            exit(0);
        }

        if ((iceil < 0) || (iceil > dim)) {
            ccout << "pfupperalpha: iceil error: " << iceil << endl;
            exit(0);
        }
        
        if ((pi < 0) || (pi > 16384)) {
            ccout << "pfupperalpha: power index error: " << pi << endl;
            exit(0);
        }

        double value;
        #pragma omp critical 
        {
            if (opttarget == target_prob) {
                value = bkzconstants::pfuppertable[dim][iceil][pi];
            } 
            if (opttarget == target_volume) {
                value = bkzconstants::pfuppertablev[dim][iceil][pi];
            } 
        }
        if (value == 0) {
            PruningFunction PF;
            bkzfloat prob = pow((bkzfloat)10,(bkzfloat)-pi);
            if (opttarget == target_prob) {
                bkzconstants::pfuppertable[dim][iceil][pi] = pruning_func::SetPowerPruningFunction(PF,1,dim,prob,iceil,parallel,vl,target_prob,PPFbase);
                value = bkzconstants::pfuppertable[dim][iceil][pi];
            }
            if (opttarget == target_volume) {
                bkzconstants::pfuppertablev[dim][iceil][pi] = pruning_func::SetPowerPruningFunction(PF,1,dim,prob,iceil,parallel,vl,target_volume,PPFbase);
                value = bkzconstants::pfuppertablev[dim][iceil][pi];
            }
        }

        return value;
    }

}

template <typename T> bkzfloat ENUMCostUB(std::vector<T>& c,int istart,int iend,bkzfloat radius,bkzfloat prob,int parallel=1,char input=INPUT_NONSQUARED) {
    PruningFunction PF;
    int dim = iend-istart+1;
    pruning_func::SetPruningFunction(c,PF,istart,iend,radius,prob,pf_systematic,parallel,VL0,"");
    bkzfloat ret = pruning_func::Rigid_upper_cost(PF.rd,c,radius,istart,iend-istart+1,parallel);
    return ret;
}

template <typename T> bkzfloat ENUMCostUB(LatticeBasis<T>& B,bkzfloat radius,bkzfloat prob,int parallel=1) {
    B.updateGSBasis();
    return ENUMCostUB(B.gs.c,1,B.dim,radius,prob,parallel);
}

template <typename T> bkzfloat ENUMCost(std::vector<T>& c,int istart,int iend,bkzfloat radius,bkzfloat prob,char pfoption,int parallel=1,char input=INPUT_NONSQUARED) {
    if (prob>=1.0) {
        return FullENUMCost<bkzfloat>(c,istart,iend,radius,input);
    }

    PruningFunction PF;
    int dim = iend-istart+1;
    pruning_func::SetPruningFunction(c,PF,istart,iend,radius,prob,pfoption,parallel,VL2,"");
    return pruning_func::Rigid_upper_cost(PF.rd,c,radius,istart,iend-istart+1,parallel);
}

template <typename T> bkzfloat ENUMCost(LatticeBasis<T>& B,bkzfloat radius,bkzfloat prob,char pfoption,int parallel=1) {
    B.updateGSBasis();
    return ENUMCostUB(B.gs.c,1,B.dim,radius,prob,pfoption,parallel);
}

#endif