#ifndef _inc_pruning_func_cpp
#define _inc_pruning_func_cpp


//Define the input form of pruning function
//MI=monotonically increasing
//MD=monotonically decreasing
#define orderMI 1
#define orderMD 2

//Target for optimizing 
#define target_prob 0x01
#define target_volume 0x02


//#define pfdebug

#include "pruningfuncsupport.cpp"
#include "vectorenumeration.cpp"

extern double Estimate_delta(quad_float* c,int bs);


namespace progressive_bkz {
    double total_pruningfunc_time = 0;
}



namespace pruning_func {
    
    int initpf = 0;
    std::vector<std::vector<double> > dtable;     //table of pruning function 

    //used for caching pruning functions
    int buffnum=160000;
    int** buffa;
    double** buffb;
    double** small_modify_a;
    int** buffva;   //for volume factor cache
    quad_float** buffvb;
    double** gnr_share_double;
    
    //used for cache
    std::string pfsavedir;
    
    void init_pruning_func() {
        
        if (initpf>0) return;
        
        if (FileExists("pfdata2table.dat")==true) {
            loaddoublevector(dtable,"pfdata2table.dat");
            //ccout << "read: " << dtable.size() << endl;
        } else {
            ccout << "pfdata2table.dat not found!" << endl;
            exit(0);
        }

        int i,j;
        buffa = new int*[buffnum];
        buffb = new double*[buffnum];
        for (i=0;i<buffnum;i++) {
            buffa[i] = new int[5];
            buffa[i][0] = -1;   //not defined
        }
        for (i=0;i<buffnum;i++) {
            buffb[i] = new double[601];
        }

        buffva = new int*[buffnum];
        buffvb = new quad_float*[buffnum];
        for (i=0;i<buffnum;i++) {
            buffva[i] = new int[5];
            buffva[i][0] = -1;   //not defined
        }
        for (i=0;i<buffnum;i++) {
            buffvb[i] = new quad_float[601];
        }
 
        gnr_share_double = new double*[5];
        for (i=0;i<5;i++) {
            gnr_share_double[i] = new double[latticemaxdim];
        }
        
        small_modify_a = new double*[30];
        for (i=0;i<30;i++) {
            small_modify_a[i] = new double[500];
            for (j=0;j<500;j++) small_modify_a[i][j] = -1;
        }


        pfsavedir =  makefullpath(ReadConf("bkz.conf","pfcache"));
        mkdirRecursive(pfsavedir.c_str(), 0777);
        initpf = 1;
    }
   
    int curvehash(double p0,int cdim,int cdelta) {
        int chash;
        chash = -log(p0)*1000.0;
        chash *= 612;
        chash += cdim * 17311 + cdelta * 45987;
        chash %= buffnum;
        return chash;
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
    
    bool pfreadfromcache(double* R,double p0,int dim,double delta,int opt) {
        int chash = curvehash(p0,dim,delta);
        if (buffa[chash][0]!=-1) {
            //the curve is in cache
            if ((buffa[chash][1]==(int)(-log(p0)*1000.0)) && (buffa[chash][2]==dim) && (buffa[chash][3]==delta)) {
                for (int i=1;i<=dim;i++) R[i] = buffb[chash][i];
                return true;
            } else {
            }
        }
        std::ostringstream fname;
        mkdir(pfsavedir.c_str(),0777);
        fname << pfsavedir << "/app" << min(p0,1.0) << "_" << delta << "_" << dim << "_" << opt;
        //ccout << "pfcachename=" << fname.str() << endl;
        if (FileExists(fname)==true) {
            loadarray(R,dim,fname.str());
            return true;
        }
        return false;
    }

    bool pfsavetocache(double* R,double p0,int dim,double delta,int opt) {
        //save to cache
        int chash = curvehash(p0,dim,delta);
        buffa[chash][0] = chash;
        buffa[chash][1] = (int)(-log(p0)*1000.0);
        buffa[chash][2] = dim;
        buffa[chash][3] = delta;
        int i;
        for (i=1;i<=dim;i++) buffb[chash][i] = R[i]; 


        std::ostringstream fname;
        fname << pfsavedir << "/app" << min(p0,1.0) << "_" << delta << "_" << dim << "_" << opt;
        savearray(R,dim,fname.str());
        return true;
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
        //ccout << "i=" << i << " " << prob << " " << dim << " " << delta << endl;
        //ccout << "pfsize=" << pf[i].size() << endl;
        //ccout << pf.size() << endl;

        R.resize(pf[i].size());
        for (j=0;j<(int)pf[i].size();j++) R[j] = pf[i][j];

    }
#ifdef __omit    
    RR set_pruning_function_anyprob_ver2(double* R,double p0,double delta,int dim,int opt=0) {

        double startcputime = clock();
        RR pp;      //final probability
        
        //read_from_cache
#ifndef pfdebug
 /*
        if (pfreadfromcache(R,p0,dim,delta,opt)==true) {
            return R[0];
        }
  */ 
#endif
        //generating function
        if (dim<20) {
            //for small dimensions
            if (dim<=10) {
                int i;
                for (i=1;i<=dim;i++) R[i] = 1.0;
                return to_RR(1.0);
            } else {
                RR pp;
                pp = set_pruning_function_anyprob_ver2(buffb[0],p0,delta,20,opt);
                int i;
                for (i=1;i<=dim;i++) {
                    R[i] = buffb[0][i+(20-dim)];
                }
                return pp;
            }            
        }

        if (p0>=1.0) {
            int i;
            for (i=1;i<=dim;i++) R[i] = 1.0;
            return to_RR(1.0);
        }
        
        double pu,pl;
        double du = 0,dl = 0;
        double a;
        int dimu,diml;
        int i;
        
        getprobUL(p0,pu,pl);
        getdeltaUL(delta,du,dl);
        getdimUL(dim,dimu,diml);
/*        
        ccout << "pu=" << pu << " pl=" << pl << endl;
        ccout << "du=" << du << " dl=" << dl << endl;
        ccout << "dimu=" << dimu << " diml=" << diml << endl;
  */      
        std::vector<double> Ruuu,Ruul,Rulu,Rull;        //{u,l}^3 is corresponding to prob,dim,delta resp.
        std::vector<double> Rluu,Rlul,Rllu,Rlll;        //{u,l}^3 is corresponding to prob,dim,delta resp.
        double d0 = delta;      //for readability
        
        
        //marge delta

        //for upper prob
        getfuncfromtable(dtable,Ruuu,pu,dimu,du);
        getfuncfromtable(dtable,Ruul,pu,dimu,dl);
        a = (d0-dl) / (du-dl);
        for (i=0;i<(int)Ruuu.size();i++) Ruuu[i] = Ruul[i] * (1.0-a) + Ruuu[i] * a;   
            if (diml==dimu) {
            Rulu.resize(Ruuu.size());    
            for (i=0;i<(int)Ruuu.size();i++) Rulu[i] = Ruuu[i];   
        } else {
            getfuncfromtable(dtable,Rulu,pu,diml,du);
            getfuncfromtable(dtable,Rull,pu,diml,dl);
            a = (d0-dl) / (du-dl);
            //ccout << "a=" << a << endl;
            for (i=0;i<(int)Rulu.size();i++) Rulu[i] = Rull[i] * (1.0-a) + Rulu[i] * a;   
            //for (i=0;i<Rulu.size();i++) ccout << "Rulu[" << i << "]=" << Rulu[i] << endl;
        }
        
        if (pl>0) {
            //the same operation to the lower prob
            getfuncfromtable(dtable,Rluu,pl,dimu,du);
            getfuncfromtable(dtable,Rlul,pl,dimu,dl);
            a = (d0-dl) / (du-dl);
            for (i=0;i<(int)Rluu.size();i++) Rluu[i] = Rlul[i] * (1.0-a) + Rluu[i] * a;   

            if (diml==dimu) {
                Rllu.resize(Rluu.size());    
                for (i=0;i<(int)Rluu.size();i++) Rllu[i] = Rluu[i];   
            } else {
                getfuncfromtable(dtable,Rllu,pl,diml,du);
                getfuncfromtable(dtable,Rlll,pl,diml,dl);
                a = (d0-dl) / (du-dl);
                for (i=0;i<(int)Rllu.size();i++) Rllu[i] = Rlll[i] * (1.0-a) + Rllu[i] * a;   
            }
        }

        //Here Ruuu and Rulu are interpolated function for upper prob
        //     Rluu and Rllu are interpolated function for lower prob
        //Next, interpolate the dimension factor
        //for upper prob
        double* RU = gnr_share_double[0]; 
        double* RL = gnr_share_double[1]; 
        double au,al = 0;
        int shift=0;    //=1 if dim is odd
        if (dim%2==1) shift = 1;
        
        if (diml==dimu) {
            for (i=0;i<dim;i++) {
                if (pl>0) RL[i+1] = Rllu[i+2];
                RU[i+1] = Ruuu[i+2];
            }
            au = Ruuu[shift];
            if (pl>0) al = Rllu[shift];
        } else {
            a = (1.0 * dim-diml) / (dimu-diml);
            double il;
            double id;
            for (i=0;i<dim;i++) {
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
#ifdef pfdebug
                //for debug
                ccout << "RU: " << RU[i+1] << endl;
                ccout << "i=" << i << endl;
                ccout << id << " "  << Rulu[2+(int)il] << endl; 
                ccout << id << " "  << Ruuu[2+(int)il] << endl;
#endif
            }
            au = a * Ruuu[shift] + (1.0-a) * Rulu[shift];
            if (pl>0) al = a * Rluu[shift] + (1.0-a) * Rllu[shift];
        }

#ifdef pfdebug
        //for debug
        for (i=1;i<=dim;i++) {
            ccout << "RL[" << i << "]=" << RL[i] << " ";
            ccout << "RU[" << i << "]=" << RU[i] << " " << endl;
        }
       
        if (pl>0) ccout << "probL(RL) = " << Rigid_lower_prob(RL,dim) << endl; 
        if (pl>0) ccout << "probU(RL) = " << Rigid_upper_prob(RL,dim) << endl; 
        ccout << "probL(RU) = " << Rigid_lower_prob(RU,dim) << endl; 
        ccout << "probU(RU) = " << Rigid_upper_prob(RU,dim) << endl; 
        if (pl>0) ccout << "probex(RL) =" << Approx_prob(RL,dim) << endl;
        ccout << "probex(RU) =" << Approx_prob(RU,dim) << endl;

        ccout << "au=" << au << " al=" << al << endl;
        if (pl>0) ccout << "interpolate_prob(RL)=" << al * Rigid_upper_prob(RL,dim) + (1.0-al) * Rigid_lower_prob(RL,dim) << endl;
        ccout << "interpolate_prob(RU)=" << au * Rigid_upper_prob(RU,dim) + (1.0-au) * Rigid_lower_prob(RU,dim) << endl;
#endif
        
        if (pl>0) {
            for (i=2;i<=dim;i++) RL[i] = max(RL[i],RL[i-1]);
            for (i=1;i<=dim;i++) RL[i] = min(RL[i],1.0);
            for (i=dim/2;i>=1;i--) {
                if (RL[i]<0) RL[i] = RL[i+1];
            }
        }

        for (i=2;i<=dim;i++) RU[i] = max(RU[i],RU[i-1]);
        for (i=1;i<=dim;i++) RU[i] = min(RU[i],1.0);
        for (i=dim/2;i>=1;i--) {
            if (RU[i]<0) RU[i] = RU[i+1];
        }
#ifdef pfdebug        
        for (i=1;i<=dim;i++) ccout << "RL[" << i << "]=" << RL[i] << endl;
        for (i=1;i<=dim;i++) ccout << "RU[" << i << "]=" << RU[i] << endl;
#endif
        //Assume 0<RU,RL<1
        if (pl>0) for (i=1;i<=dim;i++) RL[i] = log(RL[i]);
        for (i=1;i<=dim;i++) RU[i] = log(RU[i]);
        
#ifdef pfdebug        
        for (i=1;i<=dim;i++) ccout << "RL[" << i << "]=" << RL[i] << endl;
        for (i=1;i<=dim;i++) ccout << "RU[" << i << "]=" << RU[i] << endl;
#endif

        
        double f = 0.5;
        double df = 0.5;
        while (1) {
            if (pl>0) {
                for (i=1;i<=dim;i++) R[i] = exp(f*RU[i] + (1.0-f) * RL[i]);
            } else {
                for (i=1;i<=dim;i++) R[i] = exp((1.5-f)*RU[i]);
            }
            //modify
            for (i=2;i<=dim;i++) R[i] = max(R[i],R[i-1]);
            for (i=1;i<=dim;i++) R[i] = min(R[i],1.0);
            for (i=dim/2;i>=1;i--) {
                if (R[i]<0) R[i] = R[i+1];
            }
            
            //for (i=1;i<=dim;i++) ccout << "R[" << i << "]=" << R[i] << endl;
            a = f * au + (1.0-f) * al;
            if (pl==0) a = au;
            //pp = a * Rigid_upper_prob(R,dim) + (1.0-a) * Rigid_lower_prob(R,dim);
#ifdef pfdebug
            ccout << "UB=" <<  Rigid_upper_prob(R,dim) << endl;
            ccout << "LB=" <<  Rigid_lower_prob(R,dim) << endl;
            ccout << "f=" << f << " pred_pp=" << pp << endl;
#endif
            if (pp>p0) f-=df;
            if (pp<p0) f+=df;
            df *= 0.5;
            if (fabs(pp-p0)/p0 < 0.01) break;
            if (df<1e-6) break;
        }
        progressive_bkz::total_pruningfunc_time += (clock()-startcputime) / CLOCKS_PER_SEC; 

#ifndef pfdebug
        pfsavetocache(R,p0,dim,delta,opt);
#endif
        return pp;
    }
       
    
    RR set_pruning_function_anyprob(double* R,double p0,double delta,int dim,int opt=0,int version=1) {
        init_pruning_func();
        return set_pruning_function_anyprob_ver2(R,p0,delta,dim,opt);
    }
    
    void computeVfromR(quad_float* V,double* R,int dim) {
        RR t;
        int i,k;
        double* R2 = gnr_share_double[3];
        for (k=2;k<dim;k+=2) {
          if (R[k]!=0) {
            for (i=1;i<=k/2;i++) R2[i] = R[2*i]/R[k];
            t = EvenSimplex::EvenSimplexVolume(k/2,R2,opt_volume_prob);
            t *= VolumeBall(k,sqrt(R[k]));
            V[k/2]=to_quad_float(t);
          }
        }
    }

    void set_volume_factors(quad_float* V,double p0,double delta,int dim,char opt=0) {

        init_pruning_func();

        int cdim,cdelta;
        cdim = dim - 60;
        if (cdim<0) cdim = 0;
        cdelta = (int)((delta-1.005)*10000.0+0.5);
        if (cdelta<0) cdelta = 0;
        if (cdelta>=150) cdelta = 150;
        int chash = curvehash(min(p0,1.0),dim,cdelta);

        int i;
        if (buffva[chash][0]!=-1) {
            if ((buffva[chash][1]==(int)(-log(min(p0,1.0)*1000.0))) && (buffva[chash][2]==dim) && (buffva[chash][3]==cdelta)) {
                for (i=1;i<=dim/2+1;i++) V[i] = buffvb[chash][i];
                return;
            } else {
            }
        }

        double* R = gnr_share_double[2];
        set_pruning_function_anyprob(R,p0,delta,dim,opt,1); 

        computeVfromR(V,R,dim);


        buffva[chash][0] = chash;
        buffva[chash][1] = (int)(-log(min(p0,1.0)*1000.0));
        buffva[chash][2] = dim;
        buffva[chash][3] = cdelta;
        for (i=1;i<=dim/2;i++) buffvb[chash][i] = V[i]; 
    }
    
    RR optimize_pruning_function(double* pf,double clim,quad_float* c,int n,double tprob,int vl,int timelimit=-1,char order=orderMD,char opttarget=target_prob) {
        //pf[1,...bs] monotonically increasing if order==orderMI
        //pf[1...bs] monotonically decreasing if order==orderMD
        //c[1,...,bs] squared value
        //pf[bs] and pf[bs-1] must be 1 to compute upper/lower bounds
        //clim is the pruning radius (non-squared)
        
        double startopttime = clock();
        sampling_tools::initialize();
        
        int i;
        RR problb,probub,probex,a;
        RR cost,mincost;
        double persec =  lattice_enum::enum_speed_bench(omp_get_max_threads()) * 1000000.0;
        if (persec<0) persec = 3e+7;
        unsigned int* th = sampling_tools::rstates; //for generating random numbers        
        
        
        //reverse
        if (order==orderMD) for (i=0;i<n/2;i++) swap(pf[i+1],pf[n-i]);
        // Assume here that pf[1...n] is monotonically increasing
        
       if ((opttarget==target_prob) && (tprob>=1.0)) {
           RR ret = 2 * Rigid_upper_cost(pf,c,n,clim); 
           if (order==orderMD) for (i=0;i<n/2;i++) swap(pf[i+1],pf[n-i]);
           return ret;
       }
        
        //pf[1] = max(pf[2],max(pf[1],min(0.01,tprob)));
        for (i=2;i<=n;i++) pf[i] = max(pf[i],pf[i-1]);  //monotonically increasing
        for (i=1;i<=n;i++) pf[i] = min(1.0,max(0.0,pf[i])); //between zero and one
        pf[n-1] = pf[n] = 1.0;
        for (i=n;i>=1;i--) {
            if (pf[i]==0) pf[i] = pf[i+1];
        }

        
        int docount = 0;
        do {

            if (opttarget==target_prob) probex = to_double(Approx_prob(pf,n));
            if (opttarget==target_volume) {
                probex = to_double(Approx_final_count(pf,n,c,clim));
            }
            
            ccout << "probex=" << probex << endl;
            ccout << "tprob=" << tprob << endl;
            
            if (probex < tprob * 0.99) {
                for (i=1;i<=n;i++) pf[i] = min(pf[i]*1.01,1.0);
                for (i=n/2;i>=1;i--) {
                    if (pf[i]<=0) pf[i] = pf[i+1];
                }
                for (i=2;i<=n;i++) pf[i] = max(pf[i],pf[i-1]);  //monotonically increasing
            }

            if (docount++>30) {
                ccout << "infinite loop?";
                return to_RR(0);
            }
        } while (probex < tprob * 0.99);

        if (opttarget==target_prob) {         
            problb = Rigid_lower_prob(pf,n);
            probub = Rigid_upper_prob(pf,n);
            if (vl>=3) {
                ccout << "dim=" << n ;
                ccout << " init pruning cost: target_prob=" << tprob; 
                ccout << " rigid_lower_prob=" << problb;
                ccout << " rigid_upper_prob=" << probub << "\r";
                ccout.flush();
            }        
            //probex = to_double(Approx_prob(pf,n));
            if (probub > problb * 1.005) {
                a = (probex - problb) / (probub - problb);
            } else {
                //both values are very close to each other
                a = 0.5;
            }
        }
        cost =  mincost = 2 * Rigid_upper_cost(pf,c,n,clim);
        if (opttarget==target_prob) {         
            if (vl>=3) {
                ccout << "dim=" << n ;
                ccout << " init pruning cost: target_prob=" << tprob; 
                ccout << " rigid_lower_prob=" << problb;
                ccout << " rigid_upper_prob=" << probub;
                ccout << " Exact_prob=" << probex << endl;
                ccout << "exact ~ " << a << "*ub+" << 1-a << "*lb" << endl;
                ccout << "initial cost=" << cost << " (" << cost / persec << "sec)" << endl;
    #ifdef pfdebug
                if ((a<0) || (a>1)) {
                    ccout << "pf error?: " << endl;
                    for (i=0;i<=n;i++) ccout << pf[i] << ",";
                    ccout << endl;
                    for (i=1;i<=10;i++) ccout << "Re-exact=" << to_double(Approx_prob(pf,n)) << endl;
                    for (i=1;i<=10;i++) ccout << "Re-direct2=" << to_double(sampling_tools::PruningProbDirect(pf,n)) << endl;

                    exit(0);
                }

    #endif            
            }
        }
        
        double* bpf = gnr_share_double[4];  //backup of pruning function
        for (i=1;i<=n;i++) bpf[i] = pf[i];
        int numloop = 0;
        
        
        double startcputime = clock();
        double expecttotal; //Elapsed time for optimizing pruning function + Expected time for lattice enumeration
        double minexpect;
        expecttotal = minexpect = to_double(cost / persec);
        double tick = 0.07 - 0.001 * n;
        do {
            double cputime = (clock() - startcputime) / CLOCKS_PER_SEC;
            
            tick = sampling_tools::unique(-2.0/n,2.0/n,th);
            for (i=1;i<=n;i++) pf[i] = bpf[i];
            
            //perturbation
            int center;
            int band = min(5,n/20);

            do {
                center = 1 + sampling_tools::normal_dist(n*0.5,n*1.4,th);
            } while ((center<-band) || (center>n+band));
  
            double ss;
            ss = sampling_tools::unique(100.0/n/n,0.3,th);
            for (i=center-band;i<=center+band;i++) {
                if ((1<=i) && (i<=n)) {
                    pf[i] += tick * exp(-(i-center)*(i-center)*ss);
                }
            }
            
            if (opttarget==target_prob) {
                pf[1] = max(pf[2],max(pf[1],min(0.01,tprob)));
            } else {
                pf[1] = max(pf[2],max(pf[1],1e-12));
            }
            for (i=2;i<=n;i++) pf[i] = max(pf[i],pf[i-1]);  //monotonically increasing
            for (i=1;i<=n;i++) pf[i] = min(1.0,max(0.0,pf[i])); //between zero and one
            if (opttarget==target_volume) for (i=1;i<=n;i+=2) pf[i] = pf[i+1];
            pf[n-1] = pf[n] = 1.0;


            
            if (opttarget==target_prob) probex = a * Rigid_upper_prob(pf,n) + (1-a) * Rigid_lower_prob(pf,n);
            if (opttarget==target_volume) {
                probex = to_double(Approx_final_count(pf,n,c,clim));
            }


            int ilp=0;
            while ((probex < tprob) & (ilp++<20)) {
                int center = 1+rand()%n;
                double diff = sampling_tools::unique(-1.0/n,1.0/n,th);;
                for (i=center-band;i<=center+band;i++) {
                    if ((1<=i) && (i<=n)) {
                        pf[i] += diff * exp(-(i-center)*(i-center)*0.1);
                    }
                }
                if (opttarget==target_prob) {
                    pf[1] = max(pf[2],max(pf[1],min(0.01,tprob)));
                } else {
                    pf[1] = max(pf[2],max(pf[1],1e-12));
                }
                for (i=2;i<=n;i++) pf[i] = max(pf[i],pf[i-1]);  //monotonically increasing
                for (i=1;i<=n;i++) pf[i] = min(1.0,max(0.0,pf[i])); //between zero and one
                pf[n-1] = pf[n] = 1.0;
                if (opttarget==target_volume) for (i=1;i<=n;i+=2) pf[i] = pf[i+1];

                
                
                if (opttarget==target_prob) probex = a * Rigid_upper_prob(pf,n) + (1-a) * Rigid_lower_prob(pf,n);
                if (opttarget==target_volume) {
                    probex = to_double(Approx_final_count(pf,n,c,clim));
                }
                if (probex>tprob) {
                    cost =  0.5*(Rigid_upper_cost(pf,c,n,clim) + Rigid_lower_cost(pf,c,n,clim));
                    if (cost < mincost * 0.9999) break;
                }
/*                
                diff =  0.0001 * (tprob - probex) / (probex - probex2);
                for (i=c-band;i<=c+band;i++) {
                    if ((1<=i) && (i<=n)) {
                        pf[i] += diff * exp(-(i-c)*(i-c)*0.1);
                    }
                }
                for (i=2;i<=n;i++) pf[i] = max(pf[i],pf[i-1]);  //monotonically increasing
                for (i=1;i<=n;i++) pf[i] = min(1.0,max(0.0,pf[i])); //between zero and one
                pf[n-1] = pf[n] = 1.0;
                 
                probex = a * Rigid_upper_prob(pf,n) + (1-a) * Rigid_lower_prob(pf,n);
                //ccout << probex << endl;
 
 */
            };
        cost =  0.5*(Rigid_upper_cost(pf,c,n,clim) + Rigid_lower_cost(pf,c,n,clim));
            if (vl>=3) {
                ccout << "loop=" << numloop << " sec=" << cost/persec << " (minsec=" << mincost / persec << ") Rtime=" << timelimit - (clock() - startcputime)/CLOCKS_PER_SEC << " prob=" << probex << "\r";
                //ccout << "loop=" << numloop << " sec=" << cost/persec << " (minsec=" << mincost / persec << ") Elapsed time=" << cputime << " prob=" << probex << "\r";
                ccout.flush();
            }
            if ((cost < mincost * 0.9999) && (probex > tprob)) {
                mincost = cost;
                for (i=1;i<=n;i++) bpf[i] = pf[i];
                expecttotal =  to_double(mincost / persec) + cputime;
                if (minexpect<0) minexpect = expecttotal;
                if (expecttotal < minexpect) minexpect = expecttotal;
                //startcputime = clock();
            }
            numloop++;
            
            if (timelimit==-1) {
                if (minexpect *1.05 <  mincost / persec + cputime) break;
            } else {
                if (cputime > timelimit) break;
            } 

#ifdef pfdebug
            if (numloop%10==0) {
#else 
            if ((numloop%100==0) && (opttarget==target_prob)) {
#endif
                double startcputime2 = clock();
                problb = Rigid_lower_prob(bpf,n);
                probub = Rigid_upper_prob(bpf,n);
                probex = to_double(Approx_prob(bpf,n));
                if (vl>=3) {
                    ccout << "probex=" << probex << endl;
                    ccout << "update a:" << a << "->";
                }
                if ((probub > problb * 1.005)  || (cost/persec>10000)) {
                    a = (probex - problb) / (probub - problb);

                    if (a<0) {
                        ccout << endl;
                        ccout << endl;
                        ccout << "ex=" << probex << endl;
                        ccout << "ub=" << probub << endl;
                        ccout << "lb=" << problb << endl;
                        ccout << "a=" << a << endl;
                        ccout << "dim=" << n << endl;
                                
                        ccout << "bpf=" << endl;
                        for (i=1;i<=n;i++) ccout << bpf[i] << " ";
                        ccout << endl;
                        ccout << "pf=" << endl;
                        for (i=1;i<=n;i++) ccout << pf[i] << " ";
                        ccout << endl;
                        break;
                        //exit(0);
                        //return to_RR(0);
                    }
  
                } else {
                    a =0.5;
                }
                        
                //if (vl>=3) ccout << a << endl;
                startcputime += clock() - startcputime2;        //do not count up cputime of this section
            }
        } while (1);
        
#ifdef pfdebug
        //final stat
        problb = Rigid_lower_prob(pf,n);
        probub = Rigid_upper_prob(pf,n);
        probex = to_double(Approx_prob(pf,n));
        ccout << endl;
        ccout << "probex=" << probex << " (target=" << tprob << ")" << endl;
        ccout << "final_a:" << a << "->";
        a = (probex - problb) / (probub - problb);
        ccout << a << endl;
#endif 
        
        for (i=1;i<=n;i++) pf[i] = bpf[i];
        //for (i=1;i<=n;i+=2) pf[i] = pf[i+1];
        mincost = Rigid_upper_cost(pf,c,n,clim);

        ccout << "pf=" << endl;
        for (i=1;i<=n;i++) ccout << pf[i] << " ";
        ccout << endl;

        //reverse
        if (order==orderMD) for (i=0;i<n/2;i++) swap(pf[i+1],pf[n-i]);
        //reverse
        //for (i=0;i<n/2;i++) swap(pf[i+1],pf[n-i]);
        if (vl>=3) ccout << endl;
        
        progressive_bkz::total_pruningfunc_time += (clock()-startopttime) / CLOCKS_PER_SEC; 
/*        
        ccout << "last_volume=" << to_double(Approx_final_count(pf,n,c,clim)) << endl;
        ccout << "upper_cost=" << Rigid_upper_cost(pf,c,n,clim) << endl;
        ccout << "lower_cost=" << Rigid_lower_cost(pf,c,n,clim) << endl;
        ccout << "rigid_lower_prob=" << Rigid_lower_prob(pf,n) << endl;
        ccout << "rigid_upper_prob=" << Rigid_upper_prob(pf,n) << endl;
  */      
        return mincost;
        
    }

    void invert(double* pf,int a,int b) {
        int i;
        for (i=0;i<(b-a)/2;i++) {
            swap(pf[a+i],pf[b-i]);
        }
    }
#endif
    
}

void set_approximate_extreme_pruning(double* pf,quad_float* c,int bs,double prob,double delta) {
 }


#endif
