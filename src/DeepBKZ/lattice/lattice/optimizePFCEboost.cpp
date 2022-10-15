#ifndef _optimize_pfce_boost_cpp
#define _optimize_pfce_boost_cpp


#include "bkzconstants.hpp"
#include "pfuncwrapper.cpp"


template <typename DFLOAT> void samplingpf(std::vector<DFLOAT>& ret,std::vector<DFLOAT>& center,std::vector<DFLOAT>& sigma,std::vector<DFLOAT>& lowerbound,DFLOAT bound) {
    //fill the data to ret[1..n]
    boost::random::mt19937 rng;  
    using boost::math::normal;
    normal gaussian;
    
    int n = ret.size();
    ret[0] = 0;
    for (int i=1;i<n;i++) {
        DFLOAT s = sigma[i];
        
        DFLOAT c = center[i];
        do {
              ret[i] = sampling_tools::normal<DFLOAT>(c,s);
        } while (ret[i]<=0);
    }

    std::sort(ret.begin(),ret.end());
    for (int i=1;i<n;i++) {
        ret[i] = max(ret[i],lowerbound[i]);
    }
    for (int i=n-2;i>=1;i--) {
        if (ret[i] < ret[i+1]*0.1) ret[i] = ret[i+1]*0.1;     //To prevent big drop
    }
    if (bound>0) {
        for (int i=1;i<n;i++) ret[i] = min((bkzfloat)bound,ret[i]); 
        ret[n-1] = max((bkzfloat)bound,ret[n-1]);
    } else {

    }
}

template <typename DFLOAT> struct CoeffData {
    std::vector<DFLOAT> r; 
    DFLOAT cost;
    DFLOAT prob;
    void init(DFLOAT tcost,DFLOAT  tprob) {
        cost = tcost;
        prob = tprob;
    }
};


template <typename DFLOAT> class CoeffList {

public:

    std::vector<CoeffData<DFLOAT> >  X;
    CoeffList() {

    }

    void addelement(std::vector<DFLOAT>& pftemp,DFLOAT tprob,DFLOAT tcost) {
        CoeffData<DFLOAT> D;
        copyvec(D.r,pftemp);
        D.cost = tcost;
        D.prob = tprob;

        X.push_back(D);
    }

    bool replacemaxelement(std::vector<DFLOAT>& pftemp,DFLOAT tprob,DFLOAT tcost) {

        int maxi = 0;
        DFLOAT mcost = X[0].cost;
        for (int i=1;i<X.size();i++) {
            if (mcost < X[i].cost) {
                mcost = X[i].cost;
                maxi = i;
            }
        }
        if (mcost > tcost) {
            copyvec(X[maxi].r,pftemp);
            X[maxi].cost = tcost;
            X[maxi].prob = tprob;
            return true;
        } else {
            return false;
        }
        return false;
    }

          std::vector<std::vector<bkzfloat> > YT;
          void predict(std::vector<DFLOAT>& c,std::vector<DFLOAT>& v,DFLOAT beta,DFLOAT bound) {
            //c=center
            //v=variance
            int n = X[0].r.size();
            c.resize(n);
            v.resize(n);

            YT.resize(X.size());
            for (int i=0;i<X.size();i++) {
                YT[i].resize(2);
                YT[i][0] = X[i].cost;
                YT[i][1] = i;
            }
            std::sort(YT.begin(),YT.end());

            //weight of point Y
            std::vector<bkzfloat> wt;
            wt.resize(X.size());
            for (int i=1;i<X.size()-1;i++) {
                wt[i] = 0.5*(YT[i+1][0]+YT[i][0]) - 0.5*(YT[i][0]+YT[i-1][0]);
            }
            wt[0] = YT[1][0]-YT[0][0];
            wt[X.size()-1] = YT[X.size()-1][0]-YT[X.size()-2][0];

            //minY
            bkzfloat minY = YT[0][0];

            for (int i=1;i<n;i++) {

                bkzfloat A,B;
                bkzfloat h1,h2;

                c[i] = 0;
                for (int j=0;j<X.size();j++) c[i] += X[j].r[i];
                c[i] /= X.size();

                h1 =  h2 = 0;
                int c1=0;
                int c2=0;
                for (int j=0;j<X.size();j++) {
                    if (j < X.size()/2) {
                        h1 += X[j].r[i];
                        c1++;
                    } else {
                        h2 += X[j].r[i];
                        c2++;
                    }
                } 
                h1 /= c1;
                h2 /= c2;


                v[i] = 0;
                for (int j=0;j<X.size();j++) {
                    v[i] += (X[j].r[i] - c[i]) * (X[j].r[i] - c[i]) ;
                }
                v[i] /= X.size();
                v[i] = 2.0 * sqrt(v[i]);

                //re-compute center
                int yi = YT[0][1];
                c[i] *= exp(beta * (log(X[yi].r[i]) - log(c[i])));
                if (bound > 0) {
                    c[i] = max((bkzfloat)1e-120,min((bkzfloat)bound,c[i]));
                } else {
                    c[i] = max((bkzfloat)1e-120,c[i]);
                }
                v[i] = 0.3 * c[i];
            }        
        }

    void getcv(std::vector<DFLOAT>& c,std::vector<DFLOAT>& v) {
        //c=center
        //v=variance
        int n = X[0].r.size();
        c.resize(n);
        v.resize(n);
        for (int i=1;i<n;i++) {
            c[i] = 0;
            for (int j=0;j<X.size();j++) {
                c[i] += X[j].r[i];
            }
            c[i] /= X.size();

            v[i] = 0;
            for (int j=0;j<X.size();j++) {
                v[i] += (X[j].r[i] - c[i]) * (X[j].r[i] - c[i]) ;
            }
            v[i] /= X.size();
            v[i] = 2.0 * sqrt(v[i]);
            if (v[i]==0) v[i] = 0.1 * c[i];
        }        
    }
    bkzfloat getvariance() {
        if (X.size()==0) return 0;
        int n = X[0].r.size();
        std::vector<bkzfloat> c;
        std::vector<bkzfloat> v;
        c.resize(n);
        v.resize(n);
        bkzfloat vret;
        for (int i=1;i<n;i++) {
            c[i] = 0;
            for (int j=0;j<X.size();j++) {
                c[i] += X[j].r[i];
            }
            c[i] /= X.size();

            v[i] = 0;
            for (int j=0;j<X.size();j++) {
                v[i] += (X[j].r[i] - c[i]) * (X[j].r[i] - c[i]) ;
            }
            v[i] /= X.size();
            vret += 2.0 * sqrt(v[i]);
        }        
        return vret / n;
    }
};

#include "monotone_opt.cpp"

namespace pruning_func {

        void adjust_pruning_function(std::vector<bkzfloat>& pf) {
            int n = pf.size()-1;
            for (int i=2;i<=n;i++) pf[i] = max(pf[i],pf[i-1]);  //monotonically increasing
            for (int i=1;i<=n;i++) pf[i] = min((bkzfloat)1.0,max((bkzfloat)0.0,pf[i])); //between zero and one
            pf[n-1] = pf[n] = 1.0;
            for (int i=n;i>=1;i--) {
                if (pf[i]==0) pf[i] = pf[i+1];
            }
        }

    bkzfloat Approx_final_volume(std::vector<bkzfloat>& pf,bkzfloat clim,int istart,int iend) {
        int dim = iend-istart+1;
        static std::vector<bkzfloat> R2;
        R2.resize(1+dim/2);
        int shift = 0;
        if (dim%2==1) dim--;
        for (int i=1;i<=dim/2;i++) {
            R2[i] = pf[(istart-1)+2*i-shift]/pf[(istart-1)+dim-shift];
        }
        bkzfloat ret = EvenSimplexVolumeWrapper<bkzfloat>(R2,dim/2,opt_volume_prob,1.0) ;
        shift = 1;
        for (int i=1;i<=dim/2;i++) {
            R2[i] = pf[(istart-1)+2*i-shift]/pf[(istart-1)+dim-shift];
        }
        ret += EvenSimplexVolumeWrapper<bkzfloat>(R2,dim/2,opt_volume_prob,1.0) ;
        ret *= 0.5;
        ret *= bkzconstants::vol_unit_ball(dim);
        ret *= pow(clim,dim);
        return ret; 
        
    }

    template <typename T> bkzfloat Approx_final_count(std::vector<bkzfloat>& pf,std::vector<T>& c,bkzfloat clim,int istart,int iend) {
        bkzfloat volume = Approx_final_volume(pf,clim,istart,iend);
        for (int i=istart;i<=iend;i++) volume /= (bkzfloat)c[i];
        //cout << "count=" << volume*0.5 << endl;
        return volume * 0.5;    //0.5 is by symmetry
    }
        
        
    
template <typename T>    bkzfloat optimize_pruning_function_crossentropy(std::vector<bkzfloat>& pf,std::vector<T>& c,bkzfloat ptarget,bkzfloat clim,int timemax,int opttarget,int sharpness,int pfinit,int vl,int option=0,std::string optst="") {
        //pfinit=0 ... use input pf as the initial bounding function
        init_pruning_func();
        int n = c.size()-1;

        //Algorithm parameters
        int NC = 10;
        int NCStart = NC;
        int noupdateNC=0;
        int noupdate_limit=10;  //break after this number of loops without list update
        int maxN = max(15,n/10);
        double mmfactor = pow(0.1,1.0/noupdate_limit);

        noupdate_limit=10;
        mmfactor = pow(0.1,1.0/noupdate_limit);

        //Extract and set options
        std::map<std::string,std::string> options;
        ExtractOptions(options,optst);

        std::ofstream of;
        if (options["logfile"]!="") {
            of.open(options["logfile"],ios::app);
        }
        if (options["timemax"]!="") {
            timemax = atoi( options["timemax"].c_str());
            if (timemax < 1) timemax = 1;
        }
        if (options["updatemax"]!="") {
            noupdate_limit = atoi( options["updatemax"].c_str());
        }

        if (vl>=3) {
            ccout << "Cross-entropy based pruning function optimization" << endl;
            ccout << "dim=" << n << " ptarget=" << ptarget << endl;
            ccout << "radius=" << clim << endl;
        }
        if (ptarget<0) {
            ccout << "OPT_pfunc_CE: ptarget error " << ptarget << endl;
            exit(0);
        }
        

        //temporal values
        std::vector<bkzfloat> newsample,ppftemp,minpf;
        newsample.resize(n+1);
        ppftemp.resize(n+1);
        minpf.resize(n+1);
        bkzfloat maxbenefit = -1;
        bkzfloat mincost = -1;
        bkzfloat tprob;
        bkzfloat tcost;

        int parallel=1;
        
        if ((opttarget==target_prob) && (ptarget >= 1.0)) {
            for (int i=1;i<=n;i++) pf[i] = 1.0;
            return FullENUMCost<bkzfloat>(c,1,n,(T)clim);
        }

        std::vector<bkzfloat> center;
        std::vector<bkzfloat> sigma; 
        center.resize(n+1);
        sigma.resize(n+1);

        std::vector<bkzfloat> lb;
        lb.resize(n+1);
        if (options["samplebound"]=="true") {
            ccout << "LB=[";
            for (int i=1;i<n;i++) {
                lb[i] = ibeta_inv<bkzfloat>(0.5*i,0.5*(n-i),ptarget);
                ccout << lb[i] << " ";
            }
            ccout << endl;
        }

        double startopttime = clock();

        //initial approx
        if (opttarget==target_prob) {
            if (pfinit==1) {
                bkzfloat delta = Detect_delta(c,1,n,INPUT_NONSQUARED);
                SetApproximatePruningFunction(pf,n,delta,ptarget);
            }
        } else {
            //Todo: for target_volume
        }

        if (opttarget == target_prob) {
            tprob = Rigid_lower_prob(pf,1,n,parallel,ptarget) + Rigid_upper_prob(pf,1,n,parallel,ptarget);
            tprob *= 0.5;
        } else {
            tprob = Approx_final_count(pf,c,clim,1,n);
        }

        char strategy;
        if (tprob > ptarget) strategy = 0;  //Searching low probs.
        if (tprob <= ptarget) strategy = 1;  //Searching high probs.

        double alpha = 1.0;
        double skip = 1.0;
        if (strategy == 0) skip = 10.0;
        double mult = 1.0;  //Surrogate optimization optputs SurrogateCost>=ptarget * mult

        tcost = Rigid_upper_cost(pf,c,clim,1,n,parallel,ptarget);

        if (vl>=1) {
            ccout << "start_prob=" << tprob << " ";
            ccout << "start_cost=" << tcost << "    \r";
            ccout.flush();
        }

        //Regist the data to List
        CoeffList<bkzfloat> CL;    //for generating initial center for surrogating optimization

        for (int i=1;i<=n;i++) {
            center[i] = pf[i];
            sigma[i] = center[i] * 0.3;
        }

        int loop=0;
        double startcputime = clock();

        int phase = 0;
        bool endflag = false;
        while (endflag==false) {
            if (phase>0) {
                loop++;
                if (CL.X.size()<=3) {
                    //normal sampling
                    if (CL.X.size()>0) CL.getcv(center,sigma);
                    for (int i=1;i<=n;i++) sigma[i] *= mult;
                    samplingpf(pf,center,sigma,lb,(bkzfloat)1.2);
                } else {
                    //prediction sampling
                    double beta = max(0.0,2.0 - 0.2 * (NC-NCStart)) ;
                    if (beta!=0)  {
                        CL.predict(center,sigma,beta,(bkzfloat)1.2);
                    } else {
                        CL.getcv(center,sigma);
                    }
                    for (int i=1;i<=n;i++) sigma[i] *= mult;
                    samplingpf(pf,center,sigma,lb,(bkzfloat)1.2);
                }            
            } else {
                for (int i=1;i<=n;i++)  sigma[i] = center[i] * 0.1 * mult;
                samplingpf(pf,center,sigma,lb,(bkzfloat)1.2);
            }

            
            //processing generated curves
            int pploop=0;

            ValueTable<bkzfloat> tt;
            tt.direction = +1;
            bkzfloat alpha,range;
            bkzfloat logprob = log(ptarget);
            bkzfloat tprob;
            int localloop=0;
            while (++localloop<10) {
                if (tt.data.size()==0) {
                    alpha = 1.0;
                } else 
                if (tt.data.size()==1) {
                    if (tt.data[0][1] < ptarget) {
                        alpha = 1.1;
                    } else {
                        alpha = 0.9;
                    }
                } else {
                    alpha = tt.nextcandidate(logprob,range);
                }
                for (int i=1;i<=n;i++) {
                    if (pf[i] < 1.0) {
                        ppftemp[i] = alpha * pf[i];
                        if (ppftemp[i] > 1.0) ppftemp[i] = 1.0;
                    } else {
                        ppftemp[i] = 1.0;
                    }                        
                }
                if (opttarget == target_prob) {
                    tprob = Rigid_lower_prob(ppftemp,1,n,parallel,ptarget) + Rigid_upper_prob(ppftemp,1,n,parallel,ptarget);
                    if ((tprob < 0) || (tprob>=1)) {
                        tprob = Rigid_lower_prob(ppftemp,1,n,parallel) + Rigid_upper_prob(ppftemp,1,n,parallel);
                    }
                   tprob *= 0.5;
               } else {
                   tprob = Approx_final_count(ppftemp,c,clim,1,n);
               }
               if (tprob>=0) tt.add(alpha,log(tprob));
               if ((ptarget <= tprob) && (tprob <= (1.0+0.1*mult)*ptarget)) break;
            }

                if (tprob > ptarget) {
                    tcost = Rigid_upper_cost(ppftemp,c,clim,1,n,parallel,ptarget);
                    pploop++;
                    if (vl>=3) {
                        bkzfloat variance = CL.getvariance();
                        ccout << "p=" << tprob << " cost=" << tcost << " mult=" << mult << " loop=" << loop << " NC=" << NC << " nu=" << noupdateNC << " var=" << variance << " \r";
                        ccout.flush();
                        if ((CL.X.size()>5) && (variance < 0.0001)) endflag = true;
                    }

                    
                    
                    if (CL.X.size() <= NC) {
                       CL.addelement(ppftemp,tprob,tcost);
                       noupdateNC=0;
                   } else {
                       phase=1;
                       mult = pow(mmfactor, noupdateNC  );
                       mult = max(mult,0.01);
                       if (CL.replacemaxelement(ppftemp,tprob,tcost)==false) {
                           if ((noupdate_limit>1) &&  (noupdateNC++ > noupdate_limit)) {
                               noupdateNC = 0;
                               NC = min(NC+1,maxN);
                               if (vl>=2) {
                                   ccout << "#Set=" << NC << " mincost=" << mincost << " loop=" << loop; 
                                   ccout << " Rtime=" << abs(timemax - (clock() - startcputime)/CLOCKS_PER_SEC) << "  \r";
                                   if (vl>=3) ccout << endl;
                                   ccout.flush();
                               }
                               mult = 1.0;
                           } else {

                           } 
                       } else {
                           noupdateNC=0;
                       }
                       if (option==1) {
                           return mincost;
                       }
                   }
                    if (of.is_open()==true) {
                       of << loop << "\t" << (clock() - startcputime)/CLOCKS_PER_SEC << "\t" << mincost << endl;
                   }
                }

                if (tprob > ptarget) {
                    if (mincost==-1) {
                        mincost = tcost;
                    }
                    if (mincost < tcost) {
                    } else {
                        for (int i=1;i<=n;i++) {
                            minpf[i] = min((bkzfloat)1.0,ppftemp[i]);
                        }
                        mult=1.0;
                        if (mult>1.0) mult=1.0;
                        if (vl>=3) ccout << endl;
                    }
                    mincost = min(mincost,tcost);
                }
            if  ((timemax>0) && (timemax - (clock() - startcputime)/CLOCKS_PER_SEC<0)) break;   //timelimit
            if (CL.X.size()  >= maxN ) break;
        }

        if ((sharpness==pf_crossentropy_exact)  && (opttarget==target_prob)) {
            double alpha = 1.0;
            bkzfloat a = pruning_func::Approx_prob(minpf,1,n);            
            if (a > ptarget) {
                alpha = 1.05;
            } else {
                alpha = 0.95;
            }

            double step= 0.05;
            
            do {
                for (int i=1;i<=n;i++) {
                    ppftemp[i] = pow(minpf[i],alpha);
                }                
                bkzfloat a =  pruning_func::Approx_prob(ppftemp,1,n);            
                step *= 0.5;
                if (a > ptarget) alpha += step;
                if (a < ptarget) alpha -= step;

            } while (fabs(step)>0.001);
            for (int i=1;i<=n;i++) minpf[i] = ppftemp[i];
        }

        if (of.is_open()==true) {
           of << endl;
       }
        for (int i=1;i<=n;i++) pf[i] = minpf[i];
        return mincost;
    }

    template <typename T> bkzfloat optimize_pruning_functionCE(std::vector<bkzfloat>& pf,bkzfloat clim,std::vector<T>& c,bkzfloat tprob,int sharpness,int pfinit,int vl,int timemax=10,char opttarget=target_prob) {
        //Input: valid Pruning function to modify
        //       clim is radius (non-squared)
        //       cd[1..n] is |b*i| (non-squared)
        //       tprob = target probability in (0,1]
        //       vl = verbose level
        //Output: found cost
        
        int n = c.size()-1;
        
        double startopttime = clock();
        adjust_pruning_function(pf);

       if ((opttarget==target_prob) && (tprob>=1.0)) {
           for (int i=1;i<=n;i++) pf[i] = 1.0; //between zero and one
           bkzfloat ret = FullENUMCost<bkzfloat,T>(c,1,n,(T)clim);
           return ret;
       }

        bkzfloat mincost;
        
#ifdef _allow_cachefiles
        //cache file
        std::ostringstream fname;
        init_pruning_func();

        std::string chash = vectorhash(c,1);
        fname << pfsavedir.c_str() << "/cep_n" << n << "c" << clim;

        mkdirRecursive(fname.str().c_str(), 0777);
        fname << "/p" << tprob << "_" << chash << "_" << timemax << "_" << (int)opttarget << "_" << sharpness;
        if (0) {
            loadstdvector(pf,fname.str());
            fname << "_cost";
            LoadElement(mincost,fname.str().c_str());
            if (mincost>0) return mincost;
        } else {
#endif
            mincost = optimize_pruning_function_crossentropy(pf,c,tprob,clim,timemax,opttarget,sharpness,pfinit,vl,2);
#ifdef _allow_cachefiles
        }
        savestdvector(pf,fname.str());
        fname << "_cost";
        SaveElement(mincost,fname.str().c_str());
#endif
        return mincost;
    }
}


#endif
