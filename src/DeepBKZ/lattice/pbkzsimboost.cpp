#ifndef _inc_pbkz_simulate_boost
#define _inc_pbkz_simulate_boost


#include "bkzconstants.hpp"

void PBKZParam(double& alpha,double& r,bkzfloat& prob,int beta) {
    //Base parameter of PBKZ
    r = bkzconstants::mainparam_r(beta);
    prob = bkzconstants::mainparam_prob(beta);
    alpha = bkzconstants::mainparam_alpha(beta);
}


#define modesingle 0x01
#define modeGNR 0x02
#define modeapprox 0x04

template <typename T> T ibeta_inv_wrapper(T a,T b,T x) {
#ifdef _ibeta_approx_wrapper
    if (a==0) {
    }
    if (b==1) {
        return mypower(x,(double)a);
    }
    if (x==0) {
        return 0;
    }
    if (x==1) {
        return 1;
    }

    if (x<1e-20) {
        T lower = 0;
        T upper = 1.0;
        T lv = 0;
        T uv = 1.0;
        while (1) {
            T med = (lower+upper)/2;
            T mv = ibeta(a,b,med);
            if (abs(mv-x)/x < 1e-10) {
                return med;
            }
            if ((lv<x) && (x <= mv)) {
                upper = med;
                uv = mv;
            } else 
            {
                lower = med;
                lv = mv;
            }
        }
    } else {
        return  ibeta_inv<T>(a,b,x);
    }
#else
            return  ibeta_inv<T>(a,b,x);
#endif
}

template <typename T> void ENUMCostLB_maketable(std::vector<T>& lbt,int istart,int iend,bkzfloat prob,char mode=modesingle) {
    lbt.resize(iend+1);
    int dim = iend-istart+1;
    for (int i=iend;i>=istart+2;i--) {  //+2 by omitting last parts
        bkzfloat lf;
        if (i-istart>0) {
            bkzfloat aa = 0.5*(iend-i+1);
            bkzfloat bb = 0.5*(i-istart);
            if (mode==modesingle) {
                lf = ibeta_inv_wrapper<bkzfloat>(aa,bb,prob);   //lower bound of radius
                lf = mypower(lf,0.5*(iend-i+1));
            } 
            if (mode==modeGNR) {
                lf = aa * beta(aa,bb) * prob;
            } 
        } else {
            lf = 0;
        }
        lbt[i] = (T)lf;
    }
}

template <typename T> bkzfloat ENUMCostLBwithtable(std::vector<T>& c,std::vector<T>& lbt,int istart,int iend,bkzfloat radius,char input=INPUT_NONSQUARED) {

    bkzfloat cost,lcost;
    lcost = 1.0;
    int dim = iend-istart+1;
    for (int i=iend;i>=istart+2;i--) {  //+2 by omitting last parts
        bkzfloat lf = lbt[i];
        if (input==INPUT_NONSQUARED) {
            lcost *= radius / c[i];
        } else {
            lcost *= radius / sqrt(c[i]);
        }
        bkzfloat localcost = lcost * bkzconstants::vol_unit_ball(iend-i+1) * lf;
        cost += localcost;
    }
    return 0.5*cost;        //halved by symmetry
}

template <typename T> bkzfloat ENUMCostLB(std::vector<T>& c,int istart,int iend,bkzfloat radius,bkzfloat prob,char input=INPUT_NONSQUARED,char mode=modesingle,bkzfloat numbases=-1,char istargetvolume=0) {
    bkzfloat cost,lcost;
    lcost = 1.0;
    if (prob>1.0) prob=1.0;
    int dim = iend-istart+1;
    for (int i=iend;i>=istart+2;i--) { 
        bkzfloat lf;
        if (i-istart>0) {
            bkzfloat aa = 0.5*(iend-i+1);
            bkzfloat bb = 0.5*(i-istart);
            if (istargetvolume==1) bb++;
            if (mode & modesingle) {
                //cout << aa << " " << bb << " " << prob << endl;
                lf = ibeta_inv_wrapper<bkzfloat>(aa,bb,prob);   //lower bound of radius
                //cout << "lf=" << lf << endl;
                lf = mypower(lf,0.5*(iend-i+1));
            } 
            if (mode & modeGNR) {
                if (numbases==-1) {
                    //Assume infinitely many bases
                    lf = aa * beta(aa,bb) * prob;
                } else {
                    lf = ibeta_inv_wrapper<bkzfloat>(aa,bb,prob/numbases);
                    lf = mypower(lf,0.5*(iend-i+1))*numbases;
                    //ccout << i << " " << lf<< " " << aa << " " << bb << endl;
                }
            } 
        } else {
            lf = 0;
        }
        if (input==INPUT_NONSQUARED) {
            lcost *= (bkzfloat)(radius / (bkzfloat)c[i]);
        } else {
            lcost *= (bkzfloat)(radius / sqrt((bkzfloat)c[i]));
        }
        bkzfloat localcost = lcost * bkzconstants::vol_unit_ball(iend-i+1) * lf;
        cost += localcost;
        if (mode & modeapprox) {
            if (i < iend -10) {
                if (localcost < cost * 1e-4) {
                    break;
                }
            }
        }
    } 
    return 0.5*cost;        //halved by symmetry
}

template <typename T> bkzfloat ENUMCostLBsimple(std::vector<T>& c,int istart,int iend,bkzfloat radius,bkzfloat prob,char input=INPUT_NONSQUARED,char mode=modesingle,bkzfloat numbases=-1,char istargetvolume=0) {
    bkzfloat cost,lcost;
    for (int i=iend;i>=istart+2;i--) { 
        bkzfloat lf;
        lf = ibeta_inv_wrapper<bkzfloat>(0.5*(iend-i+1), 0.5*(i-istart),prob);   //lower bound of radius
        lf = pow(lf,0.5*(iend-i+1));
        lcost *= radius / c[i];
        bkzfloat localcost = lcost * bkzconstants::vol_unit_ball(iend-i+1) * lf;
        cost += localcost;
    }
    return 0.5*cost;        //halved by symmetry
}


template <typename T> bkzfloat ENUMCostLBvol(std::vector<T>& c,int istart,int iend,bkzfloat radius,bkzfloat vol,char input=INPUT_NONSQUARED,char mode=modesingle) {
    //Finding lower cost of enumeration subroutine s.t. R^n*vol(C)/det(L) = vol
    bkzfloat cost,lcost;
    lcost = 1.0;

    //adjusting volume
    bkzfloat v = vol;
    for (int i=istart;i<=iend;i++) v*= c[i];
    v /= pow(radius,iend-istart+1);
    v /= bkzconstants::vol_unit_ball(iend-istart+1); 
    if (v>=1) v=1.0;

    int dim = iend-istart+1;
    for (int i=iend;i>=istart+2;i--) {  //+2 by omitting last parts
        bkzfloat lf;
        if (i-istart>0) {
            bkzfloat aa = 0.5*(iend-i+1);
            bkzfloat bb = 0.5*(i-istart+3);
            
            if (mode==modesingle) {
                lf = ibeta_inv_wrapper<bkzfloat>(aa,bb,v);
                lf = mypower(lf,0.5*(iend-i+1));
            } 
            if (mode==modeGNR) {
                lf = aa * beta(aa,bb) * v;
            } 
        } else {
            lf = 0;
        }
        if (input==INPUT_NONSQUARED) {
            lcost *= radius / c[i];
        } else {
            lcost *= radius / sqrt(c[i]);
        }
        bkzfloat localcost = lcost * bkzconstants::vol_unit_ball(iend-i+1) * lf;
        cost += localcost;
    }
    return 0.5*cost;        //halved by symmetry
}

void AdjustGHValue(std::vector<bkzfloat>& c,int istart,int iend,double alpha) {
    
    //Adjust cc[istart] so that it is equal to alpha*(GH[istart..iend])
    int dim = iend-istart+1;
    c[istart] = 0;
    for (int i=istart+1;i<=iend;i++) {
        c[istart] += log(c[i]);
    }
    c[istart] /= (-1.0+dim);
    c[istart] -= log(bkzconstants::vol_unit_ball(dim))/(-1.0+dim);
    c[istart] += 1.0 * dim / (dim-1.0) * log(alpha);
    c[istart] = exp(c[istart]);
}

void SimulateLLL(std::vector<bkzfloat>& cc,int dim,bkzfloat delta=1.022) {
    cc.resize(1+dim);
    bkzfloat q = exp(- log(delta) * 4.0*dim/(dim-1.0));
    int i;
    cc[1] = 1.0;
    for (int i=2;i<=dim;i++) cc[i] = cc[i-1] * sqrt(q);
}


void PBKZSimulate_first(std::vector<bkzfloat>& cc,int dim,int blockbeta,double alpha,bkzfloat prob) {

    cc.resize(1+dim);
    cc[dim] = 1.0;
    for (int i=dim-1;i>=1;i--) {
        //Simulate: apply BKZ-(beta,alpha,p) for c[i..i+beta-1]
        int bs = min(dim-i+1,blockbeta);  //blocksize
        bkzfloat lp = 0.5 * prob * pow(alpha,bs);  //# lattice point pair (v,-v) in the searching area
        bkzfloat ratio = lp * beta(lp,(1.0+bs)/bs);
        double alpha2 = alpha * ratio.convert_to<double>() ; //expected minimum length of found vector
        if (bs<50) alpha2 = max(alpha2,lattice_tools::smallghconst[bs]);
        AdjustGHValue(cc,i,i+bs-1,alpha2);
    }
}

void PBKZSimulate_second(std::vector<bkzfloat>& cc,std::vector<double>& falphatable,std::vector<double>& ualphatable,std::vector<bkzfloat>& probtable, int dim,int blockbeta,double alpha,bkzfloat prob) {
    //alphatale[...] is radius parameter in simulated GS basis
    //i.e., constant alpha s.t. |b*i| = alpha*GH(sub-basis)
    using namespace bkzconstants;

    bkzfloat lp = 0.5 * pow(alpha,blockbeta) * prob;  //# lattice point pair (v,-v) in the searching area
    bkzfloat ratio = lp * beta(lp,(1.0+blockbeta)/blockbeta);

    double alpha2 = alpha * ratio.convert_to<double>() ; //expected minimum length of found vectors
    if (blockbeta<50) {
        alpha2 = max(alpha2,lattice_tools::smallghconst[blockbeta]);
        alpha = alpha2  / ratio.convert_to<double>();
    }
    bkzfloat baseCost;
    bkzfloat radius = alpha*lattice_tools::LatticeGH(cc,dim-blockbeta+1,dim,INPUT_NONSQUARED);
    baseCost = ENUMCostLB(cc,dim-blockbeta+1,dim,radius,prob,INPUT_NONSQUARED);
    
    bkzfloat logbasecost = log(baseCost);
    
    falphatable.resize(dim+1);
    ualphatable.resize(dim+1);
    probtable.resize(dim+1);
    
    for (int i=1;i<=dim-blockbeta+1;i++) {
        falphatable[i] = alpha2;
        ualphatable[i] = alpha;
        probtable[i] = prob;
    } 

    for (int i=dim-blockbeta+2;i<=dim-1;i++) {
        falphatable[i] = 0;
        ualphatable[i] = 0;
        probtable[i] = 0;   //mark to compute
    }
    
    int j=blockbeta/2;
    int mythread =  omp_get_thread_num();
    bkzfloat* localghtable = (bkzfloat*)shared_memory::allocate1<bkzfloat>(177800+mythread,dim); 
    for (int i=dim-blockbeta+2;i<=dim-1;i+=j) {
        int bs =  min(blockbeta,dim-i+1);
        localghtable[i] = lattice_tools::LatticeGH(cc,i,i+bs-1,INPUT_NONSQUARED);
    }
    while (1) {
        for (int i=dim-blockbeta+2;i<=dim-1;i+=j) {
            if (probtable[i] == 0) {
                int bs =  min(blockbeta,dim-i+1);
                bkzfloat tprob = prob;
                bkzfloat step = prob;
                bkzfloat talpha;
                int phase = 0;
                bkzfloat localgh = localghtable[i];
                ValueTable<bkzfloat> tt;
                tt.direction = +1;
                bkzfloat range;
                do {
                    if (tt.data.size()==0) {
                        tprob = prob;
                        for (int k=i-1;k>=dim-blockbeta+2;k--) {
                            if (probtable[k] !=0) {
                                tprob = probtable[k];   //lower bound
                                break;
                            }
                        }
                    } else 
                    if (tt.data.size()==1) {
                        tprob = 0.5;
                        for (int k=i+1;k<=dim-1;k++) {
                            if (probtable[k] !=0) {
                                tprob = probtable[k];   //upper bound
                                break;
                            }
                        }
                    } else { 
                        tprob = tt.nextcandidate(logbasecost,range);            
                    }
                    if (tprob > 0.5) {
                        tprob = 0.5;
                        phase = 1;
                    }
                    talpha = pow(tprob * 0.5,-1.0/bs);
                    radius = talpha* localgh;
                    bkzfloat tcost;
                    tcost = ENUMCostLB(cc,i,i+bs-1,radius,tprob,INPUT_NONSQUARED,modesingle + modeapprox);
                    if ((tprob==0.5) && (tcost < baseCost)) break;
                    tt.add(tprob,log(tcost));
                    //cout << range / tprob << " " << tprob << " " << log(tcost) << endl;
                } while ((tt.data.size()<=2) || (range / tprob > 1e-2));
                probtable[i] = tprob;
                //cout << "i=" << i << " " << tprob << endl;
                if (tprob==0.5) {
                    for (int k=i;k<=dim-1;k++) {
                        probtable[k] = 0.5;
                    }
                }
            }
        }
        j /=2;
        if (j<=1) j = 1;
        bool nocountflag = false; 
        for (int i=dim-blockbeta+2;i<=dim-1;i++) {
            if (probtable[i] == 0) nocountflag = true;
        }        
        if (nocountflag == false) break;
    }
    for (int i=dim-blockbeta+2;i<=dim-1;i++) {
        int bs =  min(blockbeta,dim-i+1);
        bkzfloat alpha = pow(probtable[i] * 0.5,-1.0/bs);
        ualphatable[i] = alpha.convert_to<double>();
        falphatable[i] = ualphatable[i]  * bs / (1.0+bs);
        if (bs<50) falphatable[i]  = max(falphatable[i],lattice_tools::smallghconst[bs]);
        if ((bs>=50) && (falphatable[i]  < pow(2.77258872223978,1.0/bs) )) falphatable[i] = pow(2.77258872223978,1.0/bs);
    }    
}

double totalsimulatingtime = 0;

void AlphatoGS(std::vector<bkzfloat>& cc,std::vector<double>& alpha,int blockbeta) {
    lattice_tools::initialize();
    int dim = alpha.size()-1;
    cc.resize(dim+1);
    cc[dim] = 1.0;
    for (int i=dim-1;i>=1;i--) {
        int bs =  min(blockbeta,dim-i+1);
        if (bs<50) alpha[i] = max(alpha[i],lattice_tools::smallghconst[bs]);
        AdjustGHValue(cc,i,i+bs-1,alpha[i]);
    }
}


void SimulateAlphatable(std::vector<double>& falphatable,std::vector<double>& ualphatable,int dim,int beta) {
    using namespace::bkzconstants;

    bkzfloat value;
    #pragma omp critical
    {
        value = 0;
        if (bkzconstants::ealpha.size() >= dim+1) {
            if (bkzconstants::ealpha[dim].size() >= beta+1) {
                value = bkzconstants::ealpha[dim][beta][0];
            }            
        }
    }
    if (value==0) {
        std::vector<bkzfloat> cc;
        std::vector<bkzfloat> probtable;
        double r;
        double alpha;
        bkzfloat prob;
        ccout << "Simulate alpha table dim=" << dim << " beta=" << beta << "    \r";
        ccout.flush();
        PBKZParam(alpha,r,prob,beta);
        double ss = gettimeofday_sec();
        PBKZSimulate_first(cc,dim,beta,alpha,prob);
        PBKZSimulate_second(cc,falphatable,ualphatable,probtable,dim,beta,alpha,prob);

        #pragma omp critical
        {
            bkzconstants::ealpha[dim][beta].resize(1+dim);
            bkzconstants::ualpha[dim][beta].resize(1+dim);
            bkzconstants::eprob[dim][beta].resize(1+dim);
        }
        for (int i=0;i<=dim;i++) ealpha[dim][beta][i] = falphatable[i];
        ealpha[dim][beta][0] = 1;

        for (int i=0;i<=dim;i++) ualpha[dim][beta][i] = ualphatable[i];
        ualpha[dim][beta][0] = 1;

        for (int i=0;i<=dim;i++) eprob[dim][beta][i] = probtable[i];
        eprob[dim][beta][0] = 1;

    } else {
        vec_copy(falphatable,bkzconstants::ealpha[dim][beta]);
        vec_copy(ualphatable,bkzconstants::ualpha[dim][beta]);
    }
     
}

void SimulateProbtable(std::vector<bkzfloat>& probtable,int dim,int beta) {
    using namespace::bkzconstants;
    bkzfloat value;
    #pragma omp critical
    {
        value = 0;
        if (bkzconstants::ealpha.size() >= dim+1) {
            if (bkzconstants::ealpha[dim].size() >= beta+1) {
                value = bkzconstants::ealpha[dim][beta][0];
            }            
        }
    }
    if (value<=0) {
        std::vector<bkzfloat> cc;
        std::vector<double> falphatable;
        std::vector<double> ualphatable;
        double r;
        double alpha;
        bkzfloat prob;

        ccout << "Simulate prob table dim=" << dim << " beta=" << beta << "    \r";
        ccout.flush();
        PBKZParam(alpha,r,prob,beta);
        double ss = gettimeofday_sec();
        PBKZSimulate_first(cc,dim,beta,alpha,prob);
        PBKZSimulate_second(cc,falphatable,ualphatable,probtable,dim,beta,alpha,prob);

        #pragma omp critical
        {
            bkzconstants::ealpha[dim][beta].resize(1+dim);
            bkzconstants::ualpha[dim][beta].resize(1+dim);
            bkzconstants::eprob[dim][beta].resize(1+dim);
        }
        for (int i=0;i<=dim;i++) ealpha[dim][beta][i] = falphatable[i];
        ealpha[dim][beta][0] = 1;

        for (int i=0;i<=dim;i++) ualpha[dim][beta][i] = ualphatable[i];
        ualpha[dim][beta][0] = 1;

        for (int i=0;i<=dim;i++) eprob[dim][beta][i] = probtable[i];
        eprob[dim][beta][0] = 1;
    } else {
        vec_copy(probtable,bkzconstants::eprob[dim][beta]);
    }
}
        
void PBKZSimulate(std::vector<bkzfloat>& cc,int dim,int beta) {
    //simulating |b*i| of BKZ-(beta,alpha,prob)  (return c[1..dim] is non-squared value)

    bkzfloat value;
    #pragma omp critical
    {
            value = bkzconstants::pbkzsimlengthtable[dim][beta][0];
            //cout << "value=" << value << endl;
    }
    
    if (value==0) {
        #pragma omp critical
        {
            ccout << "Simulate basis after BKZ-" << beta << " in dim-" << dim << "                   \r";
            ccout.flush();
        }
        cc.resize(dim+1);
        bkzconstants::loadcachetable();

        double ss = gettimeofday_sec();
        std::vector<double> falphatable,ualphatable;
        SimulateAlphatable(falphatable,ualphatable,dim,beta);
        AlphatoGS(cc,falphatable,beta);
        totalsimulatingtime += (gettimeofday_sec() - ss);

        bkzconstants::pbkzsimlengthtable[dim][beta].resize(dim+1);
        for (int i=1;i<=dim;i++) bkzconstants::pbkzsimlengthtable[dim][beta][i] = cc[i];

        bkzconstants::pbkzsimlengthtable[dim][beta][0] = 1;
        #pragma omp critical
        {
            ccout << "Simulate basis after BKZ-" << beta << " in dim-" << dim << " ... finished                       \r";
            ccout.flush();
        }
    } else {
        cc.resize(dim+1);
        vec_copy(cc,bkzconstants::pbkzsimlengthtable[dim][beta]);
    }
}

template <typename T> void GSASimulate(std::vector<bkzfloat>& cc,int dim,T r) {
    cc.resize(1+dim);
    cc[1] = 1.0;
    for (int i=2;i<=dim;i++) cc[i] = cc[i-1] * sqrt(r);
}

void GSASimulate_delta(std::vector<bkzfloat>& cc,int dim,double delta) {

    double r = delta_to_r(delta,dim);
    cc.resize(1+dim);
    cc[1] = 1.0;
    for (int i=2;i<=dim;i++) cc[i] = cc[i-1] * sqrt(r);
}


void NormalizeSimulateGSdet(std::vector<bkzfloat>& cc,bkzfloat det) {
    //Adjust |b*i| so that given determinant
    bkzfloat cdet = 1;
    for (int i=1;i<cc.size();i++) cdet *= cc[i];
    int dim = cc.size()-1;
    cdet = pow(det/cdet,1.0/dim);
    for (int i=1;i<cc.size();i++) cc[i] *= cdet;
}     

namespace bkzconstants{

    bkzfloat simlogfec(int dim,int beta) {
        //cccout << "get_simlogfec: dim=" << dim << " beta=" << beta << endl;
        bkzconstants::loadcachetable();
        initialize();
        
        if (dim>2000) {
            ccout << "simulate log-fec: dimension error: " << dim << endl;
            exit(0);
        }
        
        if (dim < beta) {
            ccout << "simulate log-fec: blocksize error: " << dim << " " << beta << endl;
            exit(0);
        }

        bkzfloat value;
        #pragma omp critical
        {
            value = simlogfec_table[dim][beta];
            //ccout << "value=" << value << endl;
        }
        if (value<=0) {
            #pragma omp critical
            {
                ccout << "Computing simlogfec dim=" << dim << " beta=" << beta << "                         \r";
                ccout.flush();
            }

            std::vector<bkzfloat> cc;
            PBKZSimulate(cc,dim,beta);
            
            bkzfloat fec = FullENUMCost<bkzfloat>(cc,1,dim,lattice_tools::LatticeGH(cc,INPUT_NONSQUARED),INPUT_NONSQUARED);
            simlogfec_table[dim][beta] = log(fec);
            #pragma omp critical
            {
                ccout << "Computing simlogfec dim=" << dim << " beta=" << beta;
                ccout << " ... finished                         \r" ;
                ccout.flush();
            }
        }
        savecachetable();
        return simlogfec_table[dim][beta];   
        
    }
}

void gammaHKZSimulate(std::vector<bkzfloat>& cc,int dim,double gamma) {
    cc.resize(1+dim);
    cc[dim] = 1.0;
    lattice_tools::initialize();
    for (int i=dim-1;i>=1;i--) {
        //Simulate: apply BKZ-(beta,alpha,p) for c[i..i+beta-1]
        int bs = dim-i+1;  //blocksize
        double alpha = 1.0;
        if (bs<50) alpha = max(alpha,lattice_tools::smallghconst[bs]);
        if (i==1) alpha = gamma;
        AdjustGHValue(cc,i,i+bs-1,alpha);
    }
}

void CompleteBKZSimulate(std::vector<bkzfloat>& cc,int dim,int beta) {
    cc.resize(1+dim);
    cc[dim] = 1.0;
    lattice_tools::initialize();
    for (int i=dim-1;i>=1;i--) {
        //Simulate: apply BKZ-(beta,alpha,p) for c[i..i+beta-1]
        int bs = min(beta,dim-i+1);  //blocksize
        double alpha = 1.0;
        if (bs<50) alpha = max(alpha,lattice_tools::smallghconst[bs]);
        AdjustGHValue(cc,i,i+bs-1,alpha);
    }
    
}


void PBKZSimulate_by_fec(std::vector<bkzfloat>& cc,int dim,double fec) {
    //Output |b*i| s.t. FEC(B)=given fec

    bkzconstants::loadcachetable();
    int beta = 9;
    if (bkzconstants::simlogfec(dim,beta) < fec) {
        bkzfloat delta = 1.022;
        bkzfloat ddelta = 0.003;
        do {
            delta -= ddelta;
            SimulateLLL(cc,dim,delta);
            auto currentfec = log(FullENUMCost<bkzfloat>(cc,1,dim,lattice_tools::LatticeGH(cc,INPUT_NONSQUARED),INPUT_NONSQUARED));
            if (currentfec < fec) delta+=ddelta;
            ddelta *= 0.5;
        } while (ddelta > 1e-5);
        return;
    }
    
    while (beta<dim) {
        beta++;
        if (bkzconstants::simlogfec(dim,beta) < fec) {
            //Desired basis is between beta-1 and beta
            std::vector<bkzfloat> cc1;
            PBKZSimulate(cc,dim,beta);
            PBKZSimulate(cc1,dim,beta-1);
            double fecu = bkzconstants::simlogfec(dim,beta-1).convert_to<double>();
            double fecl = bkzconstants::simlogfec(dim,beta).convert_to<double>();
            double cgamma = (fec-fecl)/(fecu-fecl);
            for (int i=1;i<=dim;i++) {
                cc[i] = exp(cgamma*log(cc1[i]) + (1.0-cgamma)*log(cc[i]));
            }
            return;
        }
    }
    
}


#endif

