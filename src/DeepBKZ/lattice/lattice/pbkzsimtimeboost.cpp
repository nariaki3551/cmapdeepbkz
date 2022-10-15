#ifndef _inc_pbkzsimtimeboost_cpp
#define _inc_pbkzsimtimeboost_cpp

#define opt_wantnumloop 0x01

//Constants used for progressive BKZ simulator
namespace bkzconstants{

    int get_Ca(int beta) {
        int Ca = 32;
        if (beta>=80) Ca = 16;
        if (beta>=100) Ca = 8;
        return Ca;
    }

    bkzfloat simulatebkzcost(std::vector<bkzfloat>& cc,std::vector<double>& ualphatable,int beta) {
        bkzfloat ret=0;
        int dim = cc.size()-1;
        bkzfloat partdet = 1.0;
        for (int i=1;i<dim;i++) {
            int bs = min(dim-i+1,beta);
            partdet *= cc[i];
            bkzfloat radius = ualphatable[i] * pow(partdet,1.0/i) *  bkzconstants::ghconstant(i);
            bkzfloat fec = FullENUMCost<bkzfloat>(cc,i,i+bs-1,radius,INPUT_NONSQUARED);
            ret += min(fec,(bkzfloat)(bkzconstants::target_cost_ub_modify(beta) * get_Ca(beta)));
        }
        return ret;
    }

    //Time for LLL subroutine 
    bkzfloat simulatebkz_lll(std::vector<bkzfloat>& cc,int beta) {
        bkzfloat ret=0;
        int dim = cc.size()-1;
        for (int i=1;i<dim;i++) {
            int bs = min(dim-i+1,beta);
            ret += bs*bs*dim;   //blocksize^2 * dim
        }
        return ret;
    }

            
    //Simulating number of loops
    void simulatebkzloop(std::vector<bkzfloat>& cc,std::vector<double>& falphatable,int beta) {
        int dim = cc.size() - 1;
        bkzfloat localdet=1; 

        for (int i=1;i<min(beta,dim);i++) {
            localdet *= cc[i];
        }
        for (int i=1;i<dim;i++) {
            int bs = min(dim-i+1,beta);
            if (bs ==beta) localdet *= cc[i+bs-1];
            //simulating found vector
            double alpha = falphatable[i];
            if (i >= dim-beta + 1) {
                if (dim-i+1 <= 50) {
                    alpha = min(alpha,lattice_tools::smallghconst[dim-i+1]);
                }
            }

            bkzfloat found = alpha * pow(localdet,1.0/bs) * bkzconstants::ghconstant(bs);  
            if (cc[i] > found) {
                bkzfloat diff = 1;
                diff = cc[i]/found;
                cc[i] /= diff;
                cc[i+1] *= diff;
            }
            localdet /= cc[i];
            
        }
    }

    void simtime(bkzfloat& t_lll,bkzfloat& t_enum,double& loop,int dim,int sbeta,int ebeta,int usebeta,int option=0) {
        //Simulating computing time of BKZ-(usebeta) in dim
        //starting from BKZ-(sbeta) beais, target is BKZ-(ebeta)
        if (ebeta > usebeta) {
            //wrong beta setting?
            ccout << "simtime error: end_beta > used_beta" << endl;
            ccout << sbeta << " " << ebeta << " " << usebeta << endl;
        }
        if (ebeta == usebeta) {
            ccout << "warning: target_beta=applied_beta=" << ebeta << ": it may be inifinite loop" << endl;
        }

        bkzfloat venum,vlll,vloop;
        #pragma omp critical
        {
            venum = simtime_enum[dim][sbeta][ebeta][usebeta-ebeta];
            vlll = simtime_lll[dim][sbeta][ebeta][usebeta-ebeta];
            vloop = simtime_loop[dim][sbeta][ebeta][usebeta-ebeta];
        }
        
        if ((venum==0) || (vlll==0) || (vloop==0)) {
            #pragma omp critical
            {
                ccout << "Simulating BKZ Time: " << sbeta << "->" << ebeta << " (by BKZ-" << usebeta << ")                      \r";
                ccout.flush();
            }
            std::vector<bkzfloat> cc;
            PBKZSimulate(cc,dim,sbeta);

            bkzfloat target = simlogfec(dim,ebeta);

            int sloop=0;
            bkzfloat logfec,prevlogfec;

            std::vector<double> falphatable,ualphatable;    //found alpha and used alpha
            std::vector<bkzfloat> probtable;
            SimulateAlphatable(falphatable,ualphatable,dim,usebeta);
            SimulateProbtable(probtable,dim,usebeta);

            logfec = prevlogfec = log(FullENUMCost<bkzfloat>(cc,1,dim));

            t_enum = 0;
            t_lll = 0;

            bkzfloat local_enum,local_lll;
            do {
                local_enum = simulatebkzcost(cc,ualphatable,usebeta);
                local_lll = simulatebkz_lll(cc,usebeta);
                t_enum += local_enum;
                t_lll += local_lll;
                simulatebkzloop(cc,falphatable,usebeta);
                prevlogfec = logfec;
                logfec = log(FullENUMCost<bkzfloat>(cc,1,dim));
                
                sloop++;
                if ((sloop>1+dim/100) || (logfec * 1.0000001> prevlogfec)){   //too many loops, or converged?
                    sloop = 1000 * dim;
                    break;
                }
            } while (logfec > target);
            
            bkzfloat adjust =  (target-prevlogfec) / (logfec-prevlogfec);
            if (adjust < 0) adjust = 0;
            loop = sloop -1 + adjust.convert_to<double>();

            t_enum -= local_enum;
            t_enum += local_enum * adjust;

            t_lll -= local_lll;
            t_lll += local_lll* adjust;
            if (t_enum<=0) t_enum=1e-10;
            if (t_lll<=0) t_lll=1e-10;
            if (loop<=0) t_enum=1e-10;
            simtime_enum[dim][sbeta][ebeta][usebeta-ebeta] = t_enum; 
            simtime_lll[dim][sbeta][ebeta][usebeta-ebeta] = t_lll.convert_to<double>(); 
            simtime_loop[dim][sbeta][ebeta][usebeta-ebeta] = loop;

            #pragma omp critical
            {
                ccout << "Simulating BKZ Time: " << sbeta << "->" << ebeta << " (by BKZ-" << usebeta << ")  ... finished.               \r";
                ccout.flush();
            }
        } else {
            t_enum =  simtime_enum[dim][sbeta][ebeta][usebeta-ebeta]; 
            t_lll =  simtime_lll[dim][sbeta][ebeta][usebeta-ebeta]; 
            loop = simtime_loop[dim][sbeta][ebeta][usebeta-ebeta];
        }        
        bkzconstants::savecachetable();
    } 

    void simtimeopt(bkzfloat& t_lll,bkzfloat& t_enum,int& usebeta,int dim,int sbeta,int ebeta,bkzfloat const_lll,bkzfloat const_enum) {
        if (sbeta==ebeta) {
            t_lll = 0;
            t_enum = 0;
            return;
        }
        
        //find the optimal beta to minimize total cost from BKZ-sbeta to BKZ-ebeta
        //Time is estimated by const_lll*t_lll + const_enum*t_enum
        bkzfloat cost,mincost;
        bkzfloat prevcost;
        int minbeta;
        mincost = -1;
        
        ccout << "Simulating BKZ Time: " << sbeta << "->" << ebeta << "                  \r";
        ccout.flush();
        double loop;
        for (usebeta = ebeta+1;usebeta<dim;usebeta++) {
            simtime(t_lll,t_enum,loop,dim,sbeta,ebeta,usebeta);
            cost = t_lll * const_lll + t_enum * const_enum;
            if ((mincost<0) || (cost < mincost)) {
                mincost = cost;
                minbeta = usebeta;
            }

            if (usebeta >= minbeta+3) {
                break;
            }
        }
        simtime(t_lll,t_enum,loop,dim,sbeta,ebeta,minbeta);
        usebeta = minbeta;
        ccout << "Simulating BKZ Time: " << sbeta << "->" << ebeta << "  ...finished \r";
        ccout.flush();
    }
}


struct BKZProcess {
    int sbeta;
    int ebeta;
    int usebeta;
    bkzfloat t_lll,t_enum,cost;
};


class BKZStrategy {

private:

public:
    
    std::vector<BKZProcess> process;
    bkzfloat totalcost;

    BKZStrategy& operator =(BKZStrategy& newstr) {
        totalcost = newstr.totalcost;
        process.resize(newstr.process.size());
        for (int i=0;i<process.size();i++) {
            process[i].sbeta = newstr.process[i].sbeta;
            process[i].ebeta = newstr.process[i].ebeta;
            process[i].usebeta = newstr.process[i].usebeta;
            process[i].t_lll = newstr.process[i].t_lll;
            process[i].t_enum = newstr.process[i].t_enum;
            process[i].cost = newstr.process[i].cost;
        }
        return *this;
    }
    
    void savetofile(std::string fname) {
        std::ofstream of;
        of.open(fname,ios::trunc);
        of << totalcost << endl;
        of << process.size() << endl;
        for (int i=0;i<process.size();i++) {
            of << process[i].sbeta << "\t";
            of << process[i].ebeta << "\t";
            of << process[i].usebeta << "\t";
            of << process[i].t_lll << "\t";
            of << process[i].t_enum << "\t";
            of << process[i].cost << endl;
        }        
        of.close();
    }

    void loadfromfile(std::string fname) {
        std::ifstream ifs;
        ifs.open(fname,ios::in);
        ifs >> totalcost;
        int ps;
        ifs >> ps;
        process.resize(ps);
        for (int i=0;i<ps;i++) {
            ifs >> process[i].sbeta;
            ifs >> process[i].ebeta;
            ifs >> process[i].usebeta;
            ifs >> process[i].t_lll;
            ifs >> process[i].t_enum;
            ifs >> process[i].cost;
        }        
        ifs.close();
    }
    void addstrategy(int newsbeta,int newebeta,int ubeta,bkzfloat newt_lll,bkzfloat newt_enum,bkzfloat newcost) {
        int size = process.size();
        process.resize(1+size);
        process[size].sbeta = newsbeta;
        process[size].ebeta = newebeta;
        process[size].usebeta = ubeta;
        process[size].t_lll = newt_lll;
        process[size].t_enum = newt_enum;
        process[size].cost = newcost;
    }

    void display() {
        if (process.size()==0) return;
        ccout << process[0].sbeta;
        for (int i=0;i<process.size();i++) {
            ccout << "->(" << process[i].usebeta << ")->" << process[i].ebeta;
        }
        ccout << "   cost=" << totalcost  << "                     " << endl;
    }
    
    
};  


void gen_bkzstrategy(BKZStrategy& BS,int dim,int sbeta,int ebeta,bkzfloat const_lll,bkzfloat const_enum,double loopratio,std::string stringoptions="") {


    bkzconstants::initialize();
    
    std::map<std::string,std::string> options;
    if (stringoptions!="") ExtractOptions(options,stringoptions);

    int parallel = 1;
    if (options["parallel"]!="") {
        parallel = atoi(options["parallel"].c_str());
    }

    int usebeta;
    bkzfloat t_lll,t_enum,cost;
    bkzfloat localcost;
    std::vector<BKZStrategy> BSshort;
    BSshort.resize(ebeta+1);

    bkzfloat mincost;
    int minssbeta;

    ebeta = min(ebeta,dim-1);

    //Computing necessary data
    //simulate BKZ-basis and FEC

    //Precompute time data
    
    #pragma omp parallel  num_threads(parallel)
    {
#ifdef __no__multithreads
        int mythread = 0;
#else
        int mythread = omp_get_thread_num();
#endif
        for (int beta=sbeta;beta<=min(dim-5,ebeta+2);beta++) {
            if (beta%parallel==mythread) {
                bkzconstants::simlogfec(dim,beta);   //empty call for generating table
            }
        }
    }            
    
    #pragma omp parallel for  schedule(dynamic) num_threads(parallel)
    for (int beta=sbeta;beta<=min(dim-1,ebeta+20);beta++) {
        bkzconstants::target_cost_ub(beta);
    }    

    std::vector<std::vector<int> > beta_schedule;
    for (int isbeta=sbeta;isbeta<ebeta;isbeta++) {
        for (int iebeta=isbeta+1;iebeta<=min(isbeta+20,ebeta);iebeta++) {
            for (int iubeta=iebeta+1;iubeta<=min(dim-1,min(ebeta+15,iebeta+30));iubeta++) {
                std::vector<int> bt;
                bt.resize(3);
                bt[0] = isbeta;
                bt[1] = iebeta;
                bt[2] = iubeta;
                beta_schedule.push_back(bt);
            }
        }
    }
    

    //simulate time: sbeta -> ebeta (using BKZ-ubeta)
    #pragma omp parallel for  schedule(dynamic) num_threads(parallel)
    for (int i=0;i<beta_schedule.size();i++) {
        bkzfloat t_lll,t_enum;
        double loop;
        bkzconstants::simtime(t_lll,t_enum,loop,dim,beta_schedule[i][0],beta_schedule[i][1],beta_schedule[i][2]);
        #pragma omp critical
        {
            #ifdef __no__multithreads
                    int mythread = 0;
            #else
                    int mythread = omp_get_thread_num();
            #endif
            if (mythread==0) {
                ccout << beta_schedule[i][0] << " " << beta_schedule[i][1] << " " << beta_schedule[i][2]<< " " << t_lll << " " << t_enum << "             \r";
                ccout.flush();
            }
        }
    }
         
    
    for (int beta=sbeta+1;beta<=ebeta;beta++) {

        ccout << "Simulating optimal strategy: " << sbeta << "->" << beta << "              \r";
        ccout.flush();
        //Find optimal strategy from sbeta to beta
        mincost = -1;
        for (int ssbeta=max(sbeta,beta-10);ssbeta<beta;ssbeta++) {
            //Cost of Optimal[sbeta -> ssbeta] + Direct[ssbeta -> beta]
            //beta-ssbeta is limited less than 10
            //Optimal
            if (ssbeta == sbeta) localcost = 0;
            if (ssbeta > sbeta) {
                if (BSshort[ssbeta].totalcost == 0){
                    //not computed
                    gen_bkzstrategy(BSshort[ssbeta],dim,sbeta,ssbeta,const_lll,const_enum,loopratio,stringoptions);
                }
                localcost = BSshort[ssbeta].totalcost;
            }
            //Direct:ssbeta -> beta
            bkzconstants::simtimeopt(t_lll,t_enum,usebeta,dim,ssbeta,beta,const_lll,const_enum);
            cost = t_lll * const_lll + t_enum * const_enum;
            cost *= loopratio;  //convert to time in seconds

            if ((mincost<0) || (mincost > cost + localcost)) {
                mincost = cost + localcost;
                minssbeta = ssbeta;
                BSshort[beta] = BSshort[ssbeta];
                BSshort[beta].addstrategy(ssbeta,beta,usebeta,t_lll,t_enum,cost);
                BSshort[beta].totalcost = mincost;
            }
        }
        ccout << "Simulating optimal strategy: " << sbeta << "->" << beta << "  finished \r";
        ccout.flush();
    }
    BS = BSshort[ebeta];
    bkzconstants::savecachetable(1);
}

#endif