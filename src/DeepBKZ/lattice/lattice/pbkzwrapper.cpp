#ifndef _inc_pbkzwrapper
#define _inc_pbkzwrapper

//Progressive BKZ main

#include "pbkzsimboost.cpp"
#include "bkzlibwrapper.hpp"
#include "enumwrapper.cpp"

template <typename T> inline void EraseZeroVectors(LatticeBasis<T>& B) {
    int j = 0;
    for (int i=0;i<B.L.NumRows();i++) {
        char nonzeroflag=0;
        int len = B.L[i].length();
        for (int k=0;k<len;k++) {
            if (B.L[i][k] != 0) {
                nonzeroflag = 1;
                break;
            }
        }
        if (nonzeroflag==1) {
            if (i!=j) swap(B.L[j],B.L[i]);
            j++;
        } else {
        }
    }
    B.L.SetDims(j,B.L.NumCols());
    
}


mat_ZZ BUCandidates;

template <typename T> int BasisUpdate(LatticeBasis<T>& B,mat_ZZ*U,int istart,int iend,FoundENUMVectors& EV,int vl,int option,char memoryflag=0) {
    //Insert vectors after B[istart]
    EV.duplicate();
    int nonzero = EV.countnonzero();
    if (nonzero==0) return 0;

    double ss = gettimeofday_sec();
    istart = EV.getminproj();
    
    int n = B.L.NumCols();
    int m = B.L.NumRows();

#ifdef __no__multithreads
    int th = 1;
#else
    int th = omp_get_max_threads();
#endif
    mat_ZZ** bubp;
    #pragma omp critical 
    {
        bubp = (mat_ZZ**)shared_memory::allocate1<mat_ZZ*>(256,th);
    }
#ifdef __no__multithreads
    mat_ZZ* BUBuffer = bubp[0];
#else
    mat_ZZ* BUBuffer = bubp[omp_get_num_threads()-1];
#endif
    if (BUBuffer==0) BUBuffer = new mat_ZZ;

    
    if ((nonzero > BUBuffer->NumRows()) || (n > BUBuffer->NumCols())) {
       BUBuffer->SetDims(nonzero,n);
    }
    if ((nonzero > BUCandidates.NumRows()) || (n > BUCandidates.NumCols())) {
       BUCandidates.SetDims(nonzero,n);
    }

    if (B.L.NumRows() < iend+nonzero ) {
        B.L.SetDims(iend + nonzero,n);
    }


    //escape to buffer memory
    for (int i=0;i<nonzero;i++) {
        (*BUBuffer)[i] = B.L[iend+i];
    }
    T* cback= (T*)shared_memory::allocate1<T>(750,iend+nonzero+1);
    for (int i=istart;i<=iend+nonzero;i++) cback[i] = B.gs.c[i];

    //make zero buffer
    for (int i=iend;i>=istart;i--) {
        swap(B.L[i+nonzero-1],B.L[i-1]);
    }

    //Here, B[istart...istart+nonzero-1] is free
    B.GScomputed = istart-1;
    
    //insert vectors (greedy)
    int j = 0;
    int numinsert=0;
    for (int i=0;i<EV.data.size();i++) {
        if (EV.data[i].norm > 0) {
            coeff_to_zz(BUCandidates[j++],EV.data[i].coeffs,B.L,EV.data[i].projection-1 + nonzero );
        }
    }
    int numcandidates = j;
    
    j = istart;
    T localdet = 1;
    T newdet=1;

    while (1) {
        double minproj = -1;
        int index;
        for (int i=0;i<numcandidates;i++) {
            double proj = ProjectedLength(B,1,j-1,BUCandidates[i]);
            if (proj>0.5) {
                if ((minproj == -1) || (proj < minproj)){
                    minproj = proj;
                    index = i;
                }
            }
        }

        if (minproj>0) {
            //the first inserting vector must improve the basis
            localdet *= cback[j];
            newdet *= minproj;
            if (newdet > localdet * 0.999) {
                break;
            }
            B.L[j-1] = BUCandidates[index];
            B.updateGSBasis(j,j);
            numinsert++;
            j++;
        } else {
            break;
        }
    }

    for (int i=j;i<istart+nonzero;i++) {
        clear(B.L[i-1]);
    }
    double delta = 0.99;
    int final_iend = local_LLL(B,0,delta,istart,iend+nonzero);
    
    if (final_iend != iend) {
        B.updateGSBasis();
        B.gs.displaymu();
        B.gs.displayc();
        exit(0);
    }


    B.GScomputed = final_iend;
    for (int i=0;i<nonzero;i++) {
        B.L[iend+i] = (*BUBuffer)[i] ;
    }
    for (int i=0;i<nonzero;i++) {
        B.gs.c[iend+i+1] = cback[iend+i+1];
        B.gs.cd[iend+i+1] = B.gs.c[iend+i+1]*B.gs.c[iend+i+1];
    }
    //Todo: this re-computation of GS can be replaced by shifting 
    if (memoryflag==0) B.L.SetDims(m ,n);
    B.updateGSBasis(istart,iend);
    return numinsert;
}


#include "pbkz_smalllattice.cpp"
#include "pfuncwrapper.cpp"

template <typename PFLOAT> int EnumthenUpdate(LatticeBasis<PFLOAT>& B,mat_ZZ*U,int istart,int iend,double radius,bkzfloat prob,int pfoption,int vl=0,std::string stringoptions="") {
        //Note: after update, B.gs.c and B.gs.cd have correct values for 1..dim
        //      but mu after index iend+1 are in old values
    
        double ss = gettimeofday_sec(); //timer
        std::map<std::string,std::string> options;
        if (stringoptions!="") ExtractOptions(options,stringoptions);

        B.updateGSBasis(1,iend);        //Assume that GS-basis is computed 
        for (int j=istart;j<=iend;j++) SizeReduce(B,U,j,1,nogsupdate);    

        
        
        int parallel = 1;
        if (options["parallel"]!="") {
            parallel = atoi(options["parallel"].c_str());
        }

        if ((vl>=1) && (pfoption!=0)) {
            ccout << "set pruning function [" << istart << ":" << iend << "]    \r";
            ccout.flush();
        }
        radius = min(radius,(double)B.gs.c[istart]);

        long int elim = 10000000;
        if (options["elim"]!="") {
            elim = atoi(options["elim"].c_str());
        }
        
        PruningFunction PF;
        bkzfloat realprob = pruning_func::SetPruningFunction(B,PF, istart,iend,radius , prob, pfoption,vl);
        vec_ZZ zerovec;
        
        FoundENUMVectors EV;
        EV.clearmemory(16);
        
        int finishmode = finish_mode_nowaste;
        int enumoption = enum_mode_find_short_projection;
        
        int slide = iend-istart;
        
        int dim = iend-istart+1;
        
        if ( dim <= 40) {
            parallel = 1;
        }
        
        if ((vl>=1) && (pfoption!=0)) {
            ccout << "process enum [" << istart << ":" << iend << "]    \r";
            ccout.flush();
        }
        
        if (parallel==1) {
            if (dim <= 140) {
#ifndef __lite__compile
                lattice_enum::ENUMCoreboost_double(EV,slide,istart,iend-istart+1,elim,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption);
#else
                lattice_enum::ENUMCoreboost<double,double>(EV,slide,istart,iend-istart+1,elim,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption);
#endif
            } else {
                lattice_enum::ENUMCoreboost<long double,long double>(EV,slide,istart,iend-istart+1,elim,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption);
            }
        } else {
            if (dim <= 140) {
#ifndef __lite__compile
                lattice_enum::ENUMCoreboostmt_double(EV,slide,istart,iend-istart+1,elim,parallel,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption,finishmode);
#else
    #ifndef __no__multithreads
                lattice_enum::ENUMCoreboostmt<double,double>(EV,slide,istart,iend-istart+1,elim,parallel,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption,finishmode);
    #else
                lattice_enum::ENUMCoreboost<double,double>(EV,slide,istart,iend-istart+1,elim,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption);
    #endif
#endif
            } else {
    #ifndef __no__multithreads
                lattice_enum::ENUMCoreboostmt<long double,long double>(EV,slide,istart,iend-istart+1,elim,parallel,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption,finishmode);
    #else
                lattice_enum::ENUMCoreboost<long double,long double>(EV,slide,istart,iend-istart+1,elim,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption);
    #endif
            }
        }

        if ((vl>=1) && (pfoption!=0)) {
            ccout << "basis update [" << istart << ":" << iend << "]    \r";
            ccout.flush();
        }

        int ret=0;
        if (iend-istart+1 <= 60) {
            ret =  BasisUpdateGSsmall<double,PFLOAT>(B,U,istart,iend,EV,VL0,0,1);
        } else 
        if (iend-istart+1 <= 90) {
           ret =  BasisUpdateGSsmall<long double,PFLOAT>(B,U,istart,iend,EV,VL0,0,1);
        } else {
           ret =  BasisUpdateGSsmall<float15,PFLOAT>(B,U,istart,iend,EV,VL0,0,1);
        }
        return ret;
}
#include "bkzpreprocessboost.cpp"
#include "pbkzsimtimeboost.cpp"

template <typename PFLOAT> void PrunedBKZ(LatticeBasis<PFLOAT>& B,mat_ZZ*U,int beta,double alpha,double prob,double tour,int vl=0,std::string stringoptions="") {


        std::map<std::string,std::string> options;
        ExtractOptions(options,stringoptions);

        int istart=1;
        int iend=B.dim;
        if (options["istart"]!="") {
            istart = atoi( options["istart"].c_str());
            if (istart < 1) istart = 1;
            if (istart > iend) istart = iend;
        }
        if (options["iend"]!="") {
            iend = atoi( options["iend"].c_str());
            if (iend > B.dim) iend = B.dim;
        }

        if (options["doLLL"]!="") {
            BigLLL(B.L,0,0.999,VL1);
        }

        B.updateGSBasis();
        
        
        int n = B.dim;
        double ss = gettimeofday_sec();
        int tswap = 0;

        int t = 0;
        int first_update_index=0;
        while (1) {
            int totalinsert = 0;
            int iistart = max(istart,first_update_index-beta+1);
            first_update_index = 0;
            for (int i = istart;i<iend-1;i++) {
                if (vl>=2) ccout << "time=" <<   gettimeofday_sec() - ss << " tour=" << t+1 << "/" << tour << " index=" << i << " |b1|=" << LengthOf(B.L[0]) << endl;
                int bstart = i;
                int bend = min(bstart+beta-1,iend);
                double radius = alpha * lattice_tools::LatticeGH(B.gs.c,bstart,bend);
                if (radius > B.gs.c[bstart]) radius = B.gs.c[bstart];                
                //Returned value is number of inserted vectors

                int ninsert = EnumthenUpdate<PFLOAT>(B,U,bstart,bend,radius,prob,0,vl,"parallel=1");
                if ((ninsert!=0) && (first_update_index==0)) first_update_index = i;
                totalinsert += ninsert;
                for (int i=bstart;i<=bend;i++) {
                    SizeReduce(B,U,i);
                }
                //local_LLL(B,0,0.999,bstart,bend,stnormal,VL0,swap);
                //local_LLL(B,0,0.999,1,iend,stnormal,VL0,swap);
            }
            int swap;
            local_LLL(B,0,0.999,1,iend,stnormal,VL0,swap);
            tswap += swap;

            if (vl>=1) {
                ccout << "time=" <<   gettimeofday_sec() - ss << " tour=" << t+1;
                if (tour>0) ccout << "/" << tour;
            }
            ccout << " |b1|=" << LengthOf(B.L[0]) << endl;
            //ccout << "TI=" << totalinsert << " TS=" << tswap << endl;
            t++;
            if ((tour>0) && (t>=tour)) break;
            if (totalinsert==0) break;
        }
}


template <typename PFLOAT> int DetectBeta(LatticeBasis<PFLOAT>& B,int istart,int iend,int checkupper = -1) {

    //Find PBKZ level: FEC(beta) > FEC(B) > FEC(beta+1)
    
    int dim = iend-istart+1;
    PFLOAT radius = lattice_tools::LatticeGH(B.gs.c,istart,iend);
    bkzfloat fec = FullENUMCost<bkzfloat,PFLOAT>(B.gs.c,istart,iend,radius);
    bkzfloat logfec = log(fec);
    for (int beta=10;beta<dim;beta++) {
        bkzfloat simfec = bkzconstants::simlogfec(dim,beta); 
        if (logfec > simfec) return beta-1; //level of log-fec between is sim[beta] and sim[beta-1] 
        if (checkupper > 0) {
            if (beta >= checkupper) return beta;
        }
    }
    return dim;
}

//This subroutine is called by progressive BKZ routine
template <typename PFLOAT> int processoneblock(LatticeBasis<PFLOAT>& B,mat_ZZ*U,int beta,int bstart,int bend,double falpha,double ualpha,bkzfloat uprob,int parallel,
        int vl, double& totalpptime,double& totalpftime,double& totalenumtime,double& totalenumcputime ) {
    //preprocess
    int pfoption=0;
    int bs = bend-bstart+1; //blocksize 
    int updateflag=0;

    bool enumflag=true;
    
    if ((beta >= 50) && (bs >= 30)) {
        int Ca = 4;
        bkzfloat costlimit_lb = bkzconstants::target_cost_lb(beta) * Ca;
        bkzfloat costlb;
        double lss = gettimeofday_sec(); //timer
        //if (bstart==18) debug_output=1;
        if (bend-bstart<=80) {
            //debug
            costlb = ProgressiveBKZPreprocess<PFLOAT,double>(B,U,bstart,bend,costlimit_lb,falpha,ualpha,uprob,vl);
        } else {
            costlb = ProgressiveBKZPreprocess<PFLOAT,long double>(B,U,bstart,bend,costlimit_lb,falpha,ualpha,uprob,vl);
        }

        totalpptime += (gettimeofday_sec()-lss);
        if (costlb < 1e+6) parallel=1;  //to prevent slowdown from overhead
        if (costlb > costlimit_lb) enumflag = false;
        
    }
    if (bend-bstart+1<=50) parallel=1;
    //main process
    if (vl>=2) ccout << "ghcheck:" << B.gs.c[bstart] << " " << lattice_tools::LatticeGH(B.gs.c,bstart,bend)  << " " << B.gs.c[bstart]/lattice_tools::LatticeGH(B.gs.c,bstart,bend) << " " << falpha << endl;
    if (B.gs.c[bstart] < falpha * lattice_tools::LatticeGH(B.gs.c,bstart,bend)) enumflag = false; 
    
    if (enumflag==true) {
        double lss = gettimeofday_sec(); //timer
        auto gh = lattice_tools::LatticeGH(B.gs.c,bstart,bend);
        double radius = ualpha * (double)gh;
        auto prev = B.gs.c[bstart];
        updateflag += EnumthenUpdate<PFLOAT>(B,U,bstart,bend,radius,uprob,pfoption,vl,"parallel=" + to_stdstring(parallel));

        totalenumtime += lattice_enum::current_etime;
        totalenumcputime += lattice_enum::current_cputime;

        //Post process
        for (int i=bstart;i<=bend;i++) {
            SizeReduce(B,U,i);
        }
    }
    return updateflag;
}

template <typename PFLOAT> void UpdateSubWindow(mat_ZZ& WU,LatticeBasis<PFLOAT>& B,mat_ZZ* U,int wstart) {
    
    int dim = WU.NumRows();

    mat_ZZ Btemp;
    Btemp.SetDims(dim,B.L.NumCols());
    
    for (int i=0;i<dim;i++) {
        for (int j=0;j<dim;j++) {
            if (WU[i][j] == 1) {
                Btemp[i] += B.L[j+wstart-1];
            } else 
            if (WU[i][j] == -1) {
                Btemp[i] -= B.L[j+wstart-1];
            } else 
            if (WU[i][j] != 0) {
                Btemp[i] += WU[i][j] * B.L[j+wstart-1];
            }
        }
    }
    for (int i=0;i<dim;i++) {
        B.L[i+wstart-1] = Btemp[i];
    }

    B.GScomputed = min(B.GScomputed,wstart-1);
    B.updateGSBasis(1,wstart+dim-1);;

    if (U!=0) {
        Btemp.SetDims(dim,(*U).NumCols());
        for (int i=0;i<dim;i++) {
            clear(Btemp[i]);
        }
        for (int i=0;i<dim;i++) {
            for (int j=0;j<dim;j++) {
                if (WU[i][j] != 0) {
                    Btemp[i] += WU[i][j] * (*U)[j+wstart-1];
                }
            }
        }
        for (int i=0;i<dim;i++) {
            (*U)[i+wstart-1] = Btemp[i];
        }
    }
}

template <typename PFLOAT> void StrategiedBKZ(LatticeBasis<PFLOAT>& B,mat_ZZ*U,BKZStrategy& BS,int vl=0,std::string stringoptions="") {

    B.updateGSBasis();
    bkzconstants::loadcachetable();
    
    //ExtractOptions
    std::map<std::string,std::string> options;
    ExtractOptions(options,stringoptions);
    bool flag_betashiftupdate = true;
    if (options["betashiftupdate"]=="false") {
        flag_betashiftupdate = false;
    }

    bool flag_ignoreflat = true;
    if (options["ignoreflat"]=="false") {
        flag_ignoreflat = false;
    }
    
    int endflag=0;
    int beta = 10;
    int iistart = 1;
    int iiend = B.dim;
    if (options["istart"]!="") {
        iistart = atoi( options["istart"].c_str());
        if (iistart < 1) iistart = 1;
        if (iistart > iiend) iistart = iiend;
    }
    if (options["iend"]!="") {
        iiend = atoi( options["iend"].c_str());
        if (iiend > B.dim) iiend = B.dim;
    }

    std::ofstream of;
    if (options["logfile"]!="") {
        of.open(options["logfile"],ios::app);
    }


    int parallel = 1;
    if (options["parallel"]!="") {
        parallel = atoi( options["parallel"].c_str());
        if (parallel < 1) parallel = 1;
    }

    PFLOAT target_length = 0;
    if (options["targetlen"]!="") {
        target_length = boost::lexical_cast<PFLOAT>(options["targetlen"]);
    }
    
    char ignoreflat = false;
    if (options["ignoreflat"]!="") {
        ignoreflat = true;
    }
    
    bool h_reduceSR = false;
    if (options["heuristic"]!="") {
        std::vector<std::string> items;
        boost::split(items, options["heuristic"], boost::is_any_of(","));
        for (size_t i=0;i<items.size();i++) {
            //ccout << "heuristic opt=" << items[i] << endl;
            if (items[i] == "reduceSRCost") h_reduceSR = true;
        }
    }

    double start = gettimeofday_sec();
    double ssim = 0;
    int tour = 1;
    B.updateGSBasis();
    PFLOAT gh = lattice_tools::LatticeGH(B.gs.c,iistart,iiend);
    double totalenumtime = 0;
    double totalenumcputime = 0;
    double totalpptime = 0; //Time for BKZ-preprocess
    double totalpftime = 0;
    totalsimulatingtime=0;
    
    if (gh<=0) {
        ccout << "StrategiedBKZ in pbkzwrapper.cpp: gh error: " << gh << endl;
        B.gs.displayc(iistart,iiend);
        exit(0);
    }
    
    bkzfloat prevfec =  log(FullENUMCost<bkzfloat,PFLOAT>(B.gs.c,iistart,iiend,(PFLOAT)gh));
    if (prevfec<0) {
        ccout << "StrategiedBKZ in pbkzwrapper.cpp: fec error: " << prevfec << endl;
        B.gs.displayc(iistart,iiend);
        exit(0);
    }
    
    int achievedbeta = 0;
    int bsline = 0;
    int betashift = 0;
    
    while (endflag==0) {
        int istart = iistart;
        int iend = iiend;
        if (flag_ignoreflat==true) {
            //narrow the processing range
            while ((istart < iend-1) && ( abs(B.gs.c[istart] -  B.gs.c[istart+1])/B.gs.c[istart] < 0.001)  ) istart++;
            while ((istart < iend-1) && (  abs(B.gs.c[iend] -  B.gs.c[iend-1])/B.gs.c[iend] < 0.001) )  iend--;
            istart = max(iistart,istart-5);
            iend = min(iiend,iend+5);
            if ((istart!=iistart) || (iend!=iiend)) {
                if (vl>=1) ccout << "reset processing range[" << iistart << "," << iiend << "] -> [" << istart << "," << iend << "]                       " << endl;
            }
            gh = lattice_tools::LatticeGH(B.gs.c,istart,iend);
        }
        //SaveLattice(B,"bkzstartdebug");

        //Detect blocksize
        int betadetect = DetectBeta(B,istart,iend,iend-istart+1);
        if (betadetect>=iend-istart+1) return;    //already reduced

        achievedbeta = max(achievedbeta,betadetect);

        while ((bsline < BS.process.size()) && (achievedbeta >= BS.process[bsline].ebeta)) {
            //ccout << "bsline=" << bsline << " " << achievedbeta << " " << BS.process[bsline].ebeta << " " << BS.process.size() << endl;
            bsline++;
        }
        if (bsline == BS.process.size()) {
            endflag = 1;
            break;
        }
        beta = min(BS.process[bsline].usebeta , iend-istart+1);
        bkzconstants::simlogfec(iend-istart+1,beta); //call to make (alpha,prob) table
        
        
        double alpha,r;
        bkzfloat prob;
        PBKZParam(alpha,r,prob,beta);

        int updateflag=0;

        bkzfloat fec = FullENUMCost<bkzfloat,PFLOAT>(B.gs.c,istart,iend,(PFLOAT)gh);
        bkzfloat logfec = log(fec);
        bkzfloat simfec = bkzconstants::simlogfec(iend-istart+1,achievedbeta + 1);    //next target fec 
        if (simfec == 0) {
            exit(0);
        }

        bkzconstants::savecachetable();

        if ((flag_betashiftupdate==true) && (flag_ignoreflat==false))  {
            if (logfec > prevfec) {
                betashift++;
                beta = min(BS.process[bsline].usebeta + betashift , iend-istart+1);
            } else {
                betashift--;
                if (betashift<0) betashift = 0;
            }
        }
        prevfec = logfec;

        if (vl>=1) {
            ccout << "time=" <<   gettimeofday_sec() - start - (totalsimulatingtime - ssim) <<  " tour=" << tour << " applybeta=" << beta; 
            ccout << " basisbeta=" <<  betadetect << " |b1|=" << B.gs.c[istart];
            ccout << "=" << B.gs.c[istart] / gh << "GH";
            ccout << " logfec=" << logfec << "/" << simfec;
            ccout << endl;
            //B.gs.displayc();
        }
        if (of.is_open()==true) {
            of << gettimeofday_sec() - start - (totalsimulatingtime - ssim)  << "\t";
            of << tour << "\t";
            of << beta << "\t";
            of << betadetect << "\t";
            of << B.gs.c[istart] << "\t";
            of << B.gs.c[istart] / gh<< "\t";
            of << logfec << "\t";
            of << simfec << "\t";
            of << totalenumtime << "\t";
            of << totalenumcputime << "\t";
            of << totalpptime << "\t";
            of << totalsimulatingtime - ssim << endl;
        }

        if (endflag==1) break;
        tour++;

        std::vector<double> falphatable,ualphatable;
        std::vector<bkzfloat> probtable;
        SimulateAlphatable(falphatable,ualphatable,iend-istart+1,beta);
        SimulateProbtable(probtable,iend-istart+1,beta);

        //start tour
        int wstart,wend;    //window_start and window_end
        wstart = -1;
        wend = -1;
        int wmaxsize = max(100,2*beta); //max window size (heuristic)
        int wsize;

        LatticeBasis<long double> WB;   //window basis
        //LatticeBasis<float15> WB;   //window basis
        mat_ZZ WU;  //window unitary

        bool usewindow=true;
        
        int window_inner_loop = 0;  //use when h_reduceSR == true
        
        int i = istart-1;
        //for (int i = istart;i<iend-2;i++) {

        while (i<iend-2) {
            i++;
            int bstart = i;
            int bend = min(bstart+beta-1,iend);
            int bs = bend-bstart+1; //blocksize 
            double falpha =  falphatable[i-istart+1];    //alpha expected to find
            double ualpha = ualphatable[i-istart+1];; //alpha used to search
            bkzfloat uprob = probtable[i-istart+1];

            if (usewindow==false) {
                if ((vl>=1) && (beta>=45)) {
                    ccout << "process [" << bstart << "," << bend <<  "]                              \r";
                    ccout.flush();
                }            
                //Not use window
                B.updateGSBasis(istart,bend);
                updateflag += processoneblock(B,U,beta,bstart,bend,falpha,ualpha,uprob,parallel,vl,totalpptime,totalpftime,totalenumtime,totalenumcputime);
                int swapcount;
                local_LLL(B,U,0.999,1,bend,stnormal,min(vl,VL0),swapcount);
            } else {
                //ccout << bstart << "," << bend << "," << wstart << "," << wend << "," << wsize << endl;
                if ((h_reduceSR == true) && (wend>1)){
                    if (bend > wend) {
                        bend = wend;
                    }
                    if ((bstart > wend - 4) && (bend > 0)) {
                    //at the end of window
                        if (window_inner_loop++ < 40) {
                            int wbetalevel = DetectBeta(WB,1,wsize,BS.process[bsline].usebeta);
                            if ((wbetalevel < BS.process[bsline].usebeta) || (window_inner_loop < 6)) {
                                i = wstart-1;
                                if (vl>=1) {
                                    ccout << "beta of local block[" << wstart << "," << wend << "] :" << wbetalevel << "/" << BS.process[bsline].ebeta << " reset index: i=" << i << " (lp=" << window_inner_loop << ")       \r";
                                    ccout.flush();
                                }
                            } else {
                                if (bend != iend) {
                                    if (vl>=1) {
                                        ccout << "beta of local block:" << wstart << "," << wend << "] :" << wbetalevel << "/" << BS.process[bsline].ebeta << "              \r";
                                        ccout.flush();
                                    }
                                }
                                //go to the next window
                                if (wend < iend) {
                                    bend = wend+1;
                                    bstart = bend - beta;
                                }
                            }
                        } else {
                            //go to the next window
                            if (wend < iend) {
                                bend = wend+1;
                                bstart = bend - beta;
                            }
                        }
                    } 
                }
                
                if (bend > wend) {
                    if (wend>0) {
                        if (vl>=1) ccout << "update window[" << wstart << ":" << wend << "]                              \r";
                        if (vl>=2) ccout << endl;
                        
                        //ccout << "update WB gs line:" << endl;
                        //WB.gs.displayc();
                        
                        //SaveLattice(B.L,"/home/a/bkzwork/ideal/BB" + to_string(wstart) + "-" + to_string(wend));
                        UpdateSubWindow(WU,B,U,wstart); //update previous window
                        //SaveLattice(B.L,"/home/a/bkzwork/ideal/B" + to_string(wstart) + "-" + to_string(wend));
                        //SaveLattice(WB.L,"/home/a/bkzwork/ideal/WB" + to_string(wstart) + "-" + to_string(wend));

                        //Blockwised LLL started at i=wstart
                        int llsize = 160;
                        int lstart = wstart;
                        int lend;
                        int localswap,maxshift;
                        while (1) {
                            lend = min(lstart+llsize -1,wend);
                            SmallLocalLLL<long double>(B,U,(long double)0.999,lstart,lend,VL0,localswap,maxshift);
                            for (int i=maxshift;i<=lend;i++) {
                                SizeReduce(B,U,i);
                            }
                            if (localswap > 0) {
                                if ((lstart==istart) && (lend==wend)) break;
                                if ((lstart < maxshift) && (lend==wend)) break;
                                lstart = maxshift - llsize/2;
                                if (lstart < istart) lstart = istart;
                            } else {
                                if (lend == wend) break;
                                lstart += llsize;
                            }
                            lend = min(lstart+llsize -1,wend);
                        }
                    }    
                    //SaveLattice(B.L,"/home/a/bkzwork/ideal/BL" + to_string(wstart) + "-" + to_string(wend));
                    
                    //define new window
                    wstart = max(istart,bstart-10);
                    wend = min(iend,bstart + wmaxsize-1);
                    if (wend > iend-20) wend = iend;    //heuristic extend
                    
                    wsize = wend-wstart+1;
                    if (vl>=2) ccout << "new window[" << wstart << ":" << wend << "]" << endl;

                    if ((tour!=1) || (wstart!=istart)) {
                        //if else, size reduced by the previous loop
                        for (int i=wstart;i<=wend;i++) {
                            ccout << "size reduce:" << i << "/" << wend << "                    \r";
                            ccout.flush();
                            //SizeReduce(B,U,i,max(1,wstart-10));   //reduce the column [wstart-10 .. i]
                            SizeReduce(B,U,i,1,gsupdate);   //reduce the column [wstart-10 .. i]
                        }
                    }

                    WB.L.SetDims(wsize,wsize);
                    WU.SetDims(wsize,wsize);

                    //ccout << "Extract: " << wstart << "," << wend << "   " << endl;
                    ExtractSubBlock(WB,B,wstart,wend);
                    SetIdentity(WU);
                    WB.GScomputed=-1;
                    WB.updateGSBasis();
                    
                    window_inner_loop = 0;
                    
                    if (WB.dim>=140) {
                        BigLLL(WB,&WU,0.999,1,WB.dim,min(vl,1));
                    } else {
                        int swapcount;
                        local_LLL(WB,&WU,0.999,1,WB.dim,stnormal,min(vl,VL1),swapcount);
                    }
                }
                if (vl>=1) {
                    //ccout << "enum(" << bstart << "," << bend << ")";
                    ccout << "process [" << bstart << "," << bend <<  "] (" << bstart-wstart+1 << "," << bend-wstart+1 << ") in window";
                    if (h_reduceSR == true) {
                        ccout << " wloop=" << window_inner_loop;
                    }
                    ccout << "                                                 \r";
                    ccout.flush();
                }
                updateflag += processoneblock(WB,&WU,beta,bstart-wstart+1,bend-wstart+1,falpha,ualpha,uprob,parallel,vl,totalpptime,totalpftime,totalenumtime,totalenumcputime);
                if (vl>=2) ccout << "processoneblock end." << endl;
                int swapcount;
                local_LLL(WB,&WU,0.999,1,bend-wstart+1,stnormal,min(vl,VL2),swapcount);
                if (vl>=2) ccout << endl;
            }
            //display
            if (vl>=2) ccout << "time=" <<   gettimeofday_sec() - start - (totalsimulatingtime - ssim)<<  " index=" << bstart << " beta=" << beta << " alpha=" << ualpha << " prob=" << uprob << endl;
            bkzconstants::savecachetable();
        }       //end of tour
        if (wend>0)  {
            if (vl>=1) ccout << "update window[" << wstart << ":" << wend << "]    \r";
            if (vl>=2) ccout << endl;
            UpdateSubWindow(WU,B,U,wstart); //update previous window
            if (vl>=2) B.gs.displayc();

            for (int i=wstart;i<=wend;i++) {
                if (vl>=1) ccout << "size reduce:" << i << "/" << wend << "                    \r";
                ccout.flush();
                SizeReduce(B,U,i);
            }
        }

        if (updateflag==0) beta++;
        if (vl>=2) ccout << "time=" <<  gettimeofday_sec() - start - (totalsimulatingtime - ssim) <<  " tour=" << tour << " finished" << endl;
        if (options["temporal"]!="") {
            SaveLattice(B,options["temporal"]);
        }

        if (target_length>0) {
            if (B.gs.c[iistart] < target_length) {
                break;
            }
        }
    }
    for (int i=iistart;i<=iiend;i++) {
        ccout << "LastSizeReduce i=" << i << "              \r";
        ccout.flush();
        SizeReduce(B,U,i);
    }
    if (vl>=1) ccout << "progressive bkz finish. time=" <<   gettimeofday_sec() - start - (totalsimulatingtime - ssim) << "                          " << endl;
    bkzconstants::savecachetable();
}

template <typename PFLOAT> void ProgressiveBKZ(LatticeBasis<PFLOAT>& B,mat_ZZ*U,int targetlevel,int vl=0,std::string stringoptions="") {

    //ExtractOptions
    std::map<std::string,std::string> options;
    ExtractOptions(options,stringoptions);

    int betashift = 6;
    if (options["betashift"]!="") {
        betashift = atoi( options["betashift"].c_str());
        if (betashift < 1) betashift = 1;
    }

    int startbeta=15;
    if (options["startbeta"]!="") {
        startbeta = atoi( options["startbeta"].c_str());
        if (startbeta < 1) startbeta = 1;
    }

    BKZStrategy BS;
    for (int i=startbeta-betashift+1;i<=targetlevel;i++) {
        BS.addstrategy(i-1,i,i+betashift-1,0,0,0);  //Last 0's are dummy
    }
    //BS.display();
    StrategiedBKZ(B,U,BS,vl,stringoptions);    
}

#endif

