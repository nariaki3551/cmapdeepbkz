#ifndef _inc_enum_wrapper_cpp
#define _inc_enum_wrapper_cpp

//#include "bkzlibwrapper.hpp"

#define flag_unduplicate 0x01

void coeff_to_zz(vec_ZZ& v,std::vector<int>& coeffs,mat_ZZ& L,int offset) {
    v = coeffs[0]*L[offset];
    for (int i=1;i<coeffs.size();i++) {
        if (coeffs[i]!=0) {
            v += coeffs[i] * L[i+offset];
        }
    }
}


struct FoundVectorInfo {
    double norm = -1;
    int projection = -1;
    std::vector<int> coeffs;
    
    void display() {
        ccout << "[" << norm << " " << " " << projection << ": ";
        for (int i=0;i<coeffs.size();i++) ccout << coeffs[i] << " ";
        ccout << "]" << endl;
    }
    
};

bool operator < (FoundVectorInfo& a,FoundVectorInfo& b) {
    if (a.projection < b.projection) return true;
    if (a.projection == b.projection) {
        if (a.norm < b.norm) return true;
    }
    return false;
}

void clear(FoundVectorInfo& a) {
    a.norm = -1;
    a.projection = -1;
    a.coeffs.resize(0);
    
}

struct FoundENUMVectors {
    int maxhold = 10000;
    std::vector<FoundVectorInfo> data;
    int foundnum;
    long long int totalnodes;

    long long int insidenodes;  //only used when #define depthcost is on

#ifdef depthcost
    std::vector<long int> depthcostarray;    
#endif
    
    bkzfloat expectednodes;
    double etime;
    
    void clearmemory() {
        data.clear();
        data.resize(maxhold);
        foundnum = 0;
        totalnodes = 0;
    }
    void clearmemory(int hold) {
        maxhold = hold;
        clearmemory();
    }
    
    int countnonzero() {
        int j = 0;
        for (int i=0;i<data.size();i++) {
            if (data[i].norm>=0) j++;
        }        
        return j;
    }

    int getminproj() {
        int j = -1;
        for (int i=0;i<data.size();i++) {
            if (data[i].norm>=0) {
                if (j==-1) {
                    j = data[i].projection;
                }
                j = min(j,data[i].projection);
            }
        }        
        return j;
    }
    
    int shiftproj() {
        int j = getminproj();
        for (int i=0;i<data.size();i++) {
            if (data[i].norm>=0) {
                data[i].projection -=j;
            }
        }
        return j;
    }

    int getemptyslot() {
        for (int i=0;i<data.size();i++) {
            if (data[i].norm==-1) return i;
        }
        return -1;
    }

    int getminslot() {
        int ms = 0;
        for (int i=0;i<data.size();i++) {
            if (data[i].norm==-1) return i;
            if (data[ms] < data[i]) ms = i;
        }
        return ms;
    }
    
    int getminnormslot() {
        int ms = -1;
        for (int i=0;i<data.size();i++) {
            if (data[i].norm!=-1) {
                if (ms==-1) ms = i;
                if (data[i] < data[ms]) ms = i;
            } 
        }
        return ms;
    }

    void storeelement(FoundVectorInfo& newFV,int undupflag = 0) {

        int ms = getemptyslot();
        if (ms<0) {
            ms = getminslot();
        }
        

        if (( data[ms].norm==-1 ) ||
             (data[ms].projection > newFV.projection) ||
             ((data[ms].projection == newFV.projection) && ( data[ms].norm > newFV.norm )) ) {
                //replace
            data[ms].norm = newFV.norm;
            data[ms].projection = newFV.projection;
            data[ms].coeffs.resize(newFV.coeffs.size());
            for (int i=0;i<newFV.coeffs.size();i++) data[ms].coeffs[i] = newFV.coeffs[i];
        }
        if (undupflag == flag_unduplicate) {
            duplicate();
        }
    }
    
    void display() {
        for (int i=0;i<data.size();i++) {
            if (data[i].norm>=0) {
                ccout << "slot=" << i << " proj=" << data[i].projection << " norm=" << data[i].norm << "(" << sqrt(data[i].norm) << "): ";
                for (int j=0;j<data[i].coeffs.size();j++) ccout << data[i].coeffs[j] << " ";
                ccout << endl;
            }
        }        
    }

    template <typename T> mat_ZZ get_zz_rep(LatticeBasis<T>& B,int offset=1) {
        int j = countnonzero();
        mat_ZZ ret;
        ret.SetDims(j,B.L.NumCols());
        j = 0;
        for (int i=0;i<data.size();i++) {
            if (data[i].norm >= 0) {
                coeff_to_zz(ret[j],data[i].coeffs,B.L,data[i].projection-1);
                j++;
            }
        }
        return ret;
    }
    
    void sort() {
        std::sort(data.begin(),data.end());
    }

    void simplify(FoundVectorInfo& a) {
        //this breaks norm information
        int frontzero=0;
        int backzero=a.coeffs.size()-1;

        while ((backzero > 1) && (a.coeffs[backzero]==0)) backzero--;
        a.coeffs.resize(backzero+1);

        while ((frontzero < a.coeffs.size()) && (a.coeffs[frontzero]==0)) frontzero++;
        for (int i=0;i<a.coeffs.size()-frontzero;i++) a.coeffs[i] = a.coeffs[i+frontzero];
        a.coeffs.resize(a.coeffs.size()-frontzero);
        a.projection += frontzero;
    }

    void cutbylowerbound(bkzfloat radius) {
        for (int i=0;i<data.size();i++) {
            if (data[i].norm <= radius) {
                data[i].projection = -1;
            }
        }
    }
    void duplicate() {
        for (int i=0;i<data.size();i++) {
            if (data[i].norm >= 0) {
                simplify(data[i]);
            }
        }
        
        for (int i=0;i<data.size();i++) {
            if (data[i].norm < 0) {
                data[i].projection = -1;
            }
        }
        sort();
        for (int i=1;i<data.size();i++) {
            if ( data[i-1].projection== data[i].projection) {
                if (data[i-1].coeffs.size()==data[i].coeffs.size()) {
                    char eqflag=1;
                    for (int j=0;j<data[i].coeffs.size();j++) {
                        if (data[i].coeffs[j] != data[i-1].coeffs[j] ) {
                            eqflag = 0;
                            break;
                        }
                    }
                    if (eqflag == 1) {
                        clear(data[i-1]);
                    }
                }
            }
        }
    }
    void erasesingle() {
        for (int i=0;i<data.size();i++) {
            if (data[i].norm >= 0) {
                int nonzerocount=0;
                for (int j=0;j<data[i].coeffs.size();j++) {
                    if (data[i].coeffs[j] != 0) nonzerocount++;
                }
                if (nonzerocount==1)  data[i].norm = -1;
            }
        }
    }
    
    //Experimental 
    int* coeff_lower = NULL;
    int* coeff_upper = NULL;
};


#include "vectorenumeration_boost.cpp"
#include "pfuncwrapper.cpp"

template <typename T> void GSrepresentation(std::vector<T>& gsrept,vec_ZZ& vv,mat_ZZ& L,std::vector<std::vector<T> >& mu,std::vector<T>& cd) {
    int n = L.NumCols();
    if (vv.length() != n) {
        ccout << "dimension mismatch at GSrepresentation." << endl;
        return;
    }
    int m = L.NumRows();
    gsrept.resize(m+1);

    
    //vv[0...n-1]: input vector
    //L[0...n-1][0..m-1]: lattice
    //mu[1...n][1..n]: GS coeff
    //cd[1...n]: |b*i|^2
    //compute gsrept[i]=<v,b*i> 
    //This is equal to: decomposing v=\sum alpha[i].b*[i] <=> dip[i]=alpha[i].|b*i|^2
    int i,k;
    for (i=1;i<=m;i++) {
        gsrept[i] = 0;
        for (k=0;k<n;k++) gsrept[i] += boost::lexical_cast<T>(vv[k] * L[i-1][k]);        //<d[j].b[i]>
        for (k=1;k<i;k++) gsrept[i] -= (T)(mu[i][k]) * gsrept[k];     //<d[j],b*[i]>
    }
    for (i=1;i<=m;i++) gsrept[i] /= cd[i];
}

template <typename T> void Coeffrepresentation(std::vector<T>& crept,vec_ZZ& vv,mat_ZZ& L,std::vector<std::vector<T> >& mu,std::vector<T>& cd) {
    GSrepresentation(crept,vv,L,mu,cd);
    int m = crept.size()-1;
    for (int i=m;i>=1;i--) {
        crept[i] = round(crept[i]);
        for (int j=i-1;j>=1;j--) {
            crept[j] = crept[j] - crept[i] * mu[i][j];
        }
    }
}

template <typename T> void GSrepresentation(std::vector<T>& gsrept,vec_ZZ& vv,LatticeBasis<T>& B) {
    B.updateGSBasis();
    GSrepresentation(gsrept,vv,B.L,B.gs.mu,B.gs.cd);
}

template <typename T> void Coeffrepresentation(std::vector<T>& gsrept,vec_ZZ& vv,LatticeBasis<T>& B) {
    B.updateGSBasis();
    Coeffrepresentation(gsrept,vv,B.L,B.gs.mu,B.gs.cd);
}

FoundENUMVectors currentEV; 
template <typename T> mat_ZZ ENUMbackend(LatticeBasis<T>& B,vec_ZZ& v,PruningFunction& PF,bkzfloat radius,int enumoption,int vl=0,std::string stringoptions="",int modeoption=modesvp) {
    //vec_ZZ& v is used for cvp subroutine
    currentEV.clearmemory();
    int slide = 0; 
    
    int elim = 10000;
    //elim = 1;
    
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

    bkzfloat lowerradius=0;
    if (options["lowerbound"]!="") {
        lowerradius = boost::lexical_cast<bkzfloat>(options["lowerbound"]);
    }

    if (options["vl"]!="") {
        vl = atoi(options["vl"].c_str());
    }
    
    int parallel = 1;
    if (options["parallel"]!="") {
        parallel = atoi(options["parallel"].c_str());
    }

    bool sizereduceflag = false;
    if (options["sizereduce"]=="true") {
        sizereduceflag=true;
    }
    
    B.updateGSBasis();

    int finishmode = finish_mode_nowaste;
    if (options["finishmode"]=="exact") {
        finishmode = finish_mode_exact;
    }
    if ((enumoption & enum_mode_all_vectors) != 0) {
        finishmode = finish_mode_exact;
    }
    //end of parameter 

    
    if (modeoption==modesvp) {
        if (parallel<=1) {
            if (iend-istart+1 < 150) {
                //lattice_enum::ENUMCoreboost<double,double,T>(EV,slide,istart,iend-istart+1,elim,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption);
#ifndef __lite__compile
                lattice_enum::ENUMCoreboost_double<T>(currentEV,slide,istart,iend-istart+1,elim,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption);
#else
                lattice_enum::ENUMCoreboost<double,double,T>(currentEV,slide,istart,iend-istart+1,elim,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption);
#endif                
            } else {
                lattice_enum::ENUMCoreboost<long double,long double,T>(currentEV,slide,istart,iend-istart+1,elim,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption);
            }
        } else {
            if (iend-istart+1 < 150) {
#ifndef __lite__compile
                lattice_enum::ENUMCoreboostmt_double<T>(currentEV,slide,istart,iend-istart+1,elim,parallel,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption,finishmode);
#else
    #ifndef __no__multithreads
                lattice_enum::ENUMCoreboostmt<long double,long double,T>(currentEV,slide,istart,iend-istart+1,elim,parallel,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption,finishmode);
    #else
                lattice_enum::ENUMCoreboost<long double,long double,T>(currentEV,slide,istart,iend-istart+1,elim,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption);
    #endif
#endif
            } else {
#ifndef __no__multithreads
                lattice_enum::ENUMCoreboostmt<long double,long double,T>(currentEV,slide,istart,iend-istart+1,elim,parallel,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption,finishmode);
#else 
                lattice_enum::ENUMCoreboost<long double,long double,T>(currentEV,slide,istart,iend-istart+1,elim,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption);
#endif                
            }
        }
    }     

    if (modeoption==modecvp) {
        std::vector<T> gsrept;
        GSrepresentation(gsrept,v,B);

#ifdef __lite__compile
        ccout << "error: cvp_emuerations are not defined. Please remove __lite__compile option" << endl;
        exit(0);
#else 
        if (parallel<=1) {
            if (iend-istart+1 < 150) {
                lattice_enum::ENUMCVCoreboost_double<T>(currentEV,slide,istart,iend-istart+1,elim,radius*radius,PF.rd,B.gs.mu,B.gs.cd,gsrept,vl,enumoption);
            } else {
                lattice_enum::ENUMCVCoreboost<long double,long double,T>(currentEV,slide,istart,iend-istart+1,elim,radius*radius,PF.rd,B.gs.mu,B.gs.cd,gsrept,vl,enumoption);
            }
        } else {
            if (iend-istart+1 < 150) {
                lattice_enum::ENUMCVCoreboostmt_double<T>(currentEV,slide,istart,iend-istart+1,elim,parallel,radius*radius,PF.rd,B.gs.mu,B.gs.cd,gsrept,vl,enumoption,finishmode);
            } else {
                lattice_enum::ENUMCVCoreboostmt<long double,long double,T>(currentEV,slide,istart,iend-istart+1,elim,parallel,radius*radius,PF.rd,B.gs.mu,B.gs.cd,gsrept,vl,enumoption,finishmode);
            }
        }
#endif
    }     

    if (vl>=3) currentEV.display();
    if (lowerradius>0) currentEV.cutbylowerbound(lowerradius);

    mat_ZZ ret;
    ret = currentEV.get_zz_rep(B);
    
    if ((istart>1) && (sizereduceflag==true)) {
        //Size-reduction
        LatticeBasis<T> B2;
        B2.resize(istart,B.L.NumCols());
        for (int i=0;i<istart-1;i++) {
            B2[i] = B[i];
        }
        B2.updateGSBasis(1,istart-1);
        for (int i=0;i<ret.NumRows();i++) {
            B2[istart-1] = ret[i];
            B2.GScomputed = istart-1;
            B2.updateGSBasis(1,istart);
            //B2.gs.displaymu();
            SizeReduce(B2,0,istart);
            ret[i] = B2[istart-1];
        }
    }
    return ret;
}

template <typename T> mat_ZZ ENUM(LatticeBasis<T>& B,bkzfloat radius,bkzfloat tprob,int enumoption=0,int pfoption=0,int vl=0,std::string otheroptions="") {
    PruningFunction PF;
    bkzfloat prob = pruning_func::SetPruningFunction(B,PF, radius , tprob,pfoption,vl);
    vec_ZZ zerovec;
    B.updateGSBasis();
    return ENUMbackend<T>(B,zerovec,PF,radius,enumoption,vl,otheroptions,modesvp);
}

template <typename T> mat_ZZ ENUMCV(LatticeBasis<T>& B,vec_ZZ& targetv,bkzfloat radius,bkzfloat tprob,int enumoption=0,int pfoption=0,int vl=0,std::string otheroptions="") {
    PruningFunction PF;
    bkzfloat prob = pruning_func::SetPruningFunction(B,PF, radius , tprob,pfoption);
    B.updateGSBasis();
    return ENUMbackend<T>(B,targetv,PF,radius,enumoption,vl,otheroptions,modecvp);
}

#endif