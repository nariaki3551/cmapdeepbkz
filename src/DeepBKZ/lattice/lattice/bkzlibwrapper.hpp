#ifndef _inc_bkzlibwrapper
#define _inc_bkzlibwrapper


#if defined(_include_svpchallenge_generator) || defined(_include_idealchallenge_generator)
#include "latticetools.cpp"
#include "memorytools.cpp"
#endif

//Structure for approximated Gram-Schmidt lengths and coefficients
template <typename T> struct GSapprox {
    std::vector<T> c;       //c[1..n] = |b*i|
    std::vector<T> cd;      //cd[1..n] = |b*i|^2
    std::vector<T> len;      //len[1..n] = |b_i|^2
    std::vector<std::vector<T> >mu;      //mu[1..n][1..n]
    std::vector<std::vector<T> >r;      //r[1..n][1..n]: r[i,j] = mu[i,j] * cd[j]
    
    void clear() {
        c.resize(0);
        cd.resize(0);
        len.resize(0);
        mu.resize(0);
        r.resize(0);
    }

    GSapprox() {
        clear();
    }
    
    template <typename T2> GSapprox& operator =(GSapprox<T2>& newc) {
        vec_copy(c,newc.c);
        vec_copy(cd,newc.cd);
        vec_copy(len,newc.len);
        vec_copy(mu,newc.mu);
        vec_copy(r,newc.r);
        return *this;
    }

    T minbi() {
        if (c.size()<=1) return 0;
        T ret = c[1];
        for (int i=1;i<c.size();i++) {
            ret = min(ret,c[i]);
        }
        return ret;
    }
    
    T maxbi() {
        if (c.size()<=1) return 0;
        T ret = c[1];
        for (int i=1;i<c.size();i++) {
            ret = max(ret,c[i]);
        }
        return ret;
    }

    void gslentofile(std::string filename,ios::openmode option=ios::trunc) {
        std::ofstream of;
        of.open(filename.c_str(),option);
        for (int i=1;i<max(c.size(),cd.size());i++) {
            of << i << "\t";
            if (i<c.size()) of << c[i];
            of << "\t";
            if (i<cd.size()) of << cd[i];
            of << endl;
        }
        of << endl;
        of.close();
    }

    void displaymu() {
        for (int i=1;i<mu.size();i++) {
            for (int j=1;j<mu[i].size();j++) {
                cout << mu[i][j] << " ";
            }
            cout << endl;
        }
    }

    void displaymuline(int i) {
        for (int j=1;j<=i;j++) {
            cout << mu[i][j] << " ";
        }
        cout << endl;
    }

    void displaymu(int i,int j) {
        for (int k=i;k<=j;k++) displaymuline(k);
    }
    void displayc(int istart=1,int iend=-1) {
        cout << "[";
        if (iend==-1) iend = c.size()-1;
        
        for (int i=istart;i<=iend;i++) {
            cout << c[i] << " ";
        }
        cout << "]" << endl;
    }
    void displaycd(int istart=1,int iend=-1) {
        cout << "[";
        if (iend==-1) iend = cd.size()-1;
        
        for (int i=istart;i<=iend;i++) {
            cout << cd[i] << " ";
        }
        cout << "]" << endl;
    }
    bool zerocheck(int istart=1,int iend=-1) {
        for (int i=istart;i<=iend;i++) {
            if (c[i] == 0) {
                cout << "zerocheck: " << i << endl;
                displayc(istart,iend);
                return false;
            }
        }
        return true;
    }
    
};

template<typename T> bool operator <(GSapprox<T>& gs1,GSapprox<T>& gs2) {
    int n = min(gs1.c.size(),gs2.c.size());
    for (int i=1;i<n;i++) {
        if (gs1.c[i] < gs2.c[i]) return true;
        if (gs1.c[i] > gs2.c[i]) return false;
    }
    return false;
}

template <typename T> void computeapproxline(std::vector<T>& approxelement,vec_ZZ& v) {
    int len = v.length();
    approxelement.resize(len+1);
    for (int i=0;i<len;i++) {
        if (abs(v[i]) <= 1e+8) {
            approxelement[i+1] = to_double(v[i]);   //NTL_ZZ to T
        } else {
            conv(approxelement[i+1] ,v[i]);
        }
    }
}

template <typename T,typename VEC> void computeapproxline(std::vector<T>& approxelement,VEC& v) {
    int len = v.length();
    approxelement.resize(len+1);
    for (int i=0;i<len;i++) {
        if (abs(v[i]) <= 1e+7) {
            approxelement[i+1] = v[i];
        } else {
            approxelement[i+1] = boost::lexical_cast<T>(v[i]);
        }
    }
}

template <typename T> void computeGramSchmidtline(mat_ZZ& L,GSapprox<T>& gs,std::vector<std::vector<T> >& approxelement,int index,int ioffset=1,int joffset=1) {

    computeapproxline(approxelement[index],L[index-1]);
    T ipapprox;
    for (int j=1;j<=index;j++) {
        ipapprox = 0;
        int size = approxelement[index].size();
        for (int i=1;i<size;i++) {
            ipapprox += approxelement[index][i] * approxelement[j][i];
        }
        gs.r[index][j] = ipapprox;
        if (j==index) gs.len[index] = gs.r[index][j];

        for (int k=1;k<j;k++) {
            gs.r[index][j] -=  gs.r[index][k] * gs.mu[j][k];
        }
        if (j<index) {
            if (gs.cd[j] > 0) {
                gs.mu[index][j] =  gs.r[index][j] / gs.cd[j];
            } else {
                gs.mu[index][j] =0;
            }
        }
    }
    gs.cd[index] = gs.r[index][index];
    gs.mu[index][index] = 1.0;
    
    if (gs.cd[index] > 0) {
        gs.c[index] = sqrt(gs.cd[index]);
    } else {
        gs.c[index] = 0;
    }

    
/*
    //old code
    ZZ ip;
    computeapproxline(approxelement[index],L[index-1]);
    T ipapprox;
    for (int j=1;j<=index;j++) {
        ipapprox = 0;
        int size = approxelement[index].size();
        for (int i=1;i<size;i++) {
            ipapprox += approxelement[index][i] * approxelement[j][i];
        }
        gs.mu[index][j] = ipapprox;
        if (j==index) gs.len[index] = gs.mu[index][j];

        for (int k=1;k<j;k++) {
            gs.mu[index][j] -= gs.mu[j][k] * gs.mu[index][k] * gs.cd[k];
        }
        if (j<index) {
            if (gs.mu[j][j] > 0) {
                gs.mu[index][j] /=  gs.mu[j][j];
            } else {
                gs.mu[index][j] =0;
            }
        }
    }
    gs.cd[index] = gs.mu[index][index];

    if (gs.cd[index] > 0) {
        gs.c[index] = sqrt(gs.cd[index]);
    } else {
        gs.c[index] = 0;
    }
 */ 
}

template <typename GFLOAT,typename MAT,typename VEC> struct LatticeBasisBase {

    MAT L;
    GSapprox<GFLOAT> gs;
    std::vector<std::vector<GFLOAT> > approxelement;      //approximation of L[i-1][j-1] (note: index is shifted)

    int alloced_rows=0;
    int alloced_cols=0;

    int dim;
    int GScomputed; //GS basis is computed for L[1..GScomputed]
    bool isNecessaryGSUpdate;
    
    LatticeBasisBase() {
        isNecessaryGSUpdate=true;        
        GScomputed = 0;
        dim=0;
    }

    LatticeBasisBase& operator =(const MAT& newL) {
        L = newL;
        dim = newL.NumRows();
        gs.clear();
        alloced_rows=0;
        alloced_cols=0;
        isNecessaryGSUpdate=true;
        GScomputed = 0;
        return *this;
    }
    
    LatticeBasisBase& operator =(LatticeBasisBase& newL) {
        L = newL.L;
        dim = newL.dim;
        gs = newL.gs;
        isNecessaryGSUpdate=true;        
        GScomputed = 0;
        return *this;
    }

    template <typename T2> LatticeBasisBase& operator =(LatticeBasisBase<T2,MAT,VEC>& newL) {
        L = newL.L;
        dim = newL.dim;
        gs = newL.gs;
        isNecessaryGSUpdate=true;        
        GScomputed = 0;
        return *this;
    }
    
    VEC& operator [] (int i) {
        if ((0<=i) && (i<L.NumRows())) {
            isNecessaryGSUpdate=true;        
            return L[i];
        } else {
            cout << "range_error: " << i << " size=" << dim << endl;
            exit(0);
        }
    }

    
    void memoryalloc() {

        if ((alloced_rows != dim+1) || (alloced_cols != L.NumCols()+1)) {
            gs.c.resize(dim+1);
            gs.cd.resize(dim+1);
            gs.len.resize(dim+1);
            if (dim+1 != gs.mu.size()) {
                gs.mu.resize(dim+1);
                for (int i=1;i<=dim;i++) gs.mu[i].resize(dim+1);
                gs.r.resize(dim+1);
                for (int i=1;i<=dim;i++) gs.r[i].resize(dim+1);
            }

            approxelement.resize(dim+1);
            for (int i=1;i<=dim;i++) {
                approxelement[i].resize(L.NumCols()+1);
            }
            alloced_rows = dim+1;
            alloced_cols = L.NumCols()+1;
        }        
    }
    
    void updateGSBasis(int istart,int iend) {
        dim = L.NumRows();
        memoryalloc();
        istart = max(istart,GScomputed+1);
        
        if (istart < 1) istart = 1;
        if (iend > dim) iend = dim;
        for (int i=istart;i<=iend;i++) {
            computeGramSchmidtline(L,gs,approxelement,i);
        }
        if (istart <= iend ) GScomputed = iend;
    }

    void updateGSBasis_noapprox(int istart,int iend) {
        dim = L.NumRows();
        memoryalloc();
        istart = max(istart,GScomputed+1);
        if (istart < 1) istart = 1;
        if (iend > dim) iend = dim;
        for (int i=istart;i<=iend;i++) {
            computeGramSchmidtline(L,gs,i);
        }
        if (istart <= iend ) GScomputed = iend;
    }

    inline void updateGSBasis() {
        int n = L.NumRows();
        GScomputed = 0;
        updateGSBasis(1,n);
    }

    inline void updateGSBasis_noapprox() {
        int n = L.NumRows();
        GScomputed = 0;
        updateGSBasis_noapprox(1,n);
    }

    void resize(int n,int m) {
        L.SetDims(n,m);
        dim = n;    //number of rows
        isNecessaryGSUpdate=true;        
    }
    
    GFLOAT projlen(int i,int j) {
        //projective length pi_i(b_j)
        updateGSBasis(i,j);
        GFLOAT ret = 0;
        for (int k=i;k<j;k++) ret += gs.mu[j][k]*gs.mu[j][k] * gs.cd[k];
        ret += gs.cd[j];
        return sqrt(ret);
    }
    GFLOAT projlen_noapprox(int i,int j) {
        //projective length pi_i(b_j)
        updateGSBasis_noapprox(i,j);
        GFLOAT ret = 0;
        for (int k=i;k<j;k++) ret += gs.mu[j][k]*gs.mu[j][k] * gs.cd[k];
        ret += gs.cd[j];
        return sqrt(ret);
    }
    
    GFLOAT loghalfvolume() {
        GFLOAT ret = 0;
        int n = L.NumRows();
        for (int i=1;i<=n/2;i++) {
            ret += log(gs.c[i]);
        }
        return ret;
    }
    
    GFLOAT energy() {
        GFLOAT ret = 0;
        int n = L.NumRows();
        for (int i=1;i<=n;i++) {
            ret += gs.cd[i];
        }
        return ret;
    }
    
    GFLOAT loglllpotential() {
        GFLOAT ret = 0;
        int n = L.NumRows();
        for (int i=1;i<=n;i++) {
            ret += log(gs.c[i]) * (n-i+1);
        }
        return ret;
    } 

    void displayapproxline(int i) {
        cout << "approx[" << i << "]: "; 
        for (int j=1;j<approxelement[i].size();j++) {
            cout << approxelement[i][j] << " ";
        }
        cout << endl;
    }

    void displayapprox(int i,int j) {
        for (int k=i;k<=j;k++) displayapproxline(k);
    }
    
    void updateapprox(int i,int j) {
        for (int k=i;k<=j;k++) {
            computeapproxline(approxelement[k],L[k-1]);
        }
    }
    
};

#include "pbkzmatrix.cpp"

template <typename GFLOAT> 
using LatticeBasis = LatticeBasisBase<GFLOAT,mat_ZZ,vec_ZZ>;

//T2=type for element T=type for GS basis
template <typename T2,typename T> using SmallLatticeBasis = LatticeBasisBase<T,pbkzmatrix<T2>,pbkzvector<T2>>;

template <typename T> std::ostream& operator <<(std::ostream &stdout, const LatticeBasis<T>  &arg) {
    stdout << arg.L;
    return stdout;
}

template <typename T> void LLL(LatticeBasis<T>& B,double delta=0.99) {
    int n = B.dim;
    if (n<=160) {
        LLL_XD(B.L,min(0.9,delta),0,0,1);
    } else {
        LLL_RR(B.L,min(0.9,delta),0,0,1);
    }
    LLL_QP(B.L,delta,0,0,1);
    B.isNecessaryGSUpdate = true;
    B.updateGSBasis();
}

template <typename GFLOAT,typename MAT,typename VEC> void SaveLattice(LatticeBasisBase<GFLOAT,MAT,VEC>& B,std::string fname) {
    SaveLattice(B.L,fname);
}

template <typename GFLOAT,typename MAT,typename VEC>  void LoadLattice(LatticeBasisBase<GFLOAT,MAT,VEC>& B,std::string fname) {
    LoadLattice(B.L,fname);
    B.isNecessaryGSUpdate = true;
    B.updateGSBasis();
}


template <typename T> T LatticeGH(LatticeBasis<T>& B) {
    B.updateGSBasis();
    lattice_tools::initialize();
    T det=0;
    for (int i=1;i<=B.dim;i++) {
        det += log(B.gs.c[i]);
    }
    det = exp(det/B.dim) * (T)bkzconstants::ghconstant(B.dim);
    return det;
}

template <typename T> T LogLatticePotential(LatticeBasis<T>& B) {
    //LLLPoteitial
    B.updateGSBasis();
    lattice_tools::initialize();
    T det=0;
    for (int i=1;i<=B.dim;i++) {
        det += log(B.gs.c[i]) * (B.dim-i+1);
    }
    return det;
}

template <typename T> T LatticeEnergy(LatticeBasis<T>& L,int istart=1,int iend=-1) {
    L.updateGSBasis();
    if (iend==-1) iend = L.dim;
    lattice_tools::initialize();
    T energy=0;
    for (int i=istart;i<=iend;i++) energy += L.gs.cd[i];
    return energy;
}

template <typename T> T LatticeEnergy(std::vector<T>& c,int istart=1,int iend=-1) {
    if (iend==-1) iend = c.size()-1;
    lattice_tools::initialize();
    T energy=0;
    for (int i=istart;i<=iend;i++) energy += c[i];
    return energy;
}

template <typename T> bkzfloat LatticeVolume(std::vector<T>& c,int istart,int iend,char opt=INPUT_NONSQUARED) { 
    lattice_tools::initialize();
    bkzfloat ret = 1;
    for (int i=istart;i<=iend;i++) {
        ret *= c[i];
    }
    if (opt==INPUT_SQUARED) return sqrt(ret);
    return ret;
}

template <typename T> bkzfloat LatticeVolume(LatticeBasis<T>& B) { 
    lattice_tools::initialize();
    B.updateGSBasis();
    return LatticeVolume(B.gs.c,1,B.dim);
}
    
template <typename T> T RootHermiteFactor(LatticeBasis<T>& B) {
    B.updateGSBasis();
    lattice_tools::initialize();
    T rhf=0;
    for (int i=1;i<=B.dim;i++) rhf += log(B.gs.c[i]);   //determinant
    rhf = B.gs.c[1] / exp(rhf/B.dim);   
    return pow(rhf,1.0/B.dim);
}


//Macro for multi process subroutine

#define __premproc(n) pid_t pid = 1;\
    int numactiveproc=0;\
    int maxproc = n;\

#define __singlecorework maxproc = 1;pid=0;

#define __startmproc(iter_val) for (int iter_val=0;iter_val<maxproc;iter_val++) {\
    if (pid!=0) {\
        pid = fork();\
    }\
    if (pid!=0) {\
        numactiveproc++;\
        if (numactiveproc>=maxproc) {\
             int status;\
             pid = wait(&status);\
             if (ECHILD==errno) {\
                break;\
            }\
            numactiveproc--;\
        }\
    } else {\

#define __endmproc _exit(0);}}

template <typename T> bkzfloat LatticeVolumeRoot(std::vector<T>& c,int istart,int iend,int opt=INPUT_NONSQUARED) {
    lattice_tools::initialize();
    int i;
    bkzfloat det = 0;
    int dim = iend-istart+1;
    for (int i=istart;i<=iend;i++) {
        det += log((bkzfloat)c[i]);
    }
    if (opt==INPUT_SQUARED) det *= 0.5;
    det = exp(det/dim);     
    return det;
}
template <typename T> bkzfloat LatticeVolumeRoot(std::vector<T>& c,int opt=INPUT_NONSQUARED) {
    lattice_tools::initialize();
    return LatticeVolumeRoot(c,1,c.size()-1,opt);
}

template <typename T> bkzfloat LatticeVolumeRoot(LatticeBasis<T>& B) { 
    lattice_tools::initialize();
    B.updateGSBasis();
    return LatticeVolumeRoot(B.gs.c,INPUT_NONSQUARED);
}

bkzfloat LatticeGH(std::vector<bkzfloat>& c,char opt=INPUT_SQUARED) {
    lattice_tools::initialize();
    int i;
    bkzfloat det = 0;
    int dim = c.size()-1;
    for (int i=1;i<=dim;i++) {
        det += log(c[i]);
    }
    if (opt==INPUT_SQUARED) det *= 0.5;
    det = exp(det/dim);     
    return det * bkzconstants::ghconstant(dim);
}

template <typename T> void display_vector(std::vector<T>& a,int istart=0,int iend=0) {
    if (iend==0) iend = a.size()-1;
    cout << "[";
    for (int i=istart;i<=iend;i++) cout << a[i] << " ";
    cout << "]" << endl;
}

template <typename T> void output_vector(std::vector<T>& a,int istart=0,int iend=0,std::string filename="") {
    if (iend==0) iend = a.size()-1;
    std::ofstream of;
    of.open(filename.c_str(),ios::trunc);
    for (int i=istart;i<=iend;i++) {
        of << i << "\t" << a[i] << endl;
    }
    of.close();
}


 template <typename T> T ProjectedLength(LatticeBasis<T>& B,int istart,int iend, vec_ZZ& vector) {
     
    B.updateGSBasis(1,iend);
    ZZ ip;
    T ipapprox;
    T* mu = (T*)shared_memory::allocate1<T>(123,iend+1);
    T result = LengthOf<T>(vector);
    result *= result;

    for (int j=istart;j<=iend;j++) {
        ipapprox = 0;
        for (int i=0;i<vector.length();i++) {
            //inner product of two vectors
            if ((vector[i]!=0) && (B.L[j-1][i]!=0)) { 
                T part;
                part = log(to_double( log(abs(vector[i] )))) + log(to_double( log(abs( B.L[j-1][i] ))));
                char sign1=1;
                char sign2=1;
                if (vector[i] < 0) sign1 = -1;
                if (B.L[j-1][i]  < 0) sign2 = -1;
                if (sign1*sign2==1) {
                    ipapprox += exp(part);
                } else {
                    ipapprox -= exp(part);
                }
            }
        }
        mu[j] = ipapprox;

        for (int k=1;k<j;k++) {
            mu[j] -= B.gs.mu[j][k] * mu[k] * B.gs.cd[k];
        }
        mu[j] /=  B.gs.cd[j];
        result -= mu[j] * mu[j] * B.gs.cd[j];
        if (result < 0.5) return 0;
    }
    return sqrt(result);
 }

using namespace boost::multiprecision;


#ifdef _use_mpfr_pfunc

template <typename T> ZZ LatticeVolumeZZ(LatticeBasis<T>& B) {
 
     //find exact determinait in integer
    T approxdet = LatticeVolume(B);

    double logdet = boost::lexical_cast<double>(log(approxdet));
    int prec = logdet / 2.302585092;
    prec += 30;
    
    int precisionback = mpfr_float::default_precision() ;
    mpfr_float::default_precision(prec);

    LatticeBasis<mpfr_float> Bm;
    Bm = B;
    Bm.updateGSBasis();
    
    mpfr_float detmp = 1;
    for (int i=1;i<=Bm.dim;i++) detmp *= Bm.gs.c[i];
    
    detmp = round(detmp);
    
    std::string t = detmp.convert_to<std::string>();
    ZZ ret = to_ZZ(t.c_str());

    mpfr_float::default_precision(precisionback);
    return ret;
 }
#endif

template <typename T>void UpdateProjectiveLength(std::vector<std::vector<T > >& p,GSapprox<T>& gs,int istart,int iend) {
    //p[j][i] = pi_j(b_i)
    
    for (int i=istart;i<=iend;i++) {
        T subproj=0;
        for (int j=i;j>=istart;j--) {
            if (j==i) {
                subproj += gs.cd[j];
            } else {
                subproj += gs.cd[j] * gs.mu[i][j] * gs.mu[i][j];
            }
            p[j][i] = subproj;
        }
    }
}    

template <typename T> void Copysublattice(LatticeBasis<T>& dstB,LatticeBasis<T>& srcB,int istart,int iend) {
    if (istart < 1) istart = 1;
    if (iend > srcB.dim) iend = srcB.dim;
    int col = srcB.L.NumCols();
    int idim = iend-istart+1;
    dstB.L.SetDims(idim,col);
    for (int i=istart;i<=iend;i++) {
        dstB.L[i-istart] = srcB.L[i-1];
    }
}



#endif