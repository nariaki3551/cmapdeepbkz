#ifndef _inc_pbkz_small_lattice_cpp
#define _inc_pbkz_small_lattice_cpp

#include "enumwrapper.cpp"


//Gram Schmidt subroutines for lattices with small elements

template <typename T> std::ostream& operator <<(std::ostream &stdout, pbkzmatrix<T>  &arg) {
    stdout << "[";
    for (int i=0;i<arg.size();i++) {
        stdout << arg[i] << endl;
    }
    stdout << "]";
    return stdout;
}

//double=type to compute GS basis, T=type to store elements (assumed to be int/float/double)
template <typename T,typename T2> inline void computeGramSchmidtline(pbkzmatrix<T>& L,GSapprox<T2>& gs,int index,int ioffset=1,int joffset=1) {

    T ipapprox;
    int len = L[index-1].length();

    for (int j=1;j<=index;j++) {
        ipapprox = 0;
        for (int i=0;i<len;i++) {
            ipapprox += L[index-1][i] * L[j-1][i];
        }
        gs.mu[index][j] = ipapprox;
        if (j==index) gs.len[index] = gs.mu[index][j];

        for (int k=1;k<j;k++) {
            gs.mu[index][j] -= gs.mu[j][k] * gs.mu[index][k] * gs.cd[k];
        }
        if (j<index) gs.mu[index][j] /=  gs.mu[j][j];
    }
    gs.cd[index] = gs.mu[index][index];
    gs.c[index] = sqrt(gs.cd[index]);
}




 template <typename T2,typename T> T ProjectedLength(SmallLatticeBasis<T2,T>& B,int istart,int iend, pbkzvector<T>& vector) {
    //T2=type for element T=type for GS basis
    B.updateGSBasis_noapprox(1,iend);
    T ip;
    T* mu = (T*)shared_memory::allocate1<T>(123,iend+1);
    T result = squared_norm(vector);

    for (int j=istart;j<=iend;j++) {
        ip = 0;
        for (int i=0;i<vector.length();i++) {
            ip += vector[i] * B.L[j-1][i];
        }
        mu[j] = ip;

        for (int k=1;k<j;k++) {
            mu[j] -= B.gs.mu[j][k] * mu[k] * B.gs.cd[k];
        }
        mu[j] /=  B.gs.cd[j];
        result -= mu[j] * mu[j] * B.gs.cd[j];
        if (result < 0.5) return 0;
    }
    return sqrt(result);
 }
 
template <typename CFLOAT,typename T> int BasisUpdateGSsmall(LatticeBasis<T>& B,mat_ZZ*U,int istart,int iend,FoundENUMVectors& EV,int vl,int option,char memoryflag=0) {

    //Insert vectors after B[istart]

    EV.duplicate();
    EV.erasesingle();
    
    int nonzero = EV.countnonzero();
    if (nonzero==0) return 0;

    double ss = gettimeofday_sec();

    istart = EV.shiftproj();
    
    //Allocate mu-buffer 
    int cols = iend-istart+1;
    int rows = cols + nonzero;

    SmallLatticeBasis<CFLOAT,CFLOAT> LB;    //local basis
    //SmallLatticeBasis<CFLOAT,CFLOAT> LBinit,LBinit2;   //Buffers for debug
    pbkzmatrix<int> LU;
    LB.L.SetDims(rows,cols);
    LU.SetDims(rows,cols);

    //Find multiply factor
    T mings = B.gs.c[istart];
    for (int i=1;i<cols;i++) {
        mings = min(mings,B.gs.c[istart+i]);
    }    
    mings = mings/10000;

    //diagonal
    for (int i=0;i<cols;i++) {
        for (int j=0;j<i;j++) {
            T ll = B.gs.mu[i+istart][j+istart] * B.gs.c[j+istart]/mings;
            LB.L[i][j] = (CFLOAT)ll;
        }
        T ll = B.gs.c[i+istart]/mings;
        LB.L[i][i] = (CFLOAT)ll;
    }
    for (int i=cols;i<rows;i++) {
        for (int j=0;j<cols;j++) {
            LB.L[i][j] = 0;
        }
    }
    
    //LBinit = LB;  //save basis for debug
    for (int i=0;i<rows;i++) {
        for (int j=0;j<cols;j++) {
            LU[i][j] = 0;
            if (i==j) LU[i][i] = 1;
        }
    }
    std::vector<CFLOAT> cback;
    cback.resize(cols);
    for (int i=0;i<cols;i++) cback[i] = LB.L[i][i];

    //found vectors
    int k = cols;
    int numinsert=0;
    int usedinsert = 0;
    for (int i=0;i<EV.data.size();i++) {
        if (EV.data[i].norm > 0) {
            //insert
            for (int j=0;j<EV.data[i].coeffs.size();j++) {
                addvector(LB.L[k], EV.data[i].coeffs[j] , LB.L[j + EV.data[i].projection]  );
                LU[k][j + EV.data[i].projection] = EV.data[i].coeffs[j];
            }
            k++;
        }
    }
    int numcandidates = k-cols;
    
    //Greedy reordering
    int j = 0;
    CFLOAT localdet = 0;
    CFLOAT newdet=0;
    bool firstflag = true;  //true => any found vector is not inserted
    
    while (1) {
        double minproj = -1;
        int index;
        int jj = j;
        if (firstflag==true) jj = cols;   //first vector must be found one 
        
        for (int i=jj;i<rows;i++) {
            //check all projective length
            double proj = conv_to_double(ProjectedLength(LB,1,j,LB.L[i]));
            if (proj>0.5) {
                if ((minproj == -1) || (proj < minproj)){
                    minproj = proj;
                    index = i;
                }
            }
        }
        if (minproj>0) {
            if (j>=cols) {  //for debug
                ccout << "BasisUpdateGSsmall: invalid read? j=" << j << " cols=" << cols << endl;
                ccout << "minproj=" << minproj << endl;
                exit(0);
            }

            if (firstflag==true) {
                if (cback[j] <= minproj) {      //invalid read
                    //not inserted
                    index = j;
                } else {
                    firstflag=false;
                }
            }

            int nonzero=0;
            for (int i=0;i<LU[index].length();i++) {
                if (LU[index][i]!=0) nonzero++;
            }

            if (firstflag==false) {
                localdet += log(cback[j]);      //invalid read
                newdet += log(minproj);
            }
            if (j>=1) {
                if (newdet > localdet + 0.0001) {
                    break;
                }
            }
            swap(LU[j],LU[index]);
            swap(LB.L[j],LB.L[index]);
            LB.updateGSBasis_noapprox(j,j);
            if ((usedinsert==numinsert) && (index >= cols)) usedinsert++;
            numinsert++;

            j++;
        } else {
            break;
        }
    }
    //ccout << "#insert=" << numinsert << endl;
    //ccout << "#used=" << usedinsert << endl;

    CFLOAT delta = 0.99;
    int swap;
    int final_iend = local_LLL_small(LB,&LU,delta,1,rows,stnormal,VL0,swap,(CFLOAT)1000); 
        //error check
       if (final_iend != cols) {
        ccout << std::setprecision(12);
        ccout << endl;
        ccout << "BasisUpdateGSsmall local_LLL error: fiend=" << final_iend << " " << cols << endl;
        ccout << LB.L << endl;
        LB.updateGSBasis_noapprox(1,rows);
        LB.gs.displaymu();
        LB.gs.displayc();
        exit(0);
    }
    
    
    //Reconstruct original matrix
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
    if (BUBuffer==0) {
        BUBuffer = new mat_ZZ;
       bubp[omp_get_num_threads()-1] = BUBuffer;
    }
    
    if ((cols > BUBuffer->NumRows()) || (n > BUBuffer->NumCols())) {
       BUBuffer->SetDims(cols,n);
    }

    if ( (BUBuffer->NumCols() != B.L.NumCols()) || (BUBuffer->NumRows() < cols )) {
       BUBuffer->SetDims(cols,B.L.NumCols());
    }
    for (int i=0;i<cols;i++) {
        clear((*BUBuffer)[i]);
        for (int j=0;j<cols;j++) {
            if (LU[i][j]!=0) RowTransformwrap((*BUBuffer)[i],  B.L[istart + j-1], to_ZZ(-LU[i][j]));
        }
    }
    for (int i=0;i<cols;i++) {
        B.L[istart+i-1] = (*BUBuffer)[i];
    }    

    B.GScomputed = istart-1;
    B.updateGSBasis(1,iend);
    
    if (U!=0) {
        if ( (BUBuffer->NumCols() != (*U)[0].length()) || (BUBuffer->NumRows() < cols )) {
           BUBuffer->SetDims(cols,(*U)[0].length());
        }
        for (int i=0;i<cols;i++) {
            clear((*BUBuffer)[i]);
            for (int j=0;j<cols;j++) {
                if (LU[i][j]!=0) (*BUBuffer)[i] += LU[i][j] * (*U)[istart+j-1] ; 
            }
        }
        for (int i=0;i<cols;i++) {
            (*U)[istart+i-1] = (*BUBuffer)[i];
        }    
    }
    return usedinsert;
}
      


#endif