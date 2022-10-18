#ifndef _inc_randomize_wrapper_cpp
#define _inc_randomize_wrapper_cpp


template <typename T> void Randomize(LatticeBasis<T>& B,int seed=0,int bit=10,int istart=1,int iend=0) {
    //Randomize (projective) lattices
    //Randomize B[istart..iend] by multiplying random unimodular matrix
    //Bitsize of unimodular matrix is defined by bit

    mat_ZZ U;
    if (iend==0) iend = B.dim;

    if (iend-istart==B.dim-1) {
        gen_random_unimodular2(U,B.dim,seed,bit,VL0);
    } else {
        mat_ZZ U2;
        gen_random_unimodular2(U2,iend-istart+1,seed,bit,VL0);
        U.SetDims(B.dim,B.dim);
        for (int i=0;i<B.dim;i++) {
            U[i][i] = 1;
        }
        for (int i=0;i<iend-istart+1;i++) {
            for (int j=0;j<iend-istart+1;j++) {
                U[istart-1+i][istart-1+j] = U2[i][j];
            }
        }
    }
    B.L = U * B.L;
    B.isNecessaryGSUpdate = true;
    B.updateGSBasis();
}

template <typename T> vec_ZZ RandomPoint(LatticeBasis<T>& B,int seed,int dim=0) {
    //Generate uniformly random point in the box \prod[-b*i/2,b*i/2]
    B.updateGSBasis();
    vec_ZZ ret;

    ZZ det = LatticeVolumeZZ(B);
    if (dim==0) dim = B.dim;
    ret.SetLength(dim);
    SetSeed(to_ZZ(seed));
    for (int i=0;i<dim;i++) {
        ret[i] = RandomBnd(det);
    }
    ret = NearestPlane(ret,B);
    
    return ret;
    
}

template <typename T> mat_ZZ RandomPoints(LatticeBasis<T>& B,int seed,int numrows,int dim=0) {

    B.updateGSBasis();
    vec_ZZ r;
    ZZ det = LatticeVolumeZZ(B);
    if (dim==0) dim = B.dim;
    r.SetLength(dim);

    mat_ZZ ret;
    ret.SetDims(numrows,dim);
    
    SetSeed(to_ZZ(seed));
    
    for (int k=0;k<numrows;k++) {
        for (int i=0;i<dim;i++) {
            r[i] = RandomBnd(det);
        }
        ret[k] = NearestPlane(r,B);
    }
    return ret;
}

#endif