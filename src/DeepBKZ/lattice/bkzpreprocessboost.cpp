#ifndef _inc_bkzpreprocess_boost_cpp
#define _inc_bkzpreprocess_boost_cpp

//Subroutines for pbkz preprocess

namespace bkzconstants{
    bkzfloat target_cost_lb_mult(int dim,int index) {
         if (dim>5000) {
            ccout << "target_cost_lb_mult: dimension error: " << dim << endl;
            exit(0);
        }
        if (dim < index) {
            ccout << "target_cost_lb_mult: blocksize error: " << dim << " " << index << endl;
            exit(0);
        }
        double alpha,r;
        bkzfloat prob;
        PBKZParam(alpha,r,prob,dim);

        #pragma omp critical
        {
            int cs = target_cost_lb_mult_table.size();
            if (dim >= cs) {
                target_cost_lb_mult_table.resize(dim+1);
                for (int i=cs;i<=dim;i++) {
                    target_cost_lb_mult_table[i].resize(i+1);
                    for (int j=0;j<target_cost_lb_mult_table[i].size();j++) target_cost_lb_mult_table[i][j] = 0;
                }
            }
        }
        bkzfloat value;
        #pragma omp critical
        {
            value = target_cost_lb_mult_table[dim][index];
        }
        if (value==0) {
            bkzfloat lf;
            lf = ibeta_inv(0.5*index,0.5*(dim-index),prob);
            lf = mypower(lf,0.5*index);
            target_cost_lb_mult_table[dim][index] = lf;
        }
        return target_cost_lb_mult_table[dim][index];
    }

    template <typename T> bkzfloat target_cost_lb_basis(std::vector<T>& cc,int istart,int iend) {

        if (iend-istart+1 <=10) return 0;
        
        //Assume cc[1..n]=|b*i| (non-squared)
        bkzfloat cost,lcost;
        lcost = 1.0;

        double alpha,r;
        bkzfloat prob;
        PBKZParam(alpha,r,prob,iend-istart+1);

        bkzfloat radius = alpha * lattice_tools::LatticeGH(cc,istart,iend,INPUT_NONSQUARED);
        int dim = iend-istart+1;
        for (int i=iend;i>=istart+2;i--) {
            lcost *= radius / cc[i];
            bkzfloat localcost = lcost * bkzconstants::vol_unit_ball(iend-i+1) *  target_cost_lb_mult(dim,iend-i+1);
            cost += localcost;
        }
        return 0.5*cost;        //halved by symmetry
    }
    
    bkzfloat target_cost_lb(int dim) {

        if ((dim<10) || (dim>=2000)) {
            ccout << "target_cost_lb: dimension error " << dim << endl;
            exit(0);
        }  

        #pragma omp critical
        {
            int cs = target_cost_lb_table.size();
            if (dim >= cs) {
                target_cost_lb_table.resize(dim+1);
                for (int i=cs;i<=dim;i++) target_cost_lb_table[i] = 0;
            }
        }
        
        bkzfloat value;
        #pragma omp critical
        {
            value = target_cost_lb_table[dim];
        }
        if (value==0) {
            double alpha,r;
            bkzfloat prob;
            PBKZParam(alpha,r,prob,dim);
            
            //make gsa basis
            std::vector<bkzfloat> cc;
            PBKZSimulate_first(cc,dim,dim,alpha,prob);
            target_cost_lb_table[dim] = target_cost_lb_basis(cc,1,dim);
        }
        return target_cost_lb_table[dim];
    }
    
    bkzfloat target_cost_ub(int dim) {

        if ((dim<10) || (dim>=2000)) {
            ccout << "target_cost_lb: dimension error " << dim << endl;
            exit(0);
        }  

        #pragma omp critical
        {
            int cs = target_cost_ub_table.size();
            if (dim >= cs) {
                target_cost_ub_table.resize(dim+1);
                for (int i=cs;i<=dim;i++) target_cost_ub_table[i] = 0;
            }
        }
        bkzfloat value;
        #pragma omp critical
        {
            value = target_cost_ub_table[dim];
        }
        if (value==0) {
            #pragma omp critical
            {
                ccout << "Simulating target_cost_ub dim=" << dim << "             \r";
                ccout.flush();
            }
            double alpha,r;
            bkzfloat prob;
            PBKZParam(alpha,r,prob,dim);

            //make gsa basis
            std::vector<bkzfloat> cc;
            PBKZSimulate_first(cc,dim,dim,alpha,prob);
            #pragma omp critical
            {
                ccout << "Simulating target_cost_ub_first dim=" << dim << "             \r";
                ccout.flush();
            }
            bkzfloat radius = alpha * lattice_tools::LatticeGH<bkzfloat>(cc,1,dim);
            #pragma omp critical
            {
                ccout << "Simulating target_cost_ub_gh dim=" << dim << "             \r";
                ccout.flush();
            }
            target_cost_ub_table[dim] = ENUMCostUB(cc,1,dim,radius,prob);
            #pragma omp critical
            {
                ccout << "Simulating target_cost_ub dim=" << dim << "  ... finished                \r";
                ccout.flush();
            }
        }
        return target_cost_ub_table[dim];
    }
        
    bkzfloat target_cost_ub_modify(int dim) {
        if (dim%2==0) return target_cost_ub(dim);
        return sqrt(target_cost_ub(dim-1) * target_cost_ub(dim+1));
    }
}

#define __bkzpreprocessbig
    #include "bkzpreprocessboostmain.cpp"
#undef __bkzpreprocessbig

#define __bkzpreprocesssmall
    #include "bkzpreprocessboostmain.cpp"
#undef __bkzpreprocesssmall

int dbcount = 0;    //for debug

template <typename PFLOAT,typename CFLOAT> bkzfloat ProgressiveBKZPreprocess(LatticeBasis<PFLOAT>&B,mat_ZZ* U,int istart,int iend,bkzfloat costlimit,bkzfloat alphalimit,double ualpha,bkzfloat uprob,int vl=0) {
            //costlb = ProgressiveBKZPreprocess<PFLOAT,double>(B,U,bstart,bend,costlimit_lb,falpha,ualpha,uprob,vl-1);
    B.updateGSBasis(1,iend);
    PFLOAT radius = lattice_tools::LatticeGH(B.gs.c,istart,iend) * ualpha;
    bkzfloat current_cost_lb = ENUMCostLB(B.gs.c,istart,iend,(bkzfloat)radius,uprob) ;
    
    if (vl>=2) ccout << "preprocess cost=" << current_cost_lb << " / " << costlimit << endl;

    //debug
    debug_display(
        if (isnan(current_cost_lb) == true) {
            ccout << "iend=" << iend << endl;
            B.gs.displayc();
            LatticeBasis<PFLOAT> tB = B;
            tB.GScomputed = 0;
            tB.updateGSBasis();
            tB.gs.displayc();

            exit(0);
        }
    );
    
    
    if (current_cost_lb < costlimit) return current_cost_lb;
    
    bkzfloat alpha = (bkzfloat)(B.gs.c[istart] / lattice_tools::LatticeGH(B.gs.c,istart,iend));
    if ((alphalimit >0) && (alpha < alphalimit)) return current_cost_lb;
    int dim = iend-istart+1;

    LatticeBasis<CFLOAT> LB;   //local basis
    mat_ZZ LU;
    LB.L.SetDims(dim,dim);
    LU.SetDims(dim,dim);
    ExtractSubBlock(LB,B,istart,iend);
    SetIdentity(LU);
    LB.updateGSBasis();
    current_cost_lb = ProgressiveBKZPreprocessCore<PFLOAT,CFLOAT>(LB,&LU,1,iend-istart+1,costlimit,alphalimit,ualpha,uprob,vl);

    if (vl>=2) ccout << "PPCore end. lb=" << current_cost_lb << endl;
    
    
    if (current_cost_lb < 0) {
        //error
        SaveLattice(LB,"lb");
        ccout << "error: istart=" << istart << " " << iend << endl;
        LB.gs.displayc();
        exit(0);
    }
    
    //Recover original basis
    mat_ZZ Btemp;
    Btemp.SetDims(dim,B.L.NumCols());

    for (int i=0;i<dim;i++) {
        for (int j=0;j<dim;j++) {
            Btemp[i] += LU[i][j] * B.L[j+istart-1];
        }
    }
    for (int i=0;i<dim;i++) {
        B.L[i+istart-1] = Btemp[i];
    }

    if (U!=0) {
        Btemp.SetDims(dim,(*U).NumCols());
        for (int i=0;i<dim;i++) {
            clear(Btemp[i]);
        }
        for (int i=0;i<dim;i++) {
            for (int j=0;j<dim;j++) {
                Btemp[i] += LU[i][j] * (*U)[j+istart-1];
            }
        }
        for (int i=0;i<dim;i++) {
            (*U)[i+istart-1] = Btemp[i];
        }
    }
    B.GScomputed = istart-1;
    B.updateGSBasis(istart,iend);
    return current_cost_lb;
}

#endif
