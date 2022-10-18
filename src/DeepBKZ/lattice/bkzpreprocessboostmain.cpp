
#include "bkzlibwrapper.hpp"

#ifndef _bkzpreprocessboostmain_common
#define  _bkzpreprocessboostmain_common
    
//Heuristic preprocessing subroutine to find a basis of low enumeration cost
template <typename PFLOAT,typename CFLOAT> bkzfloat ProgressiveBKZPreprocessCore(LatticeBasis<CFLOAT>&B,mat_ZZ* U,int istart,int iend,bkzfloat costlimit,bkzfloat alphalimit,double ualpha,bkzfloat uprob,int vl=0) {

    bkzfloat current_cost_lb,current_cost_lb_min=-1;
    int no_update_current_cost = 0;

    //Greedy BKZ
    int endflag = 0;
    int bstart;
    int bsize;
    int lowest_bsize = 20;
    bsize = lowest_bsize;
    int bmax = bsize;
    
    int loop=0;
    char direction = +1;    //-1=backward +1=forward

    if (direction == 1) bstart = istart;
    if (direction == -1) bstart =  istart + (iend-istart+1-bsize);
    
    int SRcomputed = 0;

    double radius = ualpha * lattice_tools::LatticeGH(B.gs.c,istart,iend);
    std::vector<CFLOAT> lbt;
    lbt.resize(iend+1);
    ENUMCostLB_maketable(lbt,istart,iend,uprob);

    
    
    
    while (endflag == 0) {
        //search suitable (bstart,bend) s.t. large |b*i| and small cost
        int localbs = min(bsize,iend-bstart+1);
        
        bkzfloat ghratio = B.gs.c[bstart]/lattice_tools::LatticeGH(B.gs.c,bstart,bstart+localbs-1);
        char executeflag = 0;
        bkzfloat blockcost;
        if (ghratio > 1.0 + 1.0/localbs) {
            //check cost
             blockcost =  bkzconstants::target_cost_lb_basis(B.gs.c,bstart,bstart+localbs-1);
             if (blockcost < 2 * bkzconstants::target_cost_lb(localbs)) {
                executeflag = 1;
            }
        }
        
        if (bsize==lowest_bsize) executeflag = 1;

        if (executeflag!=0) {
            double alpha,r;
            bkzfloat prob;
            PBKZParam(alpha,r,prob,localbs);
            int result;
            for (int i=SRcomputed+1;i<=bstart+localbs-1;i++) {
                SizeReduce(B,U,i);
            }
            double radius = alpha * (double)lattice_tools::LatticeGH(B.gs.c,bstart,bstart+localbs-1);
            result = EnumthenUpdate<CFLOAT>(B,U,bstart,bstart+localbs-1,radius,prob,0,vl,"");
            SRcomputed = bstart-1;    //Output of EnumthenUpdate do not guarantee size reduction

            if (result!=0) {
                if (vl>=2) {
                    ccout << "check bstart=" << bstart << " bsize=" << localbs << "           \r";
                    ccout.flush();
                }
                bmax = max(bmax,localbs);
                if (++loop%100==0) {
                    int swap;
                    local_LLL(B,U,0.999,istart,iend,stnormal,VL0,swap);
                }
            }
        }

        bstart += direction;
        int trigger_index = iend-bsize+2;
        if (direction == -1) trigger_index = istart-1;

        if (bstart==trigger_index) {
            
            current_cost_lb = ENUMCostLBwithtable(B.gs.c,lbt,istart,iend,radius); 

            if (current_cost_lb_min == -1 ) {
                current_cost_lb_min = current_cost_lb;
            }            
            if (current_cost_lb_min < current_cost_lb) {
                current_cost_lb_min = current_cost_lb;
                no_update_current_cost= 0;
            } else {
                no_update_current_cost++;
            }
            bkzfloat alpha = B.gs.c[istart] / lattice_tools::LatticeGH(B.gs.c,istart,iend);

            //Todo: preventing infinite loop
            
            if (vl>=1) {
                ccout << "preprocess block=(" << istart << "," << iend << ") cost_lb=" << current_cost_lb << "/" << costlimit;
                if (alphalimit>0) ccout << " alpha=" << alpha << "/" << alphalimit << " ";
                ccout << "  \r";
                ccout.flush();
                if (vl>=2) ccout << endl;
            }
            if (isnan(current_cost_lb) == true) {
                ccout << endl;
                ccout << endl;
                B.gs.displayc();
                return -1;
            }
            //here current_cost_lb=nan, and alpha=inf
            
            
            if (current_cost_lb < costlimit) break;
            if ((alphalimit >0) && (alpha < alphalimit)) break;
            if (no_update_current_cost > iend-istart) break;    //here, threshold=dim is heuristic
            
            bsize--;
            if (bsize<lowest_bsize) {
                int shift = 10;
                if (bmax > 50) shift = 5;
                bsize = min(bmax+shift,(iend-istart+1)-20);
                if (vl>=3) {
                    B.gs.displayc(istart,iend);
                }
            }
            
            if (direction == 1) bstart = istart;
            if (direction == -1) bstart =  istart + (iend-istart+1-bsize);
        }
    }
    return current_cost_lb;
}


#endif



#ifdef __bkzpreprocessbig
template <typename PFLOAT,typename CFLOAT> void ExtractSubBlock(LatticeBasis<CFLOAT>&  LB,LatticeBasis<PFLOAT>& B,int istart,int iend,std::string stringoptions="") {
#endif

#ifdef __bkzpreprocesssmall
template <typename PFLOAT,typename CFLOAT> void ExtractSubBlock(SmallLatticeBasis<CFLOAT,CFLOAT>&  LB,LatticeBasis<PFLOAT>& B,int istart,int iend,std::string stringoptions="") {
#endif
    
    //ccout << "extract subblock: " << istart << " " << iend << endl;
    //Extracting sublattice B[istart..iend] to LB using Gram-Schmidt data

    int ristart = istart;
    int riend = iend;
    int vl =0;
    
    if (stringoptions!= "") {
        std::map<std::string,std::string> options;
        ExtractOptions(options,stringoptions);
        if (options["ristart"]!="") {
            ristart = atoi( options["ristart"].c_str());
            ccout << "ristart=" << ristart << endl;
        }
        if (options["riend"]!="") {
            riend = atoi( options["riend"].c_str());
            ccout << "riend=" << riend << endl;
        }
        //vl = 2; //tentative code for debug
    } else {
    }
    B.updateGSBasis(1,iend);
    if (vl>=2) {
        B.gs.displayc();
    }

    //Find factor
    PFLOAT mings = B.gs.c[ristart];
    PFLOAT maxgs = B.gs.c[ristart];
    for (int i=ristart+1;i<=riend;i++) {
        if (B.gs.c[i]!=0) {
            mings = min(mings,B.gs.c[i]);
            maxgs = max(maxgs,B.gs.c[i]);
        }
    }    
    mings = mings * (mings / maxgs)  / (iend-istart+1);    //heuristic param.
    int dim = iend - istart + 1;

    if (vl>=2) {
        ccout << "Extract: mings=" << mings << endl;
        B.gs.displayc();
        
    }

    LB.L.SetDims(dim,dim);
    
    
    //diagonal
    for (int i=0;i<dim;i++) {
        for (int j=0;j<i;j++) {
            PFLOAT ll = B.gs.mu[i+istart][j+istart] * B.gs.c[j+istart]/mings;
#ifdef __bkzpreprocesssmall
            LB.L[i][j] = (CFLOAT)ll;
#endif
#ifdef __bkzpreprocessbig
            //need to convert to NTL_ZZ
            conv_to_ZZ(LB.L[i][j] ,ll);
#endif
        }
        PFLOAT ll = B.gs.c[i+istart]/mings;
#ifdef __bkzpreprocesssmall
        LB.L[i][i] = (CFLOAT)ll;
#endif
#ifdef __bkzpreprocessbig
        conv_to_ZZ(LB.L[i][i] ,ll);
#endif
        for (int j=i+1;j<dim;j++) {
             LB.L[i][j]  = 0;
        }
    }
    //Todo: fill gs
    
}   //End of ExtractSubBlock


#ifdef __bkzpreprocessbig
    void SetIdentity(mat_ZZ&  U) {
#endif

#ifdef __bkzpreprocesssmall
    template <typename T> void SetIdentity(pbkzmatrix<T>& U) {
#endif
        for (int i=0;i<U.NumCols();i++) {
            for (int j=0;j<U.NumRows();j++) {
                if (i==j) {
                    U[i][j] = 1;
                } else {
                    U[i][j] = 0;
                }
            }
        }
}
    
    
    