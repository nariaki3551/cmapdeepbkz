
#include "enumwrapper.cpp"
//Note: this template generates {single,multi}-threaded subroutines for {ENUM,ENUM-CVP} problems.



//#ifdef SSE_ROUND
#ifdef __SSE4_1__
        #include <emmintrin.h>
        #include <smmintrin.h>
#endif

#ifdef __AVX__
        #include <immintrin.h>
#endif


#ifdef define_single_enum_double 
    template <typename IFLOAT> void ENUMCoreboost_double(FoundENUMVectors& EV,int slide,int jj,int bs,double elim,bkzfloat clim,
                                                  std::vector<bkzfloat>& normalpf,std::vector<std::vector<IFLOAT> >& mu,
                                                  std::vector<IFLOAT>& cd,int vl,int enummode,int recursive=0) {
#endif

        
#ifdef define_single_enumcvp_double 
    template <typename IFLOAT> void ENUMCVCoreboost_double(FoundENUMVectors& EV,int slide,int jj,int bs,double elim,bkzfloat clim,
                                                  std::vector<bkzfloat>& normalpf,std::vector<std::vector<IFLOAT> >& mu,
                                                  std::vector<IFLOAT>& cd,std::vector<IFLOAT>& gsrept,int vl,int enummode,int recursive=0) {
#endif

#ifdef define_multi_enum_double 
    template <typename IFLOAT> void ENUMCoreboostmt_double(FoundENUMVectors& EV,int slide,int jj,int bs,double elim,int parallel,bkzfloat clim,
                                                  std::vector<bkzfloat>& normalpf,std::vector<std::vector<IFLOAT> >& mu,
                                                  std::vector<IFLOAT>& cd,int vl,int enummode,int finishmode,int recursive=0) {
#endif

#ifdef define_multi_enumcvp_double 
    template <typename IFLOAT> void ENUMCVCoreboostmt_double(FoundENUMVectors& EV,int slide,int jj,int bs,double elim,int parallel,bkzfloat clim,
                                                  std::vector<bkzfloat>& normalpf,std::vector<std::vector<IFLOAT> >& mu,
                                                  std::vector<IFLOAT>& cd,std::vector<IFLOAT>& gsrept,int vl,int enummode,int finishmode,int recursive=0) {
#endif

#ifdef define_single_enum 
    template <typename PFLOAT,typename PSHORT,typename IFLOAT> void ENUMCoreboost(FoundENUMVectors& EV,int slide,int jj,int bs,double elim,bkzfloat clim,
                                                  std::vector<bkzfloat>& normalpf,std::vector<std::vector<IFLOAT> >& mu,
                                                  std::vector<IFLOAT>& cd,int vl,int enummode,int recursive=0) {
#endif

#ifdef define_single_enumcvp 
    template <typename PFLOAT,typename PSHORT,typename IFLOAT> void ENUMCVCoreboost(FoundENUMVectors& EV,int slide,int jj,int bs,double elim,bkzfloat clim,
                                                  std::vector<bkzfloat>& normalpf,std::vector<std::vector<IFLOAT> >& mu,
                                                  std::vector<IFLOAT>& cd,std::vector<IFLOAT>& gsrept,int vl,int enummode,int recursive=0) {
#endif

#ifdef define_multi_enum 
    template <typename PFLOAT,typename PSHORT,typename IFLOAT> void ENUMCoreboostmt(FoundENUMVectors& EV,int slide,int jj,int bs,double elim,int parallel,bkzfloat clim,
                                                  std::vector<bkzfloat>& normalpf,std::vector<std::vector<IFLOAT> >& mu,
                                                  std::vector<IFLOAT>& cd,int vl,int enummode,int finishmode,int recursive=0) {
#endif

#ifdef define_multi_enumcvp 
    template <typename PFLOAT,typename PSHORT,typename IFLOAT> void ENUMCVCoreboostmt(FoundENUMVectors& EV,int slide,int jj,int bs,double elim,int parallel,bkzfloat clim,
                                                  std::vector<bkzfloat>& normalpf,std::vector<std::vector<IFLOAT> >& mu,
                                                  std::vector<IFLOAT>& cd,std::vector<IFLOAT>& gsrept,int vl,int enummode,int finishmode,int recursive=0) {
#endif
        

#ifdef define_single_enum_coeff_double    
    template <typename IFLOAT> void ENUMCoeffCoreboost_double(FoundENUMVectors& EV,int slide,int jj,int bs,double elim,bkzfloat clim,
                                                  std::vector<bkzfloat>& normalpf,std::vector<std::vector<IFLOAT> >& mu,
                                                  std::vector<IFLOAT>& cd,int vl,int enummode,int recursive=0) {
#endif

        
#if defined(define_multi_enum) || defined(define_multi_enumcvp) || defined(define_multi_enum_double) || defined(define_multi_enumcvp_double)
        #define define_multi
#endif
        
#if defined(define_single_enumcvp) || defined(define_multi_enumcvp)|| defined(define_single_enumcvp_double) || defined(define_multi_enumcvp_double)
        #define define_cvp
#endif

#if defined(define_single_enum) || defined(define_multi_enum) || defined(define_single_enum_double) || defined(define_multi_enum_double) || defined(define_single_enum_coeff_double ) 
        #define define_svp
#endif

#if defined(define_single_enum_double) || defined(define_multi_enum_double) || defined(define_single_enumcvp_double) || defined(define_multi_enumcvp_double) || defined(define_single_enum_coeff_double ) 
        #define define_double
        #define PFLOAT double
        #define PSHORT double
#endif

#if defined(define_single_enum_coeff_double)
        #define define_coeff
#endif
        
        
        //PFLOAT is type to store Gram-Schmidt presentations
        //PSHORT is type to store coefficients(integers)

         PFLOAT* ctilda = (PFLOAT*)shared_memory::allocate1<PFLOAT>(0+100*recursive,bs+2) - jj; 
         PSHORT* vvec = (PSHORT*)shared_memory::allocate1<PSHORT>(1+100*recursive,bs+2) - jj; 
         PFLOAT* yvec = (PFLOAT*)shared_memory::allocate1<PFLOAT>(2+100*recursive,bs+2) - jj; 
         PSHORT* utildavec = (PSHORT*)shared_memory::allocate1<PSHORT>(4+100*recursive,bs+2) - jj; 
         PFLOAT* cshift = (PFLOAT*)shared_memory::allocate1<PFLOAT>(5+100*recursive,bs+2) - jj; 
         long int *Deltavec = (long int*)shared_memory::allocate1<long int>(0+100*recursive,bs+2) - jj; 
         long int *deltavec = (long int*)shared_memory::allocate1<long int>(1+100*recursive,bs+2) - jj; 

#ifdef define_multi
        //memory for ENUM1a
        PFLOAT *qctilda = (PFLOAT*)shared_memory::allocate1<PFLOAT>(20+100*recursive,bs+2) - jj;
        PSHORT *qvvec = (PSHORT*)shared_memory::allocate1<PSHORT>(21+100*recursive,bs+2) - jj;
        PFLOAT *qyvec = (PFLOAT*)shared_memory::allocate1<PFLOAT>(22+100*recursive,bs+2) - jj;
        PSHORT *qutildavec =(PSHORT*)shared_memory::allocate1<PSHORT>(24+100*recursive,bs+2) - jj;
        long int *qDeltavec = (long int*)shared_memory::allocate1<long int>(25+100*recursive,bs+2) - jj;
        long int *qdeltavec =(long int*)shared_memory::allocate1<long int>(26+100*recursive,bs+2) - jj;
        int sbk,tbk;    //backup of s and t in ENUM1a
        
        int numslot = 10 * parallel;
        
        PFLOAT** soctilda = (PFLOAT**)shared_memory::allocate2<PFLOAT>(11+100*recursive,numslot,bs+4);
        PSHORT** sovvec = (PSHORT**)shared_memory::allocate2<PSHORT>(12+100*recursive,numslot,bs+4);
        PFLOAT** soyvec = (PFLOAT**)shared_memory::allocate2<PFLOAT>(13+100*recursive,numslot,bs+4);
        PSHORT** soutildavec = (PSHORT**)shared_memory::allocate2<PSHORT>(15+100*recursive,numslot,bs+4);
        long int** soDeltavec = (long int**)shared_memory::allocate2<long int>(16+100*recursive,numslot,bs+4);
        long int** sodeltavec = (long int**)shared_memory::allocate2<long int> (17+100*recursive,numslot,bs+4);
        int* sos = (int*)shared_memory::allocate1<int>(25+100*recursive,numslot);
        for (int i=0;i<numslot;i++) sos[i] = -1;
#endif
        

#ifdef define_coeff
        //coefficient range check
        if ((EV.coeff_lower==NULL) || (EV.coeff_upper==NULL)) {
            ccout << "Coefficient range not defined!" << endl;
            return;
        }
#endif
        //int* coeff_lower = EV.coeff_lower;

        
        int kk = jj+bs-1;
         slide = bs;
         
         if (enummode &  enum_mode_all_vectors) slide=0;
         if (enummode &  enum_mode_find_shortest) {
             slide=0;
#ifdef define_svp
             if ((enummode & enum_mode_count_only) ==0 ) clim = min(clim,(bkzfloat)cd[jj]);
#endif
         }
         if (enummode &  enum_mode_find_abort) slide=0;
         
         int lmctlimit = 1000000;
         if (typeid(PFLOAT)==typeid(double)) lmctlimit = 50000000;
         if (typeid(PFLOAT)==typeid(long double)) lmctlimit = 20000000;
         
#ifdef depthcost
         long int totaldepthcost=0;
         long int** depthcostarray;
#ifdef define_multi
         depthcostarray = new long int*[parallel];
#else
         depthcostarray = new long int*[1];
#endif
#ifdef define_multi
         for (int i=0;i<parallel;i++) {
            depthcostarray[i] = new long int[jj+bs+1];
            for (int j=0;j<jj+bs+1;j++) depthcostarray[i][j] = 0;
         }
#else
         depthcostarray[0] = new long int[jj+bs+1];
         for (int j=0;j<jj+bs+1;j++) depthcostarray[0][j] = 0;
#endif
         EV.depthcostarray.resize(jj+bs+1);
         for (size_t i=0;i<EV.depthcostarray.size();i++) {
             EV.depthcostarray[i] = 0;
         }
         
#endif

         PFLOAT* pruning_func = (PFLOAT*)shared_memory::allocate1<PFLOAT>(7+100*recursive,bs+1);
         PFLOAT* pftemp = pruning_func - jj;
         for (int i=0;i<bs;i++) {
             pruning_func[i] = (PFLOAT)normalpf[kk-i] * (PFLOAT)clim;
             if (vl>=5) ccout << "pf[" << i << "]=" << normalpf[kk-i] << " " << clim << " " << pruning_func[i]  << " " << cd[jj+i] << endl;
         }

         if (enummode & enum_mode_find_short_projection)  {
             for (int i=0;i<0.25*(kk-jj);i++) {
                 //debug_display( ccout << "pruning_func[" << i << "]=" << pruning_func[i] << " " << cd[jj+i] << endl; )
                 pruning_func[i] = max(pruning_func[i],(PFLOAT)cd[jj+i]);
                 if (vl>=5) ccout << "pf'[" << i << "]=" << normalpf[kk-i] << " " << clim << " " << pruning_func[i]  << endl;
             }
         }
         //initialize memory for vector holding
         EV.clearmemory();

         //transmu
         PFLOAT** transmu = (PFLOAT**)shared_memory::allocate2<PFLOAT>(0+100*recursive,kk+1,kk+1);
         for (int i=0;i<bs;i++) {
             for (int j=i;j<bs;j++) {
                 transmu[i][j] = (PFLOAT)mu[j+jj][i+jj];
             }
         }
         //cd[i]=|b*_i|
         for (int i=jj;i<=kk;i++) {
             cshift[i] = (PFLOAT)cd[i];
              //if (vl>=5) ccout << "cshift[" << i << "]=" << cshift[i] << endl;
         }

         
#ifdef define_multi
         int jjs = kk-1;
         int numsubtree;
         //Automatic detection of jjs
         while (1) {
             std::vector<bkzfloat> subnormalpf;
             subnormalpf.resize(kk+1);
             for (int i=0;i<kk-jjs+1;i++) {
                 subnormalpf[jjs+i] = normalpf[jj+i];
             }
             
             
#ifdef define_double
            ENUMCoreboost_double<IFLOAT>(EV,0,jjs,kk-jjs+1,1,clim,subnormalpf,mu,cd,VL0,enum_mode_count_only,recursive+1);
#else
            ENUMCoreboost<PFLOAT,PSHORT>(EV,0,jjs,kk-jjs+1,1,clim,subnormalpf,mu,cd,VL0,enum_mode_count_only,recursive+1);
#endif
            int foundnum = EV.foundnum;
            //debug_display (ccout << "jjs=" << jjs << " " << foundnum << endl; );
            EV.clearmemory();
            if ((foundnum > 20*parallel) || (kk-1-jjs  > 0.5 * bs)) {
                if (vl>=2) {
                    ccout << "parallel subtree depth=" << kk-jjs << "   " << endl;
                    ccout << "parallel subtree nodes=" << foundnum << "   " <<  endl;
                }
                numsubtree = foundnum;
                break;
            }
            jjs--;
         }
#endif


#ifdef define_cvp
         PFLOAT* target_mud = (PFLOAT*)shared_memory::allocate1<PFLOAT>(50+100*recursive,bs+2) - jj;
         for (int i=jj;i<=kk;i++) target_mud[i] = (PFLOAT)gsrept[i];
         //computing bias=nearest_plane(L,t)
         //this algorithm computes vectors close to t-bias
         PFLOAT* cvpbias = (PFLOAT*)shared_memory::allocate1<PFLOAT>(51+100*recursive,bs+2) - jj;
         for (int i=kk;i>=jj;i--) {
             cvpbias[i] = round(target_mud[i]);
             target_mud[i] -= cvpbias[i];
             for (int j=jj;j<i;j++) {
                 target_mud[j] -= cvpbias[i] * (PFLOAT)mu[i][j];
             }
             //ccout << "target_mud[" << i << "]" << target_mud[i] << " " << "cvpbias[" << i << "]" << cvpbias[i] << endl;
         }
#endif

         //initialize variables for lattice enumeration
         for (int i = jj; i <= kk+1; i++) {
            ctilda[i] = utildavec[i] = yvec[i] = 0;
            Deltavec[i] = 0;
            vvec[i] = 0;

            deltavec[i] = 1;
#ifdef define_cvp
            if ( (i<=kk) && (target_mud[i] < 0)) deltavec[i] = -1;
            yvec[i] = -target_mud[i];
#endif
         }
         utildavec[jj] = 1;
         yvec[jj] = vvec[jj] = 0;
         Deltavec[jj] = 0;
         deltavec[jj] = 1;   
#ifdef define_cvp
            yvec[jj] = -target_mud[jj];
#endif
            
#ifdef define_cvp
         //initial_value=a distance when a[i]=0 (projective norm of target point)
         ctilda[kk] = target_mud[kk]*target_mud[kk]*cshift[kk];
         for (int i=kk-1;i>=jj;i--) {
              ctilda[i] = ctilda[i+1] + target_mud[i]*target_mud[i] * cshift[i];
              utildavec[i] = round(target_mud[i]);
         }
#endif
         
         bool flag_cvpzeroadd = false;

#ifdef define_cvp
         //if target is shorter than radius, zero vector is a solution
         //sometimes the nearest vector is not found in the resulting list
         PFLOAT tnorm = 0;
         for (int i=kk;i>=jj;i--) {
             tnorm += target_mud[i] * target_mud[i] * cshift[i];
         }
         if (tnorm <= clim) {
             //add it as the solution
             storetomemory(EV,utildavec,tnorm,jj,kk,enummode);
             flag_cvpzeroadd = true;
         }
         //exit(0);
#endif

         
#ifdef define_multi
        qdeltavec[jjs] = 1;   
        sbk = tbk = jjs;
        for (int i = jj; i <= kk+1; i++) {
            qctilda[i] = qutildavec[i] = qyvec[i] = 0;
            qDeltavec[i] = 0;
            qvvec[i] = 0;
            qdeltavec[i] = 1;
#ifdef define_cvp
            if ( (i<=kk) && (target_mud[i] < 0)) qdeltavec[i] = -1;
            qyvec[i] = -target_mud[i];
#endif
        }

#ifdef define_cvp
        qutildavec[jjs] =1;
        qyvec[jjs] = qvvec[jjs] = 0;
        qDeltavec[jjs] = 0;
        qyvec[jjs] = -target_mud[jjs];
        sbk = tbk = jjs;
         qctilda[kk] = target_mud[kk]*target_mud[kk]*cshift[kk];
         for (int i=kk-1;i>=jj;i--) {
              qctilda[i] = qctilda[i+1] + target_mud[i]*target_mud[i] * cshift[i];
              qutildavec[i] = round(target_mud[i]);
         }
#endif


         
         
         
         
         
        int ea1thread=-1;  //thread number working for 
        int totalsubstore=0;
        
        //flags used in subroutines
        int endflag=0;  // set endflag=kk to abort all enums
        char ea1flag=0;   //ea1flag=0 <=> No thread works for ENUM1a
                          //       =1 <=>  a thread works
                          //       =2 <=> ENUM1a is finished
        char ea1breakflag=0;    //sufficient number of vectors may be slotted in ENUM1a

#endif

        long double gmct = 0;   //global enum count (/10^6)
#ifdef depthcost
        long double gmct_inside = 0;   //global enum count (/10^6)
#endif
        long double gh = conv_to_long_double(lattice_tools::LatticeGH(cd,jj,jj+bs-1,INPUT_SQUARED));
        
         /* Display messages*/
         if (vl>=3) {
#ifdef define_svp
             ccout << "Start ENUM-SVP (initial radius=" << sqrt(clim) << ")" << endl;
#endif
#ifdef define_cvp
             ccout << "Start ENUM-CVP (initial radius=" << sqrt(clim) << ")" << endl;
#endif
             std::vector<bkzfloat> cc;
             cc.resize(cd.size());
             if (vl>=5) {
                for (int i=jj;i<min(jj+bs,(int)cd.size());i++) ccout << i << " " << cd[i] << endl;
             }
             if (vl>=5) {
                for (int i=jj;i<min(jj+bs,(int)normalpf.size());i++) ccout << i << " " << normalpf[i] << endl;
             }
             for (int i=jj;i<min(jj+bs,(int)cd.size());i++) cc[i] = sqrt((bkzfloat)cd[i]);
             bkzfloat upper = pruning_func::Rigid_upper_cost<bkzfloat,bkzfloat>(normalpf,cc,sqrt(clim),jj,bs);
             EV.expectednodes = upper/2;
             ccout << "Rigid_upper_bound=" << upper << endl;
             ccout << "LatticeGH=" << gh << endl;
         }
       /* end of displays */
       
        
        double es = gettimeofday_sec();   //timer
        double cts = getcputime();
        
 //ENUM main
#ifdef define_multi
        int vfinished=0;    //# finished trees

#pragma omp parallel  num_threads(parallel)
{
#endif
            
        //variables used in multithread
#ifdef __no__multithreads
        int mythread =  0;
#else
        int mythread =  omp_get_thread_num();
#endif
        long long int lmct = 0;  
#ifdef depthcost
        long long int lmct_inside = 0;  
#endif
        
        int s,t;
        s = t = jj;
        char cflag;
        PFLOAT t1d;
        PFLOAT dtemp;
        int findex=0;
        PFLOAT mlen;       //minimum length of found vector
        
#ifdef define_multi
        // ENUM-0
        if (mythread==0) {
            while (t <= jjs-1) {
#else
            while (t <= kk) {
#endif
                dtemp = (yvec[t]+utildavec[t]);
                dtemp = dtemp * dtemp * cshift[t];
                ctilda[t] = ctilda[t+1] + dtemp;
                cflag = 0;

#ifndef define_coeff

                if (ctilda[t] < pftemp[t]) {    //pruning by the projective lengths
#else
                if ((ctilda[t] < pftemp[t]) && ( EV.coeff_lower[t] <= utildavec[t]  ) && ( EV.coeff_upper[t] >= utildavec[t]  )  ) {    //pruning by the projective lengths and coefficients
#endif

#ifdef depthcost
                    lmct_inside++;  
#endif
                    if (t<=jj+slide) {
                        //a vector s.t. pi_{i}|v|<|b*_i| (i=jj+slide+4) is found
                        if ( (enummode & (enum_mode_find_shortest + enum_mode_find_abort + enum_mode_find_short_projection)) != 0) {
                            if ( (enummode & enum_mode_find_short_projection) != 0) { //for t=jj stored after
                                #pragma omp critical
                                {
                                    if ((enummode & enum_mode_count_only) ==0) storetomemory(EV,utildavec,ctilda[t],t,kk,enummode);
                                    EV.foundnum++;
                                }
                                if (t-jj+6 < slide) {
                                    slide = max(0,t-jj+6);
                                }
                            } else                         
                            if (ctilda[t]<cshift[t]*0.9999) {
                                //Store the found vector to local memory 
                                #pragma omp critical
                                {
                                    if ((enummode & enum_mode_count_only) ==0 )storetomemory(EV,utildavec,ctilda[t],t,kk,enummode);
                                    EV.foundnum++;
                               }
                                if (t-jj < slide) {
                                    slide = t-jj;
                                }
                                mlen = 0;
                                if (t==jj) {
                                    mlen = ctilda[t];
                                }
                                if (mlen>0) {
                                    if ((mlen < clim) && (slide==0) && ((enummode & enum_mode_all_vectors)==0)) {
                                        //update pruning function if the pruning radius is updated
#ifdef define_svp
                                        if ((enummode & enum_mode_count_only) ==0 ) {
                                            clim = mlen;
                                            for (int i=0;i<kk-jj+1;i++) {
                                                pruning_func[i] = (PFLOAT)normalpf[kk-i] * (PFLOAT)clim;
                                            }
                                            if ((enummode & enum_mode_find_short_projection)!=0) {
                                                for (int i=0;i<0.25*(kk-jj);i++) {
                                                    pruning_func[i] = max(pruning_func[i],(PFLOAT)cd[jj+i]);
                                                }
                                            }
                                        }
#endif
                                        if (enummode & enum_mode_find_abort) {
#ifdef define_multi
                                            endflag = 1;
#endif
                                            break;
                                        }
                                    }
                                }
                            }
                        }

                        if (t==jj) {
                            if (enummode ==0 ){
                                #pragma omp critical
                                {
                                    EV.foundnum++;
                                    storetomemory(EV,utildavec,ctilda[t],t,kk,enummode);
                                }
                            }


                            if ((enummode & enum_mode_find_abort) !=0 ){
                                #pragma omp critical
                                {
                                        if ( (enummode & enum_mode_find_short_projection) == 0) {
                                            EV.foundnum++;
                                            if ((enummode & enum_mode_count_only) ==0 ) storetomemory(EV,utildavec,ctilda[t],t,kk,enummode);
                                        }
                                }
                                break;
                            }

                            if ((enummode & enum_mode_all_vectors+enum_mode_count_only + enum_mode_find_short_projection) !=0) {
                                #pragma omp critical
                                {
                                    if ( (enummode & enum_mode_find_short_projection) == 0) {
                                        EV.foundnum++;
                                        if ((enummode & enum_mode_count_only) ==0 ) storetomemory(EV,utildavec,ctilda[t],t,kk,enummode);
                                    }
                                }
                                //store v + k*b_1 with small norm
                                PFLOAT rr = yvec[t]+utildavec[t];
                                PFLOAT sign = +1;
                                PFLOAT utback = utildavec[t];

    #ifdef define_svp
                                if ((rr>0) && (s>jj)) sign=-1;
    #endif
    #ifdef define_cvp
                                if (rr>0) sign=-1;
    #endif

                                int i=2;
                                while (1) {
                                    dtemp = yvec[t]+utback + sign*(int)(i/2);
                                    dtemp = dtemp * dtemp * cshift[t];
                                    ctilda[t] = ctilda[t+1]+dtemp;
                                    if (ctilda[t]< pftemp[t]) {
                                        utildavec[t] = utback +  sign*(int)(i/2);
                                        #pragma omp critical
                                        {
                                            if ((enummode & enum_mode_count_only) ==0 )storetomemory(EV,utildavec,ctilda[t],t,kk,enummode);
                                            EV.foundnum++;
                                        }
                                     } else {
                                        break;
                                     }
                                    if (s>jj) {
                                        i++;
                                        sign=-sign;
                                    } else {
                                        i += 2;
                                    }
                                }
                            }

                        }

                    }

                    if (t > jj) {
                        cflag = 1;
                        t--;
                        t1d = 0;
                        PFLOAT* ttmu;
                        PSHORT* tut;
                        ttmu = transmu[t-jj];
                        tut = utildavec + t + 1;
                        ttmu += t+1-jj;
                        int i = s - t + 1;


                #if defined(define_double) && defined(__AVX__)
                        __m256d ymm;
                        double atmp[4] __attribute__((aligned(32)));
                        i = s - t;
                        ymm = _mm256_setzero_pd();
                        while(i & 3){
                            t1d += (*tut) * (*ttmu);
                            ++tut, ++ttmu, --i;
                        }
                        while(i > 0){
                            ymm = _mm256_add_pd(ymm, _mm256_mul_pd(_mm256_loadu_pd(tut), _mm256_loadu_pd(ttmu)));
                            tut += 4, ttmu += 4, i -= 4;
                        }
                        _mm256_store_pd(atmp, ymm);
                        for(int j = 0; j < 4; ++j) t1d += atmp[j];
                #else 
                          while (--i) {
                             t1d += (*tut)*(*ttmu);
                             ttmu++;
                             tut++;
                          }
                #endif


    #ifdef define_svp
                        yvec[t] = t1d;
    #endif
    #ifdef define_cvp
                        yvec[t] = t1d- target_mud[t];
                        t1d = yvec[t];
    #endif

    #if defined(define_double) && defined(__SSE4_1__)
                    __m128d tmp;
                    tmp = _mm_load_sd(&t1d);
                    tmp = _mm_round_sd(tmp, tmp, 0);
                    _mm_store_sd(&t1d, tmp);
                    t1d = -t1d;
    #else
                    t1d = -round(t1d);
    #endif
                        utildavec[t] = vvec[t] = t1d;
                        Deltavec[t] = 0;
                        if (utildavec[t] > -yvec[t]) 
                           deltavec[t] = -1;
                        else
                           deltavec[t] = 1;
                      }
#ifdef depthcost
                depthcostarray[0][t]++;
                totaldepthcost += (s-t+1);     //cost for GS-coefficients 
#endif
                }

                if (cflag==0) {
                   t++;
    #ifdef define_svp
                   if (s>t) {
                       Deltavec[t] = -Deltavec[t];
                   } else {
                       s = t;
                   }
    #endif

    #ifdef define_cvp
                       Deltavec[t] = -Deltavec[t];
                       s = max(s,t);
    #endif

                   Deltavec[t] += (1+((Deltavec[t]*deltavec[t]) >> 31))*deltavec[t];
                   utildavec[t] = vvec[t] + Deltavec[t];
               }

                //counter
                if ((++lmct)^lmctlimit) {
                } else {
                    gmct += lmct / 1000000.0;
                    lmct %= 1000000;
                   if (gmct/1000 >elim) {
                       break;
                   }
#ifdef depthcost
                    for (int i=0;i<=s;i++) {
                           EV.depthcostarray[i] += depthcostarray[0][i];
                           depthcostarray[0][i] = 0;
                   }
#endif
#ifndef define_multi
                    if ((enummode & enum_mode_all_vectors) && (vl>=1)) {
                        ccout << "ENUM(" << s-jj+1 << "/" << kk-jj+1 << ") count=" <<gmct / 1000.0 + lmct / 1e+8 << "G";
                        ccout << "#found vector=" << EV.foundnum << " time="  << gettimeofday_sec() - es << "\r";
                        ccout.flush();
                    } else {
                        if (vl>=1) {
                            ccout << "ENUM(" << s-jj+1 << "/" << kk-jj+1 << ") count=" <<gmct / 1000.0 + lmct / 1e+8 << "G";
                            ccout << " time="  << gettimeofday_sec() - es;
                            ccout << " speed=" << (gmct + lmct / 1e+6) / (gettimeofday_sec() - es) << "M/s  " << "\r";
                            ccout.flush();
                        }
                    }
#endif
#ifdef define_multi
                 t += endflag;
#endif
                }
#ifdef depthcost
                if ((lmct_inside)^lmctlimit) {
                } else {
                     gmct_inside += lmct_inside / 1000000.0;
                     lmct_inside %= 1000000;
                }
#endif

            
            }  //end of enum
            gmct += 0.000001 * lmct;
            lmct = 0;
#ifdef depthcost
            gmct_inside += 0.000001 * lmct_inside;
            lmct_inside = 0;
#endif
            
#ifdef define_multi
        }       // end of if(mythread==0)
        //below the subroutine for multithread
        int vcount=0;
        //start of ENUM1a
        do {     
            if (ea1flag==0) {
                #pragma omp critical
                {
                    if (ea1thread==-1) {
                        ea1thread = mythread; 
                    }
                }
            }            
            if ((mythread==ea1thread) && (ea1flag==0) && (endflag==0)) {
                int current_used_slot=0;
                for (int i=0;i<numslot;i++) {
                    if (sos[i]!=-1) current_used_slot++;
                }
                if (current_used_slot > 0.9 * numslot) {
                    ea1flag=0;  //no thread works for ENUM-1a
                } else {
                    //Start of ENUM1a
                    //Searching vector candidates while current_used_slot < 0.9*numslot
                    t = tbk;
                    s = sbk;
                    ea1flag = 1;
                    int vindex=0;
                    //Start of ENUM1a-core
                    while (t <= kk) {
                        dtemp = (qyvec[t]+qutildavec[t]);
                        dtemp = dtemp * dtemp * cshift[t];
                        qctilda[t] = qctilda[t+1] + dtemp;
                        cflag = 0;

                        if (qctilda[t] <  pftemp[t]) {
                            if (t >= jjs) {
                                cflag = 1;
                                t--;
                                t1d = 0;
                                PFLOAT* ttmu;
                                PSHORT* tut;
                                ttmu = transmu[t-jj];
                                tut = qutildavec + t + 1;
                                ttmu += t+1-jj;
                                int i = s - t + 1;
                                while (--i) {
                                    t1d += (*tut)*(*ttmu);
                                    ttmu++;
                                    tut++;
                                }
#ifdef define_svp
                                qyvec[t] = t1d;
#endif
#ifdef define_cvp
                                qyvec[t] = t1d - target_mud[t];
#endif
                                t1d = -round(qyvec[t]);   

                                qutildavec[t] = qvvec[t] = t1d;
                                qDeltavec[t] = 0;


                                if (qutildavec[t] > -qyvec[t]) { 
                                    qdeltavec[t] = -1;
                                } else {
                                    qdeltavec[t] = 1;
                                }

                                if (t+1==jjs) {
                                    //a vector whose projection is under the pruning function is found
                                    char zeroflag = 0;
                                    for (int i=jjs;i<=kk;i++) {
                                        if (qutildavec[i]!=0) {
                                            zeroflag = 1;
                                            break;
                                        }
                                    }
                                            
                                    //search the empty slot
                                    if ((qctilda[t+1]!=0) && (zeroflag==1)) {
                                        //the subtree of qctilda(projective length)=0 is already searched in ENUM0
                                        while (sos[vindex]!=-1){
                                            vindex++;
                                            if (vindex==numslot) {
                                                vindex = 0;
                                                ea1breakflag=1;
                                            }
                                        }
                                        totalsubstore++;
                                        //store the state into slot
                                        for (i=jj;i<=kk+1;i++) {
                                            soctilda[vindex][i-jj+1] = qctilda[i];
                                            sovvec[vindex][i-jj+1] = qvvec[i];
                                            soyvec[vindex][i-jj+1] = qyvec[i];
                                            soutildavec[vindex][i-jj+1] = qutildavec[i];
                                            sodeltavec[vindex][i-jj+1] = qdeltavec[i];
                                            soDeltavec[vindex][i-jj+1] = qDeltavec[i];
                                        }
                                        sos[vindex] = s;    //set the slag for ENUM-1b process
                                        vindex++;
                                        if (vindex==numslot) ea1breakflag=1;
                                    }

                                    cflag=0;
                                }
                            }
                        } 

                        if (cflag==0) {
                            t++;
#ifdef define_svp

                            if (s>t) {
                                qDeltavec[t] = -qDeltavec[t];
                            } else {
                                s = t;
                            }
#endif
#ifdef define_cvp
                            qDeltavec[t] = -qDeltavec[t];
                            s = max(s,t);
#endif
                            qDeltavec[t] += (1+((qDeltavec[t]*qdeltavec[t]) >> 31))*qdeltavec[t];
                            qutildavec[t] = qvvec[t] + qDeltavec[t];
                        }
                        lmct++; //counter

                        //Does sufficient number of vectors stored?
                        if (ea1breakflag==1) {
                            current_used_slot=0;
                            for (int i=0;i<numslot;i++) {
                                if (sos[i]!=-1) current_used_slot++;
                            }
                            if (current_used_slot>=numslot*0.9) {
                                break;
                            }
                            vindex = 0;
                            ea1breakflag=0;
                        }
                        t += endflag;
                    } //End of ENUM1a-core

                    sbk = s;
                    tbk = t;
                    if (t==kk+1) {
                        //ENUM1a is finished
                        ea1flag=2;
                    } else {
                        ea1flag=0;  //no thread works for ENUM-1a
                    }
                    #pragma omp critical
                    {
                        ea1thread=-1;
                        #pragma omp atomic    
                        gmct += lmct / 1000000;
                        lmct %= 1000000;
                    }            
                }  //End of ENUM-1a
            }            

            //Start of ENUM-1b
            int vindex = mythread;

            //sos=-1 <=> slot is empty
            //sos>0 <=> slot is waiting for process
            //-10000<sos<-100 <=> slot is reserved for process
            //sos=-29999 <=> slot is processed
            do {
                if (ea1flag==2) {
                    //if ENUM-1a is finished, search an empty slot incrementaly
                    #pragma omp critical
                    {
                        vindex = 0;
                        while (vindex<numslot) {
                            if (sos[vindex]>0) {
                                sos[vindex]-=10000;
                                break;
                            }
                            vindex++;
                        }
                    }
                } else {

                    if (sos[vindex]>0) {
                        sos[vindex] -=10000;
                    }
                }

                if ((vindex<numslot) && (-10000<sos[vindex]) && (sos[vindex]<-100) && (endflag==0)) {  

                    //Start of ENUM-1b for the slot
                    t = jjs-1;
                    s = sos[vindex]+10000;
                    sos[vindex]=-29999;
#ifdef depthcost
                    long int* mycostarray = depthcostarray[mythread];
#endif
                    #pragma omp critical
                    {
                        vcount++;
                    }

                    //memory copy
                    PSHORT* putildavec = soutildavec[vindex]-jj+1; 
                    PFLOAT* pyvec = soyvec[vindex]-jj+1; 
                    PFLOAT* pvvec = sovvec[vindex]-jj+1; 
                    PFLOAT* pctilda = soctilda[vindex]-jj+1; 
                    long int* pdeltavec = sodeltavec[vindex]-jj+1;
                    long int* pDeltavec = soDeltavec[vindex]-jj+1;

                    //Start of ENUM-1b-core
                    while (t < jjs) {
                        dtemp = (pyvec[t]+putildavec[t]);
                        dtemp = dtemp * dtemp * cshift[t];
                        pctilda[t] = pctilda[t+1] + dtemp;
                        cflag = 0;


                        if (pctilda[t] < pftemp[t]) {
                            if (t<=jj+slide) {
                                //a vector s.t. pi_{i}|v|<|b*_i| (i=jj+slide+4) is found
                                if ( (enummode & (enum_mode_find_shortest + enum_mode_find_abort + enum_mode_find_short_projection)) != 0) {
                                    if ( (enummode & enum_mode_find_short_projection) != 0) {
                                        #pragma omp critical
                                        {
                                            if ((enummode & enum_mode_count_only) ==0) storetomemory(EV,putildavec,pctilda[t],t,kk,enummode);
                                            EV.foundnum++;
                                        }
                                        if (t-jj+6 < slide) {
                                            slide = max(0,t-jj+6);
                                        }
                                    } else                         
                                    if (pctilda[t]<cd[t]*cd[t]*0.9999) {
                                        //Store the found vector to local memory 
                                        #pragma omp critical
                                        {
                                           if ((enummode & enum_mode_count_only) ==0 ) storetomemory(EV,putildavec,pctilda[t],t,kk,enummode);
                                           EV.foundnum++;
                                        }
                                        if (t==jj) mlen = ctilda[t];
                                        if (mlen >= 0) {
                                            if (t-jj < slide) {
                                                slide = t-jj;
                                            }
                                        }
                                        if ((mlen > 0 ) && (mlen < clim) && (slide==0) && ((enummode & enum_mode_all_vectors)==0)) {
                                            //update pruning function if the pruning radius is updated
                                            #pragma omp critical
                                            {
#ifdef define_svp
                                            if ((enummode & enum_mode_count_only) ==0 ) {
                                                clim = mlen;
                                                for (int i=0;i<kk-jj+1;i++) {
                                                    pruning_func[i] = (PFLOAT)normalpf[kk-i] * (PFLOAT)clim;
                                                }
                                                if ((enummode & enum_mode_find_short_projection)!=0) {
                                                    for (int i=0;i<0.25*(kk-jj);i++) {
                                                        pruning_func[i] = max(pruning_func[i],(PFLOAT)cd[jj+i]);
                                                    }
                                                }
                                            }
#endif
                                                if (enummode & enum_mode_find_abort) {
                                                    endflag = kk;
                                                }
                                            }
                                        }
                                    }
                                }
                                if (t==jj) {
                                    if (enummode ==0 ){
                                        #pragma omp critical
                                        {
                                            EV.foundnum++;
                                            storetomemory(EV,utildavec,ctilda[t],t,kk,enummode);
                                        }
                                    }
                                    if ((enummode & enum_mode_find_abort) !=0 ){
                                        #pragma omp critical
                                        {
                                            if ( (enummode & enum_mode_find_short_projection) == 0) {
                                                EV.foundnum++;
                                                if ((enummode & enum_mode_count_only) ==0 ) storetomemory(EV,putildavec,pctilda[t],t,kk,enummode);
                                            }
                                        }
                                        if (pctilda[t]>0) {
                                            endflag = kk;
                                        }
                                    }
                                    if ((enummode &  enum_mode_all_vectors + enum_mode_count_only + enum_mode_find_short_projection)!=0) {
                                    //holding
                                        #pragma omp critical
                                        {
                                            if ( (enummode & enum_mode_find_short_projection) == 0) {
                                                if ((enummode & enum_mode_count_only) ==0 )storetomemory(EV,putildavec,pctilda[t],t,kk,enummode);
                                                EV.foundnum++;
                                            }
                                        }
                                        //store v + k*b_1 with small norm
                                        PFLOAT rr = pyvec[t]+putildavec[t];
                                        PFLOAT sign = +1;
                                        PFLOAT utback = putildavec[t];

                                        if (rr>0) sign=-1;

                                        int i=2;
                                        while (1) {
                                            dtemp = pyvec[t]+utback + sign*(int)(i/2);
                                            dtemp = dtemp * dtemp * cshift[t];
                                            pctilda[t] = pctilda[t+1]+dtemp;

                                            if (pctilda[t]< clim) {
                                                putildavec[t] = utback +  sign*(int)(i/2);
                                                #pragma omp critical
                                                {
                                                    if ((enummode & enum_mode_count_only) ==0 )storetomemory(EV,putildavec,pctilda[t],t,kk,enummode);
                                                    EV.foundnum++;
                                                }
                                            } else {
                                                break;
                                            }

                                                i++;
                                                sign=-sign;

                                        }
                                    }
                                }
                            }

                            if (t != jj) {
                                cflag = 1;
                                t--;
                                t1d = 0;
                                PFLOAT* ttmu;
                                PSHORT* tut;
                                ttmu = transmu[t-jj];
                                tut = putildavec + t + 1;
                                ttmu += t+1-jj;
                                int i = s - t + 1;

                                #if defined(define_double) && defined(__AVX__)
                                        __m256d ymm;
                                        double atmp[4] __attribute__((aligned(32)));
                                        i = s - t;
                                        ymm = _mm256_setzero_pd();
                                        while(i & 3){
                                            t1d += (*tut) * (*ttmu);
                                            ++tut, ++ttmu, --i;
                                        }
                                        while(i > 0){
                                            ymm = _mm256_add_pd(ymm, _mm256_mul_pd(_mm256_loadu_pd(tut), _mm256_loadu_pd(ttmu)));
                                        // ymm = _mm256_fmadd_pd(_mm256_loadu_pd(tut), _mm256_loadu_pd(ttmu), ymm);
                                            tut += 4, ttmu += 4, i -= 4;
                                        }
                                        _mm256_store_pd(atmp, ymm);
                                        for(int j = 0; j < 4; ++j) t1d += atmp[j];
                                #else 
                                          while (--i) {
                                             t1d += (*tut)*(*ttmu);
                                             ttmu++;
                                             tut++;
                                          }
                                #endif


#ifdef define_svp
                                pyvec[t] = t1d;
#endif
#ifdef define_cvp
                                pyvec[t] = t1d- target_mud[t];
                                t1d = pyvec[t];
#endif

#if defined(define_double) && defined(__SSE4_1__)
                                __m128d tmp;
                                tmp = _mm_load_sd(&t1d);
                                tmp = _mm_round_sd(tmp, tmp, 0);
                                _mm_store_sd(&t1d, tmp);
                                t1d = -t1d;
#else
                                t1d = -round(t1d);
#endif
                                putildavec[t] = pvvec[t] = t1d;
                                pDeltavec[t] = 0;

                                if (putildavec[t] > -pyvec[t]) { 
                                    pdeltavec[t] = -1;
                                } else {
                                    pdeltavec[t] = 1;
                                }
#ifdef depthcost
                                mycostarray[t]++;
                                totaldepthcost += (s-t+1);     //cost for GS-coefficients 
#endif
                            }
                        }

                        if (cflag==0) {
                            t++;
                            pDeltavec[t] = -pDeltavec[t];
                            pDeltavec[t] += (1+((pDeltavec[t]*pdeltavec[t]) >> 31))*pdeltavec[t];
                            putildavec[t] = pvvec[t] + pDeltavec[t];
                            s = max(s,t);
                        }

                        //counter 
                        if ((++lmct)^lmctlimit) {
                        } else {
                            #pragma omp critical
                            {
                                gmct += lmct / 1000000;
                                lmct %= 1000000;
                                if (vl>=1) {
                                    ccout << "count=" << gmct*0.001 << "G tree=" << vfinished << "/" << numsubtree;
                                    ccout << " time="  << gettimeofday_sec() - es;
                                    ccout << " speed=" << 1.0*gmct / (gettimeofday_sec() - es) << "M/s  "; // slide=" << (int)slide;
                                    ccout << "\r";
                                    ccout.flush();
                                }
                                if (gmct /1000 > elim) {
                                    endflag = kk;
                               }
#ifdef depthcost
                                for (int i=0;i<s;i++) {
                                        EV.depthcostarray[i] += mycostarray[i];
                                }
                               int64_t totald=0;
                               //ccout << endl;
                               //ccout << "count(" << mythread << "): "; 
                                for (int i=0;i<s;i++) {
                                        EV.depthcostarray[i] += depthcostarray[mythread][i];
                                        totald += EV.depthcostarray[i];
                                        //ccout << EV.depthcostarray[i]  << " ";
                                        depthcostarray[mythread][i] = 0;
                                }
                                //ccout << "totald=" << totald << endl;
#endif
                            }
                            t +=  endflag;
                        }

                    } // End of ENUM-1b core
                    sos[vindex]=-1;
                    #pragma omp critical
                    {
                        vfinished++;
                        gmct += lmct / 1000000;
                        lmct %= 1000000;
                    }
                    
                }

                if (ea1flag<=1) {
                    //ENUM-a1 is not finished
                    vindex += parallel;
                } else {
                    //ENUM-a1 is finished
                    vindex++;
                }
            } while (vindex<numslot);

            //if ENUM1a was finished
            if (ea1flag==2) {
                if (finishmode==finish_mode_nowaste) {
                    //counting current "waiting" slot
                    int current_used_slot=0;
                    for (int i=0;i<numslot;i++) {
                        if (sos[i]!=-1) current_used_slot++;
                    }

                    if (current_used_slot==0) {
                        endflag = kk;
                    }
                }
                if (finishmode==finish_mode_exact) {
                    //counting current "empty, waiting and under processing" slot
                    int current_used_slot=0;
                    for (int i=0;i<numslot;i++) {
                        if (sos[i]!=-1) current_used_slot++;
                    }
                    if (current_used_slot==0) {
                        endflag = kk;
                    }
                }
            } 
        } while (endflag==0);
        endflag = kk;

        //summerize # processed nodes
        #pragma omp critical
        {
            gmct += 0.000001 * lmct;
            lmct = 0;
        }

}   //end of whole parallel section
#endif

            
            
        if ((enummode & enum_mode_all_vectors) && (vl>=1)) {
            //ccout << endl;
        }
            
#ifdef define_cvp
        if (flag_cvpzeroadd == true) {
            EV.duplicate();
        }         
        addbias(EV,cvpbias);
#endif
        //Count # found vector
        es = gettimeofday_sec() - es;
        if (vl>=3) {
            ccout << "#nodes=" << gmct / 1000.0<< "Gnodes ";
            ccout << "Time=" << es << " enum_speed=" << gmct / es << "Mnodes/sec ";
            ccout << "final_radius=" << sqrt(clim) << endl;
            ccout << "final_slide=" << slide << endl;
        }
        EV.totalnodes = gmct * 1000000.0;
        EV.etime = es;
        current_totalnodes = EV.totalnodes;
        current_etime = EV.etime;
        current_cputime = getcputime() - cts;


        
#ifdef depthcost
        long int totalnodes=0;
        for (int i=0;i<jj+bs+1;i++) {
            //ccout << "depth[" << i << "]=" << depthcostarray[i] << endl;
#ifdef define_multi
            for (int j=0;j<parallel;j++) {
                EV.depthcostarray[i] += depthcostarray[j][i];
                totalnodes += depthcostarray[j][i];
            }

#else
            EV.depthcostarray[i] += depthcostarray[0][i];
            totalnodes += depthcostarray[0][i];
#endif
        }

#ifdef define_multi
            for (int j=0;j<parallel;j++) {
                delete [] depthcostarray[j];
            }
            delete [] depthcostarray;
#else
            delete [] depthcostarray[0];
            delete [] depthcostarray;
#endif

        for (int i=1;i<(jj+bs+1)/2;i++) {
            //swap(EV.depthcostarray[i],EV.depthcostarray[jj+bs+1-i]); //inverse
        }
        //ccout << "totalnodes=" << totalnodes << endl;
        //ccout << "final_totaldepthcost=" << totaldepthcost << endl; 
        EV.insidenodes = gmct_inside * 1000000.0;
#endif
    }

#undef define_multi
#undef define_cvp
#undef define_svp

#ifdef define_double
    #undef PFLOAT
    #undef PSHORT
#endif

#undef define_double
