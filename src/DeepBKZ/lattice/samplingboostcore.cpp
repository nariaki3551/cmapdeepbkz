//This file defines core subroutines of sampling_tools
//for even and odd dimensions

//SFLOAT: type for sampling
//DFLOAT: type of input data

//Approximating the ratio between cylinder intersections by random point sampling
#ifdef define_even
      template <typename SFLOAT,typename DFLOAT> DFLOAT ApproxRatio_even(std::vector<DFLOAT>& rd,int istart,int iend,int numsamples,int parallel) {
#endif
#ifdef define_odd
      template <typename SFLOAT,typename DFLOAT> DFLOAT ApproxRatio_odd(std::vector<DFLOAT>& rd,int istart,int iend,int numsamples,int parallel) {
#endif

        using namespace sampling_tools;
          
        long long int ktotal=0;
        long long int ktotal2=0;
        int endflag=0;

        int dim = iend-istart+1;
#ifdef define_even         
        int ds = dim/2; //dimension of sampling vector
#endif
#ifdef define_odd         
        int ds = (dim+1)/2; //dimension of sampling vector
#endif

        //Below, shifted parameters are defined in [1..dim] or [1..ds]
        SFLOAT* eR = (SFLOAT*)shared_memory::allocate1<SFLOAT>(1701,2+dim/2);
        SFLOAT* eL = (SFLOAT*)shared_memory::allocate1<SFLOAT>(1702,2+dim/2);
#ifdef define_even         
        for (int i=1;i<=ds;i++) eR[i] = (SFLOAT)rd[istart+2*i-1];
#endif
#ifdef define_odd
        for (int i=1;i<ds;i++) eR[i] = (SFLOAT)rd[istart+2*i-1];
        eR[ds] = 1.0;
#endif
        
        for (int i=1;i<=ds;i++) eL[i] = 0;

        SFLOAT* rdtemp = (SFLOAT*)shared_memory::allocate1<SFLOAT>(1703,1+dim);
        for (int i=1;i<=dim;i++) rdtemp[i] = (SFLOAT)rd[istart+i-1];
        
        
        SFLOAT** yy = (SFLOAT**)shared_memory::allocate2<SFLOAT>(1704,parallel,4+dim);
        SFLOAT** yyc = (SFLOAT**)shared_memory::allocate2<SFLOAT>(1705,parallel,4+dim);
        
        
        #pragma omp parallel  num_threads(parallel) 
        {

#ifdef __no__multithreads
            int mythread =  0;
#else
            int mythread =  omp_get_thread_num();
#endif
            SFLOAT* y = yy[mythread];
            SFLOAT* yc = yyc[mythread];
            
            //Initial point
            for (int i=1;i<=ds;i++) y[i] = eR[i] / 100.0;
            
            //Shuffle
            int j = 0;
            while ((j<100*ds) && (endflag==0)) {
                    sampling_tools::SamplefromDKLR<SFLOAT>(y,eL,eR,ds,ds,mythread);
                    j++;
            }
 
            //main sampling
            int klocal=0;
            int klocal2 = 0;
            while (endflag==0) {
                SamplefromCpKLR<SFLOAT>(y,yc,eL,eR,ds,mythread);

 #ifdef define_odd
                SFLOAT mult = 1.0 / sqrt(1.0-yc[dim-1]*yc[dim-1]); 
                for (j=1;j<dim;j++) {
                    yc[j] *= mult;
                }
#endif
                if (checkR(yc,rdtemp,dim-1)==true) klocal++;
                klocal2++;
                if ((klocal2 & 0x100) !=0) {
                #pragma omp critical
                    {
                        ktotal2 += klocal2;
                        klocal2 = 0;
                        ktotal += klocal;
                        klocal = 0;
                        if (ktotal2>numsamples) endflag = 1;
                    }                    
                }
            }
            #pragma omp atomic  
            ktotal += klocal;
            #pragma omp atomic  
            ktotal2 += klocal2;
        }   // end of omp parallel
        return 1.0 * ktotal / ktotal2;
}

