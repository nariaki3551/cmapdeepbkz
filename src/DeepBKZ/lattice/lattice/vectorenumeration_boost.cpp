#ifndef _inc_vector_enumeration_boost
#define _inc_vector_enumeration_boost

//short and close vector enumeration 

#define modesvp 0x01
#define modecvp 0x02

//Definitions to control the behaviours of enumeration
#define enum_mode_find_shortest 0x01
//update pruning radius during enumeration

#define enum_mode_all_vectors 0x02  
//enumerate all vectors s.t. |v|<bound and under the pruning function

#define enum_mode_find_abort 0x04
#define enum_mode_find_then_abort 0x04
//Abort the enumeration when a vector is found

#define enum_mode_find_short_projection 0x08 
//Find vectors s.t. \pi_(v) are small

#define enum_mode_except_trivial 0x10
//Exclude trivial vectors (b1,...,bn) from outputs

#define enum_mode_count_only 0x20
//Counting short vectors

#define finish_mode_nowaste 0x01
//In multithread-mode, when a non-working thread exists, then halt the subroutine

#define finish_mode_exact 0x02
//In multithread-mode, wait all working threads has been done.
//There might exist idle threads at final in computation

namespace lattice_enum {

    double current_totalnodes;
    double current_etime;
    double current_cputime;

    template <typename T,typename T2> void storetomemory(FoundENUMVectors& EV,T2* utildavec, T ctilda, int jj, int kk,int enummode,int undupflag = 0) {
        //Store the data of found vector to cache memory
        if (enummode==enum_mode_find_short_projection) {
            //If the vector is zero or b1, no vector is stored
            int zeroflag=0;
            for (int i=jj+1;i<=kk;i++) {
                if (utildavec[i]!=0) {
                    zeroflag = 1;
                }
            }
            if (zeroflag==0) return;
        }

        FoundVectorInfo FVtemp;
        FVtemp.norm = (double)ctilda;
        FVtemp.projection = jj;
        FVtemp.coeffs.resize(kk-jj+1);
        for (int i=jj;i<=kk;i++) {
            FVtemp.coeffs[i-jj] = (int)utildavec[i];
        }
        EV.storeelement(FVtemp,undupflag);
    }
    
    template <typename T> void addbias(FoundENUMVectors& EV,T* cvpbias) {
        //Add the coefficient from target vector in close vector mode
        for (int i=0;i<EV.data.size();i++) {
            int jj = EV.data[i].projection;
            for (int j=0;j<EV.data[i].coeffs.size();j++) EV.data[i].coeffs[j] += cvpbias[j+jj];
        }
    }



#define define_single_enum    
#include "vectorenumeration_boostcore.cpp"
#undef define_single_enum    

    
#ifndef __no__multithreads
#define define_multi_enum    
#include "vectorenumeration_boostcore.cpp"
#undef define_multi_enum   
#endif
    
#ifndef __lite__compile
    #define define_single_enum_double    
    #include "vectorenumeration_boostcore.cpp"
    #undef define_single_enum_double    

    #define define_multi_enum_double    
    #include "vectorenumeration_boostcore.cpp"
    #undef define_multi_enum_double   

    #define define_single_enumcvp_double    
    #include "vectorenumeration_boostcore.cpp"
    #undef define_single_enumcvp_double    

    #define define_multi_enumcvp_double   
    #include "vectorenumeration_boostcore.cpp"
    #undef define_multi_enumcvp_double    

    #define define_single_enumcvp    
    #include "vectorenumeration_boostcore.cpp"
    #undef define_single_enumcvp    

    #define define_multi_enumcvp    
    #include "vectorenumeration_boostcore.cpp"
    #undef define_multi_enumcvp    

        //Below is experimental subroutine, find vectors whose projective lengths and combination coefficients are small
    #define define_single_enum_coeff_double    
    #include "vectorenumeration_boostcore.cpp"
    #undef define_single_enum_coeff_double    

#endif



    
} // end of namespace lattice_enum

#endif
