/*
	progressive BKZ library by NICT security fundemental lab.
	https://www2.nict.go.jp/security/pbkzcode/
	contact email:pbkz-info@ml.nict.go.jp
 */

#ifndef _inc_pbkz_hpp
#define _inc_pbkz_hpp

//configures

int debug_output=0;

//#define depthcost
//Note: if you want to get some statistics, enable this option
//it counts number of nodes precisely, but slows down the enumeration


#define _include_svpchallenge_generator
//If you want to use SVP challenge generator

#define _include_idealchallenge_generator
//If you want to use Ideal-SVP challenge generator

#define _allow_cachefiles
//If you want to speed up the computation by using the cache files
//Cache file directory is defined in bkz.conf

#define _ibeta_approx_wrapper
//Note: with this option, in pbkzsimboost.cppI use an approximation procedure to compute incomplete beta function 
//          instead of boost's ibeta_inv() when x<10^(-20) because it sometimes returns "exceeds 200 iterations" error (boost 1_64_0)

#define _use_mpfr_pfunc
//Note: Use mpfr multiprecision library in the pruning function (pfuncwrapper.cpp) to manipulate very large dimensions

#include "libraries.hpp"

typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<20> > bkzfloat;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<10> > float10;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<15> > float15;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<20> > float20;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<25> > float25;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<30> > float30;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<35> > float35;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<40> > float40;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<45> > float45;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<50> > float50;

#ifndef __lite__compile
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<80> > float80;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<100> > float100;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<150> > float150;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<200> > float200;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<300> > float300;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<400> > float400;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<500> > float500;
#endif

#ifdef debug_display_info
    #define debug_display(X)    if (debug_output!=0) { X };
#else
    #define debug_display(X)    {};
#endif

        
//Fundemental tools
#include "ccout.cpp"
#include "stringtools.cpp"        
#include "timetools.cpp"        
#include "filetools.cpp"        
#include "memorytools.cpp"
#include "latticetools.cpp"        
#include "misc.cpp"

#include <lattice/gen_uni_mat.cpp>

//Challenge instance generator
//from https://www.latticechallenge.org/svp-challenge/download/generator.zip
// and https://www.latticechallenge.org/ideallattice-challenge/download/generator.zip
#ifdef _include_svpchallenge_generator
    #include "../external/tools.h"
#endif
#ifdef _include_idealchallenge_generator
    #include <NTL/ZZX.h>
    //The function set(X) conflicts with boost library
    #define set NTL::set
    #include "../external/ideal.h"
    #undef set
#endif

#if defined(_include_svpchallenge_generator) || defined(_include_idealchallenge_generator)
    #include <lattice/genlattice.cpp>
#endif

//To define verbose level (for readability in function call)
#define VL0 0
#define VL1 1
#define VL2 2
#define VL3 3
#define VL4 4
#define VL5 5
#define VL6 6
#define VL7 7


#include <lattice/bkzlibwrapper.hpp>
#include <lattice/genwrapper.cpp>
#include <lattice/bkzconstants.hpp>
#include <lattice/pfuncwrapper.cpp>
#include <lattice/enumwrapper.cpp>
#include <lattice/reductionwrapper.cpp>
#include <lattice/randomizewrapper.cpp>

#ifndef __no__bkz
    #include <lattice/pbkzsimboost.cpp>
#endif 
#include <lattice/pbkzwrapper.cpp>
#ifndef __no__bkz
    #include <lattice/pbkzsimtimeboost.cpp>
#endif



#endif
