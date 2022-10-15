//
//  Used external libraries
//

//boost library
//#include <boost/functional/hash.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#ifdef _use_mpfr_pfunc
    #include <boost/multiprecision/mpfr.hpp>
#endif
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/cpp_int.hpp> 
#include <boost/version.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/functional/hash.hpp>
#include <boost/multiprecision/random.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/math/distributions.hpp>
#include <boost/algorithm/string.hpp>


using namespace boost::multiprecision;


//system and STL
#include <cstdlib>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>     
#include <algorithm>
#include <functional>        
#include <map>
#include <errno.h>
#include <sys/file.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <time.h>
#include <random>
#include <string>
#include <vector>
#include <map>

using namespace std;

//NTL
#include <NTL/ZZ.h>
#include <NTL/matrix.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/LLL.h>

//gsl
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//misc
#include <sys/ioctl.h>

NTL_CLIENT

#ifndef __no__multithreads
    //OpenMP
    #include <omp.h>
#endif
        
