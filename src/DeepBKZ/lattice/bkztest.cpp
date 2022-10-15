#ifndef  _inc_bkztest_cpp
#define _inc_bkztest_cpp

#include "lattice/bkzlibwrapper.hpp"
#include "lattice/pfuncwrapper.cpp"
#include "lattice/bkzconstants.hpp"
#include "lattice/genwrapper.cpp"
#include "lattice/pbkzsimtimeboost.cpp"

#define ccout_separate ccout_centertitle("-","");

void ccout_centertitle(std::string filler,std::string title) {
    struct winsize w;
    ioctl(0,TIOCGWINSZ,&w);
    int vlen = w.ws_col;
    int flen = vlen - title.length();
    if (flen < 0) {
        ccout << title << endl;
        return;
    } else {
        std::string disp ;
        int ff = flen / 2;
        for (int i=0;i<ff;i++) {
            flen--;
            disp += filler[i % filler.length()];
        }
        disp += title;
        for (int i=0;i<flen;i++) {
            disp += filler[i % filler.length()];
        }
        ccout << disp << endl;
    }
}

void ccout_title(std::string title) {
    ccout_centertitle("-","");
    ccout_centertitle("-","");
    ccout_centertitle("-",title);
}

void ccout_subtitle(std::string title) {
    ccout_centertitle("-",title);
}

void announcement() {
    
    ccout_title("Progressive BKZ library test (version 202205)");
    
    cout << "Notes:" << endl;;
    cout << "1. It takes about 30-60 minutes." << endl;
    cout << "2. It creates cache files in directory described in bkz.conf" << endl;
    cout << "3. It requires challenge-600 and LWE_40_005.txt from Lattice Challenge" << endl;
    cout << "   and LWE Challenge for complete test" << endl;
    cout << endl;
    cout << "The latest version is published at https://www2.nict.go.jp/security/pbkzcode/" << endl;
    cout << "Contact mailing list: pbkz-info@nict.go.jp" << endl;
}

bool basic_test() {

    ccout_title("PBKZlib test 1: basic subroutines");

    ccout_subtitle("Read setting file: ");
    if (FileExists("bkz.conf")==false) {
        cout << "bkz.conf not exist in the current directory." << endl;
        return false;
    }

    std::vector<std::string> cache;
    cache.push_back("constantscache");
    cache.push_back("pfcache");
    cache.push_back("simcache");
    cache.push_back("svpccache");

    for (int i=0;i<cache.size();i++) {
        std::string filename;
        filename = ReadConf("bkz.conf",cache[i]);
        cout << "dir. of " << cache[i] << "=" << filename << endl;
        if (filename == "") {
            cout << "directory is empty" << endl;
            return false;
        }
    }    
    
    bkzconstants::initialize();
    cout << "full-path of cache file directory: " << bkzconstants::cachedir << endl;

    cout << "size of double=" << sizeof(double) << endl;
    cout << "size of long double=" << sizeof(long double) << endl;
    cout << "maximum number of threads: " << omp_get_num_procs() << endl;

    ccout_separate
    ccout_subtitle( "test for memory allocation");
    for (int i=0;i<10;i++) {
        int id = i;
        int size = 100;
        double* a = (double*)shared_memory::allocate1<double>(id,size);
        for (int j=0;j<size;j++) a[j] = j+ id + size + 101;
        int* b = (int*)shared_memory::allocate1<int>(id,size);
        for (int j=0;j<size;j++) b[j] = j + id + size + 102;
    }

    ccout_separate
    ccout_subtitle("Generating SVP challenge basis: ");
    LatticeBasis<double> dB;
    for (int dim=50;dim<=200;dim+=50) {
        int seed = 0;
        dB = svpc::getbasis(dim,seed);
        cout << "dim=" << dim << " seed=" << seed << endl;
        cout << "B[1,1]=" << dB.L[0][0] << endl;
        if (dB.L[0][0] == 0) return false;
    }
    
    ccout_separate
    ccout_subtitle( "Save and Load: ");

    std::string filename;
    filename = "svpc_" + to_stdstring(dB.dim);
    cout << "save lattice basis to " + filename << endl;
    SaveLattice(dB,filename);
    LatticeBasis<long double> lB;

    cout << "load lattice basis to " + filename << endl;
    cout << "matrix element check" << endl;
    LoadLattice(lB,filename);
    for (int i=0;i<dB.dim;i++) {
        for (int j=0;j<dB.dim;j++) {
            if (dB[i][j] - lB[i][j] !=0) {
                cout << "matrix error: " << i << "," << j << endl;
                cout << dB[i] << " : " << lB[i] << endl;
                return false;
            }
        }
    }
    
    int prec_back = cout.precision();
    std::cout << std::setprecision(12);
    ccout_separate
    ccout_subtitle("test for constant table");

    cout << "volume and surface of n-unit balls" << endl;
    for (int i=10;i<=100;i+=10) {
        cout << "dim=" << i << ": volume=" << bkzconstants::vol_unit_ball(i) << "   surface=" << bkzconstants::surface_unit_ball(i) << endl;
    }

    cout << "gaussian heuristic constants: V_n(1)^(-1/n)" << endl;
    for (int i=2;i<=1000;i*=2) {
        cout << "dim=" << i << " GH=" << bkzconstants::ghconstant(i) << endl;
    }
    
    for (int i=10;i<=100;i+=10) {
        cout << "factorial(" << i << ")=" << bkzconstants::factorial(i) << endl;
    }

    cout << "Binary coefficients:" << endl;
    for (int i=10;i<=100;i+=10) {
        for (int j=5;j<=i/2;j+=5) {
            cout << "(" << i << " choose " << j << ")=" << bkzconstants::binary_coeff(i,j) << endl;
        }
    }
    
    ccout_separate
    ccout_subtitle( "test for sampling tools");
    std::cout << std::setprecision(8);

    cout << "random: ";
    for (int i=0;i<10;i++) {
        cout << sampling_tools::random<double>() << " "; //generating random number in [0,1]
    }
    cout << endl;

    cout << "normal: ";
    for (int i=0;i<10;i++) {
        cout << sampling_tools::normal<double>() << " "; //generating normal distribution
    }
    cout << endl;

    int num = 0;
    for (int i=0;i<100000;i++) {
        if (sampling_tools::random<double>() <  0.1 ) num++;
    }
    cout << "Pr[random < 0.1] =" << 1.0 * num / 100000 << endl;

    num = 0;
    for (int i=0;i<100000;i++) {
        if (sampling_tools::normal<double>() >  0.5 ) num++;
    }
    cout << "Pr[normal > 0.5] =" << 1.0 * num / 100000 << endl;
    
    
    std::cout << std::setprecision(prec_back);
    return true;

}
    
bool gs_test() {
    
    bkzconstants::initialize();
    std::cout << std::setprecision(8);

    ccout_title("PBKZlib test 2: Gram-Schmidt subroutines");
    
    LatticeBasis<bkzfloat> fB;
    int dim = 100;

    ccout_separate
    ccout_subtitle("Generate random unimodular matrix dim=" + to_stdstring(dim));
    mat_ZZ U;
    gen_random_unimodular2(U,dim,0,10,VL0);
    cout << "U[0]=" << U[0] << endl;
    
    fB = U;
    fB.updateGSBasis();
    
    cout << "list of |b*i|: ";
    fB.gs.displayc();

    cout << "Last line of mu-matrix:";
    fB.gs.displaymuline(dim);
    
    ccout_separate
    ccout_subtitle("Check size reduction");

    for (int i=1;i<=dim;i++) SizeReduce(fB,0,i);
    cout << "Last line of mu-matrix after size reduce:";
    fB.gs.displaymuline(dim);
    for (int i=1;i<=dim-1;i++) {
        if (abs(fB.gs.mu[dim][i]) > 0.5) {
            cout << "GS-coefficient error: mu[" << i << "]=" << fB.gs.mu[dim][i] << endl;
        }
    }
    return true;
}

bool pftable_test() {
    ccout_title("BKZlib test3: precomputed table of pruning coefficients");

    PruningFunction PF;
    std::vector<bkzfloat> pt;
    pt.push_back(0.1);
    pt.push_back(1e-4);
    pt.push_back(1e-12);
    
    int parallel = omp_get_num_procs();

    for (int dim=60;dim<=140;dim+=40) {
        for (int j=0;j<pt.size();j++) {
            bkzfloat ptarget = pt[j];
            double delta = 1.01;
            ccout_separate
            cout << "test: dim=" << dim << " ptarget=" << ptarget << " delta=" << delta << endl;
            pruning_func::fromtable(PF,ptarget,dim, delta);
            cout << "Rigid_Lower_Prob=" << pruning_func::Rigid_lower_prob(PF) << " ";
            cout << "Rigid_Upper_Prob=" << pruning_func::Rigid_upper_prob(PF) << " ";
            bkzfloat psucc = pruning_func::Approx_prob(PF,parallel,100);
            cout << "Approx_Prob=" << psucc << endl;
            if (psucc < 0.5 * ptarget) return false;
            if (psucc > 2.0 * ptarget) return false;
        }
    }
    return true;
}


bool lll_test() {

    ccout_title("BKZlib test4: LLL and pruning coefficients");
    int prec_back = cout.precision();
    std::cout << std::setprecision(8);

    
    int dim = 110;
    int seed = 10;
    LatticeBasis<double> dB;

    ccout_separate
    ccout_subtitle("Generate  LLL basis for svp challenge dim=" + to_stdstring(dim) + " seed=" + to_stdstring(seed));
    dB = svpc::getlllbasis(dim,seed);
    dB.updateGSBasis();

    cout << "first vector = " << dB.L[0] << endl;
    
    cout << "|b*i|: ";
    dB.gs.displayc();

    double gh = LatticeGH(dB);
    cout << "GaussianHeuristic(L)=" << gh << endl;

    bkzfloat fec = FullENUMCost<bkzfloat>(dB,gh);
    cout << "FullENUMCost(L)=" << fec << endl;

    bkzfloat prob=0.0003;
    bkzfloat radius = 1.1 * gh;
    PruningFunction PF;

    ccout_separate
    cout << "Pruning coefficient for L (radius=1.1GH,p=" << prob << ", even dimension, from table interpolation)" << endl;
    pruning_func::SetPruningFunction(dB,PF, radius , prob,pf_fromtable);  
    PF.display();
    cout << "Rigid_Lower_Prob=" << pruning_func::Rigid_lower_prob(PF) << endl;
    cout << "Rigid_Upper_Prob=" << pruning_func::Rigid_upper_prob(PF) << endl;
    cout << "Approx_Prob=" << pruning_func::Approx_prob(PF) << endl;
    cout << "Rigid_Upper_Cost=" << pruning_func::Rigid_upper_cost(PF,dB.gs.c,radius) << endl;
    
    ccout_separate
    cout << "Pruning coefficient for L (radius=1.1GH,p=" << prob << ", even dimension, from systematic upper bound)" << endl;
    pruning_func::SetPruningFunction(dB,PF, radius , prob, pf_systematic);  
    PF.display();
    cout << "Rigid_Lower_Prob=" << pruning_func::Rigid_lower_prob(PF) << endl;
    cout << "Rigid_Upper_Prob=" << pruning_func::Rigid_upper_prob(PF) << endl;
    cout << "Approx_Prob=" << pruning_func::Approx_prob(PF) << endl;
    cout << "Rigid_Upper_Cost=" << pruning_func::Rigid_upper_cost(PF,dB.gs.c,radius) << endl;
    
    ccout_separate
    cout << "Pruning coefficient for L (radius=1.1GH,p=" << prob << ", even dimension, from systematic upper bound + optimization)" << endl;
    pruning_func::SetPruningFunction(dB,PF, radius , prob,pf_crossentropy,1,VL2);  
    PF.display();
    cout << "Rigid_Lower_Prob=" << pruning_func::Rigid_lower_prob(PF) << endl;
    cout << "Rigid_Upper_Prob=" << pruning_func::Rigid_upper_prob(PF) << endl;
    cout << "Approx_Prob=" << pruning_func::Approx_prob(PF) << endl;
    cout << "Rigid_Upper_Cost=" << pruning_func::Rigid_upper_cost(PF,dB.gs.c,radius) << endl;

    ccout_separate
    cout << "Pruning coefficient for L (radius=1.1GH,p=" << prob << ", even dimension, from systematic upper bound + optimization + post process)" << endl;
    pruning_func::SetPruningFunction(dB,PF, radius , prob,pf_crossentropy_exact,1,VL2);  
    PF.display();
    cout << "Rigid_Lower_Prob=" << pruning_func::Rigid_lower_prob(PF) << endl;
    cout << "Rigid_Upper_Prob=" << pruning_func::Rigid_upper_prob(PF) << endl;
    cout << "Approx_Prob=" << pruning_func::Approx_prob(PF) << endl;
    cout << "Rigid_Upper_Cost=" << pruning_func::Rigid_upper_cost(PF,dB.gs.c,radius) << endl;
    cout << endl;
    cout << endl;


    dim++;
    cout << "Change dimension to " << dim << endl;
    ccout_separate
    cout << "Pruning coefficient for L (radius=1.1GH,p=" << prob << ", odd dimension, from table interpolation)" << endl;
    pruning_func::SetPruningFunction(dB,PF, radius , prob, pf_fromtable);  
    PF.display();
    cout << "Rigid_Lower_Prob=" << pruning_func::Rigid_lower_prob(PF) << endl;
    cout << "Rigid_Upper_Prob=" << pruning_func::Rigid_upper_prob(PF) << endl;
    cout << "Approx_Prob=" << pruning_func::Approx_prob(PF) << endl;
    cout << "Rigid_Upper_Cost=" << pruning_func::Rigid_upper_cost(PF,dB.gs.c,radius) << endl;
    
    ccout_separate
    cout << "Pruning coefficient for L (radius=1.1GH,p=" << prob << ", odd dimension, from systematic upper bound)" << endl;
    pruning_func::SetPruningFunction(dB,PF, radius , prob, pf_systematic);  
    PF.display();
    cout << "Rigid_Lower_Prob=" << pruning_func::Rigid_lower_prob(PF) << endl;
    cout << "Rigid_Upper_Prob=" << pruning_func::Rigid_upper_prob(PF) << endl;
    cout << "Approx_Prob=" << pruning_func::Approx_prob(PF) << endl;
    cout << "Rigid_Upper_Cost=" << pruning_func::Rigid_upper_cost(PF,dB.gs.c,radius) << endl;
    
    ccout_separate
    cout << "Pruning coefficient for L (radius=1.1GH,p=" << prob << ", odd dimension, from systematic upper bound + optimization)" << endl;
    pruning_func::SetPruningFunction(dB,PF, radius , prob,pf_crossentropy);  
    PF.display();
    cout << "Rigid_Lower_Prob=" << pruning_func::Rigid_lower_prob(PF) << endl;
    cout << "Rigid_Upper_Prob=" << pruning_func::Rigid_upper_prob(PF) << endl;
    cout << "Approx_Prob=" << pruning_func::Approx_prob(PF) << endl;
    cout << "Rigid_Upper_Cost=" << pruning_func::Rigid_upper_cost(PF,dB.gs.c,radius) << endl;

    ccout_separate
    cout << "Pruning coefficient for L (radius=1.1GH,p=" << prob << ", odd dimension, from systematic upper bound + optimization + post process)" << endl;
    pruning_func::SetPruningFunction(dB,PF, radius , prob,pf_crossentropy_exact);  
    PF.display();
    cout << "Rigid_Lower_Prob=" << pruning_func::Rigid_lower_prob(PF) << endl;
    cout << "Rigid_Upper_Prob=" << pruning_func::Rigid_upper_prob(PF) << endl;
    cout << "Approx_Prob=" << pruning_func::Approx_prob(PF) << endl;
    cout << "Rigid_Upper_Cost=" << pruning_func::Rigid_upper_cost(PF,dB.gs.c,radius) << endl;
    cout << endl;
    cout << endl;
    std::cout << std::setprecision(prec_back);
    return true;
}    
   

bool simplex_prec_test() {


    ccout_title("BKZlib test5: Precision test for pruning functions (double, long double, float*)");
    cout << "Check subroutine for computing success probability" << endl;
    cout << "Input coefficient R_i=1 for all i, expected output is prob=1" << endl;
    int prec_back = cout.precision();
    std::cout << std::setprecision(12);
    int dim = 10;

    do {
        std::vector<bkzfloat> rd;
        rd.resize(dim+1);
        for (int i=1;i<=dim;i++) rd[i] = 1.0;
        bkzfloat pp = EvenSimplexVolumeWrapperMain<bkzfloat,double>(rd,dim,opt_volume_prob);
        cout << " dim=" << dim << " result(double)=" << pp << endl;
        if ( abs(pp-1) > 1e-3) break;
    } while (dim++ < 1000);

    
    ccout_separate
    do {
        std::vector<bkzfloat> rd;
        rd.resize(dim+1);
        for (int i=1;i<=dim;i++) rd[i] = 1.0;
        bkzfloat pp = EvenSimplexVolumeWrapperMain<bkzfloat,long double>(rd,dim,opt_volume_prob);
        cout << " dim=" << dim << " result(long double)=" << pp << endl;
        if ( abs(pp-1) > 1e-3) break;
    } while (dim++ < 1000);
    
    ccout_separate
    do {
        std::vector<bkzfloat> rd;
        rd.resize(dim+1);
        for (int i=1;i<=dim;i++) rd[i] = 1.0;
        bkzfloat pp = EvenSimplexVolumeWrapperMain<bkzfloat,float10>(rd,dim,opt_volume_prob);
        cout << " dim=" << dim << " result(float10)=" << pp << endl;
        if ( abs(pp-1) > 1e-3) break;
    } while (dim++ < 1000);
 
/*
    //Note: below is used for high-precision experiments
    ccout_separate
    do {
        std::vector<bkzfloat> rd;
        rd.resize(dim+1);
        for (int i=1;i<=dim;i++) rd[i] = 1.0;
        bkzfloat pp = EvenSimplexVolumeWrapperMain<bkzfloat,float30>(rd,dim,opt_volume_prob);
        cout << " dim=" << dim << " result(float30)=" << pp << endl;
        if ( abs(pp-1) > 1e-3) break;
    } while (dim++ < 200);
 

    ccout_separate
    do {
        std::vector<bkzfloat> rd;
        rd.resize(dim+1);
        for (int i=1;i<=dim;i++) rd[i] = 1.0;
        bkzfloat pp = EvenSimplexVolumeWrapperMain<bkzfloat,float50>(rd,dim,opt_volume_prob);
        cout << " dim=" << dim << " result(float50)=" << pp << endl;
        if ( abs(pp-1) > 1e-3) break;
    } while (dim++ < 500);

    ccout_separate
    do {
        std::vector<bkzfloat> rd;
        rd.resize(dim+1);
        for (int i=1;i<=dim;i++) rd[i] = 1.0;
        bkzfloat pp = EvenSimplexVolumeWrapperMain<bkzfloat,float80>(rd,dim,opt_volume_prob);
        cout << " dim=" << dim << " result(float80)=" << pp << endl;
        if ( abs(pp-1) > 1e-3) break;
    } while (dim++ < 500);

    ccout_separate
    do {
        std::vector<bkzfloat> rd;
        rd.resize(dim+1);
        for (int i=1;i<=dim;i++) rd[i] = 1.0;
        bkzfloat pp = EvenSimplexVolumeWrapperMain<bkzfloat,float100>(rd,dim,opt_volume_prob);
        cout << " dim=" << dim << " result(float100)=" << pp << endl;
        if ( abs(pp-1) > 1e-3) break;
    } while (dim++ < 5000);

    ccout_separate
    do {
        std::vector<bkzfloat> rd;
        rd.resize(dim+1);
        for (int i=1;i<=dim;i++) rd[i] = 1.0;
        bkzfloat pp = EvenSimplexVolumeWrapperMain<bkzfloat,float150>(rd,dim,opt_volume_prob);
        cout << " dim=" << dim << " result(float150)=" << pp << endl;
        if ( abs(pp-1) > 1e-3) break;
        dim+=5;
    } while (dim < 5000);
    
    ccout_separate
    do {
        std::vector<bkzfloat> rd;
        rd.resize(dim+1);
        for (int i=1;i<=dim;i++) rd[i] = 1.0;
        bkzfloat pp = EvenSimplexVolumeWrapperMain<bkzfloat,float200>(rd,dim,opt_volume_prob);
        cout << " dim=" << dim << " result(float200)=" << pp << endl;
        if ( abs(pp-1) > 1e-3) break;
        dim+=5;
    } while (dim < 5000);

    ccout_separate
    do {
        std::vector<bkzfloat> rd;
        rd.resize(dim+1);
        for (int i=1;i<=dim;i++) rd[i] = 1.0;
        bkzfloat pp = EvenSimplexVolumeWrapperMain<bkzfloat,float300>(rd,dim,opt_volume_prob);
        cout << " dim=" << dim << " result(float300)=" << pp << endl;
        if ( abs(pp-1) > 1e-3) break;
        dim+=5;
    } while (dim < 5000);
    exit(0);
 */  

    ccout_subtitle("A bounding function of p approx 0.1");
    
    dim = 10;
    ccout_separate
    do {
        std::vector<bkzfloat> rd,rd2;
        rd2.resize(dim+1);
        pruning_func::SetApproximatePruningFunction(rd2,dim,(bkzfloat)1.02,(bkzfloat)0.1);
        rd.resize(1+dim/2);
        for (int i=1;i<=dim/2;i++) rd[i] = rd2[2*i-1];
        bkzfloat ppexact = EvenSimplexVolumeWrapperMain<bkzfloat,mpfr_float>(rd,dim/2,opt_volume_prob);
        bkzfloat pp = EvenSimplexVolumeWrapperMain<bkzfloat,double>(rd,dim/2,opt_volume_prob);
        cout << " dim=" << dim/2 << " exact=" << ppexact << " result(double)=" << pp << " error=" << abs(pp-ppexact)/pp << endl;
        if ( abs(pp-ppexact)/ppexact > 1e-3) break;
        dim+=2;
    } while (1);

    ccout_separate
    do {
        std::vector<bkzfloat> rd,rd2;
        rd2.resize(dim+1);
        pruning_func::SetApproximatePruningFunction(rd2,dim,(bkzfloat)1.02,(bkzfloat)0.1);
        rd.resize(1+dim/2);
        for (int i=1;i<=dim/2;i++) rd[i] = rd2[2*i-1];
        bkzfloat ppexact = EvenSimplexVolumeWrapperMain<bkzfloat,mpfr_float>(rd,dim/2,opt_volume_prob);
        bkzfloat pp = EvenSimplexVolumeWrapperMain<bkzfloat,long double>(rd,dim/2,opt_volume_prob);
        cout << " dim=" << dim/2 << " exact=" << ppexact << " result(long double)=" << pp << " error=" << abs(pp-ppexact)/pp << endl;
        if ( abs(pp-ppexact)/ppexact > 1e-3) break;
        dim+=2;
    } while (1);

    ccout_separate
    do {
        std::vector<bkzfloat> rd,rd2;
        rd2.resize(dim+1);
        pruning_func::SetApproximatePruningFunction(rd2,dim,(bkzfloat)1.02,(bkzfloat)0.1);
        rd.resize(1+dim/2);
        for (int i=1;i<=dim/2;i++) rd[i] = rd2[2*i-1];
        bkzfloat ppexact = EvenSimplexVolumeWrapperMain<bkzfloat,mpfr_float>(rd,dim/2,opt_volume_prob);
        bkzfloat pp = EvenSimplexVolumeWrapperMain<bkzfloat,float10>(rd,dim/2,opt_volume_prob);
        cout << " dim=" << dim/2 << " exact=" << ppexact << " result(float10)=" << pp << " error=" << abs(pp-ppexact)/pp << endl;
        if ( abs(pp-ppexact)/ppexact > 1e-3) break;
        
        dim+=2;
    } while (1);


    ccout_separate
    do {
        std::vector<bkzfloat> rd,rd2;
        rd2.resize(dim+1);
        pruning_func::SetApproximatePruningFunction(rd2,dim,(bkzfloat)1.02,(bkzfloat)0.1);
        rd.resize(1+dim/2);
        for (int i=1;i<=dim/2;i++) rd[i] = rd2[2*i-1];
        bkzfloat ppexact = EvenSimplexVolumeWrapperMain<bkzfloat,mpfr_float>(rd,dim/2,opt_volume_prob);
        bkzfloat pp = EvenSimplexVolumeWrapperMain<bkzfloat,float10>(rd,dim/2,opt_volume_prob);
        cout << " dim=" << dim/2 << " exact=" << ppexact << " result(float10)=" << pp << " error=" << abs(pp-ppexact)/pp << endl;
        if ( abs(pp-ppexact)/ppexact > 1e-3) break;
        
        dim+=2;
    } while (1);

    
/*
    //Note: below is used for high-precision experiments
   ccout_separate
    do {
        std::vector<bkzfloat> rd,rd2;
        rd2.resize(dim+1);
        pruning_func::SetApproximatePruningFunction(rd2,dim,(bkzfloat)1.02,(bkzfloat)0.1);
        rd.resize(1+dim/2);
        for (int i=1;i<=dim/2;i++) rd[i] = rd2[2*i-1];
        bkzfloat ppexact = EvenSimplexVolumeWrapperMain<bkzfloat,mpfr_float>(rd,dim/2,opt_volume_prob);
        bkzfloat pp = EvenSimplexVolumeWrapperMain<bkzfloat,float30>(rd,dim/2,opt_volume_prob);
        cout << " dim=" << dim/2 << " exact=" << ppexact << " result(float30)=" << pp << " error=" << abs(pp-ppexact)/pp << endl;
        if ( abs(pp-ppexact)/ppexact > 1e-3) break;
        dim+=2;
    } while (1);
 


    ccout_separate
    do {
        std::vector<bkzfloat> rd,rd2;
        rd2.resize(dim+1);
        pruning_func::SetApproximatePruningFunction(rd2,dim,(bkzfloat)1.02,(bkzfloat)0.1);
        //display_vector(rd2);
        rd.resize(1+dim/2);
        for (int i=1;i<=dim/2;i++) rd[i] = rd2[2*i-1];
        bkzfloat ppexact = EvenSimplexVolumeWrapperMain<bkzfloat,mpfr_float>(rd,dim/2,opt_volume_prob);
        bkzfloat pp = EvenSimplexVolumeWrapperMain<bkzfloat,float50>(rd,dim/2,opt_volume_prob);
        //bkzfloat ppexact =pp;
        cout << " dim=" << dim/2 << " exact=" << ppexact << " result(float50)=" << pp << " error=" << abs(pp-ppexact)/pp << endl;
        if ( abs(pp-ppexact)/ppexact > 1e-3) break;
        dim+=10;
    } while (1);

    ccout_separate
    do {
        std::vector<bkzfloat> rd,rd2;
        rd2.resize(dim+1);
        pruning_func::SetApproximatePruningFunction(rd2,dim,(bkzfloat)1.02,(bkzfloat)0.1);
        rd.resize(1+dim/2);
        for (int i=1;i<=dim/2;i++) rd[i] = rd2[2*i-1];
        bkzfloat ppexact = EvenSimplexVolumeWrapperMain<bkzfloat,mpfr_float>(rd,dim/2,opt_volume_prob);
        bkzfloat pp = EvenSimplexVolumeWrapperMain<bkzfloat,float100>(rd,dim/2,opt_volume_prob);
        cout << " dim=" << dim/2 << " exact=" << ppexact << " result(float100)=" << pp << " error=" << abs(pp-ppexact)/pp << endl;
        if ( abs(pp-ppexact)/ppexact > 1e-3) break;
        dim+=10;
    } while (1);
 
    ccout_separate
    do {
        std::vector<bkzfloat> rd,rd2;
        rd2.resize(dim+1);
        pruning_func::SetApproximatePruningFunction(rd2,dim,(bkzfloat)1.02,(bkzfloat)0.1);
        rd.resize(1+dim/2);
        for (int i=1;i<=dim/2;i++) rd[i] = rd2[2*i-1];
        bkzfloat ppexact = EvenSimplexVolumeWrapperMain<bkzfloat,mpfr_float>(rd,dim/2,opt_volume_prob);
        bkzfloat pp = EvenSimplexVolumeWrapperMain<bkzfloat,float150>(rd,dim/2,opt_volume_prob);
        cout << " dim=" << dim/2 << " exact=" << ppexact << " result(float150)=" << pp << " error=" << abs(pp-ppexact)/pp << endl;
        if ( abs(pp-ppexact)/ppexact > 1e-3) break;
        dim+=20;
    } while (1);

    ccout_separate
    do {
        std::vector<bkzfloat> rd,rd2;
        rd2.resize(dim+1);
        pruning_func::SetApproximatePruningFunction(rd2,dim,(bkzfloat)1.02,(bkzfloat)0.1);
        rd.resize(1+dim/2);
        for (int i=1;i<=dim/2;i++) rd[i] = rd2[2*i-1];
        bkzfloat ppexact = EvenSimplexVolumeWrapperMain<bkzfloat,mpfr_float>(rd,dim/2,opt_volume_prob);
        bkzfloat pp = EvenSimplexVolumeWrapperMain<bkzfloat,float200>(rd,dim/2,opt_volume_prob);
        cout << " dim=" << dim/2 << " exact=" << ppexact << " result(float200)=" << pp << " error=" << abs(pp-ppexact)/pp << endl;
        if ( abs(pp-ppexact)/ppexact > 1e-3) break;
        dim+=20;
    } while (1);

    ccout_separate
    do {
        std::vector<bkzfloat> rd,rd2;
        rd2.resize(dim+1);
        pruning_func::SetApproximatePruningFunction(rd2,dim,(bkzfloat)1.02,(bkzfloat)0.1);
        rd.resize(1+dim/2);
        for (int i=1;i<=dim/2;i++) rd[i] = rd2[2*i-1];
        bkzfloat ppexact = EvenSimplexVolumeWrapperMain<bkzfloat,mpfr_float>(rd,dim/2,opt_volume_prob);
        bkzfloat pp = EvenSimplexVolumeWrapperMain<bkzfloat,float300>(rd,dim/2,opt_volume_prob);
        cout << " dim=" << dim/2 << " exact=" << ppexact << " result(float300)=" << pp << " error=" << abs(pp-ppexact)/pp << endl;
        if ( abs(pp-ppexact)/ppexact > 1e-3) break;
        dim+=20;
    } while (1);
 */  
    

    ccout_separate
    ccout_subtitle("Precision test for pruning functions (double, long double, bkzfloat)");
    cout << "Check subroutine for computing success probability" << endl;
    cout << "Input coefficient is generated for probability = 10^(-24)" << endl;
    std::cout << std::setprecision(12);
    int usedprecision = 4;
    for (int dim=60;dim<=200;dim+=10) {
        std::vector<bkzfloat> rd,rd2;
        rd2.resize(1+dim);
        //pruning_func::SetApproximatePruningFunction(rd2,dim,(bkzfloat)1.02,(bkzfloat)1e-24);
        pruning_func::SetApproximatePruningFunction(rd2,dim,(bkzfloat)1.02,(bkzfloat)0.1);
        rd.resize(1+dim/2);
        for (int i=1;i<=dim/2;i++) rd[i] = rd2[2*i-1];

        bkzfloat ppexact = EvenSimplexVolumeWrapperMain<bkzfloat,mpfr_float>(rd,dim/2,opt_volume_prob);
        bkzfloat pp = EvenSimplexVolumeWrapperMain<bkzfloat,mpfr_float>(rd,dim/2,opt_volume_prob,1,usedprecision);
         cout << "prec=" << usedprecision << " dim=" << dim/2 << " result=" << pp << " exact=" << ppexact  << " error=" << abs(pp-ppexact)/pp << endl;
        if ( abs(pp-ppexact)/ppexact > 1e-3) {
            dim-=6;
            if (dim<10) dim = 10;
            usedprecision++;
        } else {
            //cout << "prec=" << usedprecision << " dim=" << dim << " result=" << pp << endl;
        }
    }
    std::cout << std::setprecision(prec_back);
    return true;
}

bool enum_test1() {
    bkzconstants::initialize();
    int prec_back = cout.precision();
    std::cout << std::setprecision(8);
    ccout_title("BKZlib test6A: enumeration subroutine (SVP, without pruning)");

    mat_ZZ VV;
    LatticeBasis<double> B;
    int dim = 44;

    B = svpc::getlllbasis(dim,0);
    double radius = 1.1 * LatticeGH(B);
    bkzfloat prob = 1.0;
    cout << "dim=" << dim << " radius=" << radius << endl;
    cout << "full enumeration (prec=double, mode=enum_mode_all_vectors)" << endl;
    VV = ENUM(B,radius,prob,enum_mode_all_vectors,0,VL0,"");
    cout << "#found=" << VV.NumRows() << endl;
    int numvectors = VV.NumRows();
    for (int i=VV.NumRows()-1;i>=1;i--) {
        if (LengthOf<double>(VV[i]) < LengthOf<double>(VV[i-1])) swap(VV[i],VV[i-1]);
    }
    double lambda1 = LengthOf(VV[0]);
    cout << "shortest length=" << lambda1 << endl;

    cout << "full enumeration (prec=double, mode=enum_mode_find_shortest)" << endl;
    VV = ENUM(B,radius,prob,enum_mode_find_shortest,0,VL0,"");
    cout << "#found=" << VV.NumRows() << endl;
    for (int i=VV.NumRows()-1;i>=1;i--) {
        if (LengthOf<double>(VV[i]) < LengthOf<double>(VV[i-1])) swap(VV[i],VV[i-1]);
    }

    if (LengthOf(VV[0]) != lambda1) {
        cout << "shortest length error: " << LengthOf(VV[0]) << " " << lambda1 << endl;
        cout << VV << endl;
        return false;
    }
    cout << "found shortest length=" << LengthOf(VV[0]) << endl;
    
    cout << "full enumeration (prec=double, mode=enum_mode_find_short_projection)" << endl;
    VV = ENUM(B,radius,prob,enum_mode_find_short_projection,0,VL0,"");
    cout << "#found=" << VV.NumRows() << endl;

    cout << "full enumeration (prec=double, mode=num_mode_find_abort)" << endl;
    VV = ENUM(B,radius,prob,enum_mode_find_abort,0,VL0,"");
    cout << "found_length=" << LengthOf(VV[0]) << " (shortest=" << lambda1 <<")" << endl;
 
    int parallel = omp_get_num_procs();
    cout << "max_parallel=" << parallel << endl;

    cout << "full enumeration (prec=double, mode=enum_mode_all_vectors, parallell)" << endl;
    VV = ENUM(B,radius,prob,enum_mode_all_vectors,0,VL0,"parallel=" + to_stdstring(parallel) + " finishmode=exact");
    cout << "#found=" << VV.NumRows() << endl;

    if (VV.NumRows() != numvectors) {
        cout << "Vector count error: " << VV.NumRows() << " " << numvectors << endl;
        cout << "radius=" << radius << endl;
        for (int i=0;i<VV.NumRows();i++) {
            cout << "Length[" << i+1 << "]=" << LengthOf(VV[i]) << endl;
        }
        cout << VV << endl;
        return false;
    }
  
    ccout_separate
    cout << "full enumeration for projected sublattice" << endl;
    dim = 30;
    int subdim = 20;
    int k = 2;
    cout << "dim=" << dim << endl;
    cout << "norm=sqrt(" << k << ")" << endl;

    LatticeBasis<long double> LS,LC; //lattice for SVP and CVP
    mat_ZZ diag;
    diag.SetDims(dim,dim);
    for (int i=0;i<dim;i++) {
        for (int j=0;j<dim;j++) {
            diag[i][j] = 0;
            if (i==j) {
                diag[i][i] = 1;
            }
        }
    }
    LS = diag;
    LS.updateGSBasis();
    double norm = 0;
    norm = sqrt((double)k) + 0.0001;

    int numfound = 0;
    int subnumfound = 0;
    
    for (int i=1;i<=k;i++) {
        numfound += boost::math::binomial_coefficient<double>(dim,i) * pow(2,i-1);
    }
    for (int i=1;i<=k;i++) {
        subnumfound += boost::math::binomial_coefficient<double>(subdim,i) * pow(2,i-1);
    }
    cout << "expected number of vectors (whole basis)=" << numfound << endl;
    cout << "expected number of vectors (sub " << subdim << "-dim basis)=" << subnumfound << endl;

    int istart = 1;
    int iend = subdim;

    cout << "full enumeration (single core, whole lattice basis)" << endl;
    VV = ENUM(LS,norm,1.0,enum_mode_all_vectors,0,VL0);
    cout << "# found vectors=" << VV.NumRows() << endl;
    if (VV.NumRows() != numfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }

    cout << "full enumeration (single core, whole lattice basis, parallel)" << endl;
    VV = ENUM(LS,norm+ 0.0001,1.0,enum_mode_all_vectors,0,VL0,"parallel=" + to_stdstring(parallel));
    cout << "# found vectors=" << VV.NumRows() << endl;
    if (VV.NumRows() != numfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }


    istart = 1; iend = subdim;

    cout << "full enumeration (single core, sub lattice basis)" << endl;
    VV = ENUM(LS,norm+ 0.0001,1.0,enum_mode_all_vectors,0,VL0,"istart="+to_stdstring(istart) + " iend=" + to_stdstring(iend));
    cout << "# found vectors=" << VV.NumRows() << endl;
    if (VV.NumRows() != subnumfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }
    istart = dim - subdim + 1;     iend = dim;
    VV = ENUM(LS,norm+ 0.0001,1.0,enum_mode_all_vectors,0,VL0,"istart="+to_stdstring(istart) + " iend=" + to_stdstring(iend));
    cout << "# found vectors=" << VV.NumRows() << endl;
    if (VV.NumRows() != subnumfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }

    
    istart = 1;    iend = subdim;
    cout << "full enumeration (parallel, sub lattice basis)" << endl;
    VV = ENUM(LS,norm+ 0.0001,1.0,enum_mode_all_vectors,0,VL0,"istart="+to_stdstring(istart) + " iend=" + to_stdstring(iend) + " parallel=" + to_stdstring(parallel));
    cout << "# found vectors=" << VV.NumRows() << endl;
    if (VV.NumRows() != subnumfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }
    istart = dim - subdim + 1;     iend = dim;
    VV = ENUM(LS,norm+ 0.0001,1.0,enum_mode_all_vectors,0,VL0,"istart="+to_stdstring(istart) + " iend=" + to_stdstring(iend) + " parallel=" + to_stdstring(parallel));
    cout << "# found vectors=" << VV.NumRows() << endl;
    if (VV.NumRows() != subnumfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }

    
    Randomize(LS,0,3,istart,subdim);
    ccout_separate
    cout << "full enumeration (single core, whole lattice basis, after randomization)" << endl;
    VV = ENUM(LS,norm,1.0,enum_mode_all_vectors,0,VL0);
    cout << "# found vectors=" << VV.NumRows() << endl;
    if (VV.NumRows() != numfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }

    ccout_separate
    cout << "full enumeration (parallel, whole lattice basis, after randomization)" << endl;
    VV = ENUM(LS,norm+ 0.0001,1.0,enum_mode_all_vectors,0,VL0,"parallel=" + to_stdstring(parallel));
    cout << "# found vectors=" << VV.NumRows() << endl;
    if (VV.NumRows() != numfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }

    istart = 1; iend = subdim;
    ccout_separate
    cout << "full enumeration (single core, sub lattice basis, after randomization)" << endl;
    VV = ENUM(LS,norm+ 0.0001,1.0,enum_mode_all_vectors,0,VL0,"istart="+to_stdstring(istart) + " iend=" + to_stdstring(iend));
    cout << "# found vectors=" << VV.NumRows() << endl;
    if (VV.NumRows() != subnumfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }
    istart = dim - subdim + 1;     iend = dim;
    VV = ENUM(LS,norm+ 0.0001,1.0,enum_mode_all_vectors,0,VL0,"istart="+to_stdstring(istart) + " iend=" + to_stdstring(iend));
    cout << "# found vectors=" << VV.NumRows() << endl;
    if (VV.NumRows() != subnumfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }

    
    istart = 1;    iend = subdim;
    ccout_separate
    cout << "full enumeration (parallel, sub lattice basis, after randomization)" << endl;
    VV = ENUM(LS,norm+ 0.0001,1.0,enum_mode_all_vectors,0,VL0,"istart="+to_stdstring(istart) + " iend=" + to_stdstring(iend) + " parallel=" + to_stdstring(parallel));
    cout << "# found vectors=" << VV.NumRows() << endl;
    if (VV.NumRows() != subnumfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }
    istart = dim - subdim + 1;     iend = dim;
    VV = ENUM(LS,norm+ 0.0001,1.0,enum_mode_all_vectors,0,VL0,"istart="+to_stdstring(istart) + " iend=" + to_stdstring(iend) + " parallel=" + to_stdstring(parallel));
    cout << "# found vectors=" << VV.NumRows() << endl;
    if (VV.NumRows() != subnumfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }
    
    std::cout << std::setprecision(prec_back);

    return true;
}

bool enum_test2() {

    bkzconstants::initialize();
    int prec_back = cout.precision();
    std::cout << std::setprecision(8);
    ccout_title("BKZlib test6B: enumeration subroutine (CVP, without pruning)");
    cout << "test for enumeration subroutine (CVP, without pruning)" << endl;

    mat_ZZ VV;
    LatticeBasis<double> LC;

    int dim = 30;
    int subdim = 20;
    int k = 1;
    
    cout << "dim=" << dim << " subdim=" << subdim << endl;
    
    mat_ZZ diag;
    diag.SetDims(dim,dim);
    for (int i=0;i<dim;i++) {
        for (int j=0;j<dim;j++) {
            diag[i][j] = 0;
            if (i==j) diag[i][i] = 10;
        }
    }
    LC = diag;
    LC.updateGSBasis();

    vec_ZZ tv;
    tv.SetLength(dim);
    for (int i=0;i<tv.length();i++) {
        tv[i] = 4.0 + i*10;
        if (i%2==0) tv[i] = 6.0+i*10;
    }
    cout << "target vector=" << tv << endl;

    double norm = 0.0001;
    double subnorm = 0.0001;
    for (int i=0;i<dim;i++) {
        if (i<k) {
                norm += 6.0*6.0;
                subnorm += 6.0*6.0;
        } else {
                norm += 4.0*4.0;
                if (i < subdim) subnorm += 4.0*4.0;
        }
    }

    norm = sqrt(norm);
    subnorm = sqrt(subnorm);
    
    
    cout << "norm=" << norm << " subnorm=" << subnorm << endl;

    int numfound = 0;
    int subnumfound = 0;
    for (int i=0;i<=k;i++) {
        numfound += boost::math::binomial_coefficient<double>(dim,i);
        subnumfound += boost::math::binomial_coefficient<double>(subdim,i);
    }
    
    int istart,iend;
    
    cout << "expected number of close vectors (whole basis)=" << numfound << endl;
    cout << "expected number of close vectors (sub " << subdim << "-dim basis)=" << subnumfound << endl;
    int parallel = omp_get_num_procs();
    
    ccout_separate
    cout << "full enumeration BDD (single core, whole lattice basis)" << endl;
    VV = ENUMCV(LC,tv,norm,1.0,enum_mode_all_vectors,0,VL0);
    cout << "# found vectors=" << VV.NumRows() << endl;
    //cout << VV << endl;

    if (VV.NumRows() != numfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }
   
    ccout_separate
    cout << "full enumeration BDD (parallel, whole lattice basis)" << endl;
    VV = ENUMCV(LC,tv,norm,1.0,enum_mode_all_vectors,0,VL0," parallel=" + to_stdstring(parallel));
    cout << "# found vectors=" << VV.NumRows() << endl;
    if (VV.NumRows() != numfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }

    istart = 1; iend = subdim;
    ccout_separate
    cout << "full enumeration BDD(single core, sub lattice basis)" << endl;
    VV = ENUMCV(LC,tv,subnorm,1.0,enum_mode_all_vectors,0,VL0,"istart="+to_stdstring(istart) + " iend=" + to_stdstring(iend));
    cout << "# found vectors=" << VV.NumRows() << endl;
    if (VV.NumRows() != subnumfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }
    istart = dim - subdim + 1;     iend = dim;
    VV = ENUMCV(LC,tv,subnorm,1.0,enum_mode_all_vectors,0,VL0,"istart="+to_stdstring(istart) + " iend=" + to_stdstring(iend));
    cout << "# found vectors=" << VV.NumRows() << endl;
    if (VV.NumRows() != subnumfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }

    istart = 1; iend = subdim;
    ccout_separate
    cout << "full enumeration BDD(parallel, sub lattice basis)" << endl;
    VV = ENUMCV(LC,tv,subnorm,1.0,enum_mode_all_vectors,0,VL0,"istart="+to_stdstring(istart) + " iend=" + to_stdstring(iend)+ " parallel=" + to_stdstring(parallel));
    cout << "# found vectors=" << VV.NumRows() << endl;
    if (VV.NumRows() != subnumfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }
    istart = dim - subdim + 1;     iend = dim;
    VV = ENUMCV(LC,tv,subnorm,1.0,enum_mode_all_vectors,0,VL0,"istart="+to_stdstring(istart) + " iend=" + to_stdstring(iend)+" parallel=" + to_stdstring(parallel));
    cout << "# found vectors=" << VV.NumRows() << endl;
    if (VV.NumRows() != subnumfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }
    
    Randomize(LC,0,3,istart,istart+5);
    
    ccout_separate
    cout << "full enumeration BDD (single core, whole lattice basis, after randomization)" << endl;
    VV = ENUMCV(LC,tv,norm,1.0,enum_mode_all_vectors,0,VL0);
    cout << "# found vectors=" << VV.NumRows() << endl;
    if (VV.NumRows() != numfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }
   
    ccout_separate
    cout << "full enumeration BDD (parallel, whole lattice basis, after randomization)" << endl;
    VV = ENUMCV(LC,tv,norm,1.0,enum_mode_all_vectors,0,VL0," parallel=" + to_stdstring(parallel));
    cout << "# found vectors=" << VV.NumRows() << endl;
    if (VV.NumRows() != numfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }

    istart = 1; iend = subdim;
    ccout_separate
    cout << "full enumeration BDD (single core, sub lattice basis, after randomization)" << endl;
    VV = ENUMCV(LC,tv,subnorm,1.0,enum_mode_all_vectors,0,VL0,"istart="+to_stdstring(istart) + " iend=" + to_stdstring(iend));
    cout << "# found vectors=" << VV.NumRows() << endl;
    if (VV.NumRows() != subnumfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }
    istart = dim - subdim + 1;     iend = dim;
    VV = ENUMCV(LC,tv,subnorm,1.0,enum_mode_all_vectors,0,VL0,"istart="+to_stdstring(istart) + " iend=" + to_stdstring(iend));
    cout << "# found vectors=" << VV.NumRows() << endl;
    if (VV.NumRows() != subnumfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }

    istart = 1; iend = subdim;
    ccout_separate
    cout << "full enumeration BDD (parallel, sub lattice basis, after randomization)" << endl;
    VV = ENUMCV(LC,tv,subnorm,1.0,enum_mode_all_vectors,0,VL0,"istart="+to_stdstring(istart) + " iend=" + to_stdstring(iend)+ " parallel=" + to_stdstring(parallel));
    cout << "# found vectors=" << VV.NumRows() << endl;
    if (VV.NumRows() != subnumfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }
    istart = dim - subdim + 1;     iend = dim;
    VV = ENUMCV(LC,tv,subnorm,1.0,enum_mode_all_vectors,0,VL0,"istart="+to_stdstring(istart) + " iend=" + to_stdstring(iend)+" parallel=" + to_stdstring(parallel));
    cout << "# found vectors=" << VV.NumRows() << endl;
    if (VV.NumRows() != subnumfound) {
        cout << "enum error: " << VV.NumRows() << " " << numfound << endl;
        return false;
    }
    std::cout << std::setprecision(prec_back);
    return true;
}


bool enum_test3() {
    bkzconstants::initialize();
    int prec_back = cout.precision();
    std::cout << std::setprecision(8);
    ccout_title("BKZlib test6C: enumeration subroutine (SVP, with pruning)");

    ccout_separate
    ccout_subtitle("Checking mild probability");

    LatticeBasis<double> dB;

    for (int dim = 48;dim<=60;dim+=4) {
        int seed = 0;
        dB = svpc::getlllbasis(dim,seed);
        double radius = 1.05 * LatticeGH(dB);
        bkzfloat prob = pow(1.05,-dim);
        ccout_separate
        cout << "find a vector shorter than 1.05GH: dim=" << dim <<  endl;
        cout << "GH=" << LatticeGH(dB) << " radius=" << radius << " p_succ=" << prob << endl;
        mat_ZZ VV;
        VV = ENUM(dB,radius,prob,enum_mode_all_vectors,0,VL0,"");
        cout << "#found=" << VV.NumRows() << endl;
        for (int i=0;i<VV.NumRows();i++) {
            cout << "vector=" << VV[i] << endl;
            cout << "length=" << LengthOf(VV[i]) << "=" << LengthOf(VV[i])/LatticeGH(dB) << "GH" << endl;
        }
    }

    ccout_separate
    ccout_subtitle("Checking low probability");

    int parallel = omp_get_num_procs()/2;
    LatticeBasis<long double> lB;
    for (int dim = 60;dim<=120;dim+=20) {
        int seed = 0;
        lB = svpc::getlllbasis(dim,seed);
        lB.updateGSBasis();
        double radius = 3.5 * LatticeGH(lB);
        bkzfloat prob = pow(3.5,-dim);
        ccout_separate
        cout << "find a vector shorter than 3.5GH: dim=" << dim <<  endl;
        cout << "GH=" << LatticeGH(lB) << " radius=" << radius << " p_succ=" << prob << endl;
        mat_ZZ VV;
        VV = ENUM(lB,radius,prob,enum_mode_find_abort,0,VL2,"parallel=" + to_stdstring(parallel) + " finishmode=exact");
        cout << "#found=" << VV.NumRows() << "                                                          " << endl;
        for (int i=0;i<min(3,VV.NumRows());i++) {
            cout << "vector=" << VV[i] << endl;
            cout << "length=" << LengthOf(VV[i]) << "=" << LengthOf(VV[i])/LatticeGH(lB) << "GH" << endl;
        }
    }
    
    std::cout << std::setprecision(prec_back);

    return true;
}

bool preprocess_test() {

    bkzconstants::initialize();

    int prec_back = cout.precision();
    std::cout << std::setprecision(8);

    int dim = 100;
    int seed = 0;

    ccout_title("BKZlib test7: Preprocessing subroutine for progressive BKZ dim=" + to_stdstring(dim));
    LatticeBasis<long double> B;
    B = svpc::getlllbasis(dim,seed);
    
    int subbeta = 80;
    int bstart = 10;

    int Ca = bkzconstants::get_Ca(subbeta);
    bkzfloat costlimit_lb = bkzconstants::target_cost_lb(subbeta) * Ca;
    cout << "costlimit_lb=" << costlimit_lb << endl;
    double falpha = 1.05;
    double ualpha = falpha * (1.0+subbeta)/subbeta;
    bkzfloat uprob = pow(falpha,-subbeta)/2;
    
    B.updateGSBasis();
    cout << "|b*i| before: ";
    B.gs.displayc(1);
    ProgressiveBKZPreprocess<long double,long double>(B,0,bstart,bstart+subbeta-1,costlimit_lb,falpha,ualpha,uprob,VL1);
    cout << "|b*i| after: ";
    B.gs.displayc(1);
    std::cout << std::setprecision(prec_back);
    
    return true;

}



bool reduction_test() {

    bkzconstants::initialize();
    int prec_back = cout.precision();
    std::cout << std::setprecision(8);
    ccout_title("BKZlib test8: progressive BKZ");

    LatticeBasis<double> B;
    vec_ZZ phi;
    int index = 172;
    int seed=0;


    ccout_separate
    ccout_subtitle("generate ideal lattice challenge instance index=" +to_stdstring(index) + " seed=" +to_stdstring(seed) );

    B = isvpc::getlllbasis(index,seed,phi);
    B.updateGSBasis();
    cout << "dimension=" << B.dim << endl;
    cout << " phi=" << phi << endl;
    cout << "|b*i|: ";
    B.gs.displayc();


    int istart,iend;
    double alpha = 1.1;
    bkzfloat prob = 0.1;
    
    ccout_separate
    ccout << "Enum-then-update subroutine, blocksize=40" << endl;
    for (int i=0;i<10;i++) {
        istart = 1 + 3*i;
        iend = istart + 40 -1;
        double radius = alpha * lattice_tools::LatticeGH(B.gs.c,istart,iend);
        EnumthenUpdate<double>(B,0,istart,iend,radius,prob,0,VL0,"");
        local_LLL(B,0,0.999);
    }
    cout << "|b*i|: ";
    B.gs.displayc();
     
    ccout_separate
    ccout << "Enum-then-update subroutine, blocksize=60" << endl;
    prob = 0.01;
    for (int i=0;i<10;i++) {
        istart = 1 + 2*i;
        iend = istart + 60-1;
        double radius = alpha * lattice_tools::LatticeGH(B.gs.c,istart,iend);
        EnumthenUpdate<double>(B,0,istart,iend,radius,prob,0,VL0,"");
        local_LLL(B,0,0.999);
    }
    cout << "|b*i|: ";
    B.gs.displayc();
    cout << "|b1|=" << LengthOf(B.L[0]) / LatticeGH(B) << endl;
    
    //index = 206;    //dim=102
    index = 119;    //dim=96
    //index = 141;    //dim=92
    
    seed = 0;
    LatticeBasis<long double> lB;
    lB = isvpc::getlllbasis(index,seed,phi);
    ccout_separate
    ccout_separate
    ccout_subtitle("Progressive BKZ for ideal lattice challenge instance index=" +to_stdstring(index) + " dim=" +to_stdstring(lB.dim) +" seed=" +to_stdstring(seed) );
    cout << "(It may take a few minutes|.)" << endl;
    cout << "After LLL: |b1|=" << LengthOf(lB.L[0]) / LatticeGH(lB) << endl;
    debug_output = 0;
    int parallel = omp_get_num_procs();
    double gh = LatticeGH(lB);

    parallel = omp_get_num_procs() /2;
    for (int i=0;i<20;i++) {
        bkzfloat prob;
        ccout_separate
        cout << "loop=" << i+1 << endl;
        ProgressiveBKZ<long double>(lB,0,50+i,VL1,"");
        bkzconstants::savecachetable();
        prob = pow(1.05,-lB.dim);
        double radius = 1.05 * LatticeGH(lB);
        bkzfloat bcost;
        do {
            prob *= 0.98;
            bcost = ENUMCostLB(lB.gs.c,1,lB.dim,radius,prob);
            cout << "prob=" << prob << " CostLB=" << bcost<< "     \r";
            cout.flush();
        } while (bcost > 1e+8);
        cout << endl;
        EnumthenUpdate<long double>(lB,0,1,lB.dim,radius,prob,0,VL1,"parallel=" + to_stdstring(parallel));
        if (LengthOf(lB.L[0]) / gh < 1.05) break; 
    } 
    cout << "found vector=" << lB.L[0] << endl;
    cout << "|b1|=" << LengthOf(lB.L[0]) / gh << "*GH" << endl;

    std::cout << std::setprecision(prec_back);
    
    return true;
}

bool simulator_test() {
    bkzfloat coeff_lll = 1.61028e-07;
    bkzfloat coeff_enum = 5.07206e-11;  //constants in BKZ simulator, not exactly ENUM speed
    double Loopratio = 1.42616;

    bkzconstants::initialize();
    std::cout << std::setprecision(8);
    ccout_title("BKZlib test9: PBKZ simulator");
    int dim;
    int sbeta ;
    int ebeta ;
    int usebeta;
/*

    ccout_separate
    ccout_subtitle("Progressive BKZ parameter table");
    for (int beta=10;beta<100;beta+=10) {
        double talpha,r;
        bkzfloat tprob;
        PBKZParam(talpha,r,tprob,beta);
        cout << "beta=" << beta << " r=" << r << " radius_alpha=" << talpha << " p_succ" << tprob << endl;
    }


    for (int beta=10;beta<=100;beta+=5) {
        cout << "target_cost_ub(" << beta << ")=" <<  bkzconstants::target_cost_ub_modify(beta) << "                           " << endl;
    }

    bkzconstants::savecachetable();

    ccout_separate
    ccout_subtitle("Target log-FEC table (calls subroutine to simulate |b*i| after BKZ-reduction)");
    dim = 200;
    for (int beta=20;beta<=60;beta+=10) {
        std::vector<bkzfloat> cc;
        PBKZSimulate(cc,dim,beta);
        cout << "simulate 200-dim after PBKZ-" << beta << endl;
        display_vector(cc,1);
        cout << "logfec[" << dim << "," << beta << "]=" << bkzconstants::simlogfec(dim,beta) << endl;
    }    
 
    ccout_separate
    ccout_subtitle("Timing simulation" );

    bkzfloat t_lll,t_enum;
    dim = 160;
    sbeta = 20;
    ebeta = 25;
    usebeta = 26;
    double loop;

    cout << "Simulate BKZ(" << sbeta << ") -> BKZ(" << ebeta << ") using PBKZ-" << usebeta << endl;
    bkzconstants::simtime(t_lll,t_enum,loop,dim,sbeta,ebeta,usebeta,opt_wantnumloop); 
    cout << "LLL-factor=" << t_lll << " ENUM-factor=" << t_enum << " expected # loops=" << loop * Loopratio<< endl;
    bkzconstants::savecachetable();

    dim = 200;  sbeta=30;   ebeta=60; usebeta=61;

    cout << "Simulate BKZ(" << sbeta << ") -> BKZ(" << ebeta << ") using PBKZ-" << usebeta << endl;
    bkzconstants::simtime(t_lll,t_enum,loop,dim,sbeta,ebeta,usebeta,opt_wantnumloop); 
    cout << "LLL-factor=" << t_lll << " ENUM-factor=" << t_enum << " expected # loops=" << loop * Loopratio << endl;
    bkzconstants::savecachetable();
    
    dim = 200;  sbeta=40;   ebeta=45;
    cout << "Simulate minimum cost of BKZ(" << sbeta << ") -> BKZ(" << ebeta << ")" << endl;
    cout << "Total cost is based on optimization cost estimated by lll_coeff*lll_factor+enum_coeff*enum_factor" << endl;


    cout << "lll_coeff=" << coeff_lll << " enum_coeff=" << coeff_enum << endl;
    bkzconstants::simtimeopt(t_lll,t_enum,usebeta,dim,sbeta,ebeta,coeff_lll,coeff_enum);
    cout << "LLL-factor=" << t_lll << " ENUM-factor=" << t_enum << " usebeta=" << usebeta << endl;
    bkzconstants::savecachetable();
    ccout_separate
 */
    dim = 200;  sbeta=20;   ebeta=45;
    BKZStrategy BS;
    
    ccout_separate
    ccout_subtitle("Find BKZ strategy: BKZ-" +to_stdstring(sbeta) + " -> " + to_stdstring(ebeta));
    cout << "(It may take a few minutes.)" << endl;
    gen_bkzstrategy(BS,dim,sbeta,ebeta,coeff_lll,coeff_enum,Loopratio);

    cout << "Optimized strategy: ";
    BS.display();
    bkzconstants::savecachetable();

    LatticeBasis<float10> B;
    ccout_separate
    B = svpc::getlllbasis(dim,10);
    cout << "progressive BKZ: LLL -> (25) -> 20                 " << endl;
    ProgressiveBKZ<float10>(B,0,20,VL1,"startbeta=25");

    cout << "Strategied BKZ" << endl;
    StrategiedBKZ<float10>(B,0,BS,VL1,"");
    
    return true;
    
}

bool speed_test() {

    bkzconstants::initialize();
    int prec_back = cout.precision();
    std::cout << std::setprecision(8);
    ccout_title("BKZlib test10: Speed test (multithreads)");
    cout << "max threads=" << omp_get_num_procs() << endl;
    cout << "used threads=" << omp_get_num_procs()/2 << endl;
            
    int dim = 100;
    LatticeBasis<double> B;
    B = svpc::getlllbasis(dim,0);
    double radius = LatticeGH(B);
    PruningFunction PF;
    pruning_func::SetPruningFunction(B,PF, radius , 0.001,0);  
    double elim;    //limiting giga-nodes
    FoundENUMVectors EV; 
    int parallel = omp_get_num_procs()/2;

    ccout_separate
    cout << "Benchmark for ENUM(single core,double precision, dim=" << dim << ")" << endl;
    elim = 0.4;
    double estart = gettimeofday_sec(); //timer
    lattice_enum::ENUMCoreboost_double(EV,0,1,dim,elim,radius*radius,PF.rd,B.gs.mu,B.gs.cd,VL0,enum_mode_count_only+enum_mode_find_shortest);
    double time = gettimeofday_sec() - estart;
    cout << "processed_nodes=" << EV.totalnodes << endl;
    cout << "elapsed_time=" << EV.etime << endl;
    cout << "speed=" << EV.totalnodes / EV.etime << "nodes/sec" << endl;

    ccout_separate
    cout << "Benchmark for ENUM(multicore(" << parallel << "),double precision, dim=" << dim << ")" << endl;
    elim = 0.2 * parallel;
    estart = gettimeofday_sec(); //timer
    lattice_enum::ENUMCoreboostmt_double(EV,0,1,dim,elim,parallel,radius*radius,PF.rd,B.gs.mu,B.gs.cd,VL0,enum_mode_count_only+enum_mode_find_shortest,finish_mode_exact);
    time = gettimeofday_sec() - estart;
    cout << "processed_nodes=" << EV.totalnodes << endl;
    cout << "elapsed_time=" << EV.etime << endl;
    cout << "speed=" << EV.totalnodes / EV.etime << "nodes/sec" << endl;


    LatticeBasis<long double> lB;
    lB = svpc::getlllbasis(dim,0);
    lB.updateGSBasis();
    ccout_separate
    cout << "Benchmark for ENUM(single core,long double precision, dim=" << dim << ")" << endl;
    elim = 0.2;
    estart = gettimeofday_sec(); //timer
    lattice_enum::ENUMCoreboost<long double,long double,long double>(EV,0,1,dim,elim,radius*radius,PF.rd,lB.gs.mu,lB.gs.cd,VL0,enum_mode_count_only+enum_mode_find_shortest);
    time = gettimeofday_sec() - estart;
    cout << "processed_nodes=" << EV.totalnodes << endl;
    cout << "elapsed_time=" << EV.etime << endl;
    cout << "speed=" << EV.totalnodes / EV.etime << "/sec" << endl;

    ccout_separate
    cout << "Benchmark for ENUM(multicore(" << parallel << "),long double precision, dim=" << dim << ")" << endl;
    elim = 0.1 * parallel;
    estart = gettimeofday_sec(); //timer
    lattice_enum::ENUMCoreboostmt<long double,long double,long double>(EV,0,1,dim,elim,parallel,radius*radius,PF.rd,lB.gs.mu,lB.gs.cd,VL0,enum_mode_count_only+enum_mode_find_shortest,finish_mode_exact);
    time = gettimeofday_sec() - estart;
    cout << "processed_nodes=" << EV.totalnodes << endl;
    cout << "elapsed_time=" << EV.etime << endl;
    cout << "speed=" << EV.totalnodes / EV.etime << "/sec" << endl;
    
    LatticeBasis<float15> fB;
    fB = svpc::getlllbasis(dim,0);
    fB.updateGSBasis();
    ccout_separate
    cout << "Benchmark for ENUM(single core,float15, dim=" << dim << ")" << endl;
    elim = 0.005;
    estart = gettimeofday_sec(); //timer
    lattice_enum::ENUMCoreboost<float15,float15,float15>(EV,0,1,dim,elim,radius*radius,PF.rd,fB.gs.mu,fB.gs.cd,VL0,enum_mode_count_only+enum_mode_find_shortest);
    time = gettimeofday_sec() - estart;
    cout << "processed_nodes=" << EV.totalnodes << endl;
    cout << "elapsed_ime=" << EV.etime << endl;
    cout << "speed=" << EV.totalnodes / EV.etime << "/sec" << endl;

    ccout_separate
    cout << "Benchmark for ENUM(multicore(" << parallel << "),float15, dim=" << dim << ")" << endl;
    elim = 0.003 * parallel;
    estart = gettimeofday_sec(); //timer
    lattice_enum::ENUMCoreboostmt<float15,float15,float15>(EV,0,1,dim,elim,parallel,radius*radius,PF.rd,fB.gs.mu,fB.gs.cd,VL0,enum_mode_count_only+enum_mode_find_shortest,finish_mode_exact);
    time = gettimeofday_sec() - estart;
    cout << "processed_nodes=" << EV.totalnodes << endl;
    cout << "elapsed_time=" << EV.etime << endl;
    cout << "speed=" << EV.totalnodes / EV.etime << "/sec" << endl;
    std::cout << std::setprecision(prec_back);
    return true;
    
}

bool simulate_svp_challenge() {

    int dim = 90;
    int seed=3;
    ccout_title("BKZlib test11: Time simulation and solve SVP challenge dim=" + to_stdstring(dim));
    
    bkzfloat coeff_lll = 1.61028e-07;
    bkzfloat coeff_enum = 5.07206e-11;  //constants in BKZ simulator, not exactly ENUM speed

        //alpha_enum=1.43128e-10       
        //alpha_LLL=4.34891e-07
    
    
    double Loopratio = 1.42616;

    bkzfloat enum_speed = 35000000; //benchmark at long double
    BKZStrategy BS;
    OptConvex<int,bkzfloat> gtable(dim*0.5,dim*0.9);

    BKZStrategy BSmin;
    bkzfloat Mmin;
    bkzfloat costmin=-1;
    int betamin;

    while (1) {
        int brange;
        int beta = gtable.next(brange);
        if (brange <= 1) break;    //optimal value found

        gen_bkzstrategy(BS,dim,20,beta,coeff_lll,coeff_enum,Loopratio);
        cout << "Time to find BKZ-" << beta << " basis in sec: " << BS.totalcost << "               " << endl;

        if ((costmin>0) && (BS.totalcost>costmin)) {
            gtable.add(beta,BS.totalcost);
        } else {
            std::vector<bkzfloat> cc;
            PBKZSimulate(cc,dim,beta);

            bkzfloat localmincost;

            bkzfloat M=1;   //number of used bases
            bkzfloat prob = pow(1.05,-dim) * 2 / M;
            bkzfloat radius = 1.05 * lattice_tools::LatticeGH(cc);
            bkzfloat cost,cost2;

            //ccout << "prob=" << prob << " radius=" << radius << endl;
            PruningFunction PF;
            OptConvex<bkzfloat,bkzfloat> costtable(1,1);

            cout << "beta=" << beta << endl;

            //M=1
            M=1;
            pruning_func::SetPruningFunction(cc,PF,1,dim,radius,prob,pf_crossentropy_exact,1,VL0);
            cost = pruning_func::Rigid_upper_cost(PF,cc,radius)/enum_speed +  BS.totalcost;
            costtable.add(M,cost);
            cout << "M=1 cost=" << cost << endl;

            M=2;
            prob = pow(1.05,-dim) * 2 / M;
            pruning_func::SetPruningFunction(cc,PF,1,dim,radius,prob,pf_crossentropy_exact,1,VL0);
            cost2 = pruning_func::Rigid_upper_cost(PF,cc,radius)/enum_speed*2 +  BS.totalcost*2;
            costtable.add(M,cost2);
            cout << "M=2 cost=" << cost2 << endl;

            if (cost<cost2) {
                //minimum is at M=1
                gtable.add(beta,cost);
                localmincost = cost;
            } else {

                do {
                    cost = cost2;
                    M*=2;
                    prob = pow(1.05,-dim) * 2 / M;
                    pruning_func::SetPruningFunction(cc,PF,1,dim,radius,prob,pf_crossentropy_exact,1,VL0);
                    cost2 = pruning_func::Rigid_upper_cost(PF,cc,radius)/enum_speed*M +  BS.totalcost*M;
                    costtable.add(M,cost2);
                    cout << "M=" << M << " cost=" << cost2 << endl;
                } while (cost > cost2);

                bkzfloat Mrange;
                do {
                    M = round(costtable.next(Mrange));
                    prob = pow(1.05,-dim) * 2 / M;
                    pruning_func::SetPruningFunction(cc,PF,1,dim,radius,prob,pf_crossentropy_exact,1,VL0);
                    bkzfloat ec = pruning_func::Rigid_upper_cost(PF,cc,radius);
                    cost2 = ec /enum_speed*M +  BS.totalcost*M;
                    costtable.add(M,cost2);
                    cout << "M=" << M << " tcost=" << cost2 << " enumcost=" << ec << endl;
                } while (Mrange>=2);
                localmincost = costtable.getmin();
                gtable.add(beta,localmincost);
            }
            
            if ((costmin==-1) || (localmincost < costmin)) {
                costmin = localmincost;
                BSmin = BS;
                betamin = beta;
                Mmin = costtable.getminkey();
            }
            bkzconstants::savecachetable();
        }        
    }

    cout << "global_min=" << gtable.getmin() << endl;
    cout << "Preprocess blocksize beta=" << betamin << endl;
    cout << "Strategy:" << endl;
    BSmin.display();
    //BSmin.savetofile("strategy_svp_challenge" + to_stdstring(dim));
    cout << "#bases=" << Mmin << endl;

    //Solve
    LatticeBasis<long double> B;
    B = svpc::getlllbasis(dim,seed);
    ProgressiveBKZ<long double>(B,0,15,VL1,"startbeta=20");

    bkzfloat prob = pow(1.05,-dim) * 2 / Mmin;
    bkzfloat radius = 1.05 * LatticeGH(B);
    do {
        StrategiedBKZ<long double>(B,0,BSmin,VL1,"");
        //SaveLattice(B,"temp");
        EnumthenUpdate<long double>(B,0,1,dim,radius.convert_to<double>(),prob,pf_crossentropy_exact,VL3,"");
        cout << "Found Length=" << B.gs.c[1] << "=" << B.gs.c[1]/LatticeGH(B) << "*GH(L)                " << endl;
        if (B.gs.c[1]<radius) break;
        //randomize
        seed++;
        Randomize(B,seed,10,1,dim);
        ::BigLLL(B.L,0,0.999,VL1);
        
    } while (1);
    
    cout << "find vector: " << B.L[0] << endl;
    cout << "Length=" << B.gs.c[1] << "=" << B.gs.c[1]/LatticeGH(B) << "*GH(L)" << endl;

    return true;
}
bool solve_lattice_challenge() {

    ccout_title("BKZlib test12: Solve LatticeChallenge");

    std::string fname = "challenge-600";
    if (FileExists(fname)==false) {
        cout << "file(" << fname << ") does not exist" << endl;
        cout << "skip this test" << endl;
        return false;
    }

    ccout_separate
    ccout_subtitle("Read file: " + fname);

    LatticeBasis<float20> B;
     std::ifstream ifs;
        ifs.open("challenge-600");
        int n;
        ifs >> n;
        int q;
        ifs >> q;
        int r;
        ifs >> r;
        cout << "n=" << n << " q=" << q << " r=" << r << endl;
        ifs >> B.L;
        //Reverse
        for (int i=0;i<n/2;i++) swap(B.L[i],B.L[n-i-1]);
        //cut
        int rank = 200;
        cout << "Extract last " << rank << " vectors" << endl;
        B.L.SetDims(rank,n);
        B.updateGSBasis();
        cout << "Table of |b*i|: ";
        B.gs.displayc();

        ccout_separate
        ccout_subtitle("Apply LLL");

        int swapcount=0;    //dummy
        local_LLL(B,0,0.999,1,rank,stnormal,VL1,swapcount);
        B.gs.displayc();

        ccout_separate
        ccout_subtitle("Apply PBKZ");
        ProgressiveBKZ(B,0,125,VL1,"startbeta=40 sbetashift=6 parallel=6 ignoreflat targetlen=" + to_stdstring(r-1));
        cout << "b1=" << B.L[0] << endl;
        cout << "|b1|=" << LengthOf(B.L[0]) << endl; 

        if (LengthOf(B.L[0]) < r) return true;
        return false;
}

void ExtractLWEChallengeFile(mat_ZZ& A,vec_ZZ& b,double& alpha,int& q,std::string filename) {

    ifstream ff;
    ff.open(filename.c_str(),ios_base::in);
    int n,m;
    ff >> n;
    ff >> m;
    ff >> q;
    //cout << "filename=" << filename.c_str() << endl;
    cout << "n=" << n << " m=" << m << " q=" << q << endl;
    ff >> alpha;
    cout << "alpha=" << alpha << endl;
    ff >> b;
    cout << "dim(b)=" << b.length() << endl;
    ff >> A;
    cout << "size(A)=" << A.NumRows() << "," << A.NumCols() << endl;
}


bool solve_lwe_challenge() {

    ccout_title("BKZlib test13: Solve LWE Challenge");
    std::string lwefile = "LWE_40_005.txt";
    if (FileExists(lwefile)==false) {
        cout << "file(" << lwefile << ") does not exist" << endl;
        cout << "skip this test" << endl;
        return true;
    }

    ccout_separate
    ccout_subtitle("Read file: " + lwefile);

    //Read file
    mat_ZZ lweA;
    vec_ZZ lweb;
    double lwealpha;
    int lweq;
    
    ExtractLWEChallengeFile(lweA,lweb,lwealpha,lweq,lwefile);

    int n = lweA.NumCols();
    int m = n * 2.5;
    cout << "subdim=" << m << endl;
    double sigma = lwealpha * lweq;
    cout << "sigma=" << sigma << endl;
    bkzfloat maxnorm = 1.5 * sqrt((double)m) * sigma;
    cout << "maxnorm of e=" << maxnorm << endl;
    
    vec_ZZ target;
    target.SetLength(m);
    for (int i=0;i<m;i++) target[i] = lweb[i];
    
    LatticeBasis<long double> A,B;
    A.resize(m+n,m);
    for (int i=0;i<m;i++) A.L[i][i] = lweq;
    for (int i=0;i<n;i++) {
        for (int j=0;j<m;j++) A.L[i+m][j] = lweA[j][i];
    }

    ccout_separate
    ccout_subtitle("Find basis vectors for A");
    ::BigLLL(A.L,0,0.999,VL1);
    
    ccout_separate
    ccout_subtitle("Apply PBKZ");
    B.resize(m,m);
    for (int i=0;i<m;i++) B.L[i] = A.L[i];
    ProgressiveBKZ(B,0,45,VL1,"ignoreflat");
    
    ccout_separate
    ccout_subtitle("Find close vector");
    mat_ZZ VV;
    VV = ENUMCV(B,target,maxnorm,0.5,enum_mode_all_vectors,0,VL3);
    vec_ZZ candidate;
    int can=-1;    //candidate index
    for (int i=0;i<VV.NumRows();i++) {
        candidate = VV[i] - target;
        cout << "error candidate[" << i+1 << "]: " << candidate << endl;
        cout << "norm[" << i+1 << "]=" << LengthOf(candidate) << endl;
        if (LengthOf(candidate) < maxnorm) {
            can = i;
            break;
        }
    }
    if (can==-1) return false;

     //From error to secret by gaussian elimination
    ccout_separate
    ccout_subtitle("Extract secret vector from error candidate by gaussian elimitation");
    for (int i=0;i<n;i++) {
        //normalize i-th row
        ZZ a;
        int shift=1;
        do {
            a = lweA[i][i];
            if (a==0) {
                swap(lweA[i],lweA[i+shift]);
                swap(VV[can][i],VV[can][i+shift]);
                shift++;
            }
        } while (a==0);
        ZZ ai = InvMod(a,to_ZZ(lweq));
        for (int j=0;j<n;j++) {
            lweA[i][j] *= ai;
            lweA[i][j] %= to_ZZ(lweq);
            if (lweA[i][j]<0) lweA[i][j] += lweq;
        }
        VV[can][i] *= ai;
        VV[can][i] %= to_ZZ(lweq);
        if (VV[can][i]<0) VV[can][i] += lweq;

        for (int j=0;j<n;j++) {
            if (j!=i) {
                VV[can][j] -= lweA[j][i] * VV[can][i];
                lweA[j] -= lweA[j][i] * lweA[i];
            }
            for (int k=0;k<n;k++) {
                lweA[j][k] %= to_ZZ(lweq);
                if (lweA[j][k]<0) lweA[j][k] += lweq;
            }
            VV[can][j] %= to_ZZ(lweq);
            if (VV[can][j]<0) VV[can][j] += lweq;
        }
    }
    
    cout << "final solution s=[";
    for (int i=0;i<n;i++) cout << VV[can][i] << " ";
    cout << "]" << endl;
    
    return true;
}

bool simulator_param_estim() {
    //Note: this routine is not called from bkztest()
    cout << "Estimating rational coefficients A1 and A2 in PBKZ simulator" << endl;
    cout << "It may take a few minuets" << endl;
    
    int dim = 150;
    int seed = 0;
    int bkz = 80;
    std::string tflog = "simulator_param_estim.dat";

    cout << "Step 1: generating " << dim << "-dim. LLL basis" << endl;
            
    LatticeBasis<long double> B;
    B = svpc::getlllbasis(dim,0);
    cout << "Step 2: Lattice reduction by PBKZ-" << bkz << " " << endl;

    if (FileExists(tflog)==false) {
        ProgressiveBKZ<long double>(B,0,bkz,VL1,"logfile="+tflog + " betashiftupdate=false parallel=1 temporal=temp200.txt");
    }
 
    std::vector<std::vector<double> > logmatrix;
    loadtable(logmatrix,tflog);

    cout << "logsize=" << logmatrix.size() << endl;
    std::vector<double> real_time_enum,real_time_LLL;
    std::vector<double> sim_time_enum,sim_time_LLL;

    for (int i=0;i<logmatrix.size()-1;i++) {
        //Each record is generated at the start of BKZ-loop
        cout << "input fec=" << logmatrix[i][6] << " ";
        cout << "start time=" << logmatrix[i][0] << " end time=" << logmatrix[i+1][0] << " ";
        cout << "apply bkz=" << logmatrix[i][2]  << " ";


        double LLLtime = logmatrix[i+1][0] - logmatrix[i][0] - (logmatrix[i+1][8] - logmatrix[i][8]) - (logmatrix[i+1][10] - logmatrix[i][10]); //total LLL time + something
        cout << "LLL time=" << LLLtime << endl;
        
        
        //cout << endl;

        real_time_enum.push_back(logmatrix[i+1][9] - logmatrix[i][9]);  //enum cpu time
        real_time_LLL.push_back(LLLtime);

        
        //Below simulate one loop of BKZ
        std::vector<double> falphatable,ualphatable;
        int dim = B.dim;
        int usebeta = logmatrix[i][2];
        SimulateAlphatable(falphatable,ualphatable,dim,usebeta);
        
        std::vector<bkzfloat> cc;
        PBKZSimulate_by_fec(cc,dim,logmatrix[i][6]);    //Todo: output the simulated basis having target fec
        //cout << "check_fec=" << log(FullENUMCost<bkzfloat>(cc,1,dim,lattice_tools::LatticeGH(cc,INPUT_NONSQUARED),INPUT_NONSQUARED)) << endl;
        
        bkzfloat local_enum = bkzconstants::simulatebkzcost(cc,ualphatable,usebeta);
        bkzfloat local_lll = bkzconstants::simulatebkz_lll(cc,usebeta);
        cout << "local_enum=" << local_enum << "    " << endl;
        cout << "local_lll=" << local_lll << "    " << endl;
        sim_time_enum.push_back(local_enum.convert_to<double>());
        sim_time_LLL.push_back(local_lll.convert_to<double>());
        bkzconstants::savecachetable();
    }

    cout << "Step 3: Find rational constants by regression " << endl;
    
    double zy=0;
    double zz=0;
    for (int b=0;b<=real_time_LLL.size();b++) {
        if ((real_time_enum[b] !=0) && (sim_time_enum[b]!=0)) {
            zy += real_time_enum[b] * sim_time_enum[b];
            zz +=  sim_time_enum[b] * sim_time_enum[b];
        }
        
    }    
    cout << "alpha_enum=" << zy/zz << endl;
    double enumalpha = zy/zz;
    
    zy = zz = 0;
    for (int b=20;b<=bkz;b++) {
        if ((real_time_enum[b] !=0) && (sim_time_enum[b]!=0)) {
            zy += real_time_LLL[b] * sim_time_LLL[b];
            zz +=  sim_time_LLL[b] * sim_time_LLL[b];
        }
    }    
    cout << "alpha_LLL=" << zy/zz << endl;
    double LLLalpha = zy/zz;

    std::ofstream of;
    of.open(tflog+".graph.txt",ios::trunc);
    for (int b=0;b<real_time_LLL.size();b++) {
        if ((real_time_enum[b] !=0) && (sim_time_enum[b]!=0)) {
            if ((real_time_enum[b] !=0) && (sim_time_enum[b]!=0)) {
                of << b << "\t" << real_time_LLL[b] << "\t" << real_time_enum[b] << "\t" << sim_time_LLL[b] * LLLalpha << "\t" << sim_time_enum[b] * enumalpha << endl; 
            }
        }
    }
    of.close();
    
    return true;
    
}


void bkztest() {
    
    double ss = gettimeofday_sec(); //timer

    announcement();
    if (basic_test()==false) exit(0);
    if (gs_test()==false) exit(0);
    if (pftable_test()==false) exit(0);
    if (lll_test()==false) exit(0);
    if (simplex_prec_test()==false) exit(0);
    if (enum_test1()==false) exit(0);
    if (enum_test2()==false) exit(0);
    if (enum_test3()==false) exit(0);
    if (preprocess_test()==false) exit(0);
    if (reduction_test()==false) exit(0);
    if (simulator_test()==false) exit(0);
    if (speed_test()==false) exit(0);
    if (simulate_svp_challenge()==false) exit(0);
    if (solve_lattice_challenge()==false) exit(0);
    if (solve_lwe_challenge()==false) exit(0);  

    cout << endl;
    cout << "All tests have been finished." << endl;
    cout << "Total test time=" <<gettimeofday_sec() - ss << endl;
  
    exit(0);
}

#endif
