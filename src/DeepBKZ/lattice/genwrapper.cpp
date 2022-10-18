#ifndef _inc_lattice_gen_wrapper
#define _inc_lattice_gen_wrapper


//A toolkit to experiment LatticeChallenge problems
//Generated bases are saved as cache files

#define opt_random 0x00
#define opt_ideal 0x01

extern void BigLLL(mat_ZZ& L,mat_ZZ*U,double delta,int vl);

namespace svpc {
    int init = 0;
    std::string cachedir;   //cache dir
    
    void initialize() {
        if (init>0) return; 
#ifdef _allow_cachefiles
        cachedir =  makefullpath(ReadConf("bkz.conf","svpccache"));
        mkdirRecursive(cachedir.c_str(), 0777);
#endif       
        init  = 1;
    }

    std::string getfname(int n,int seed,double bit,int level) {
        //cache file name
        initialize();
        std::ostringstream fname;
#ifdef _allow_cachefiles
        fname << cachedir;
        mkdir(fname.str().c_str(),0777);
        
        //dimension
        fname << "/dim" << n;
        mkdir(fname.str().c_str(),0777);
        //type
        fname << "/random";
        mkdir(fname.str().c_str(),0777);
        //seed, bit, level
        fname << "/seed" << seed << "_bit" << bit;
        if (level==0) {
            fname << "_gen.txt";
        } else 
        if (level==1) {
            fname << "_lll.txt";
        } else        
        if ((2<=level) && (level<=n)) {
            fname << "_bkz" << level << ".txt";
        } else 
        if (level==n+1) {
            fname << "_vr.txt";
        } else        
        if (level>n+2) {
            fname << "_sv.txt";
        }        
#endif
        return fname.str();
    }
    
    mat_ZZ getbasis(int n,int seed,double bit=10) {
        initialize();
        std::ostringstream fname;
        fname << getfname(n,seed,bit,0);

        mat_ZZ L;
        if (FileExists(fname)==true) {
            LoadLattice(L,fname.str());
        } else {
            gen_svpchallenge(L,n,to_ZZ(seed),bit);
#ifdef _allow_cachefiles
            SaveLattice(L,fname.str());
#endif
        }
        return L;
    }

    mat_ZZ getlllbasis(int n,int seed,double bit=10) {
        initialize();
        std::ostringstream fname;
        fname << getfname(n,seed,bit,1);
        mat_ZZ L;
        
        if (FileExists(fname)==true) {
            LoadLattice(L,fname.str());
        } else {
                L = getbasis(n,seed,bit);
                ::BigLLL(L,0,0.999,VL1);
#ifdef _allow_cachefiles
            //ccout << "save L to: " << fname.str() << endl;
            SaveLattice(L,fname.str());
#endif
        }
        return L;
    }
}

namespace isvpc {
    //Generating instances for ideal lattice challenges
    int init = 0;
    std::string cachedir;   //cache dir
    
    void initialize() {
        if (init>0) return; 
#ifdef _allow_cachefiles
        cachedir =  makefullpath(ReadConf("bkz.conf","svpccache"));
        mkdirRecursive(cachedir.c_str(), 0777);
#endif       
        init  = 1;
    }

    std::string getfname(int n,int seed,double bit,int level) {
        initialize();
        std::ostringstream fname;
#ifdef _allow_cachefiles
        fname << cachedir;
        mkdir(fname.str().c_str(),0777);
        //dimension
        fname << "/dim" << n;
        mkdir(fname.str().c_str(),0777);
        //type
        fname << "/ideal";
        mkdir(fname.str().c_str(),0777);
        //seed, bit, level
        fname << "/seed" << seed << "_bit" << bit;
        if (level==0) {
            fname << "_gen.txt";
        } else 
        if (level==1) {
            fname << "_lll.txt";
        } else        
        if ((2<=level) && (level<=n)) {
            fname << "_bkz" << level << ".txt";
        } else 
        if (level==n+1) {
            fname << "_vr.txt";
        } else        
        if (level>n+2) {
            fname << "_sv.txt";
        }        
#endif
        return fname.str();
    }

    mat_ZZ getbasis(int n,int seed,vec_ZZ& phi,double bit=10) {
        initialize();
        
        std::ostringstream fname;
        fname << getfname(n,seed,bit,0);

        mat_ZZ L;
        if (FileExists(fname)==true) {
            LoadLattice(L,fname.str());
            fname << ".phi";
            LoadElement(phi,fname.str());
        } else {
            vec_ZZ phi;
            //n stands for index
            gen_idealsvpchallenge(L,n,to_ZZ(seed),phi); 
#ifdef _allow_cachefiles
            SaveLattice(L,fname.str());
            fname << ".phi";
            SaveElement(phi,fname.str());
#endif
        }
        return L;
    }

    mat_ZZ getlllbasis(int n,int seed,vec_ZZ& phi,double bit=10) {

        initialize();
        std::ostringstream fname;
        fname << getfname(n,seed,bit,1);
        mat_ZZ L;
        
        L = getbasis(n,seed,phi,seed); // to recover phi
        //ccout << "L.dim=" << L.NumRows() << endl;
        if (FileExists(fname)==true) {
            LoadLattice(L,fname.str());
        } else {
            L = getbasis(n,seed,phi,seed);
            ::BigLLL(L,0,0.999,VL1);
#ifdef _allow_cachefiles
            SaveLattice(L,fname.str());
#endif
        }
        return L;
    }
}


#endif