#ifndef _inc_latticetools_cpp
#define _inc_latticetools_cpp

#include "bkzconstants.hpp"

#define INPUT_SQUARED 0x40
#define INPUT_NONSQUARED 0x20

//Basic toolkits for lattice

namespace lattice_tools {
    
    int init=0;
    
    double* smallghconst;
    double* simhkzconst;
    void initialize() {
        if (init!=0) return;

        int i,j;
        smallghconst = new double[51];
        simhkzconst = new double[51];

        //Copied from Chen-Nguyen's paper
        double CNconstant[] = {0.593208,0.582161,0.561454,0.544344,0.522066,0.502545,0.479628,
    0.459819,0.438675,0.413708,0.392483,0.370717,0.344447,0.322574,
    0.297318,0.273761,0.249247,0.225483,0.199940,0.173832,0.147417,
    0.123425,0.100035,0.074487,0.043089,0.020321,-0.013844,-0.042863,
    -0.068204,-0.093892,-0.124345,-0.151097,-0.183912,-0.214122,
    -0.241654,-0.274612,-0.302966,-0.330965,-0.367514,-0.391956,
    -0.426507,-0.457813,-0.488113,-0.518525,-0.554184,-0.585479,
    -0.617705,-0.646749,-0.671864,-0.687300};    

        for (i=0;i<25;i++) swap(CNconstant[i],CNconstant[49-i]);
        for (i=0;i<=49;i++) {
            CNconstant[i] = exp(CNconstant[i]);
        }
        double det;
        for (i=1;i<=49;i++) {
            det = 0;
            for (j=0;j<=i;j++) det += log(CNconstant[j]);
            det = exp(1.0 * det/(i+1)) * (double)bkzconstants::ghconstant(i+1); 
            smallghconst[i+1] = CNconstant[i] / det;
        }
        
        for (i=0;i<25;i++) swap(CNconstant[i],CNconstant[49-i]);
        det=0;
        for (j=0;j<50;j++) det += log(CNconstant[j]);
        det = exp(det / 1.0 / 50);   //det^(1/n)
        for (i=1;i<=50;i++) simhkzconst[i] = CNconstant[i-1] / det;     //constant |b1|=C*det(L)^(1/n) in n-dim last basis
        init = 1;
    }
    
    void finish() {
        //Todo
    }
    
    template <typename T> T LatticeGH(std::vector<T>& c,int istart,int iend,char opt=INPUT_NONSQUARED) { 
        initialize();
        T ret = 0;
        for (int i=istart;i<=iend;i++) {
            ret += log(c[i]);
        }
        ret /= (iend-istart+1);
        if (opt==INPUT_SQUARED) ret *= 0.5;
        ret = exp(ret);
        lattice_tools::initialize();
        return ret *= (T)bkzconstants::ghconstant(iend-istart+1);
    }

    //Gaussian Heuristic length of lattice
    template <typename T> T LatticeGH(std::vector<T>& c,char opt=INPUT_NONSQUARED) { 
        return LatticeGH(c,1,c.size()-1,opt);
    }

}

template <typename T> std::string vectorhash(std::vector<T>& v,int start=0) {
    std::size_t h1 = 0;
    std::string sh;
    for (int i=start;i<v.size();i++) {
        h1 += (i+1)*std::hash<std::string>{}(to_stdstring(v[i]));
    }
    sh = to_hexstring(h1);        
    while (sh.length()<16) sh = "0" + sh;
    return sh;
}

namespace bkzconstants {

    void  mainparam_update(int k) {
        //Base parameter (alpha,r,p) at k-dimension of progressive BKZ
        lattice_tools::initialize();
        int cs = mainparam_table.size();
        if (k >= cs) {
            #pragma omp critical        
            {
                mainparam_table.resize(k+1);
                for (int i=cs;i<=k;i++) {
                    mainparam_table[i].resize(3);
                    for (int j=0;j<3;j++) {
                        mainparam_table[i][j] = 0;
                    }
                }
            }
        }
        if (mainparam_table[k][0]==0) {
            double alpha;
            double r;
            bkzfloat prob;
            int beta = k;

            if (beta >= 84) {
                //Heuristic parameter from experiments
                r = -40.60241266/beta/(beta-1) - 1.603898*log(beta)/(beta+113.3386);
                r = exp(r);
                double abar = 1.0 / (double)bkzconstants::ghconstant(beta) * exp(-0.25*(beta-1.0)*log(r));
                alpha = abar * ((1.0+beta)/beta);
                prob = pow((bkzfloat)alpha,-beta)*2;
            } else {
                r = -1.50065628562105e-06 * beta*beta +  0.000444247990862759 * beta +  0.932324451217323;  //fitting
                double abar = 1.0 / (double)bkzconstants::ghconstant(beta) * exp(-0.25*(beta-1.0)*log(r));
                double minalpha = pow(1.38629436111989,1.0/beta);
                minalpha = pow(2.77258872223978,1.0/k);
                abar = max(abar,minalpha);
                alpha = abar * ((1.0+beta)/beta);
                prob = pow((bkzfloat)alpha,-beta)*2;
                if (beta<=55) prob=min(0.15+0.03*(55-beta),1.0);
            }
            mainparam_table[k][0] = r;
            mainparam_table[k][1] = alpha;
            mainparam_table[k][2] = prob;
        }
    }

    double mainparam_r(int beta) {
        mainparam_update(beta);
        return mainparam_table[beta][0].convert_to<double>();
    }

    double mainparam_alpha(int beta) {
        mainparam_update(beta);
        return mainparam_table[beta][1].convert_to<double>();
    }

    bkzfloat mainparam_prob(int beta) {
        mainparam_update(beta);
        return mainparam_table[beta][2];
    }

}

#endif

