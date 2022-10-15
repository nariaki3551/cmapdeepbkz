#ifndef _inc_pbkzsupport_cpp
#define _inc_pbkzsupport_cpp

#include "pbkzsharemem.cpp"

// Estimate (Hermite factor)^(1/n) 
double Estimate_delta(quad_float* c,int bs) {
   
    //Estimating delta by least square
    double r,xs,ys,xxs;
    int i;
    xs = bs*(bs-1.0)/2.0;
    xxs = bs*(bs-1.0)*(2.0*bs-1.0)/6.0;
    ys = 0;
    for (i=0;i<bs;i++) ys += 0.5*log(to_double(c[i]));
    r = 0;
    for (i=0;i<bs;i++) r += 0.5*i*log(to_double(c[i]));
    r = (1.0*bs*r - xs*ys) / (1.0*bs*xxs-xs*xs);
    r = exp(2.0*r);
    double ret = exp(-log(r)*(bs-1.0)/(4.0*bs)); 
    if ((ret<=0.9) || (ret>=1.1) ) {
        if (bs>10) {
            ccout << "delta-error: " << endl;
            for (i=0;i<bs;i++) {
                ccout << "c[" << i << "]=" << c[i] << endl;
            }
        }
    }
    return  ret;
}


double InitRadius(quad_float* c,int jj,int bs,double radius,char mode) {
    //Input sequence of |b*i|^2, indexes of local block, parameters
    //Output: (initial radius)^2
  double clim = 0, ret = 0;
  
    if (mode=='G') {
        //return=(alpha*GH)^2
        clim = lattice_tools::LatticeGH(c+jj-1,bs,INPUT_SQUARED);      // computing from c[jj...jj+bs-1]
        clim *= radius;
        clim *= clim;
        clim = min(clim,to_double(c[jj]));
        ret = clim;
    }
    else if (mode=='0') ret = radius * radius;

    return ret;
}


extern RR optimize_pruning_functionCEpre(double* pf,double clim,quad_float* c,int n,double tprob,int vl,int timemax,char order,char opttarget);


quad_float ENUMCost(quad_float* c,int jj,int bs,double p0,double radius,char mode,double delta=0,char input=INPUT_SQUARED,char optimize_level=0){
    
}
 


int compute_approx(mat_ZZ& B,int m,int n,quad_float** B1,quad_float** mu,quad_float* b,quad_float* c,mat_ZZ* U) {
}

void compute_approxRR(mat_ZZ& B,int m,int n,mat_RR& B1,mat_RR& mu,vec_RR& b,vec_RR& c,mat_ZZ* U) {
}

void convert_to_double(vector<vector<double> >& d,vector<vector<int> >& l,int offset) {
}

template <typename T> double l2norm(vector<T>& d,char sq=0) {
    int i;
    double ret=0;
    for (i=0;i<(int)d.size();i++) {
        ret += (double)d[i]*(double)d[i];
    }
    if (sq==0) {
        return sqrt(ret);
    } else {
        return ret;
    }
}

template <typename T> double l2ip(vector<T>& d,vector<T>& d2) {
    int i;
    double ret=0;
    for (i=0;i<(int)d.size();i++) {
        ret += (double)d[i]*(double)d2[i];
    }
    return ret;
}

void computeorthogonal(std::vector<double>& a,std::vector<double>& b) {

    if (a.size()!=b.size()) return;
    int i;
    double mu = 0;
    double norm = 0;
    for (i=0;i<(int)a.size();i++) {
        mu += a[i] * b[i];
        norm += b[i] * b[i];
    }
    mu /= norm;
    for (i=0;i<(int)a.size();i++) {
        a[i] -= mu * b[i];
    }
}

void Insert(mat_ZZ& L,vec_ZZ& v,int index) {
    int n,m;
    int i,j;
    m = L.NumCols();
    n = L.NumRows();

    L.SetDims(n+1,m);
    for (i=n-1;i>=index;i--) {
        L[i+1] = L[i];
    }
    L[index] = v;
    LLL_QP(L,0.999,0,0,1);
    
    j = 0;
    while (LengthOf(L[j])==0) j++;
    for (i=0;i<=n-j;i++) {
        L[i] = L[i+j];
    }
    L.SetDims(n,m);
    return;
}
void gen_hkzlattice(quad_float* c,int n,double alpha) {
    
    //generating |b*i|=alpha*GH(L') for all i
    //output c[i]=|b*i|
    int i,j;


    double CNconstant[] = {0.593208,0.582161,0.561454,0.544344,0.522066,0.502545,0.479628,
0.459819,0.438675,0.413708,0.392483,0.370717,0.344447,0.322574,
0.297318,0.273761,0.249247,0.225483,0.199940,0.173832,0.147417,
0.123425,0.100035,0.074487,0.043089,0.020321,-0.013844,-0.042863,
-0.068204,-0.093892,-0.124345,-0.151097,-0.183912,-0.214122,
-0.241654,-0.274612,-0.302966,-0.330965,-0.367514,-0.391956,
-0.426507,-0.457813,-0.488113,-0.518525,-0.554184,-0.585479,
-0.617705,-0.646749,-0.671864,-0.687300};    
    
    for (i=1;i<=min(50,n);i++) {
        c[i] = exp(CNconstant[50-i]);
    }
    
    for (i=51;i<=n;i++) {
        c[i] = 0;
        for (j=1;j<i;j++) {
            c[i] += 1.0*log(c[j]) / (i-1.0);
        }
        c[i] += 1.0*log(alpha)*i/(i-1.0);
        c[i] -= to_quad_float(1.0*log(VolumeBall(i,1.0))/(i-1.0));
        c[i] = exp(c[i]);
    }

    for (i=1;i<=n/2;i++) swap(c[i],c[n-i+1]);
    for (i=2;i<=n;i++) c[i] /= c[1];
    c[1] = 1.0;
    
}

RR minimum_cost(int n,double alpha) {

    int i;
    quad_float* c;
    c = new quad_float[n+1];
    quad_float ucost;
    
    gen_hkzlattice(c,n,1.00);

    for (i=1;i<=n;i++) c[i] *= c[i];    //c[i]=|b*_i|^2
    
    return lattice_tools::FullENUMCost(c,n);
    //ucost = 2.0*ENUMCost(c,1,n,prob,alpha,'G');
    //return ucost;
}

double compute_target_r(int beta) {
    double r;
    lattice_tools::initialize();
    if (beta<=90) {
        r = -18.2139/(beta+318.978);
        r = exp(r);
        if (beta<=35) {
            r /= exp(2.0*log(lattice_tools::smallghconst[beta])/beta);
        }    
    } else {
        r = max(-1.06889/(beta-31.0345)*log(0.417419*beta-25.4889), -18.2139/(beta+318.978));
        r = exp(r);
    }
    return r;
}

#include "pbkzsimulator.cpp"



#endif
