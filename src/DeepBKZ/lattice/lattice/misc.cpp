#ifndef _inc_latticemisc_cpp
#define _inc_latticemisc_cpp



template <typename T,typename T2> void vec_copy(std::vector<T>& a,std::vector<T2>& b) {
    a.resize(b.size());
    for (int i=0;i<b.size();i++) a[i] = b[i];
}

template <typename T,typename T2> void vec_copy(std::vector<std::vector<T> >& a,std::vector<std::vector<T2> >& b) {
    a.resize(b.size());
    for (int i=0;i<b.size();i++) vec_copy(a[i],b[i]);
}

int to_int(double b) {return (int)b;};

template <typename T> int to_int(T b) {return boost::lexical_cast<int>(b);};

void erasepoint(std::string& xstr) {
    for (int i=0;i<xstr.length();i++) {
        if (xstr[i]=='.') {
            xstr = xstr.substr(0,i);
            return;
        }
    }
    return;
}

template <typename T> std::string to_intstring(T& a) {
    std::string xstr = a.str(0,std::ios_base::fixed);
    erasepoint(xstr);
    return xstr;
}

template <> std::string to_intstring(double& d) {
    float10 a = d;
    std::string xstr = a.str(0,std::ios_base::fixed);
    erasepoint(xstr);
    return xstr;
}

template <> std::string to_intstring(long double& d) {
    float10 a = d;
    std::string xstr = a.str(0,std::ios_base::fixed);
    erasepoint(xstr);
    return xstr;
}

void RowTransformwrap(vec_ZZ& A, vec_ZZ& B, const ZZ& MU1) {
    //A=A-MU*B
    //Note: some of operations are not thread-safe
    if (MU1==0) return;
    if (MU1==1) {
        #pragma omp critical
        {
            A -= B;
        }
        return;
    }
    if (MU1==-1) {
        #pragma omp critical
        {
            A += B;
        }
        return;
    }
    #pragma omp critical
    {
        A -= B*MU1;
    }
    return;
    
}
#ifdef _use_mpfr_pfunc
    void conv(mpfr_float& a,ZZ& b) {    a = (mpfr_float)to_stdstring(b); }
#endif

#ifndef __lite__compile
void conv(float100& a,ZZ& b){    a = (float100)to_stdstring(b); }
void conv(float40& a,ZZ& b) {    a = (float40)to_stdstring(b); }
void conv(float30& a,ZZ& b) {    a = (float30)to_stdstring(b); }
#endif

void conv(float200& a,ZZ& b) {    a = (float200)to_stdstring(b); }
void conv(float150& a,ZZ& b) {    a = (float150)to_stdstring(b); }
void conv(float80& a,ZZ& b) {     a = (float80)to_stdstring(b); }
void conv(float50& a,ZZ& b) {    a = (float50)to_stdstring(b); }
void conv(float20& a,ZZ& b) {    a = (float20)to_stdstring(b); }
void conv(float15& a,ZZ& b) {    a = (float15)to_stdstring(b); }
void conv(float10& a,ZZ& b) {    a = (float10)to_stdstring(b); }

double conv_to_double(const double a) { return a; }
double conv_to_double(const long double a) { return a; }
double conv_to_double(const float10& a) { return a.convert_to<double>(); }
double conv_to_double(const float15& a) { return a.convert_to<double>(); }
double conv_to_double(const float20& a) { return a.convert_to<double>(); }
double conv_to_double(const float30& a) { return a.convert_to<double>(); }
double conv_to_double(const float40& a) { return a.convert_to<double>(); }
double conv_to_double(const float50& a) { return a.convert_to<double>(); }
double conv_to_double(const float80& a) { return a.convert_to<double>(); }
double conv_to_double(const float100& a) { return a.convert_to<double>(); }
double conv_to_double(const float150& a) { return a.convert_to<double>(); }
double conv_to_double(const float200& a) { return a.convert_to<double>(); }

long double conv_to_long_double(const double a) { return a; }
long double conv_to_long_double(const long double a) { return a; }
long double conv_to_long_double(const float10& a) { return a.convert_to<long double>(); }
long double conv_to_long_double(const float15& a) { return a.convert_to<long double>(); }
long double conv_to_long_double(const float20& a) { return a.convert_to<long double>(); }
long double conv_to_long_double(const float30& a) { return a.convert_to<long double>(); }
long double conv_to_long_double(const float40& a) { return a.convert_to<long double>(); }
long double conv_to_long_double(const float50& a) { return a.convert_to<long double>(); }
long double conv_to_long_double(const float80& a) { return a.convert_to<long double>(); }
long double conv_to_long_double(const float100& a) { return a.convert_to<long double>(); }
long double conv_to_long_double(const float150& a) { return a.convert_to<long double>(); }
long double conv_to_long_double(const float200& a) { return a.convert_to<long double>(); }

double to_int(const float10& a) { return a.convert_to<int>(); }
double to_int(const float15& a) { return a.convert_to<int>(); }
double to_int(const float20& a) { return a.convert_to<int>(); }
double to_int(const float30& a) { return a.convert_to<int>(); }
double to_int(const float40& a) { return a.convert_to<int>(); }
double to_int(const float50& a) { return a.convert_to<int>(); }
double to_int(const float80& a) { return a.convert_to<int>(); }
double to_int(const float100& a) { return a.convert_to<int>(); }
double to_int(const float150& a) { return a.convert_to<int>(); }
double to_int(const float200& a) { return a.convert_to<int>(); }

void conv(ZZ& a,bkzfloat& b) {
    conv(a,to_stdstring(b).c_str());
}

void conv(quad_float& a,ZZ& b) {
    a = to_quad_float(b);
}
void conv(double& a,ZZ& b) { a = to_double(b); }
void conv(long double& a,ZZ& b) { a = boost::lexical_cast<long double>(b); }

template <typename T> void copyvec(std::vector<T>& a,std::vector<T>& b) {
    a.resize(b.size());
    for (int i=0;i<b.size();i++) a[i] = b[i];
}

template <typename T> inline void addvector(T* a,int c,T*b,int n) {
    for (int i=0;i<n;i++) a[i] += b[i] * c;
}

template <typename T> inline void display(T* a,int n) {
    ccout << "[";
    for (int i=0;i<n;i++) ccout << a[i] << " ";
    ccout << "]" << endl;
}

template <typename T> inline T squared_norm(vec_ZZ a) {
    T ret = 0;
    ZZ ip;
    InnerProduct(ip,a,a);
    conv(ret,ip);
    return ret;
}


unsigned char ctzbuff[10000];

template <typename T> void conv_to_ZZ(ZZ& a,T& b) {
    //this function is not thread safe!
    if (abs(b) < 2e+9) {
        conv(a,(long int)b);
        return;
    }
    

    double logdet256 = boost::lexical_cast<double>(log(abs(b))) / 5.54517744448;   //5.54 is ln(256)
    int n =2+logdet256;
    if (n>=10000) {
        ccout << "memory error: please change unsigned char ctzbuff[10000]; in misc.cpp" << endl;
    }
    if (n<0) {
        ccout << "value error b=" << b << endl;
        ccout << "log=" << logdet256 << endl;
        exit(0);
    }
    T tb = b;
    if (b<0) tb = -tb;
    T p256 = pow((T)256,n-1);
    for (int i=n-1;i>=0;i--) {
        int byte =(int)(boost::lexical_cast<double>(tb / p256));
        ctzbuff[i] = byte;
        tb -= p256 * byte; 
        p256 /= 256;
    }
    ZZFromBytes(a,ctzbuff,n);    
    if (b < 0) a = -a;
    return;
}

namespace global_error_control {
    
    std::vector<std::map<string,int> > errortables;
    // errortables[i] is a table of error reports in thread i
    
    void put_error(std::string s,int code) {
        
    }
    
    
}




#endif