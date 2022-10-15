template<typename T>
std::istream& operator>>(const std::string& in, T& v)
{
    static std::istringstream* iss = 0;
    if (iss) delete iss;
    iss = new std::istringstream(in);
    *iss>>v;
    return *iss;
}

#define PARSE_MAIN_ARGS \
    std::ostringstream message; \
    message << argv[0];  \
    for (int i=0; i<argc; ++i)

#define MATCH_MAIN_ARGID(arg,dest) \
    if (i==0)  {\
	message << " [" << arg << " " << dest << "]"; \
    } else { \
	if (std::string(argv[i]) == std::string(arg))  {\
	    std::string(argv[i+1]) >> dest; \
	    i++; continue; \
	} \
    }

#define DETECT_MAIN_ARGID(arg,dest,value) \
    if (i==0)  {\
	message << "[" << arg << "]"; \
    } else { \
	if (std::string(argv[i]) == std::string(arg))  {\
	    dest = value; \
	    continue; \
	} \
    }

#define SYNTAX() \
    if (i==0) continue; \
    std::cerr << "Syntax: " << message.str() << std::endl; \
    exit(1);

//--------------------- tools to compute the Gaussian heuristic
// Computes \log(\Gamma(i/2))
// Caches results for faster calculation
RR LogGammaIs2(long i)
{
    static vec_RR lgam;
    static int imax=-1;
    if (imax==-1) { // do first time only
        lgam.SetLength(1000);
        lgam(2)=0.;
        lgam(1)=log(to_RR(M_PI))/to_RR(2);
        imax=2;
    }
    if (imax>=i) return lgam(i);
    for (int j=imax+1;j<=i;j++) {
        lgam(j)= lgam(j-2) + log(j/2.-1.);
    }
    imax=i;
    return lgam(i);
}

// Computes the radius of the ball of volume 1 in n dimensions
RR Gauss_rn(long n)
{
    RR R_PI; string("3.14159265358979323846264338328") >> R_PI;
    return exp((LogGammaIs2(n+2)-to_RR(n/2.)*log(R_PI))/to_RR(n));
}

//-----------uniform random vector on the unit sphere
// Sample a real number with Gaussian distribution
void Gauss_random(RR & r)
{
    RR u_1,u_2,v_1,v_2,un,s,t;

    conv(un,1);

    do {
        random(u_1);
        random(u_2);
        add(t,u_1,u_1);
        sub(v_1,t,un); /* v_1 = 2u_1-1 */
        add(t,u_2,u_2);
        sub(v_2,t,un); /* v_2 = 2u_2-1 */
        sqr(t,v_1);
        sqr(s,v_2);
        add(s,t,s); /* s = v_12+v_22 */
    }
    while (s >= un);

    t = -2*log(s)/s;
    SqrRoot(s,t);
    mul(r,v_1,s);
}

// Choose a uniformly random vector from the n dimensional sphere
// First chooses a Gaussian vectors and then normalizes
void Sphere_random(vec_ZZ& out, long n, const RR& norme)
{
    RR norm_v,un,norm2;
    long i;
    vec_RR v;
    v.SetLength(n);
    out.SetLength(n);

    for (i=1;i<=n;i++) {
        Gauss_random(v(i));
    }
    norm_v = SqrRoot(v*v);
    v *= norme/norm_v;
    for (int i=1; i<=n; ++i)
	RoundToZZ(out(i),v(i));
}


//------ Random lattice challenge

void generate_random_HNF(vec_ZZ& out,long n,long bit, ZZ seed)
{
    SetSeed(seed);
    out.SetLength(n); clear(out);
    ZZ p; GenPrime(p,bit*n);
    out(1) = p;
    for (int i=2; i<=n; i++)
    {
	RandomBnd(out(i),p);
    }
}

#include <NTL/RR.h>

void generate_chlgsvp_HNF(vec_ZZ& out,long n,long bit,const ZZ& seed)
{
   SetSeed(seed);
   //generate a solution vector
   ZZ p1; GenPrime(p1,bit);
   vec_ZZ sol; sol.SetLength(n);
   Sphere_random(sol,n,power2_RR(bit)*sqrt(static_cast<double>(n)));
   cerr << "expected SVP: " << sol << endl;

   RR lambda1; lambda1 = sqrt(to_RR(sol*sol));
   RR rnvol = (lambda1 * 1.01)/Gauss_rn(n);
   double dbitvol = to_double(to_RR(n) * log(rnvol)/log(2.));

   //generate the orthogonal lattice
   out.SetLength(n+1);

   //generate the volume
   ZZ p; int k=1;
   do { GenPrime(p,(k++) + long(dbitvol)); }
   while((log(p)/log(2))<dbitvol);
   cerr << "bitvol:" << (log(p)/log(2)) << " /required: "<< (dbitvol) << endl;

   //generate the hermite normal form
   out(1)=p;
   ZZ somme; somme = 0;
   for (int i=2; i<=n; i++)
   {
       RandomBnd(out(i),p);
       somme = (somme + out(i)*sol(i-1))%p;
   }
   //out(n+1) = -somme / sol(n) mod p
   InvMod(out(n+1),p-(sol(n)%p),p);
   out(n+1) = (out(n+1)*somme)%p;
}

