#ifndef _inc_genlattice_cpp
#define _inc_genlattice_cpp

//Subroutines below call the subroutines in ideal.h and tools.h in (ideal) SVP challenge generator

#ifdef _include_svpchallenge_generator
//SVP challenge problem generator 
//modified from generate_random.cpp in SVP Challenge web page

void gen_svpchallenge(mat_ZZ& B,int n,ZZ seed,int bit=10) {
    vec_ZZ v; 
    generate_random_HNF(v,n,bit,seed);
    B.SetDims(n,n); clear(B);
    B(1,1) = v(1);
    for (int i=2; i<=n; i++)
    {
	B(i,1)=v(i);
	B(i,i)=1;
    }
}
#endif

#ifdef _include_idealchallenge_generator
//Ideal lattice challenge problem generator 
//modified from generate_ideal.cpp in Ideal lattice Challenge web page
void gen_idealsvpchallenge(mat_ZZ& B,int index,ZZ seed,vec_ZZ& phivec) {

    ZZX phi=find_cyclotomic(index);  
    long n=deg(phi);
    ZZ det=find_determinant(index,10*n,seed);
    ZZ alpha=find_unity_root(index,det,phi);
    B.SetDims(n,n); 
    clear(B);
    B(1,1) = det;
    for (long i=2; i<=n; i++)
      {
	B(i,1)=det-PowerMod(alpha,i-1,det);
	B(i,i)=1;
      }
    
    phivec.SetLength(n+1);
    for (int i=0;i<=n;i++) phivec[i] = coeff(phi,i);
}
#endif


#endif