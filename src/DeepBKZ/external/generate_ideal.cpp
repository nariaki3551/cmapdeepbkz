#include <sstream>
#include <iostream>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/LLL.h>

NTL_CLIENT;

#include "tools.h"
#include "ideal.h"


int main(int argc,char** argv)
{
  long index=256;
  ZZ seed;

  clear(seed);
    
    PARSE_MAIN_ARGS {
      MATCH_MAIN_ARGID("--index",index);
      MATCH_MAIN_ARGID("--seed",seed);
      SYNTAX();
    }
    
    ZZX phi=find_cyclotomic(index);  
    long n=deg(phi);
    ZZ det=find_determinant(index,10*n,seed);
    ZZ alpha=find_unity_root(index,det,phi);
    
    mat_ZZ B; 
    B.SetDims(n,n); 
    clear(B);
    B(1,1) = det;

    for (long i=2; i<=n; i++)
      {
	B(i,1)=det-PowerMod(alpha,i-1,det);
	B(i,i)=1;
      }
    
    cout << B << endl;
    cout << n << " " << index << " " << seed<< endl;
    cout << phi << endl;
}
