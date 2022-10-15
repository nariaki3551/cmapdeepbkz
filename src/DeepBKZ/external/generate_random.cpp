#include <sstream>
#include <iostream>
#include <NTL/LLL.h>

NTL_CLIENT;

#include "tools.h"

int main(int argc,char** argv)
{
    long n = 80;
    long bit = 10;
    ZZ seed; seed = 0;

    PARSE_MAIN_ARGS {
	MATCH_MAIN_ARGID("--dim",n);
	MATCH_MAIN_ARGID("--seed",seed);
//	MATCH_MAIN_ARGID("--bit",bit);
	SYNTAX();
    }

    vec_ZZ v; generate_random_HNF(v,n,bit,seed);
    mat_ZZ B; B.SetDims(n,n); clear(B);
    B(1,1) = v(1);
    for (int i=2; i<=n; i++)
    {
	B(i,1)=v(i);
	B(i,i)=1;
    }
    cout << B << endl;
}
