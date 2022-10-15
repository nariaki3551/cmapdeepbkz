
####################
Progressive BKZ Library
Version 202205
Released date 2022/05/17

Contact email address:
pbkz-info@ml.nict.go.jp

Contact postal address: 
4-2-1, Nukui-Kitamachi, Koganei, Tokyo, 184-8795, Japan.

The latest version is available at https://www2.nict.go.jp/security/pbkzcode/.

####################
    INTRODUCTION
####################

This progressive BKZ library is an implementation of the algorithm proposed by
Y. Aono, Y. Wang, T. Hayashi and T. Takagi, 
in ÅgImproved Progressive BKZ Algorithms and their Precise Cost Estimation by Sharp SimulatorÅh,
published in Eurocrypt 2016.
The full version is available at https://eprint.iacr.org/2016/146.
The experimental results presented in this paper also can be verified using this progressive BKZ library.

Author and copyright: 
Security Fundamentals Laboratory, Cybersecurity Research Institute in
National Institute of Information and Communications Technology.

The core part of enumeration subroutine (vectorenumeration_boostcore.cpp) is
modified from the BKZ subroutine in NTL library.
We follow the copyright notices in "doc/copying.txt" in the NTL library 10.5.0:

----------------------------------------------------------------------
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
----------------------------------------------------------------------

This is an open-source library distributed under the terms of the GNU Lesser General Public License version 2.1,
without any warranty.

We notice that the change of license by following NTL.
The previous version of our library (v1.0 and v1.1) were distributed under the GNU General Public License version 3.
Also, this version is distributed under the GNU Lesser General Public License version 2.1.


####################
    REQUIREMENTS
####################

Please ensure that you have already installed the following libraries,
which are required to support. 
We recommend you to use versions no lower than the given ones.
Please see their homepages to get more information about the latest version.

gcc-9.4.0 (https://gcc.gnu.org/)
NTL-11.5.1 (https://www.shoup.net/ntl/)
boost-1.75.0 (https://www.boost.org/)
gmp-6.2.1 (https://gmplib.org/)
gsl-2.5 (https://www.gnu.org/software/gsl/)
mpfr-3.1.4 (https://www.mpfr.org/)

Noting that you may be able to use NTL older than 9.4 to submit your records to SVP Challenge as mentioned in https://www.latticechallenge.org/svp-challenge/.

####################
    INSTALLATION
####################

###
 A
###

Download and extract the file pbkzlib-xxx.zip.

###
 B
###

Our library calls some assistant subroutines from open sources.

Note: install_files.sh automatically executes the installing process,
however, if it fails, please copy them manually.

B1. Bases generator from SVP Challenge:
https://www.latticechallenge.org/svp-challenge/download/generator.zip
In the generator folder, Ågtools.hÅh is necessary.

B2. Bases generator from Ideal Lattice Challenge:
https://www.latticechallenge.org/ideallattice-challenge/download/generator.zip
In the generator folder, Ågideal.hÅh is necessary.

Please copy them to "/external" directory.
Totally the 2 supplements above are required.

To prevent the compiling error, you should modify the 137-th line of tools.h as
Sphere_random(sol,n,power2_RR(bit)*sqrt(n));
 -> Sphere_random(sol,n,power2_RR(bit)*sqrt((double)n));

Besides, the following 2 files from Lattice Challenge and LWE Challenge are optionally required.
Put them into the directory same as main.cpp.
If you don't put it, the test program will skip several (not necessary) steps.

B3. The 600-dimensional challenge instance in Lattice Challenge
https://www.latticechallenge.org/challenges/challenge-600.bz2
extract it.

B4. The (n,alpha)=(40,0.005) instance of LWE Challenge
https://www.latticechallenge.org/lwe_challenge/challenges/LWE_40_005.txt

###
 C
###

run %make in the extracted directory.
Note that according to the compilers,
some warning messages may be imported during the compiling process.
Please ignore them if the file Åga.outÅh is output successfully. 

###
 D
###

Configure the bkz.conf, to set a path for storing the cache files.
In default setting, it makes cache directories at the same place as the execution file as follows:

svpccache=./bkzcache/svpc
constantscache=./bkzcache/const
pfcache=./bkzcache/pf
simcache=./bkzcache/sim

They set the directories to cache the generated and LLL bases of (Ideal) SVP Challenge file, 
some constants used in our library, pruning functions, simulator results.

###
 E
###

After you set the configure file, run %./a.out.
The testing program will be executed.
Note that it takes about 30-60 minutes.

All the test programs are described in lattice/bkztest.cpp,
and you can find examples of how to call our library.

######################
    A SHORT EXAMPLE
######################

To use our library, please include the header: #include <lattice/pbkz.hpp>.


Our library manipulates integer lattice bases by the type LatticeBasis<T> where T is the type to represent
the Gram-Schmidt coefficients.
LatticeBasis<T> has the member "mat_ZZ L" to represent the lattice basis.

Ex1. The Data structure for integer lattices

For your basis L of mat_ZZ type, you can substitute it as 

   LatticeBasis<double> B;
   B = L;

To compute its Gram-Schmidt coefficients, call the members

   B.updateGSBasis();    //compute Gram-Schmidt
   B.gs.displaymu();    //print mu_{i,j}

Ex2. LLL reduction

To find an LLL reduced basis, you can try 

   local_LLL(B,U,delta);

where "mat_ZZ *U" returns the corresponding unimodular matrix such that (output B)=U*(input B).
Also, "double delta" stands for the constant in the LLL algorithm.

Note that local_LLL may not work if the matrix elements are huge 
such as an HNF matrix.
We can call the heuristically modified subroutine BigLLL(B,U,delta,1,B.dim);


Ex3. Progressive BKZ reduction

To find a BKZ reduced basis in our sense, you can call

   ProgressiveBKZ<double>(B,U,beta,vl,stringoption);

 - You have to set the template type <double> by following "LatticeBasis<double> B;"

 - "int beta" stands for the target PBKZ blocksize.
   Note that the output of routine is not exactly classical BKZ-beta reduced basis.
   For example, output of ProgressiveBKZ<double>(B,U,100,0,""); is weaker than classical BKZ-100.

 - "int vl" stands for the amount of displaying information.
   For example, if you set vl=1, the subroutine outputs reasonably small amount of display,
   but if vl=3, it outputs much amount of information for debugging.

 - "std::string stringoption" stands for other options given by text formats.
   A possible text is as follows.
   
   1. stringoption = "istart=10 iend=60 logfile=testlog.log parallel=4";
   means that it reduces the projected sublattice {b_10,...,b_60} and output the computing log to the file testlog.log
   and use 4 threads in enumeration.

   2. stringoption = "targetlen=100 ignoreflat";
   means that it aborts the computation if |b1|<100 
   and it shrinks the computing index to "non-flat range" during the computation.
  This option may be useful for reducing q-ary lattices.

More examples are described in the testing program lattice/bkztest.cpp,
or our website https://www2.nict.go.jp/security/pbkzcode/index.html.


####################
   REPORTING BUGS
####################

If you find a concrete bug, please report it to pbkz-info@ml.nict.go.jp.

All bug reports should include:

       The version number of progressive BKZ, and where you obtained it.
       The hardware, the operating system.
       The version numbers of libraries listed in REQUIREMENTS above, if there is any difference.
       A description of the bug behavior.
       A short program which can reproduce the bug.

It is useful if you can check whether the same problem occurs after you update your libraries.
Thank you.

Any errors or omissions in the manual can also be reported to 
the same address as pbkz-info@ml.nict.go.jp.


########################
   MISCELLANEOUS INFO.
########################

Subroutine to generate Unimodular matrices (lattice/gen_uni_mat.cpp) 
is recoded from the open source code Sage-7.1, where the function is Ågrandom_unimodular_matrix()Åh.

################
   HISTORY
################

v202205 (2022/05/18)
 - Bug fixes

v201808 (2018/08/06)
 - Small bug fixes

v201803 (2018/03/12)

v1.1 (2016/6/29)

v1.0 (2016/5/2) 
