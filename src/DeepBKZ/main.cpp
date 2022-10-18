/*
   progressive BKZ library by NICT security fundemental lab.
 */

#include <lattice/pbkz.hpp>

#include <DeepBKZ/lattice.h>
#include <DeepBKZ/DeepLLL.h>
#include <DeepBKZ/DeepBKZ.h>

int main(int argc, char** argv) {
	int n = 100;
	lattice L; L.SetDims(n, n);

	/* Input basis */
	std::ifstream File("Bases/BKZ20_SVP100_seed0.txt");
	File >> L.basis;
	// cout << L.basis << endl;

	// cout << "DeepLLL-start" << endl;
	L.SetGSO();
	// cout << L.B << endl;

	for (int b=50; b<=80; b=b+10) {
		L.DeepBKZ(1, n, b, 0.99, n, 4);
	}
}
