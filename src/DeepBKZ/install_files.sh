#! /bin/bash
if [-e external]; then
  echo "directory external found"
else
	mkdir external
fi

cd external

if [-e external/generator_svp.zip]; then
  echo "generator_svp.zip found"
else
  wget http://www.latticechallenge.org/svp-challenge/download/generator.zip -O generator_svp.zip --no-check-certificate
fi

unzip generator_svp.zip

if [-e external/generator_ideal.zip]; then
  echo "generator_ideal.zip found"
else
	wget http://www.latticechallenge.org/ideallattice-challenge/download/generator.zip -O  generator_ideal.zip --no-check-certificate
fi

unzip -o generator_ideal.zip

cd ..

if [-e challenge-600.bz2]; then
  echo "challenge-600.bz2 found"
else
	wget http://www.latticechallenge.org/challenges/challenge-600.bz2 --no-check-certificate
fi

bzip2 -d challenge-600.bz2

if [-e LWE_40_005.txt]; then
  echo "LWE_40_005.txt found"
else
	wget https://www.latticechallenge.org/lwe_challenge/challenges/LWE_40_005.txt  --no-check-certificate
fi


