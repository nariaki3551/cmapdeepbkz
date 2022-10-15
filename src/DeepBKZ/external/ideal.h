ZZX find_cyclotomic(long index)
{
  long i;
  ZZX phi;
  phi=1;
  
  for (i = 1; i <= index-1; i++)
    if (index % i == 0)
      phi *= find_cyclotomic(i);
  
    return  (ZZX(i, 1) - 1)/phi;
}

ZZ find_determinant(long index,long length, ZZ seed)
{
  SetSeed(seed);
  ZZ det=RandomLen_ZZ(length);
  det-=rem(det-1,index);
  do {
    det+=index;
  } while(ProbPrime(det, 20)!=1);
  return det;
}

ZZ find_unity_root(long index,ZZ det,ZZX phi)
{
  ZZ_p::init(det);
  
  ZZ e=(det-1)/index;
  ZZ gen;  set(gen);
  ZZ alpha;
  
  do {
    gen++;
    alpha=PowerMod(gen,e,det);
  } while(IsZero(eval(to_ZZ_pX(phi),to_ZZ_p(alpha)))==0);
  return alpha;
}
