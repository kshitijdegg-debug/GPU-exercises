#include <stdio.h>
#include <stdint.h>
#include <iostream>

// Need to find (-q^-1) mod R for Montmogery modular multiplication

/* 
The Newton–Raphson method is a root-finding algorithm for real functions, but it generalizes beautifully to modular arithmetic over powers of 2.
Our R = 2^64 is a power of 2 and Newton-Raphson is fast and convergent for power-of-two moduli.
The key reason is how modular arithmetic behaves under powers of two.

If q is odd, it’s invertible modulo 2^k.

The Newton iteration ensures that if x_n is correct modulo 2^m, then x_(n+1) is correct modulo 2^(2m).

So it doubles precision each time — that’s why after just a few steps, you get the full 64-bit inverse.

Example for intuition:
Start with x_0=1 which satisfies x_0=1 mod (2^1)

Then each step refines it to be correct mod 2^2, 2^4, 2^8,...,2^k.

Newton Raphson's method - 
Solution for f(x)=0 can be found by iteratively - 
x_(n+1) = x_n - f(x_n)/f'(x_n)

Here, f(x) = 1/x + q;  (solving (1/x+q=0) gives x=-q^(-1))

Then, x_(n+1) = 2*x_n - q*(x_n)^2

As R=2^64, we get correct -q^(-1) within 6 iterations, after x0=1.
*/

uint64_t montgomery_inverse(uint64_t q, __uint128_t R)
{
	int k=0;
	while(R!=1)
	{
		R/=2;
		k++;
	}

	uint64_t x=1;

	for(int i=0;i<k;i++)
	{
		x = 2*x - q*x*x;
	}

	return x;
}

int main()
{
	uint64_t q=(1ULL<<61)-1;
	__uint128_t R = ((__uint128_t)1<<64);	

	uint64_t R_hi = (uint64_t)(R >> 64);
    uint64_t R_lo = (uint64_t)(R & 0xFFFFFFFFFFFFFFFFULL);
    printf("q: %llu, R_hi: %llu, R_lo: %llu\n",q, R_hi,R_lo);

    uint64_t inv = montgomery_inverse(q, R);
    printf("q^-1(mod R) = %llu\n", inv);

    uint64_t minus_one_mod_R = R-1;

    uint64_t minus_one_into_q_inv_mod_R = (inv*minus_one_mod_R)%R;

    printf("(-q^-1)mod R= %llu\n",minus_one_into_q_inv_mod_R);

    uint64_t check = (minus_one_into_q_inv_mod_R*q)%R;
    uint64_t check2 = (inv*q)%R;

    printf("check: %llu, check2: %llu\n",check,check2);

	return 0;
}