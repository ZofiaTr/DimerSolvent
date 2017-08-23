#include "SBDTypeRandom.hpp" 

#include <math.h>

#define N										624
#define M										397
#define MATRIX_A								0x9908b0dfUL   /* constant vector a */
#define UMASK									0x80000000UL /* most significant w-r bits */
#define LMASK									0x7fffffffUL /* least significant r bits */
#define MIXBITS(u,v)							( ((u) & UMASK) | ((v) & LMASK) )
#define TWIST(u,v)								((MIXBITS(u,v) >> 1) ^ ((v)&1UL ? MATRIX_A : 0UL))
#define RANDOM_PI								3.1415926535897932


/// \class SBDTypeRandom
/// The class SBDTypeRandom generator implements the Mersenne Twister generator from the paper:
/// M. Matsumoto and T. Nishimura, "Mersenne Twister: A 623-dimensionally equidistributed 
/// uniform pseudorandom number generator", ACM Trans. on Modeling and Computer Simulation Vol. 8, No. 1, January pp.3-30 (1998).
///
/// Note that the default constructor of SBDTypeRandom, SBDTypeRandom::SBDTypeRandom(), has a fixed seed. This is useful when debugging code, 
/// since it allows for reproducible results. For a "production run", though, a SBDTypeRandom object should be seeded with a random number, 
/// for example using SAMSON::getTime():
/// \code
///
/// SBRandom r(SAMSON::getTime());
/// SBUUID UUID=r.randUUID(); // create a random UUID
///
/// \endcode
///
/// \shortname{SBRandom}
/// \see SBUUID


SBDTypeRandom::SBDTypeRandom() { seed(17081704UL); }
SBDTypeRandom::SBDTypeRandom(unsigned long s) { seed(s); }
SBDTypeRandom::~SBDTypeRandom() { }

unsigned long SBDTypeRandom::randUnsignedLong() {

    unsigned long y;

    if (--left == 0) nextState();
    y = *next++;

    // Tempering 

    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;

}

long SBDTypeRandom::randLong() {

    unsigned long y;

    if (--left == 0) nextState();
    y = *next++;

    // Tempering 

    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return (long)(y>>1);

}

double SBDTypeRandom::randDouble1() {

    unsigned long y;

    if (--left == 0) nextState();
    y = *next++;

    // Tempering
    
	y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return (double)y * (1.0/4294967295.0); 

}

double SBDTypeRandom::randDouble2() {

    unsigned long y;

    if (--left == 0) nextState();
    y = *next++;

	// Tempering

    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return (double)y * (1.0/4294967296.0); 

}

double SBDTypeRandom::randDouble3() {

    unsigned long y;

    if (--left == 0) nextState();
    y = *next++;

    // Tempering
    
	y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return ((double)y + 0.5) * (1.0/4294967296.0); 

}

double SBDTypeRandom::randRes53() { 

    unsigned long a=randUnsignedLong()>>5, b=randUnsignedLong()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 

} 

void SBDTypeRandom::seed(unsigned long s) {
    
	// initializes state[N] with a custom seed
	
	for (unsigned int i=0;i<N;i++) state[i]=0;

	state[0]= (unsigned long)s & 0xffffffffUL;

    for (unsigned int j=1; j<N; j++) {

        state[j] = (1812433253UL * (state[j-1] ^ (state[j-1] >> 30)) + j); 
        state[j] &= 0xffffffffUL;  /* for >32 bit machines */

	}

	left = 1; 

}

void SBDTypeRandom::nextState() {

	unsigned long *p=state;
    int j;

    left = N;
    next = state;
    
    for (j=N-M+1; --j; p++) *p = p[M] ^ TWIST(p[0], p[1]);
    for (j=M; --j; p++) *p = p[M-N] ^ TWIST(p[0], p[1]);

    *p = p[M-N] ^ TWIST(p[0], state[0]);

}
