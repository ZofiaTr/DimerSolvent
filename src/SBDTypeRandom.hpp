#pragma once



/// The class SBDTypeRandom implements a random number generator. 


class SBDTypeRandom {

public:

	/// \name Constructors and destructors
	//@{

	SBDTypeRandom();																									///< A random generator with fixed seed
	SBDTypeRandom(unsigned long s);																						///< A random generator with custom seed
	virtual ~SBDTypeRandom();

	//@}
	/// \name Random integer generators
	//@{

	unsigned long								randUnsignedLong();														///< Generates a random number on [0,0xffffffff]-interval 
	long										randLong();																///< Generates a random number on [0,0x7fffffff]-interval

	//@}
	/// \name Random real generators
	//@{

	double										randDouble1();															///< Generates a random number on [0,1]-real-interval 
	double										randDouble2();															///< Generates a random number on [0,1)-real-interval
	double										randDouble3();															///< Generates a random number on (0,1)-real-interval
	double										randRes53();															///< Generates a random number on [0,1) with 53-bit resolution

	//@}

	/// \name Seed
	//@{

	void										seed(unsigned long);													///< Reseed the random generator

	//@}

private: 

	void														nextState();

	unsigned long												state[624];																///< the array for the state vector
	int															left;
	unsigned long*												next;

};


typedef SBDTypeRandom											SBRandom;


