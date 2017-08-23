#pragma once


#if !defined(C_PROC)
	#define C_PROC 2700000000.0
#endif

#if defined(__linux)

	#if defined(__i386__) && !defined(__x86_64__) //32 bit

		#if !defined(__INTEL_COMPILER)
			typedef long long __int64; 
		#endif

		static __inline__ __int64 SBCTimeClock(void) {
		__int64 x;
		__asm__ volatile ("rdtsc\n\t" : "=A" (x));
		return x;	
		}

	#else // 64 bit

		#include <stdint.h>
		typedef uint64_t __int64;

		static __inline__ __int64 SBCTimeClock(void) {
			uint32_t lo, hi;
			__asm__ __volatile__ (      // serialize
								  "xorl %%eax,%%eax \n        cpuid"
								  ::: "%rax", "%rbx", "%rcx", "%rdx");
			/* We cannot use "=A", since this would use %rax on x86_64 */
			__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
			return (uint64_t)hi << 32 | lo;
		}

	#endif

#endif

#if defined(__APPLE_CPP__) || defined(__APPLE_CC__) || defined(__MACOS_CLASSIC__) 

	//typedef unsigned long long __int64; 
	#include <stdint.h>
	typedef uint64_t __int64;

	static __inline__ __int64 SBCTimeClock(void) {
		uint32_t lo, hi;
		__asm__ __volatile__ (      // serialize
							  "xorl %%eax,%%eax \n        cpuid"
							  ::: "%rax", "%rbx", "%rcx", "%rdx");
		/* We cannot use "=A", since this would use %rax on x86_64 */
		__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
		return (uint64_t)hi << 32 | lo;
	}

#endif //MAC

#if defined(_WIN32)

	inline
	__declspec(naked)
	__int64 SBCTimeClock() {
		__asm {
			rdtsc
				ret
		}
	}

#endif

#define SBCTime __int64

#define SBTime SBCTime

