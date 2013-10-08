/*
 * 
 * High-Resolution Time
 * 
 * This stuff is not in the header to avoid namespace pollution (from windows.h)
 * 
 * Also, the MacOSX code is completely untested.
 * 
 */

#include "hrtime.h"

// platform identification. more than one may get defined.
#if defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(__CYGWIN__)
	// using windows
	#define HRTIME_WIN32 
#endif
#if defined(__APPLE__)
	// using macosx
	#define HRTIME_MACOSX
#endif
#if defined(__GNUC__)
	// using gcc under linux or other unixes
	// FIXME other compilers?
	#define HRTIME_UNIX
#endif

#if defined(HRTIME_WIN32)
	// Use code that works on MS Visual Studio. This also works in MinGW and Cygwin.
	#include <windows.h>
	
	//typedef LARGE_INTEGER hrtime_t;
	
	void getHighResTime(hrtime_t *t) {
		QueryPerformanceCounter((LARGE_INTEGER *) t->m_data);
	}
	
	int highResTimeEquals(const hrtime_t *t0, const hrtime_t *t1) {
		return ((LARGE_INTEGER *)(t0->m_data))->QuadPart == ((LARGE_INTEGER *)(t1->m_data))->QuadPart;
	}
	
	double highResTimeToSec(const hrtime_t *begin, const hrtime_t *end) {
		LARGE_INTEGER freq;
		QueryPerformanceFrequency(&freq);
		return (((LARGE_INTEGER *)(end->m_data))->QuadPart - ((LARGE_INTEGER *)(begin->m_data))->QuadPart) / (double) freq.QuadPart;
	}
	
#elif defined(HRTIME_MACOSX) 
	// Assume we're running on MacOSX
	// This code uses calls from the CoreServices framework, so to get this to work you need to
	// add the "-framework CoreServices" parameter g++ in the linking stage. This code was adapted from:
	// http://developer.apple.com/qa/qa2004/qa1398.html
	#include <CoreServices/CoreServices.h>
	#include <mach/mach.h>
	#include <mach/mach_time.h>
	
	// typedef uint64_t hrtime_t;
	
	void getHighResTime(hrtime_t *t) {
		*((uint64_t *)(t->m_data)) = mach_absolute_time();
	}
	
	int highResTimeEquals(const hrtime_t t0, const hrtime_t t1) {
		return *((uint64_t *)(t0->m_data)) == *((uint64_t *)(t1->m_data));
	}
	
	double highResTimeToSec(const hrtime_t *begin, const hrtime_t *end) {
		uint64_t elapsed = *((uint64_t *)(end->m_data)) - *((uint64_t *)(begin->m_data));
		Nanoseconds elapsedNano = AbsoluteToNanoseconds(*(AbsoluteTime *) &elapsed);
		return (*(uint64_t *) &elapsedNano) * (double)(1e-9);
	}
	
#elif defined(HRTIME_UNIX) 
	// Assume we have POSIX call clock_gettime()
	// This works on Cygwin, but the win32 code works better.
	// On some Linux systems, you may need to link in the realtime library (librt.a or librt.so) in 
	// order to use this code.  You can do this by including -lrt on the gcc/g++ command line.
	#include <time.h>
	
	// typedef struct timespec hrtime_t;
	
	void getHighResTime(hrtime_t *t) {
		clock_gettime(CLOCK_MONOTONIC, (struct timespec *) t->m_data);
	}
	
	int highResTimeEquals(const hrtime_t *t0, const hrtime_t *t1) {
		const struct timespec *t0_impl = (const struct timespec *) t0->m_data;
		const struct timespec *t1_impl = (const struct timespec *) t1->m_data;
		return (t0_impl->tv_sec == t1_impl->tv_sec) && (t0_impl->tv_nsec == t1_impl->tv_nsec);
	}
	
	double highResTimeToSec(const hrtime_t *begin, const hrtime_t *end) {
		const struct timespec *begin_impl = (const struct timespec *) begin->m_data;
		const struct timespec *end_impl = (const struct timespec *) end->m_data;
		return (end_impl->tv_sec - begin_impl->tv_sec) + (double)(1e-9) * (end_impl->tv_nsec - begin_impl->tv_nsec);
	}
	
#endif
