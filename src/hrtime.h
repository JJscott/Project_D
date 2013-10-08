/*
 * 
 * High-Resolution Time
 * 
 * Adapted from:
 * http://homepage.cs.uiowa.edu/~cwyman/classes/spring07-22C251/code/HighResolutionTimer.h
 * 
 * @author Ben Allen
 * 
 */

#ifndef HRTIME_H
#define HRTIME_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	unsigned char m_data[16];
} hrtime_t;

void getHighResTime(hrtime_t *t);

int highResTimeEquals(const hrtime_t *t0, const hrtime_t *t1);

double highResTimeToSec(const hrtime_t *begin, const hrtime_t *end);

#ifdef __cplusplus
}
#endif

#endif // HRTIME_H
