#ifndef TIMING_H
#define TIMING_H

#include <sys/time.h>

// return time in second
inline double cpuSecond() {
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}

#endif
