#include "xtime.hpp"
#ifdef  __linux__
#include <sys/time.h>
#include <time.h>

TTime::TTime()
{
	gettimeofday(&tstart, NULL);
	ts =tstart.tv_sec*1000 +(tstart.tv_usec/1000.0);
}


double TTime::start()
{
	gettimeofday(&tstart, NULL);
	ts = tstart.tv_sec*1000 + (tstart.tv_usec/1000.0);
	return ts;
}


double TTime::stop()
{
  gettimeofday(&tstop, NULL);
  tp = tstop.tv_sec*1000 + (tstop.tv_usec/1000.0);
  te = (tstop.tv_sec - tstart.tv_sec)*1000.0;
  te  +=  (tstop.tv_usec - tstart.tv_usec)/1000.0;
  return (te);
}


double TTime::elapsed()
{
  return (te);
}

#else
//Windows
#include <windows.h>

TTime::TTime()
{
  // get ticks per second
  QueryPerformanceFrequency(&frequency);
  // start timer
  QueryPerformanceCounter(&t1);
}

double TTime::start()
{
  QueryPerformanceCounter(&t1);
  return  t1.QuadPart*1000/frequency.QuadPart;
}

double TTime::stop()
{
  // stop timer
  QueryPerformanceCounter(&t2);
  // compute and print the elapsed time in millisec
  elapsedTime = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart;
  return  elapsedTime;
}

double TTime::elapsed()
{
  return  elapsedTime;
}

#endif

/*int main(void)
{
  char buffer[30];
  struct timeval tv;

  time_t curtime;

  gettimeofday(&tv, NULL);
  curtime=tv.tv_sec;

  strftime(buffer,30,"%m-%d-%Y  %T.",localtime(&curtime));
  printf("%s%ld\n",buffer,tv.tv_usec);

  return 0;

}
*/
