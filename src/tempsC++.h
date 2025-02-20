//
// Class Chrono: This file contains functions that can be used to get some
//               timing information about your program.
//
// VERSION 1.1
// 5 - 11 - 2004
//
// AUTHOR : Francois Guertin
//

#ifndef TEMPS_H
#define TEMPS_H

#include <sys/times.h>
#include <sys/time.h>
#include <unistd.h>
#include <limits.h>
#include <iostream>


//==============================================================
class Chrono
{
public:

  //constructor and destructor
   Chrono();

  virtual ~Chrono() {};

  //Resets cumulated time to 0 and indicate that the timer is stopped.
  virtual void    Reset() = 0;

  //start the chrono
  virtual void    Start() = 0;

  //stop the chrono. The time elapsed between the stop and the next start
  //will not be accumulated.
  //the stop and the next start will not be accumulated.
  virtual void    Stop() = 0;

  //This function returns the current elapsed time.  This function can
  //be used after the clock is started and after is it stopped.
  virtual double  GetTime() = 0;

  
protected:
  double clockTicks;           //number of clock ticks per second 
  
  double clockStart;	       // Last time we consulted the clock 
  double clockTotal;	       // Total time elapsed              
  bool fStart;		       // Is the clock started or not     
};

//
//  This function initializes the structure so it can be used later.
//  It must be called before using any of the other functions,
//  otherwise they will not produce valid results.

//
inline
Chrono::Chrono()
  : clockTicks(sysconf(_SC_CLK_TCK)),
    clockStart(0.0),
    clockTotal(0.0),
    fStart(false)
{}

//==============================================================


class ChronoCPU : public Chrono
{
public:

  //constructor and destructor
   ChronoCPU();

  //Resets cumulated time to 0 and indicate that the timer is stopped.
  void    Reset();

  //start the chrono
  void    Start();

  //stop the chrono. The time elapsed between the stop and the next start
  //will not be accumulated.
  //the stop and the next start will not be accumulated.
  void    Stop();

  //This function returns the current elapsed time.  This function can
  //be used after the clock is started and after is it stopped.
  double  GetTime();

};

//  This function initializes the structure so it can be used later.
//  It must be called before using any of the other functions,
//  otherwise they will not produce valid results.
inline ChronoCPU::ChronoCPU()
  : Chrono()
{}


// Resets cumulated time to 0 and indicate that the timer is stopped.
inline void ChronoCPU::Reset()
{
  clockTotal = 0.0;
  fStart = false;
}


// This function indicates that the clock should be started now.     
inline void ChronoCPU::Start()
{
  tms tmsStart;
  times(&tmsStart);
  clockStart = tmsStart.tms_utime + tmsStart.tms_stime;
  fStart = true;
}


// This stops the clock.  This means that the time elapsed between
// the stop and the next start will not be accumulated.
inline void ChronoCPU::Stop()
{
  if (fStart)
  {
    tms tmsCur;
    times(&tmsCur);
    clockTotal += tmsCur.tms_utime + tmsCur.tms_stime - clockStart;
    
    fStart = false;
  }
  else
    std::cout <<"start() must be called before using stop()" <<std::endl;
}


// This function returns the current elapsed time.  This function can
// be used after the clock is started and after is it stopped.
inline double ChronoCPU::GetTime()
{
  if(fStart)
    {
      tms tmsCur;
      times(&tmsCur);
      clockTotal += tmsCur.tms_utime + tmsCur.tms_stime - clockStart; 
      clockStart = tmsCur.tms_utime + tmsCur.tms_stime; 
    }
  return clockTotal / clockTicks ;
}


//==============================================================


class  ChronoReal : public Chrono
{
public:

  //constructor and destructor
   ChronoReal();

  //Resets cumulated time to 0 and indicate that the timer is stopped.
  void    Reset();

  //start the chrono
  void    Start();

  //stop the chrono. The time elapsed between the stop and the next start
  //will not be accumulated.
  //the stop and the next start will not be accumulated.
  void    Stop();

  //This function returns the current elapsed time.  This function can
  //be used after the clock is started and after is it stopped.
  double  GetTime();

};


//  This function initializes the structure so it can be used later.
//  It must be called before using any of the other functions,
//  otherwise they will not produce valid results.
inline ChronoReal::ChronoReal()
  : Chrono()
{}


// Resets cumulated time to 0 and indicate that the timer is stopped.
inline void ChronoReal::Reset()
{
  clockTotal = 0.0;
  fStart = false;
}


// This function indicates that the clock should be started now.     
inline void ChronoReal::Start()
{
  timeval start;
  gettimeofday(&start, NULL);
  clockStart = start.tv_sec + start.tv_usec/1.0e6;
  fStart = true;
}


// This stops the clock.  This means that the time elapsed between
// the stop and the next start will not be accumulated.
inline void ChronoReal::Stop()
{
  if (fStart)
  {
    timeval cur;
    gettimeofday(&cur, NULL);
    clockTotal += cur.tv_sec + cur.tv_usec/1.0e6 - clockStart;
    
    fStart = false;
  }
  else
    std::cout <<"start() must be called before using stop()" <<std::endl;
}


// This function returns the current elapsed time.  This function can
// be used after the clock is started and after is it stopped.
inline double ChronoReal::GetTime()
{
  if(fStart)
    {
      timeval cur;
      gettimeofday(&cur, NULL);
      clockTotal += cur.tv_sec + cur.tv_usec/1.0e6 - clockStart; 
      clockStart = cur.tv_sec + cur.tv_usec/1.0e6; 
    }
  return clockTotal;
}

#endif
