#ifndef TIMER_H
#define TIMER_H

#include <stdio.h>
#include <iostream>
#include <string>
#include <map>
#include <stack>

#include <sys/time.h>

// To time a function, just put:
//
//  TIMER functionTimer(__FUNCTION__);
//
// at the beginning of the function. It will deactivate at the end of the
// function when the destructor is called. If you want to stop it by hand,
// call the stop() function.

class TIMER
{
public:
  // start the timer by default -- if a tick is called later,
  // it will just stomp it
  TIMER(std::string blockName); 
  ~TIMER();

  void stop();
  const double elapsed() { return _elapsed; };

  static double timing(timeval& begin = _tick, timeval& end  = _tock) {
    double beginTime = (double)begin.tv_sec + 1e-6 * begin.tv_usec;
    double endTime = (double)end.tv_sec + 1e-6 * end.tv_usec;
    return endTime - beginTime;
  };
  static int hours(int seconds) { return seconds / (60 * 60); };
  static int minutes(int seconds) {
   int mod = seconds % (60 * 60);
   return mod / 60;
  };
  static int seconds(int seconds) {
    int mod = seconds % (60 * 60);
    return mod % 60;
  };

  static void printTimings();
  static void printTimingsPerFrame(const int frames);

private:
  // begin and end of current block being timed
  static timeval _tick;
  static timeval _tock;

  // hash table of all timings
  static std::map<std::string, double> _timings;

  // call stack
  static std::stack<std::string> _callStack;

  bool _stopped;
  double _elapsed;
};

#endif
