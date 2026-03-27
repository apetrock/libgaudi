#ifndef TIMER_H
#define TIMER_H

#include <stdio.h>
#include <iostream>
#include <string>
#include <map>
#include <stack>
#include <chrono>

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
  using TimePoint = std::chrono::steady_clock::time_point;

  // start the timer by default -- if a tick is called later,
  // it will just stomp it
  TIMER(std::string blockName); 
  ~TIMER();

  void stop();
  const double elapsed() { return _elapsed; };

  static double timing(const TimePoint &begin = _tick,
                       const TimePoint &end = _tock) {
    return std::chrono::duration<double>(end - begin).count();
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
  static TimePoint _tick;
  static TimePoint _tock;

  // hash table of all timings
  static std::map<std::string, double> _timings;

  // call stack
  static std::stack<std::string> _callStack;

  bool _stopped;
  double _elapsed;
};

#endif
