#ifndef TIMER_HPP
#define TIMER_HPP

/*
  Timer class. Starts timer when created.
  Stops when told (and computes duration).

  Notes:
  - duration data is considered static after stop() has been called.
*/
#pragma once
#include <ctime>
#include <ratio>
#include <chrono>

class Timer
{
private:

  std::chrono::high_resolution_clock::time_point m_start;
  std::chrono::high_resolution_clock::time_point m_finish;
  std::chrono::duration<double>                  m_duration;

  bool m_done;

public:

  Timer();
  ~Timer(){};
  void stop(); // stop timer
  double duration();
};

#endif // TIMER_HPP
