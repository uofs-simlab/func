/*
  Timer class. Starts timer when created.
  Stops when told (and computes duration).

  Notes:
  - duration data is considered static after stop() has been called.
  - time is measured in seconds
*/
#pragma once
#include <chrono>

namespace func {
class Timer
{
  std::chrono::high_resolution_clock::time_point m_start;
  std::chrono::high_resolution_clock::time_point m_finish;
  std::chrono::duration<double>                  m_duration;

  bool m_done;

public:

  Timer() : m_duration(0), m_done(false)
  {
    m_start = std::chrono::high_resolution_clock::now();
  }
  ~Timer(){};

  // stop timer
  void stop()
  {
    m_finish = std::chrono::high_resolution_clock::now();

    if(!m_done){
      // store (finish - start) as time in seconds
      m_duration = std::chrono::duration_cast<std::chrono::duration<double>>(
          m_finish - m_start);
      m_done = true;
    }
  }

  double duration(){ return m_duration.count(); };
};
}
