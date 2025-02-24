#pragma once
#include <chrono>

namespace func {

/**
  \brief Starts a timer when created. Stops when `stop()` is called and returns
   the duration in seconds with `duration()`. 

  \note The timer's duration is static after calling `stop()`.
  \note This measures the elapsed time in seconds.
*/
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
