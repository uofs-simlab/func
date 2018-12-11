/*
  Implementation for Timer class.
 */

#include "Timer.hpp"

Timer::Timer() : m_duration(0), m_done(false)
{
  m_start = std::chrono::high_resolution_clock::now();
}

void Timer::stop()
{
  m_finish = std::chrono::high_resolution_clock::now();

  if (!m_done) {
    m_duration = std::chrono::
      duration_cast< std::chrono::duration<double> >
      ( m_finish - m_start );
    m_done = true;
}
}

double Timer::duration()
{
  return m_duration.count();
}
