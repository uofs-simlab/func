#pragma once
#include <string> // to_string()
#include <memory>
#include <mutex>

// Store a histogram of each location a function has been sampled at
class ArgumentRecord
{
  unsigned int *mp_histogram;
  unique_ptr<std::mutex[]> mp_histogram_mutex;

  unsigned int m_histSize;
  // min and max should always be the same as the table that contains the ArgumentRecord
  double m_min;
  double m_max;

  // statistic plus a mutex for locking
  double m_peak;
  std::mutex m_peak_mutex;

  public:
    ArgumentRecord(unsigned int histSize, double min, double max);

    // place x in the histogram
    void record_arg(double x);

    // make a string representation of the histogram
    std::string to_string();

    // print out the various fields in this class
    void print_details(std::ostream& out);
    void print_details_json(std::ostream& out);
    ~ArgumentRecord();
};
