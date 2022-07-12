#ifndef MAKE_UNIQUE_HPP
#define MAKE_UNIQUE_HPP
/*
 * implementation of make_unique for use with c++11
 * This is currently never used
 */
namespace std {
  template<typename T, typename... Args>
  std::unique_ptr<T> make_unique(Args&&... args)
  {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
  }
}
#endif // MAKE_UNIQUE_HPP
