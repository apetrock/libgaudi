#ifndef CREATE_FUNC_H
#define CREATE_FUNC_H

#include <memory>

#define DEFINE_CREATE_FUNC(classname)                                          \
  typedef std::shared_ptr<classname> ptr;                                      \
  template <typename... Args> static ptr create(Args &&...args) {              \
    return std::make_shared<classname>(std::forward<Args>(args)...);           \
  }

#endif // CREATE_FUNC_H