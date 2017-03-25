#ifndef SHAPE_FACTORY_H
#define SHAPE_FACTORY_H
#include "../elements/class_element_2D.h"

// http://stackoverflow.com/questions/28256550/factory-pattern-using-variadic-template
// possible to construct
template <typename T, typename... Args>
std::enable_if_t<
    std::is_constructible<T, Args...>::value,
    std::unique_ptr<T>
>
create(Args&&... args) {
    return std::make_unique<T>(std::forward<Args>(args)...);
}

// impossible to construct
template <typename T, typename... Args>
std::enable_if_t<
    !std::is_constructible<T, Args...>::value,
    std::unique_ptr<T>
>
create(Args&&... ) {
    return nullptr;
}

enum class element_type {
  _2d
};

template<class... Args>
static std::unique_ptr<ELEMENT> create_elt(element_type elt_typ, Args&&... args )
{
  switch (elt_typ) {
  case element_type::_2d:
    return create<ELEMENT_2D>(std::forward<Args>(args)...);
  default:
    return nullptr;
  }
}
#endif
