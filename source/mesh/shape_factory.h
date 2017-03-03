#ifndef SHAPE_FACTORY_H
#define SHAPE_FACTORY_H


enum class element_type {
  element_tri
};

enum class edge_type {

};

  // taken from http://stackoverflow.com/questions/28256550/factory-pattern-using-variadic-template
  // Factory method that uses SFINAE to determine if inputs for a given element type are an error or not.
  // If the type is invalid, the create function will simply return a null_ptr
template <typename T, typename... Args>
  std::enable_if_t<std::is_constructible<T,Args...>::value,
                   std::unique_ptr<T>
		  >
    create(Args&&... args) {
  return std::make_unique<T>(std::forward<Args>(args)...);
  }

template <typename T, typename... Args>
  std::enable_if_t<!std::is_constructible<T,Args...>::value,
                   std::unique_ptr<T>
		  >
    create(Args&&... ) {
    return nullptr;
  }

template<class... Args>
static ELEMENT* create_elt(int elt_typ, Args&&... args )
{
  switch (elt_typ) {
  case element_type::element_tri:
    return create<ELEMENT_TRI>(std::forward<Args>(args)...);
  default:
    return nullptr;
  }
}

template<class... Args>
static EDGE* create_edg(element_type edg_typ, Args&&... args )
{
  switch (edg_typ) {
  default:
    return nullptr;
  }
}

#endif
