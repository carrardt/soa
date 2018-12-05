#pragma once

#include <cstdlib> // for size_t
#include <type_traits>

namespace soatl
{

template<typename _field_id> struct FieldDescriptor
{
  using value_type = void;
  using Id = _field_id;
  static const char* name() { return "<uknown>"; }
};
template<typename T> using FieldId = FieldDescriptor<T>; // for backward compatiblity

template<typename k, typename... ids> struct find_index_of_id {};
template<typename k, typename f, typename... ids>
struct find_index_of_id<k,f,ids...> { static constexpr size_t index = std::is_same<k,f>::value ? 0 : (1+find_index_of_id<k,ids...>::index) ; };
template<typename k>
struct find_index_of_id<k> { static constexpr size_t index = 0; };

}


