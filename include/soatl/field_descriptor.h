#pragma once

#include <cstdlib> // for size_t

namespace soatl
{

template<size_t _id> struct FieldId { static constexpr size_t id =_id; };
template<size_t _index> struct FieldIndex { static constexpr size_t index = _index; };
template<size_t _id> struct FieldDescriptor {};

template<size_t k, size_t... ids> struct find_index_of_id {};
template<size_t k, size_t f, size_t... ids>
struct find_index_of_id<k,f,ids...> { static constexpr size_t index = (k==f) ? 0 : (1+find_index_of_id<k,ids...>::index) ; };
template<size_t k>
struct find_index_of_id<k> { static constexpr size_t index = 0; };

}


