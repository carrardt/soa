#pragma once

#include <cstdlib> // for size_t
#include "field_descriptor.h"

template<int N> struct _FieldArraysTupleHelper;
template<int N> struct _FieldArraysTupleHelper
{
	template<typename... FieldDescriptors>
	static inline void allocate( std::tuple<FieldDescriptors...>& arrays, size_t s)
	{
		using ValueType = typename std::decay< decltype( * std::get<N-1>( arrays ) ) >::type;
		std::get<N-1>( arrays ) = new ValueType [ s ];
		_FieldArraysTupleHelper<N-1>::allocate( arrays, s );
	}
};
template<> struct _FieldArraysTupleHelper<0>
{
	template<typename... FieldDescriptors>
	static inline void allocate( std::tuple<FieldDescriptors...>& arrays, size_t s) {}
};


template<typename... FieldDescriptors>
struct FieldArrays
{
	using ArrayTuple = std::tuple< typename FieldDescriptors::value_type* ... > ;
	using TupleHelper = _FieldArraysTupleHelper< sizeof...(FieldDescriptors) >;

	template<typename _T,int _Id>
	inline typename FieldDataDescriptor<_T,_Id>::value_type * get(FieldDataDescriptor<_T,_Id>)
	{
		static constexpr int index = FindFieldIndex<_Id,FieldDescriptors...>::index;
		return std::get<index>(m_field_arrays);
	}

	template<typename _T,int _Id>
	inline const typename FieldDataDescriptor<_T,_Id>::value_type * get(FieldDataDescriptor<_T,_Id>) const
	{
		static constexpr int index = FindFieldIndex<_Id,FieldDescriptors...>::index;
		return std::get<index>(m_field_arrays);
	}

	inline void allocate(size_t s)
	{
		TupleHelper::allocate( m_field_arrays , s );
	}

private:
	
	ArrayTuple m_field_arrays;
};

template<typename... FieldDescriptors>
inline
FieldArrays<FieldDescriptors...> make_field_arrays(FieldDescriptors...)
{
	return FieldArrays<FieldDescriptors...>();
}


