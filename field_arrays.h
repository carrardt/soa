#pragma once

#include <cstdlib> // for size_t
#include "field_descriptor.h"
#include <memory> // for std::alocator_traits

#include "assert.h"

#include "soatl_constants.h"

namespace soatl {

template<int N> struct _FieldArraysTupleHelper;
template<int N> struct _FieldArraysTupleHelper
{
	template<typename... FieldDescriptors>
	static inline void init( std::tuple<FieldDescriptors...>& arrays)
	{
		std::get<N-1>( arrays ) = nullptr;
		_FieldArraysTupleHelper<N-1>::init( arrays );
	}

	template<typename... FieldDescriptors>
	static inline void reallocate( std::tuple<FieldDescriptors...>& arrays, size_t s)
	{
		using ValueType = typename std::decay< decltype( * std::get<N-1>( arrays ) ) >::type;
		auto ptr = std::get<N-1>( arrays );
		if( ptr != nullptr ) { delete [] ptr; }
		if( s > 0 ) { std::get<N-1>( arrays ) = new ValueType [ s ]; }
		else { std::get<N-1>( arrays ) = nullptr;  }
		_FieldArraysTupleHelper<N-1>::reallocate( arrays, s );
	}
};

template<> struct _FieldArraysTupleHelper<0>
{
	template<typename... FieldDescriptors>
	static inline void reallocate( std::tuple<FieldDescriptors...>& arrays, size_t s) {}

	template<typename... FieldDescriptors>
	static inline void init( std::tuple<FieldDescriptors...>& arrays) {}
};

template<size_t _ChunkSize, typename... FieldDescriptors>
struct FieldArrays
{
	static constexpr size_t ChunkSize = (_ChunkSize<1) ? 1 : _ChunkSize;
	using ArrayTuple = std::tuple< typename FieldDescriptors::value_type* ... > ;
	using TupleHelper = _FieldArraysTupleHelper< sizeof...(FieldDescriptors) >;

	inline FieldArrays()
	{
		TupleHelper::init(m_field_arrays);
	}

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

	inline void adjustCapacity()
	{
		reallocate( ( ( m_size + ChunkSize - 1 ) / ChunkSize ) * ChunkSize );
	}

	inline void resize(size_t s)
	{
		m_size = s;
		if( m_size > m_capacity || ( m_size <= (m_capacity/2) && m_capacity >= 2*ChunkSize ) )
		{
			adjustCapacity();
		}
	}

	inline size_t size() const { return m_size; }
	inline size_t capacity() const { return m_capacity; }

	static inline constexpr size_t alignment() { return 1; }
	static inline constexpr size_t chunksize() { return ChunkSize; }

private:
	inline void reallocate(size_t s)
	{
		assert( ( s % ChunkSize ) == 0 );
		m_capacity = s;
		TupleHelper::reallocate( m_field_arrays , m_capacity );
	}

	ArrayTuple m_field_arrays;
	size_t m_size = 0;
	size_t m_capacity = 0;
};

template<typename... FieldDescriptors>
inline
FieldArrays<DEFAULT_CHUNK_SIZE,FieldDescriptors...> make_field_arrays(FieldDescriptors...)
{
	return FieldArrays<DEFAULT_CHUNK_SIZE,FieldDescriptors...>();
}

template<size_t C,typename... FieldDescriptors>
inline
FieldArrays<C,FieldDescriptors...> make_field_arrays(cst::chunk<C>, FieldDescriptors...)
{
	return FieldArrays<C,FieldDescriptors...>();
}

} // namespace soatl

