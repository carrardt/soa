#pragma once

#include <cstdlib> // for size_t
#include <memory> // for std::alocator_traits

#include "soatl/field_descriptor.h"
#include "soatl/constants.h"

namespace soatl {

// this class hould be private, i.e. not instantiated by user but only through specific operations
template<size_t _Alignment, size_t _ChunkSize, typename... FieldDescriptors>
struct FieldPointers
{
	static constexpr size_t ChunkSize = _ChunkSize;
	static constexpr size_t Alignment = _Alignment;
	static constexpr size_t TupleSize = sizeof...(FieldDescriptors) ;

	using ArrayTuple = std::tuple< typename FieldDescriptors::value_type* ... > ;

	inline FieldPointers(size_t s) : m_size(s)
	{
		init( std::integral_constant<size_t,TupleSize>() );
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

	inline size_t size() const { return m_size; }
	inline size_t capacity() const { return ( (size()+chunksize()-1) / chunksize() ) * chunksize(); }

	static inline constexpr size_t alignment() { return Alignment; }
	static inline constexpr size_t chunksize() { return ChunkSize; }

private:

        template<size_t N>
        inline void init( std::integral_constant<size_t,N> )
        {
                std::get<N-1>( m_field_arrays ) = nullptr;
                init( std::integral_constant<size_t,N-1>() );
        }
        inline void init( std::integral_constant<size_t,0> ) {}

	ArrayTuple m_field_arrays;
	size_t m_size = 0;
};

template<size_t A, size_t C, typename... FieldDescriptors>
static inline 
FieldPointers<A,C,FieldDescriptors...>
make_field_pointers( size_t N, cst::align<A>, cst::chunk<C>, const FieldDescriptors& ... )
{
	return FieldPointers<A,C,FieldDescriptors...>( N );
}


}


