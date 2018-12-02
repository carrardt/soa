#pragma once

#include <cstdlib> // for size_t
#include <memory> // for std::alocator_traits

#include "soatl/field_descriptor.h"
#include "soatl/constants.h"

namespace soatl {

// this class hould be private, i.e. not instantiated by user but only through specific operations
template<size_t _Alignment, size_t _ChunkSize, size_t... FieldIds>
struct FieldPointers
{
	static constexpr size_t ChunkSize = _ChunkSize;
	static constexpr size_t Alignment = _Alignment;
	static constexpr size_t TupleSize = sizeof...(FieldIds) ;

	using FieldIdsTuple = std::tuple< FieldId<FieldIds> ... > ;
	using ArrayTuple = std::tuple< typename FieldDescriptor<FieldIds>::value_type* ... > ;

	inline FieldPointers(size_t s) : m_size(s)
	{
		init( std::integral_constant<size_t,TupleSize>() );
	}

	template<size_t _id>
	inline typename FieldDescriptor<_id>::value_type* __restrict__ operator [] ( FieldId<_id> ) const
	{
		using ValueType = typename FieldDescriptor<_id>::value_type;
		static constexpr int index = find_index_of_id<_id,FieldIds...>::index;
		return (ValueType* __restrict__) __builtin_assume_aligned( std::get<index>(m_field_arrays) , Alignment );
	}

	template<size_t _id>
	inline void set_pointer( FieldId<_id> , typename FieldDescriptor<_id>::value_type* ptr )
	{
		static constexpr int index = find_index_of_id<_id,FieldIds...>::index;
#		ifndef NDEBUG
			size_t addr = reinterpret_cast<size_t>( ptr );
			assert( (addr%Alignment) == 0 );
#		endif
		std::get<index>(m_field_arrays) = ptr;
	}

	inline size_t size() const { return m_size; }
	inline size_t chunk_ceil() const { return ( (size()+chunksize()-1) / chunksize() ) * chunksize(); }
	inline size_t capacity() const { return chunk_ceil(); }

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

template<size_t A, size_t C, size_t... ids>
static inline 
FieldPointers<A,C,ids...>
make_field_pointers( size_t N, cst::align<A>, cst::chunk<C>, const FieldId<ids>& ... )
{
	return FieldPointers<A,C,ids...>( N );
}


}


