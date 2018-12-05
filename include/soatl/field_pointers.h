#pragma once

#include <cstdlib> // for size_t
#include <memory> // for std::alocator_traits

#include "soatl/field_descriptor.h"
#include "soatl/constants.h"

namespace soatl {

// this class hould be private, i.e. not instantiated by user but only through specific operations
template<size_t _Alignment, size_t _ChunkSize, typename... ids>
struct FieldPointers
{
	static constexpr size_t ChunkSize = _ChunkSize;
	static constexpr size_t Alignment = _Alignment;
	static constexpr size_t TupleSize = sizeof...(ids) ;

	using FieldIdsTuple = std::tuple< FieldId<ids> ... > ;
	using ArrayTuple = std::tuple< typename FieldDescriptor<ids>::value_type* ... > ;

	inline FieldPointers(size_t s) : m_size(s)
	{
		init( std::integral_constant<size_t,TupleSize>() );
	}

	template<typename _id>
	inline typename FieldDescriptor<_id>::value_type* __restrict__ operator [] ( FieldId<_id> ) const
	{
		using ValueType = typename FieldDescriptor<_id>::value_type;
		static constexpr int index = find_index_of_id<_id,ids...>::index;
		return (ValueType* __restrict__) __builtin_assume_aligned( std::get<index>(m_field_arrays) , Alignment );
	}

	template<typename _id>
	inline void set_pointer( FieldId<_id> , typename FieldDescriptor<_id>::value_type* ptr )
	{
		static constexpr int index = find_index_of_id<_id,ids...>::index;
#		ifndef NDEBUG
			size_t addr = reinterpret_cast<size_t>( ptr );
			assert( (addr%Alignment) == 0 );
#		endif
		std::get<index>(m_field_arrays) = ptr;
	}

	template<typename ArraySet, typename... other_ids>
	inline void set_pointers( ArraySet& arrays , const FieldId<other_ids>& ...  )
	{
		static_assert( ( arrays.chunksize() % chunksize() ) == 0 , "Cannot copy pointer from an array with a smaller chunksize" );
		TEMPLATE_LIST_BEGIN
			set_pointer( FieldId<other_ids>() , arrays[FieldId<other_ids>()] )
		TEMPLATE_LIST_END
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
	const size_t m_size;
};



// helper methods to create FieldPointers instances
template<size_t A, size_t C, typename... ids>
static inline 
FieldPointers<A,C,ids...>
make_field_pointers( size_t N, cst::align<A>, cst::chunk<C>, const FieldId<ids>& ... )
{
	return FieldPointers<A,C,ids...>( N );
}

template<size_t A, size_t C, typename... ids>
static inline 
FieldPointers<A,C,ids...>
make_field_pointers( size_t N, cst::align<A>, cst::chunk<C>, const std::tuple<FieldId<ids>...>& )
{
	return FieldPointers<A,C,ids...>( N );
}


template<typename ArraySetT, typename... ids>
static inline 
FieldPointers<ArraySetT::Alignment,ArraySetT::ChunkSize,ids...>
make_field_pointers( const ArraySetT& arrays, const FieldId<ids>& ... )
{
	return FieldPointers<ArraySetT::Alignment,ArraySetT::ChunkSize,ids...>( arrays.size() );
}


template<typename ArraySetT>
static inline 
auto
make_field_pointers( const ArraySetT& arrays )
{
	return make_field_pointers( arrays.size(), cst::align<ArraySetT::Alignment>(), cst::chunk<ArraySetT::ChunkSize>(), ArraySetT::FieldIdsTuple() );
}

}


