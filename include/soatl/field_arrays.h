#pragma once

#include <cstdlib> // for size_t
#include <memory> // for std::alocator_traits
#include <algorithm> // for std::max

#include "assert.h"

#include "soatl/field_descriptor.h"
#include "soatl/constants.h"
#include "soatl/memory.h"
#include "soatl/simd.h"

/*
The  function posix_memalign() allocates size bytes and places the address of the allocated memory in *memptr.  The address of the allocated memory will be a multiple of alignment, which must
be a power of two and a multiple of sizeof(void *).  If size is 0, then the value placed in *memptr is either NULL, or a unique pointer value that can later be successfully passed to free(3).
*/

namespace soatl {

// TODO: add alignment
template<size_t _Alignment, size_t _ChunkSize, typename... ids >
struct FieldArrays
{
	static constexpr bool assert_alignment_is_power_of_2 = AssertPowerOf2<_Alignment>::value;
  static constexpr size_t AlignmentLog2 = Log2<_Alignment>::value;
  static constexpr size_t Alignment = (1ul<<AlignmentLog2);
  static constexpr size_t AlignmentLowMask = Alignment - 1;
  static constexpr size_t AlignmentHighMask = ~AlignmentLowMask;
	static constexpr size_t ChunkSize = (_ChunkSize<1) ? 1 : _ChunkSize;
	static constexpr size_t TupleSize = sizeof...(ids);

	using FieldIdsTuple = std::tuple< FieldId<ids> ... > ;
	using ArrayTuple = std::tuple< typename FieldDescriptor<ids>::value_type* ... > ;
	using AllocStrategy = DefaultAllocationStrategy;

	inline FieldArrays()
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

	inline void resize(size_t s)
	{
		if( s != m_size )
		{
			size_t new_capacity = AllocStrategy::update_capacity(s,capacity(),chunksize());
			if( new_capacity != m_capacity )
			{
				reallocate( new_capacity );
			}
			m_size = s;
		}
	}

	inline size_t size() const { return m_size; }
	inline size_t chunk_ceil() const { return ( (size()+chunksize()-1) / chunksize() ) * chunksize(); }
	inline size_t capacity() const { return m_capacity; }

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

	template<size_t N>
	inline void reallocate_pointer( std::integral_constant<size_t,N> , size_t s )
	{
		using ValueType = typename std::decay< decltype( * std::get<N-1>( m_field_arrays ) ) >::type;
		ValueType * old_ptr = std::get<N-1>( m_field_arrays );
		ValueType * new_ptr = nullptr;
		if( s > 0 )
		{
			// here, we use posix_memalign (and not realloc) to benefit from aligned memory allocation provided by system library
			size_t a = std::max( alignment() , sizeof(void*) ); // this is required by posix_memalign.
			void* memptr = nullptr;
			int r = posix_memalign( &memptr, a, s*sizeof(ValueType) );
			assert( r == 0 );
			new_ptr = reinterpret_cast<ValueType*>( memptr );
		}
		if( old_ptr!=nullptr && new_ptr!=nullptr )
		{
			// we don't use realloc, so we have to copy elements
			size_t cs = std::min(s,m_size);
			for(size_t i=0;i<cs;i++) { new_ptr[i] = old_ptr[i]; }
		}
		std::get<N-1>( m_field_arrays ) = new_ptr;
		if( old_ptr != nullptr ) { free(old_ptr); }
		reallocate_pointer( std::integral_constant<size_t,N-1>() , s ); // recursion to next array
	}
	inline void reallocate_pointer( std::integral_constant<size_t,0> , size_t s ) {}

	inline void reallocate(size_t s)
	{
		assert( ( s % ChunkSize ) == 0 );
		reallocate_pointer( std::integral_constant<size_t,TupleSize>() , s );
		m_capacity = s;
	}

	ArrayTuple m_field_arrays;

	size_t m_size = 0;
	size_t m_capacity = 0;
};

template<typename... ids>
inline
FieldArrays<DEFAULT_ALIGNMENT,DEFAULT_CHUNK_SIZE,ids...> make_field_arrays(const FieldId<ids>& ...)
{
	return FieldArrays<DEFAULT_ALIGNMENT,DEFAULT_CHUNK_SIZE,ids...>();
}

template<size_t A, size_t C,typename... ids>
inline
FieldArrays<A,C,ids...> make_field_arrays(cst::align<A>, cst::chunk<C>, const FieldId<ids>& ...)
{
	return FieldArrays<A,C,ids...>();
}

} // namespace soatl

