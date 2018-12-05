#pragma once

#include <cstdint> // for size_t
#include <cstdlib> // for size_t
#include <memory>
#include <tuple>

#include <assert.h>

#include "soatl/constants.h"
#include "soatl/copy.h"
#include "soatl/memory.h"
#include "soatl/simd.h"

namespace soatl {

template<size_t A, size_t TI, typename... ids>
struct PackedFieldArraysHelper
{
	static constexpr size_t Alignment = A;
	static constexpr size_t AlignmentLowMask = Alignment - 1;
	static constexpr size_t AlignmentHighMask = ~AlignmentLowMask;
	using ElementType = typename std::tuple_element< TI-1 , std::tuple< typename FieldDescriptor<ids>::value_type ... > >::type ;

	static inline size_t field_offset( size_t capacity )
	{
		return PackedFieldArraysHelper<Alignment,TI-1,ids...>::field_offset( capacity )
		       + ( ( capacity * sizeof(ElementType) ) + AlignmentLowMask ) & AlignmentHighMask;
	}
};

template<size_t A, typename... ids>
struct PackedFieldArraysHelper<A,0,ids...>
{
	static inline constexpr size_t field_offset(size_t) { return 0; }

	template<size_t capacity>
	static inline constexpr size_t field_offset( std::integral_constant<size_t,capacity> capacityConstant ) { return 0; }
};


template< size_t _Alignment, size_t _ChunkSize, typename... ids>
struct PackedFieldArrays
{
	static constexpr bool assert_alignment_is_power_of_2 = AssertPowerOf2<_Alignment>::value;
	static constexpr size_t AlignmentLog2 = Log2<_Alignment>::value;
	static constexpr size_t Alignment = (1ul<<AlignmentLog2);
	static constexpr size_t AlignmentLowMask = Alignment - 1;
	static constexpr size_t AlignmentHighMask = ~AlignmentLowMask;
	static constexpr size_t ChunkSize = (_ChunkSize<1) ? 1 : _ChunkSize;
	static constexpr int TupleSize = sizeof...(ids);

	using FieldIdsTuple = std::tuple< FieldId<ids> ... > ;
	using AllocStrategy = DefaultAllocationStrategy;

	static constexpr size_t alignment() { return Alignment; }
	static constexpr size_t chunksize() { return ChunkSize; }

	template<size_t index>
	inline typename std::tuple_element<index,std::tuple< typename FieldDescriptor<ids>::value_type ... > >::type * __restrict__ operator [] ( cst::at<index> ) const
	{
		using ValueType = typename std::tuple_element<index, std::tuple< typename FieldDescriptor<ids>::value_type ... > >::type ;
		uint8_t* aptr = static_cast<uint8_t*>( m_storage_ptr ) + PackedFieldArraysHelper<Alignment,index,ids...>::field_offset(capacity()) ;
		return (ValueType* __restrict__) __builtin_assume_aligned( aptr , Alignment );
	}

	template<typename _id>
	inline typename FieldDescriptor<_id>::value_type * __restrict__ operator [] ( FieldId<_id> ) const
	{
		static constexpr size_t index = find_index_of_id<_id,ids...>::index;
		using ValueType = typename std::tuple_element<index, std::tuple< typename FieldDescriptor<ids>::value_type ... > >::type ;
		uint8_t* aptr = static_cast<uint8_t*>( m_storage_ptr ) + PackedFieldArraysHelper<Alignment,index,ids...>::field_offset(capacity()) ;
		return (ValueType* __restrict__) __builtin_assume_aligned( aptr , Alignment );
	}

	// resize container
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

	inline void* data() const { return m_storage_ptr; }
	inline size_t size() const { return m_size; }
	inline size_t capacity() const { return m_capacity; }
	inline size_t chunk_ceil() const { return ( (size()+chunksize()-1) / chunksize() ) * chunksize(); }
	inline size_t data_size() const { return allocation_size( capacity() ); }

	inline ~PackedFieldArrays()
	{
		resize(0);
	}

private:

	static inline size_t allocation_size(size_t capacity)
	{
		using ValueType = typename std::tuple_element<TupleSize-1,std::tuple< typename FieldDescriptor<ids>::value_type ... > >::type ;
		return PackedFieldArraysHelper<Alignment,TupleSize-1,ids...>::field_offset(capacity) + capacity * sizeof(ValueType);
	}

	inline void reallocate(size_t s)
	{
		assert( ( s % ChunkSize ) == 0 );

		size_t total_space = allocation_size( s );
		assert( ( s==0 && total_space==0 ) || ( s!=0 && total_space!=0 ) );

		void * new_ptr = nullptr;
		if( total_space > 0 )
		{
			// here, we use posix_memalign (and not realloc) to benefit from aligned memory allocation provided by system library
			size_t a = std::max( alignment() , sizeof(void*) ); // this is required by posix_memalign.
			int r = posix_memalign( &new_ptr, a, total_space );
			assert( r == 0 );
		}

		size_t cs = std::min(s,m_size);
		assert( ( m_storage_ptr!=nullptr && new_ptr!=nullptr ) || ( cs == 0 ) );
		
		// copy here
		PackedFieldArrays tmp;
		tmp.m_storage_ptr = new_ptr;
		tmp.m_size = cs;
		tmp.m_capacity = s;
		soatl::copy( tmp, *this, 0, cs, FieldId<ids>()... );
		tmp.m_storage_ptr = nullptr;
		tmp.m_size = 0;
		tmp.m_capacity = 0;

		if( m_storage_ptr!=nullptr ) { free(m_storage_ptr); }
		m_storage_ptr = new_ptr;
		m_capacity = s;
		assert( m_storage_ptr!=nullptr || m_capacity==0 );
	}

	void* m_storage_ptr = nullptr; // start of allocation (only usefull for deletion)
	size_t m_size = 0;	
	size_t m_capacity = 0;
};

template<typename... ids>
inline
PackedFieldArrays<DEFAULT_ALIGNMENT,DEFAULT_CHUNK_SIZE,ids...>
make_packed_field_arrays(const FieldId<ids>& ...)
{
	return PackedFieldArrays<DEFAULT_ALIGNMENT,DEFAULT_CHUNK_SIZE,ids...>();
}

template<size_t A, size_t C, typename... ids>
inline
PackedFieldArrays<A,C,ids...>
make_packed_field_arrays( cst::align<A>, cst::chunk<C>, const FieldId<ids>& ...)
{
	return PackedFieldArrays<A,C,ids...>();
}


} // namespace soatl

