#pragma once

#include <cstdint> // for size_t
#include <cstdlib> // for size_t
#include <memory>

#include <assert.h>

#include "soatl/constants.h"

namespace soatl {

template<size_t A, size_t TI, typename... FieldDescriptors>
struct PackedFieldArraysHelper
{
	static constexpr size_t Alignment = A;
	static constexpr size_t AlignmentLowMask = Alignment - 1;
	static constexpr size_t AlignmentHighMask = ~AlignmentLowMask;

	static inline size_t field_offset( size_t capacity )
	{
		using ElementType = typename std::decay< decltype( std::get<TI-1>( std::tuple<FieldDescriptors...>() )  ) >::type;
		size_t offset = PackedFieldArraysHelper<Alignment,TI-1,FieldDescriptors...>::field_offset( capacity );
		offset += ( ( capacity * sizeof(ElementType) ) + AlignmentLowMask ) & AlignmentHighMask;
		return offset;
	}
};

template<size_t A,typename... FieldDescriptors>
struct PackedFieldArraysHelper<A,0,FieldDescriptors...>
{
	static inline constexpr size_t field_offset(size_t) { return 0; }
};


template< size_t _Alignment, size_t _ChunkSize, typename... FieldDescriptors>
struct PackedFieldArrays
{
	static constexpr size_t AlignmentLog2 = Log2<_Alignment>::value;
	static constexpr size_t Alignment = (1ul<<AlignmentLog2);
	static constexpr size_t AlignmentLowMask = Alignment - 1;
	static constexpr size_t AlignmentHighMask = ~AlignmentLowMask;
	static constexpr size_t ChunkSize = (_ChunkSize<1) ? 1 : _ChunkSize;
	static constexpr int TupleSize = sizeof...(FieldDescriptors);

	static constexpr size_t alignment() { return Alignment; }
	static constexpr size_t chunksize() { return ChunkSize; }

	template<typename _T,int _Id>
	inline typename FieldDataDescriptor<_T,_Id>::value_type* get(FieldDataDescriptor<_T,_Id>) const
	{
		using ValueType = typename FieldDataDescriptor<_T,_Id>::value_type;
		static constexpr size_t index = FindFieldIndex<_Id,FieldDescriptors...>::index;
		uint8_t* aptr = static_cast<uint8_t*>( m_aligned_ptr );
		ValueType* typed_pointer = reinterpret_cast<ValueType*>( aptr + PackedFieldArraysHelper<Alignment,index,FieldDescriptors...>::field_offset(capacity()) );
		return typed_pointer;
	}

	inline void adjustCapacity()
	{
		adjustCapacity( size() );
	}

	inline void resize(size_t s)
	{
		if( s > m_capacity || ( s <= (m_capacity/2) && m_capacity >= 2*ChunkSize ) )
		{
			adjustCapacity( s );
		}
		m_size = s;
	}

	inline void* data() const { return m_aligned_ptr; }
	inline size_t size() const { return m_size; }
	inline size_t capacity() const { return m_capacity; }

	inline ~PackedFieldArrays()
	{
		resize(0);
	}

private:

	static inline size_t unaligned_allocation_size(size_t capacity)
	{
		using ElementType = typename std::decay< decltype( std::get<TupleSize-1>( std::tuple<FieldDescriptors...>() )  ) >::type;
		return PackedFieldArraysHelper<Alignment,TupleSize-1,FieldDescriptors...>::field_offset(capacity) + capacity * sizeof(ElementType);
	}

	static inline size_t aligned_allocation_size(size_t capacity)
	{
		size_t alloc_size = unaligned_allocation_size( capacity );
		return ( alloc_size + AlignmentLowMask ) & AlignmentHighMask;
	}

	inline void adjustCapacity(size_t s)
	{
		reallocate( ( ( s + ChunkSize - 1 ) / ChunkSize ) * ChunkSize );
	}

	inline void reallocate(size_t s)
	{
		assert( ( s % ChunkSize ) == 0 );
		size_t needed_space = unaligned_allocation_size( s );
		size_t total_space = aligned_allocation_size( s );
		m_storage_ptr = realloc( m_storage_ptr , total_space );
		size_t aligned_addr = reinterpret_cast<size_t>( m_storage_ptr );
		aligned_addr = ( aligned_addr + AlignmentLowMask ) & AlignmentHighMask;
		m_aligned_ptr = reinterpret_cast<void*>( aligned_addr );
		m_capacity = s;
		assert( m_aligned_ptr!=nullptr || m_capacity==0 );
	}

	void* m_aligned_ptr = nullptr; // start of data (real pointer to use, aligned)
	void* m_storage_ptr = nullptr; // start of allocation (only usefull for deletion)

	size_t m_size = 0;	
	size_t m_capacity = 0;
};

template<typename... FieldDescriptors>
inline
PackedFieldArrays<DEFAULT_ALIGNMENT,DEFAULT_CHUNK_SIZE,FieldDescriptors...>
make_packed_field_arrays(const FieldDescriptors&... fdt)
{
	return PackedFieldArrays<DEFAULT_ALIGNMENT,DEFAULT_CHUNK_SIZE,FieldDescriptors...>();
}

template<size_t A, size_t C, typename... FieldDescriptors>
inline
PackedFieldArrays<A,C,FieldDescriptors...>
make_packed_field_arrays( cst::align<A>, cst::chunk<C>, const FieldDescriptors&... fdt)
{
	return PackedFieldArrays<A,C,FieldDescriptors...>();
}


} // namespace soatl

