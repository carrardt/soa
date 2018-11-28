#pragma once

#include <cstdint> // for size_t
#include <cstdlib> // for size_t
#include <memory>

#include <assert.h>

#include "soatl/constants.h"

namespace soatl {

template< size_t _AlignmentLog2, size_t _ChunkSize, typename... FieldDescriptors>
struct PackedFieldArrays
{
	static constexpr size_t Alignment = (_AlignmentLog2<3) ? 8 : (1ul<<_AlignmentLog2);
	static constexpr size_t AlignmentLowMask = ( 1ul << _AlignmentLog2 ) - 1;
	static constexpr size_t AlignmentHighMask = ~AlignmentLowMask;
	static constexpr size_t ChunkSize = (_ChunkSize<1) ? 1 : _ChunkSize;
	static constexpr int TupleSize = sizeof...(FieldDescriptors);

	static constexpr size_t alignment() { return Alignment; }
	static constexpr size_t chunksize() { return ChunkSize; }

	template<typename _T,int _Id>
	inline typename FieldDataDescriptor<_T,_Id>::value_type* get(FieldDataDescriptor<_T,_Id>)
	{
		using ValueType = typename FieldDataDescriptor<_T,_Id>::value_type;
		static constexpr int index = FindFieldIndex<_Id,FieldDescriptors...>::index;
		uint8_t* aptr = static_cast<uint8_t*>( m_aligned_ptr );
		ValueType* typed_pointer = reinterpret_cast<ValueType*>( aptr + field_offset( std::integral_constant<int,index>() ) );
		return typed_pointer;
	}

	template<int TI>
	inline size_t field_offset(std::integral_constant<int,TI>) const
	{
		using ElementType = typename std::decay< decltype( std::get<TI-1>( std::tuple<FieldDescriptors...>() )  ) >::type;
		size_t offset = field_offset( std::integral_constant<int,TI-1>() );
		offset += ( ( m_capacity * sizeof(ElementType) ) + AlignmentLowMask ) & AlignmentHighMask;
		return offset;
	}
	inline size_t field_offset(std::integral_constant<int,0>) const
	{
		return 0;
	}

	inline size_t unaligned_allocation_size() const
	{
		using ElementType = typename std::decay< decltype( std::get<TupleSize-1>( std::tuple<FieldDescriptors...>() )  ) >::type;
		return field_offset( std::integral_constant<int,TupleSize-1>() ) + m_capacity * sizeof(ElementType);
	}

	inline size_t aligned_allocation_size() const
	{
		size_t alloc_size = unaligned_allocation_size();
		return ( alloc_size + AlignmentLowMask ) & AlignmentHighMask;
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

	inline ~PackedFieldArrays()
	{
		resize(0);
	}

	inline size_t size() const { return m_size; }
	inline size_t capacity() const { return m_capacity; }

private:
	inline void reallocate(size_t s)
	{
		assert( ( s % ChunkSize ) == 0 );
		m_capacity = s;
		size_t needed_space = unaligned_allocation_size();
		size_t total_space = aligned_allocation_size();
		m_storage_ptr = realloc( m_storage_ptr , total_space );
		size_t aligned_addr = reinterpret_cast<size_t>( m_storage_ptr );
		aligned_addr = ( aligned_addr + AlignmentLowMask ) & AlignmentHighMask;
		m_aligned_ptr = reinterpret_cast<void*>( aligned_addr );
		assert( m_aligned_ptr!=nullptr || m_capacity==0 );
	}

	size_t m_size = 0;	
	size_t m_capacity = 0;

	void* m_storage_ptr = nullptr; // start of allocation (only usefull for deletion)
	void* m_aligned_ptr = nullptr; // start of data (real pointer to use, aligned)
};

template<typename... FieldDescriptors>
inline
PackedFieldArrays<DEFAULT_ALIGNMENT_LOG2,DEFAULT_CHUNK_SIZE,FieldDescriptors...>
make_packed_field_arrays(const FieldDescriptors&... fdt)
{
	return PackedFieldArrays<DEFAULT_ALIGNMENT_LOG2,DEFAULT_CHUNK_SIZE,FieldDescriptors...>();
}

template<size_t A, size_t C, typename... FieldDescriptors>
inline
PackedFieldArrays<_priv::Log2<A>::value,C,FieldDescriptors...>
make_packed_field_arrays( cst::align<A>, cst::chunk<C>, const FieldDescriptors&... fdt)
{
	return PackedFieldArrays<_priv::Log2<A>::value,C,FieldDescriptors...>();
}


} // namespace soatl

