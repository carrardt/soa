#pragma once

#include <cstdint> // for size_t
#include <cstdlib> // for size_t
#include <memory>

#include <assert.h>

#include "soatl/constants.h"
#include "soatl/copy.h"
#include "soatl/variadic_template_utils.h"

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

	using FieldDescriptorTuple = std::tuple<FieldDescriptors...> ;

	static constexpr size_t alignment() { return Alignment; }
	static constexpr size_t chunksize() { return ChunkSize; }

	template<size_t index>
	inline typename std::tuple_element<index,std::tuple<FieldDescriptors...> >::type::value_type * get( cst::at<index> ) const
	{
		using ValueType = typename std::tuple_element<index,std::tuple<FieldDescriptors...> >::type::value_type ;
		uint8_t* aptr = static_cast<uint8_t*>( m_storage_ptr );
		ValueType* typed_pointer = reinterpret_cast<ValueType*>( aptr + PackedFieldArraysHelper<Alignment,index,FieldDescriptors...>::field_offset(capacity()) );
		return typed_pointer;
	}

	template<typename _T,int _Id>
	inline typename FieldDataDescriptor<_T,_Id>::value_type* get( FieldDataDescriptor<_T,_Id> ) const
	{
		static constexpr size_t index = FindFieldIndex<_Id,FieldDescriptors...>::index;
		return get( cst::at<index>() );
	}

	inline void adjustCapacity()
	{
		adjustCapacity( size() );
	}

	inline void resize(size_t s)
	{
		if( s > m_capacity || s <= (m_capacity-ChunkSize) || s==0 )
		{
			adjustCapacity( s );
		}
		m_size = s;
	}

	inline void* data() const { return m_storage_ptr; }
	inline size_t size() const { return m_size; }
	inline size_t capacity() const { return m_capacity; }

	inline ~PackedFieldArrays()
	{
		resize(0);
	}

private:

	static inline size_t allocation_size(size_t capacity)
	{
		using ElementType = typename std::decay< decltype( std::get<TupleSize-1>( std::tuple<FieldDescriptors...>() )  ) >::type;
		return PackedFieldArraysHelper<Alignment,TupleSize-1,FieldDescriptors...>::field_offset(capacity) + capacity * sizeof(ElementType);
	}

	inline void adjustCapacity(size_t s)
	{
		size_t newCapacity = ( ( s + ChunkSize - 1 ) / ChunkSize ) * ChunkSize ;
		if( newCapacity != m_capacity ) { reallocate( newCapacity ); }
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
		soatl::copy( tmp, *this, 0, cs, FieldDescriptors()... );
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

