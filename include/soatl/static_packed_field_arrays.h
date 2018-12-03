#pragma once

#include <cstdint> // for size_t
#include <cstdlib> // for size_t
#include <memory>
#include <tuple>

#include <assert.h>

#include "soatl/constants.h"
#include "soatl/copy.h"
#include "soatl/memory.h"

namespace soatl {


template<size_t A, size_t TI, size_t capacity, size_t... ids>
struct StaticPackedFieldArraysHelper
{
	static constexpr size_t Alignment = A;
	static constexpr size_t AlignmentLowMask = Alignment - 1;
	static constexpr size_t AlignmentHighMask = ~AlignmentLowMask;
	using ElementType = typename std::tuple_element< TI-1 , std::tuple< typename FieldDescriptor<ids>::value_type ... > >::type ;
	static constexpr size_t offset = StaticPackedFieldArraysHelper<Alignment,TI-1,capacity,ids...>::offset
				        + ( ( capacity * sizeof(ElementType) ) + AlignmentLowMask ) & AlignmentHighMask;
};

template<size_t A, size_t capacity, size_t... ids>
struct StaticPackedFieldArraysHelper<A,0,capacity,ids...>
{
	static constexpr size_t offset = 0;
};

template< size_t _Alignment, size_t _ChunkSize, size_t _NChunks, size_t... ids>
struct alignas(_Alignment) StaticPackedFieldArrays
{
	static constexpr bool assert_alignment_is_power_of_2 = AssertPowerOf2<_Alignment>::value;
	static constexpr size_t AlignmentLog2 = Log2<_Alignment>::value;
	static constexpr size_t Alignment = (1ul<<AlignmentLog2);
	static constexpr size_t AlignmentLowMask = Alignment - 1;
	static constexpr size_t AlignmentHighMask = ~AlignmentLowMask;
	static constexpr size_t ChunkSize = (_ChunkSize<1) ? 1 : _ChunkSize;
	static constexpr size_t Size = _NChunks * ChunkSize;
	static constexpr int TupleSize = sizeof...(ids);

	
	using FieldIdsTuple = std::tuple< FieldId<ids> ... > ;

	using LastValueType = typename std::tuple_element<TupleSize-1,std::tuple< typename FieldDescriptor<ids>::value_type ... > >::type ;
	static constexpr size_t AllocationSize = StaticPackedFieldArraysHelper<Alignment,TupleSize-1,Size,ids...>::offset + Size * sizeof(LastValueType);

	template<size_t _id>
	inline typename FieldDescriptor<_id>::value_type * __restrict__ operator [] ( FieldId<_id> ) 
	{
		static constexpr size_t index = find_index_of_id<_id,ids...>::index;
		using ValueType = typename std::tuple_element<index, std::tuple< typename FieldDescriptor<ids>::value_type ... > >::type ;
		uint8_t* aptr = m_storage + PackedFieldArraysHelper<Alignment,index,ids...>::field_offset(capacity()) ;
		return (ValueType* __restrict__) __builtin_assume_aligned( aptr , Alignment );
	}

	template<size_t _id>
	inline const typename FieldDescriptor<_id>::value_type * __restrict__ operator [] ( FieldId<_id> ) const
	{
		static constexpr size_t index = find_index_of_id<_id,ids...>::index;
		using ValueType = typename std::tuple_element<index, std::tuple< typename FieldDescriptor<ids>::value_type ... > >::type ;
		const uint8_t* aptr = m_storage + PackedFieldArraysHelper<Alignment,index,ids...>::field_offset(capacity()) ;
		return (const ValueType* __restrict__) __builtin_assume_aligned( aptr , Alignment );
	}

	static inline constexpr size_t alignment() { return Alignment; }
	static inline constexpr size_t chunksize() { return ChunkSize; }
	static inline constexpr size_t size() { return Size; }
	static inline constexpr size_t capacity() { return Size; }
	static inline constexpr size_t chunk_ceil() { return Size; }
	static inline constexpr size_t data_size() { return AllocationSize; }
	inline uint8_t* data() { return m_storage; }
	inline const uint8_t* data() const { return m_storage; }

private:
	uint8_t m_storage[ AllocationSize ];
};

template<size_t A, size_t C, size_t N, size_t... ids>
inline
StaticPackedFieldArrays<A,C,N,ids...>
make_static_packed_field_arrays( cst::align<A>, cst::chunk<C>, cst::count<N>, const FieldId<ids>& ...)
{
	return StaticPackedFieldArrays<A,C,((N+C-1)/C)*C,ids...>();
}

template<size_t N, size_t... ids>
inline
StaticPackedFieldArrays<DEFAULT_ALIGNMENT,DEFAULT_CHUNK_SIZE,N,ids...>
make_static_packed_field_arrays( cst::count<N>, const FieldId<ids>& ...)
{
	return StaticPackedFieldArrays<DEFAULT_ALIGNMENT,DEFAULT_CHUNK_SIZE,((N+DEFAULT_CHUNK_SIZE-1)/DEFAULT_CHUNK_SIZE)*DEFAULT_CHUNK_SIZE,ids...>();
}

} // namespace soatl

