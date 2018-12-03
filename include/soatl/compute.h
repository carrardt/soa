#pragma once

#include "soatl/field_descriptor.h"

namespace soatl
{

template<typename OperatorT, typename... T>
static inline void apply( OperatorT f, size_t first, size_t N, T* ... arraypack )
{
	size_t last = first+N;
	for(size_t i=first;i<last;i++)
	{
		f( arraypack[i] ... );
	}
}

/*
template<typename OperatorT, typename FieldArraysT, size_t... ids>
static inline void apply( OperatorT f, size_t first, size_t N, FieldArraysT& arrays, const FieldId<ids> & fids ... )
{
	apply(f,first,N, arrays[fids] ... );
}
*/

template<typename OperatorT, typename... T>
static inline void apply_simd( OperatorT f, size_t first, size_t N, T* __restrict__ ... arraypack )
{
	size_t last = first+N;
	#pragma omp simd
	for(size_t i=first;i<last;i++)
	{
		f( arraypack[i] ... );
	}
}


template<typename OperatorT, typename FieldArraysT, size_t... ids>
static inline void apply_simd( OperatorT f, size_t first, size_t N, FieldArraysT& arrays, const FieldId<ids>& ... fids )
{
	// debug mode : check that alignment and (non-)aliasing is suitable for simd vectorization
#	ifndef NDEBUG
	{
		size_t alignment = arrays.alignment();
		assert( alignment >= 16 ); // FIXME: should be replaced by some value defined by the system and compiler
		const void* ptr;
		size_t addr;
		TEMPLATE_LIST_BEGIN
			ptr = arrays[FieldId<ids>()] ,
			addr = reinterpret_cast<size_t>(ptr) ,
			assert( ( addr % alignment ) == 0 ) 
		TEMPLATE_LIST_END
		// TODO: check pointer aliasing
	}
#	endif

	apply_simd(f,first,N, arrays[fids] ... );
}


template<typename OperatorT, typename FieldArraysT, size_t... ids>
static inline void apply_simd( OperatorT f, FieldArraysT& arrays, const FieldId<ids>& ... fids )
{
	apply_simd( f, 0, arrays.size(), arrays[fids] ... );
}


} // namespace soatl


