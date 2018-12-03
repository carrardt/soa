#pragma once

#include "soatl/field_descriptor.h"
#include "soatl/variadic_template_utils.h"
#include "soatl/simd.h"


namespace soatl
{

// Non-SIMD versions

template<typename OperatorT, typename... T>
static inline void apply( OperatorT f, size_t N, T* __restrict__ ... arraypack )
{
#	ifndef NDEBUG
	check_pointers_aliasing( N , arraypack ... );
#	endif

	for(size_t i=0;i<N;i++)
	{
		f( arraypack[i] ... );
	}
}

template<typename OperatorT, typename... T>
static inline void apply( OperatorT f, size_t first, size_t N, T* __restrict__ ... arraypack )
{
	apply( f, N, arraypack+first ... );
}

template<typename OperatorT, typename FieldArraysT, size_t... ids>
static inline void apply( OperatorT f, size_t first, size_t N, FieldArraysT& arrays, const FieldId<ids> & ... fids )
{
	apply( f, N, arrays[fids]+first ... );
}


// SIMD versions

template<typename OperatorT, typename... T>
static inline void apply_simd( OperatorT f, size_t N, T* __restrict__ ... arraypack )
{
#	ifndef NDEBUG
	check_simd_pointers( N , arraypack ... );
#	endif

#	pragma omp simd
	for(size_t i=0;i<N;i++)
	{
		f( arraypack[i] ... );
	}
}

template<typename OperatorT, size_t VECSIZE, typename... T>
static inline void apply_simd( OperatorT f, size_t N, cst::chunk<VECSIZE>, T* __restrict__ ... arraypack )
{
#	ifndef NDEBUG
	check_simd_pointers( N , arraypack ... );
#	endif

	for(size_t i=0;i<N;i+=VECSIZE)
	{
#		pragma omp simd
		for(size_t j=0;j<VECSIZE;j++)
		{
			f( arraypack[i+j] ... );
		}
	}
}

template<typename OperatorT, typename... T>
static inline void apply_simd( OperatorT f, size_t N, cst::chunk<1>, T* __restrict__ ... arraypack )
{
	apply_simd(f,N,arraypack...);
}

template<typename OperatorT, typename... T>
static inline void apply_simd( OperatorT f, size_t first, size_t N, T* __restrict__ ... arraypack )
{
	apply_simd( f, N, arraypack+first ... );
}

template<typename OperatorT, typename FieldArraysT, size_t... ids>
static inline void apply_simd( OperatorT f, size_t first, size_t N, FieldArraysT& arrays, const FieldId<ids>& ... fids )
{
	apply_simd( f, N, cst::chunk<FieldArraysT::ChunkSize>(), arrays[fids]+first ... );
}

template<typename OperatorT, typename FieldArraysT, size_t... ids>
static inline void apply_simd( OperatorT f, size_t N, FieldArraysT& arrays, const FieldId<ids>& ... fids )
{
	apply_simd( f, N, cst::chunk<FieldArraysT::ChunkSize>(), arrays[fids] ... );
}

template<typename OperatorT, typename FieldArraysT, size_t... ids>
static inline void apply_simd( OperatorT f, FieldArraysT& arrays, const FieldId<ids>& ... fids )
{
	apply_simd( f, arrays.size(), cst::chunk<FieldArraysT::ChunkSize>(), arrays[fids] ... );
}


} // namespace soatl


