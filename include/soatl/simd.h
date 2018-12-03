#pragma once

#include <cstdlib> // for size_t
#include <cstdint> // for int64_t, int32_t, etc.
#include <algorithm> // for std::sort
#include <assert.h>
#include "soatl/variadic_template_utils.h"

namespace soatl
{

template<typename T>
struct SimdRequirements
{
	static constexpr size_t alignment = 64;
	static constexpr size_t chunksize = 16;
};
#define SET_ARCH_SIMD_REQUIREMENT(_T,_a,_c) \
	template<> struct SimdRequirements<_T> { \
	static constexpr size_t alignment=_a; \
	static constexpr size_t chunksize=_c; }

#if defined(__AVX512F__)

SET_ARCH_SIMD_REQUIREMENT(double   ,64,8);
SET_ARCH_SIMD_REQUIREMENT(float	   ,64,16);
SET_ARCH_SIMD_REQUIREMENT(int64_t  ,64,8);
SET_ARCH_SIMD_REQUIREMENT(uint64_t ,64,8);
SET_ARCH_SIMD_REQUIREMENT(int32_t  ,64,16);
SET_ARCH_SIMD_REQUIREMENT(uint32_t ,64,16);
SET_ARCH_SIMD_REQUIREMENT(int16_t  ,64,16);
SET_ARCH_SIMD_REQUIREMENT(uint16_t ,64,16);
SET_ARCH_SIMD_REQUIREMENT(int8_t   ,64,16);
SET_ARCH_SIMD_REQUIREMENT(uint8_t  ,64,16);

#elif defined(__AVX2__) || defined(__AVX__)

SET_ARCH_SIMD_REQUIREMENT(double   ,32,4);
SET_ARCH_SIMD_REQUIREMENT(float	   ,32,8);
SET_ARCH_SIMD_REQUIREMENT(int64_t  ,32,4);
SET_ARCH_SIMD_REQUIREMENT(uint64_t ,32,4);
SET_ARCH_SIMD_REQUIREMENT(int32_t  ,32,8);
SET_ARCH_SIMD_REQUIREMENT(uint32_t ,32,8);
SET_ARCH_SIMD_REQUIREMENT(int16_t  ,32,8);
SET_ARCH_SIMD_REQUIREMENT(uint16_t ,32,8);
SET_ARCH_SIMD_REQUIREMENT(int8_t   ,32,8);
SET_ARCH_SIMD_REQUIREMENT(uint8_t  ,32,8);

#elif defined(__SSE4_2__) || defined(__SSE4_1__) || defined(__SSE3__) || defined(__SSE2__) || defined(__SSE__)

SET_ARCH_SIMD_REQUIREMENT(double   ,16,2);
SET_ARCH_SIMD_REQUIREMENT(float	   ,16,4);
SET_ARCH_SIMD_REQUIREMENT(int64_t  ,16,2);
SET_ARCH_SIMD_REQUIREMENT(uint64_t ,16,2);
SET_ARCH_SIMD_REQUIREMENT(int32_t  ,16,4);
SET_ARCH_SIMD_REQUIREMENT(uint32_t ,16,4);
SET_ARCH_SIMD_REQUIREMENT(int16_t  ,16,4);
SET_ARCH_SIMD_REQUIREMENT(uint16_t ,16,4);
SET_ARCH_SIMD_REQUIREMENT(int8_t   ,16,4);
SET_ARCH_SIMD_REQUIREMENT(uint8_t  ,16,4);
static inline const char* simd_arch() { return "SSE"; }

#else
static inline const char* simd_arch() { return "<unknown>"; }
#endif

// prefered default alignment and chunk size on this machine
static constexpr size_t DEFAULT_ALIGNMENT = SimdRequirements<double>::alignment;
static constexpr size_t DEFAULT_CHUNK_SIZE = SimdRequirements<float>::chunksize;

// check correct non-aliasing of pointers
template<typename... T>
static inline void check_pointers_aliasing( size_t N, T* ... arraypack )
{
	static constexpr size_t na = 2 * sizeof...(T);
	size_t addr[ na ];
	size_t i = 0;
	TEMPLATE_LIST_BEGIN
		addr[i] = reinterpret_cast<size_t>(arraypack) ,
		addr[i+1] = addr[i] + N * sizeof(T) ,
		i += 2
	TEMPLATE_LIST_END
	std::sort( addr , addr+na );
	for(i=1;i<na;i++)
	{
		assert( addr[i] >= addr[i-1] && "Arrays' memory areas must not overlap" );
	}
}

// check correct alignment and non-aliasing of pointers
template<typename... T>
static inline void check_simd_pointers( size_t N, T* ... arraypack )
{
	check_pointers_aliasing(N,arraypack...);
	size_t addr;
	TEMPLATE_LIST_BEGIN
		addr = reinterpret_cast<size_t>(arraypack) ,
		assert( ( addr % SimdRequirements<T>::alignment ) == 0 && "Address does not have correct alignment" )
	TEMPLATE_LIST_END
}

} // namespace soatl

