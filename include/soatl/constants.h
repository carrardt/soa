#pragma once

#include <cstdlib> // for size_t

namespace soatl
{

namespace cst
{
	template<size_t A> struct align { static inline void f(){ static_assert(A>0,"alignment cannot be 0"); } };
	template<size_t C> struct chunk { static inline void f(){ static_assert(C>0,"chunk size cannot be 0"); } };
	template<size_t> struct at {};
	template<size_t> struct count {};

	cst::align<1> unaligned;
	cst::align<1> aligned_1;
	cst::align<2> aligned_2;
	cst::align<4> aligned_4;
	cst::align<8> aligned_8;
	cst::align<16> aligned_16;
	cst::align<32> aligned_32;
	cst::align<64> aligned_64;

	cst::chunk<1> no_chunk;
	cst::chunk<1> chunk_1;
	cst::chunk<2> chunk_2;
	cst::chunk<4> chunk_4;
	cst::chunk<8> chunk_8;
	cst::chunk<16> chunk_16;
}

template<size_t N> struct Log2 { static constexpr size_t value = Log2<N/2>::value+1; };
template<> struct Log2<1> { static constexpr size_t value = 0; };
template<> struct Log2<0> { static constexpr size_t value = 0; };

template<size_t n, bool v = n==(1ul<<Log2<n>::value) > struct AssertPowerOf2 {};
template<size_t n> struct AssertPowerOf2<n,true> { static constexpr bool value = true; };

} // namespace soatl

