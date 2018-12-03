#pragma once

namespace soatl
{

static constexpr size_t DEFAULT_ALIGNMENT = 64;
static constexpr size_t DEFAULT_CHUNK_SIZE = 16;

namespace cst
{
	template<size_t> struct align {};
	template<size_t> struct chunk {};
	template<size_t> struct at {};
	template<size_t> struct count {};
}

template<size_t N> struct Log2 { static constexpr size_t value = Log2<N/2>::value+1; };
template<> struct Log2<1> { static constexpr size_t value = 0; };
template<> struct Log2<0> { static constexpr size_t value = 0; };

template<size_t n, bool v = n==(1ul<<Log2<n>::value) > struct AssertPowerOf2 {};
template<size_t n> struct AssertPowerOf2<n,true> { static constexpr bool value = true; };

} // namespace soatl

