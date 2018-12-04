#include <string>
#include <iostream>
#include <random>
#include <cmath>
#include <chrono>

#include "soatl/field_descriptor.h"
#include "soatl/field_arrays.h"
#include "soatl/packed_field_arrays.h"
#include "soatl/variadic_template_utils.h"
#include "soatl/field_pointers.h"
#include "soatl/static_packed_field_arrays.h"
#include "soatl/compute.h"

#include "declare_fields.h"

#ifndef TEST_ALIGNMENT
#define TEST_ALIGNMENT 64
#endif

#ifndef TEST_CHUNK_SIZE
#define TEST_CHUNK_SIZE 16
#endif

#ifndef TEST_USE_SIMD
#define TEST_USE_SIMD 1
#endif

#ifndef TEST_DOUBLE_PRECISION
#define TEST_DOUBLE_PRECISION 1
#endif

#if TEST_DOUBLE_PRECISION
auto field_e = particle_e;
auto field_rx = particle_rx;
auto field_ry = particle_ry;
auto field_rz = particle_rz;
#else
auto field_e = particle_e_f;
auto field_rx = particle_rx_f;
auto field_ry = particle_ry_f;
auto field_rz = particle_rz_f;
#endif

std::default_random_engine rng;

template<typename ArraysT, size_t idDist, size_t idRx, size_t idRy, size_t idRz>
inline double benchmark(ArraysT& arrays, size_t N, soatl::FieldId<idDist> dist, soatl::FieldId<idRx> rx, soatl::FieldId<idRy> ry, soatl::FieldId<idRz> rz)
{
  static constexpr size_t nCycles = 100;
  using DistT = typename soatl::FieldDescriptor<idDist>::value_type;
  using PosT = typename soatl::FieldDescriptor<idRx>::value_type;

  std::uniform_real_distribution<> rdist(0.0,1.0);
  double result = 0.0;

  std::chrono::nanoseconds timens(0);

  for(size_t cycle=0;cycle<nCycles;cycle++)
  {
	  arrays.resize(N);
	  for(size_t i=0;i<N;i++)
	  {
		  arrays[rx][i] = rdist(rng);
		  arrays[ry][i] = rdist(rng);
		  arrays[rz][i] = rdist(rng);
    }

    PosT ax = arrays[rx][0];
    PosT ay = arrays[ry][0];
    PosT az = arrays[rz][0];

    auto t1 = std::chrono::high_resolution_clock::now();
#   if TEST_USE_SIMD 
	    soatl::apply_simd( [ax,ay,az](DistT& d, PosT x, PosT y, PosT z)
		    {
			    x = x - ax;
			    y = y - ay;
			    z = z - az;
			    d = std::sqrt( x*x + y*y + z*z );
		    }
		    , arrays, dist, rx, ry, rz );
#   else
	    soatl::apply( [ax,ay,az](DistT& d, PosT x, PosT y, PosT z)
		    {
			    x = x - ax;
			    y = y - ay;
			    z = z - az;
			    d = std::sqrt( x*x + y*y + z*z );
		    }
		    , arrays, dist, rx, ry, rz );
#   endif
    auto t2 = std::chrono::high_resolution_clock::now();
    timens += t2-t1;
    		  
    // try to optimize this out, compiler !
    for(size_t i=0;i<N;i++)
    {
      if(rdist(rng)>0.999) { result += arrays[dist][i]; }
    }
	}

  std::cout<<"time = "<<timens.count()/nCycles<<std::endl;

  return result;
}

enum ArraysImplementation
{
	FIELD_ARRAYS,
	PACKED_FIELD_ARRAYS,
	STATIC_PACKED_FIELD_ARRAYS
};

int main(int argc, char* argv[])
{
	static constexpr size_t S=300;
	int seed = 0;
	size_t N = 10000;
	ArraysImplementation arraysImpl = FIELD_ARRAYS;

  if(argc>=2)
	{
	  std::string option = argv[1];
	  if( option == "fa" ) { arraysImpl = FIELD_ARRAYS; }
	  else if( option == "pfa" ) { arraysImpl = PACKED_FIELD_ARRAYS; }
	  else if( option == "spfa" ) { arraysImpl = STATIC_PACKED_FIELD_ARRAYS; }
	}

	if(argc>=3)
	{
	  N = atoi(argv[2]);
	}

	if(argc>=4)
	{
	  seed = atoi(argv[3]);
	}

  if(arraysImpl == STATIC_PACKED_FIELD_ARRAYS) { N = S; }

  std::cout<<"SIMD arch="<<soatl::simd_arch();
# if TEST_DOUBLE_PRECISION
    std::cout<<" a="<< soatl::SimdRequirements<double>::alignment <<" c=" << soatl::SimdRequirements<double>::chunksize;
# else
    std::cout<<" a="<< soatl::SimdRequirements<float>::alignment <<" c=" << soatl::SimdRequirements<float>::chunksize;
# endif

  std::cout<<", type=";
	switch( arraysImpl )
	{
	  case FIELD_ARRAYS: std::cout<<"FieldArrays"; break;
	  case PACKED_FIELD_ARRAYS: std::cout<<"PackedFieldArrays"; break;
	  case STATIC_PACKED_FIELD_ARRAYS: std::cout<<"StaticPackedFieldArrays"; break;
  }

# if TEST_DOUBLE_PRECISION
  std::cout<<", double";
# else
  std::cout<<", float";
# endif

# if TEST_USE_SIMD 
  std::cout<<", vec";
# else
  std::cout<<", scal";
# endif

  std::cout<<", N="<<N<<", seed="<<seed<<std::endl;

	rng.seed( seed );

  double result = 0.0;

	switch( arraysImpl )
	{
	  case FIELD_ARRAYS:
	    {
	      auto arrays = soatl::make_field_arrays( soatl::cst::align<TEST_ALIGNMENT>(), soatl::cst::chunk<TEST_CHUNK_SIZE>(), field_rx, field_ry, field_rz, field_e);
	      result = benchmark(arrays,N,field_e,field_rx,field_ry,field_rz);
	    }
	    break;

	  case PACKED_FIELD_ARRAYS:
	    {
	      auto arrays = soatl::make_packed_field_arrays( soatl::cst::align<TEST_ALIGNMENT>(), soatl::cst::chunk<TEST_CHUNK_SIZE>(), field_rx, field_ry, field_rz, field_e);
	      result = benchmark(arrays,N,field_e,field_rx,field_ry,field_rz);
	    }
	    break;

	  case STATIC_PACKED_FIELD_ARRAYS:
	    {
	      auto arrays = soatl::make_static_packed_field_arrays( soatl::cst::align<TEST_ALIGNMENT>(), soatl::cst::chunk<TEST_CHUNK_SIZE>(), soatl::cst::count<S>(),
	                    field_rx, field_ry, field_rz, field_e);
	      N = arrays.size();
	      result = benchmark(arrays,N,field_e,field_rx,field_ry,field_rz);
	    }
	    break;
	}

  std::cout<<"result = "<<result<<std::endl;

	return 0;
}


