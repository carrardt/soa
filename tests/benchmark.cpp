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

std::default_random_engine rng;

template<typename ArraysT, size_t idDist, size_t idRx, size_t idRy, size_t idRz>
inline double benchmark(ArraysT& arrays, size_t N, bool use_simd, soatl::FieldId<idDist> dist, soatl::FieldId<idRx> rx, soatl::FieldId<idRy> ry, soatl::FieldId<idRz> rz)
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
    if( use_simd )
    {
	    soatl::apply_simd( [ax,ay,az](DistT& d, PosT x, PosT y, PosT z)
		    {
			    x = x - ax;
			    y = y - ay;
			    z = z - az;
			    d = std::sqrt( x*x + y*y + z*z );
		    }
		    , arrays, dist, rx, ry, rz );
    }
    else
    {
	    soatl::apply( [ax,ay,az](DistT& d, PosT x, PosT y, PosT z)
		    {
			    x = x - ax;
			    y = y - ay;
			    z = z - az;
			    d = std::sqrt( x*x + y*y + z*z );
		    }
		    , arrays, dist, rx, ry, rz );
    }
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
	bool double_precision = true;
	bool use_simd = false;

  if(argc>=2) 
	{
	  std::string option = argv[1];
	  if( option == "fa" ) { arraysImpl = FIELD_ARRAYS; }
	  else if( option == "pfa" ) { arraysImpl = PACKED_FIELD_ARRAYS; }
	  else if( option == "spfa" ) { arraysImpl = STATIC_PACKED_FIELD_ARRAYS; }
	}

  if(argc>=3) 
	{
	  std::string option = argv[2];
	  if( option=="s" ) { double_precision = false; }
	}

  if(argc>=4) 
	{
	  std::string option = argv[3];
	  if( option=="vec" ) { use_simd = true; }
	}

	if(argc>=5)
	{
	  N = atoi(argv[4]);
	}
	if(argc>=6)
	{
	  seed = atoi(argv[5]);
	}

  if(arraysImpl == STATIC_PACKED_FIELD_ARRAYS) { N = S; }

  std::cout<<"SIMD arch="<<soatl::simd_arch();
  if(double_precision)
  {
    std::cout<<" a="<< soatl::SimdRequirements<double>::alignment <<" c=" << soatl::SimdRequirements<double>::chunksize;
  }
  else
  {
    std::cout<<" a="<< soatl::SimdRequirements<float>::alignment <<" c=" << soatl::SimdRequirements<float>::chunksize;
  }
  std::cout<<", type=";
	switch( arraysImpl )
	{
	  case FIELD_ARRAYS: std::cout<<"FieldArrays"; break;
	  case PACKED_FIELD_ARRAYS: std::cout<<"PackedFieldArrays"; break;
	  case STATIC_PACKED_FIELD_ARRAYS: std::cout<<"StaticPackedFieldArrays"; break;
  }
  std::cout<<", dp="<<double_precision<<", vec="<<use_simd<<", N="<<N<<", seed="<<seed<<std::endl;
  
  
	rng.seed( seed );

  double result = 0.0;
	switch( arraysImpl )
	{
	  case FIELD_ARRAYS:
	    if(double_precision)
	    {
	      auto arrays = soatl::make_field_arrays( soatl::cst::align<TEST_ALIGNMENT>(), soatl::cst::chunk<TEST_CHUNK_SIZE>(), particle_rx, particle_ry, particle_rz, particle_e);
	      result = benchmark(arrays,N,use_simd,particle_e,particle_rx, particle_ry, particle_rz);
	    }
	    else
	    {
	      auto arrays = soatl::make_field_arrays( soatl::cst::align<TEST_ALIGNMENT>(), soatl::cst::chunk<TEST_CHUNK_SIZE>(), particle_rx_f, particle_ry_f, particle_rz_f, particle_e_f);
	      result = benchmark(arrays,N,use_simd,particle_e_f,particle_rx_f, particle_ry_f, particle_rz_f);
	    }
	    break;

	  case PACKED_FIELD_ARRAYS:
	    if(double_precision)
	    {
	      auto arrays = soatl::make_packed_field_arrays( soatl::cst::align<TEST_ALIGNMENT>(), soatl::cst::chunk<TEST_CHUNK_SIZE>(), particle_rx, particle_ry, particle_rz, particle_e);
	      result = benchmark(arrays,N,use_simd,particle_e,particle_rx, particle_ry, particle_rz);
	    }
	    else
	    {
	      auto arrays = soatl::make_packed_field_arrays( soatl::cst::align<TEST_ALIGNMENT>(), soatl::cst::chunk<TEST_CHUNK_SIZE>(), particle_rx_f, particle_ry_f, particle_rz_f, particle_e_f);
	      result = benchmark(arrays,N,use_simd,particle_e_f,particle_rx_f, particle_ry_f, particle_rz_f);
	    }
	    break;

	  case STATIC_PACKED_FIELD_ARRAYS:
	    if(double_precision)
	    {
	      auto arrays = soatl::make_static_packed_field_arrays( soatl::cst::align<TEST_ALIGNMENT>(), soatl::cst::chunk<TEST_CHUNK_SIZE>(), soatl::cst::count<S>(),
	                    particle_rx, particle_ry, particle_rz, particle_e);
	      result = benchmark(arrays,S,use_simd,particle_e,particle_rx, particle_ry, particle_rz);
	    }
	    else
	    {
	      auto arrays = soatl::make_static_packed_field_arrays( soatl::cst::align<TEST_ALIGNMENT>(), soatl::cst::chunk<TEST_CHUNK_SIZE>(),  soatl::cst::count<S>(),
	                    particle_rx_f, particle_ry_f, particle_rz_f, particle_e_f);
	      result = benchmark(arrays,S,use_simd,particle_e_f,particle_rx_f, particle_ry_f, particle_rz_f);
	    }
	    break;
	}

  std::cout<<"result = "<<result<<std::endl;

	return 0;
}


