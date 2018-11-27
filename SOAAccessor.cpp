#include <string>
#include <tuple>
#include <iostream>
#include <random>
#include <cmath>

#include "field_descriptor.h"
#include "field_arrays.h"
#include "packed_field_arrays.h"

enum ParticleField
{
	PARTICLE_RX,
	PARTICLE_RY,
	PARTICLE_RZ,
	PARTICLE_E,
	PARTICLE_ATOM_TYPE,
	PARITCLE_MOLECULE_ID,

	PARTICLE_DIST,

	PARTICLE_FIELD_COUNT
};

//#pragma omp declare simd
static inline void compute_distance( double& dist, double x, double y, double z, double x2, double y2, double z2 )
{
	x = x - x2;
	y = y - y2;
	z = z - z2;
	dist = std::exp(x2+y2+z2) / std::sqrt( x*x + y*y + z*z );
}

template<typename OperatorT, typename... ArrayP>
static inline void apply_to_arrays( OperatorT f, size_t first, size_t N, ArrayP __restrict__ ... arraypack )
{
	size_t last = first+N;
	#pragma omp simd
	for(size_t i=first;i<last;i++)
	{
		f( arraypack[i] ... );
	}
}

template<typename OperatorT, typename... FieldDescriptors >
static inline void apply_to_field_arrays( OperatorT f, size_t first, size_t N, FieldArrays<FieldDescriptors...> field_arrays)
{
	size_t last = first+N;
	#pragma omp simd
	for(size_t i=first;i<last;i++)
	{
		f( field_arrays.get(FieldDescriptors())[i] ... );
	}
}

int main(int argc, char* argv[])
{
	static constexpr size_t CHUNK_SIZE = 16;

	int seed = 0;
	size_t N = 10000;

	if(argc>=2) { N=atoi(argv[1]); }
	if(argc>=3) { seed=atoi(argv[2]); }

	FieldDataDescriptor<double,PARTICLE_RX> rx("rx");
	FieldDataDescriptor<double,PARTICLE_RY> ry("ry");
	FieldDataDescriptor<double,PARTICLE_RZ> rz("rz");
	FieldDataDescriptor<double,PARTICLE_E> e("e");
	FieldDataDescriptor<unsigned char,PARTICLE_ATOM_TYPE> atype("type");
	FieldDataDescriptor<int,PARITCLE_MOLECULE_ID> mid("molecule");
	FieldDataDescriptor<double,PARTICLE_DIST> dist("distance");

	auto cell_arrays1 = make_field_arrays( rx,ry,rz,e,dist );
	auto cell_arrays2 = soatl::make_packed_field_arrays( atype,rx,mid,ry,rz );
	// rebind operator ?
	// zip operator ?
	// embed several field_arrays ?

	// zip arrays

	cell_arrays1.allocate(N);
	cell_arrays2.resize(N);

	std::cout<<"atype " << (void*) cell_arrays2.get(atype) << " / " << (void*)( cell_arrays2.get(atype) + cell_arrays2.capacity() ) << std::endl;
	std::cout<<"rx " << (void*) cell_arrays2.get(rx) << " / " << (void*)( cell_arrays2.get(rx) + cell_arrays2.capacity() ) << std::endl;
	std::cout<<"mid " << (void*) cell_arrays2.get(mid) << " / " << (void*)( cell_arrays2.get(mid) + cell_arrays2.capacity() ) << std::endl;
	std::cout<<"ry " << (void*) cell_arrays2.get(ry) << " / " << (void*)( cell_arrays2.get(ry) + cell_arrays2.capacity() ) << std::endl;
	std::cout<<"rz " << (void*) cell_arrays2.get(rz) << " / " << (void*)( cell_arrays2.get(rz) + cell_arrays2.capacity() ) << std::endl;


	double* __restrict__ rx_ptr = static_cast<double*>( __builtin_assume_aligned( cell_arrays1.get(rx) , 64 ) );
	double* __restrict__ ry_ptr = static_cast<double*>( __builtin_assume_aligned( cell_arrays1.get(ry) , 64 ) );
	double* __restrict__ rz_ptr = static_cast<double*>( __builtin_assume_aligned( cell_arrays1.get(rz) , 64 ) );
	double* __restrict__ rx2_ptr = static_cast<double*>( __builtin_assume_aligned( cell_arrays2.get(rx) , 64 ) );
	double* __restrict__ ry2_ptr = static_cast<double*>( __builtin_assume_aligned( cell_arrays2.get(ry) , 64 ) );
	double* __restrict__ rz2_ptr = static_cast<double*>( __builtin_assume_aligned( cell_arrays2.get(rz) , 64 ) );
	double* __restrict__ e_ptr = static_cast<double*>( __builtin_assume_aligned( cell_arrays1.get(e) , 64 ) );
	unsigned char* __restrict__ at_ptr = static_cast<unsigned char*>( __builtin_assume_aligned( cell_arrays2.get(atype) , 64 ) );
	int* __restrict__ m_ptr = static_cast<int*>( __builtin_assume_aligned( cell_arrays2.get(mid) , 64 ) );
	double* __restrict__ dist_ptr = static_cast<double*>( __builtin_assume_aligned( cell_arrays1.get(dist) , 64 ) );

	std::mt19937 rng(seed);
	std::uniform_real_distribution<> rdist(0.0,1.0);

	for(size_t i=0;i<N;i++)
	{
		rx_ptr[i] = rdist(rng);
		ry_ptr[i] = rdist(rng);
		rz_ptr[i] = rdist(rng);
		rx2_ptr[i] = rdist(rng);
		ry2_ptr[i] = rdist(rng);
		rz2_ptr[i] = rdist(rng);
		e_ptr[i] = rdist(rng);
		at_ptr[i] = static_cast<unsigned int>( rdist(rng)*50 );
		m_ptr[i] = static_cast<unsigned int>( at_ptr[i] + rdist(rng)*500 );
		dist_ptr[i] = static_cast<unsigned int>( at_ptr[i] + rdist(rng)*500 );
	}
	
	apply_to_arrays( compute_distance , 0, N, dist_ptr, rx_ptr, ry_ptr, rz_ptr, rx2_ptr, ry2_ptr, rz2_ptr );
	//apply_to_field_arrays( compute_distance , 0, N, cell_arrays1 );

	for(size_t j=0;j<N;j++)
	{
		std::cout<<"dist["<<j<<"]="<<dist_ptr[j]<<std::endl;
	}
}


