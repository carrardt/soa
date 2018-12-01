#include <string>
#include <tuple>
#include <iostream>
#include <random>
#include <cmath>

#include "soatl/field_descriptor.h"
#include "soatl/field_arrays.h"
#include "soatl/packed_field_arrays.h"
#include "soatl/variadic_template_utils.h"
#include "soatl/field_pointers.h"

#include "declare_fields.h"

std::default_random_engine rng;

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


int main(int argc, char* argv[])
{
	int seed = 0;
	size_t N = 10000;

	if(argc>=2) { N=atoi(argv[1]); }
	if(argc>=3) { seed=atoi(argv[2]); }

	rng.seed( seed );

	auto rx = particle_rx;
	auto ry = particle_ry;
	auto rz = particle_rz;
	auto e = particle_e;
	auto atype = particle_atype;
	auto mid = particle_mid;
	auto tmp1 = particle_tmp1;
	auto tmp2 = particle_tmp2;
	auto dist = particle_dist;

	auto cell_arrays1 = soatl::make_field_arrays( rx,ry,rz,e,dist );
	auto cell_arrays2 = soatl::make_packed_field_arrays( soatl::cst::align<64>(), soatl::cst::chunk<8>(), atype,rx,mid,ry,rz );

	cell_arrays1.resize(N);
	cell_arrays2.resize(N);

	double* __restrict__ rx_ptr = static_cast<double*>( __builtin_assume_aligned( cell_arrays1[rx] , 64 ) );
	double* __restrict__ ry_ptr = static_cast<double*>( __builtin_assume_aligned( cell_arrays1[ry] , 64 ) );
	double* __restrict__ rz_ptr = static_cast<double*>( __builtin_assume_aligned( cell_arrays1[rz] , 64 ) );
	double* __restrict__ rx2_ptr = static_cast<double*>( __builtin_assume_aligned( cell_arrays2[rx] , 64 ) );
	double* __restrict__ ry2_ptr = static_cast<double*>( __builtin_assume_aligned( cell_arrays2[ry] , 64 ) );
	double* __restrict__ rz2_ptr = static_cast<double*>( __builtin_assume_aligned( cell_arrays2[rz] , 64 ) );
	double* __restrict__ e_ptr = static_cast<double*>( __builtin_assume_aligned( cell_arrays1[e] , 64 ) );
	unsigned char* __restrict__ at_ptr = static_cast<unsigned char*>( __builtin_assume_aligned( cell_arrays2[atype] , 64 ) );
	int* __restrict__ m_ptr = static_cast<int*>( __builtin_assume_aligned( cell_arrays2[mid] , 64 ) );
	double* __restrict__ dist_ptr = static_cast<double*>( __builtin_assume_aligned( cell_arrays1[dist] , 64 ) );

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

	for(size_t j=0;j<N;j++)
	{
		std::cout<<"dist["<<j<<"]="<<dist_ptr[j]<<std::endl;
	}


	auto zip = soatl::make_field_pointers( N, soatl::cst::align<64>(), soatl::cst::chunk<8>(), rx, ry, rz, dist );

	zip[rx] = cell_arrays2[rx];
	zip[ry] = cell_arrays2[ry];
	zip[rz] = cell_arrays2[rz];
	zip[dist] = cell_arrays1[dist];

	soatl::copy( zip , cell_arrays1 , 0, N, rx, ry, rz );

	double ax = rx2_ptr[0];
	double ay = ry2_ptr[0];
	double az = rz2_ptr[0];

	apply_to_arrays( [ax,ay,az](float& d, double x, double y, double z)
			 {
			     x -= ax;
			     y -= ay;
			     z -= az;
			     d = std::sqrt(x*x+y*y+z*z);
			 }
			 , 0, N,
			 zip[dist], zip[rx], zip[ry], zip[rz]
			 );
}


