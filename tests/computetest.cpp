#include <string>
#include <tuple>
#include <iostream>
#include <random>
#include <cmath>
#include <typeinfo>

#include "soatl/field_descriptor.h"
#include "soatl/field_arrays.h"
#include "soatl/packed_field_arrays.h"
#include "soatl/variadic_template_utils.h"
#include "soatl/field_pointers.h"

#include "declare_fields.h"

std::default_random_engine rng;

static inline void compute_distance( float& dist, double x, double y, double z, double x2, double y2, double z2 )
{
	x = x - x2;
	y = y - y2;
	z = z - z2;
	dist = /*std::exp(x2+y2+z2) / std::sqrt*/ ( x*x + y*y + z*z );
}

static inline void compute_distance_d( double& dist, double x, double y, double z, double x2, double y2, double z2 )
{
	x = x - x2;
	y = y - y2;
	z = z - z2;
	dist = std::exp(x2+y2+z2) / std::sqrt ( x*x + y*y + z*z );
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

template<typename T, size_t id>
static inline void print_field_info(const T& arrays, soatl::FieldId<id> f)
{
	auto ptr = arrays[f];
	std::cout<<soatl::FieldDescriptor<id>::name()<<" : array type is "<<typeid(ptr).name()<<std::endl;
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

	auto cell_arrays1 = soatl::make_field_arrays( soatl::cst::align<64>(), soatl::cst::chunk<8>(), rx,ry,rz,e,dist );
	auto cell_arrays2 = soatl::make_packed_field_arrays( soatl::cst::align<64>(), soatl::cst::chunk<8>(), atype,rx,mid,ry,rz );

	std::cout<<"resize arrays to "<<N<<std::endl;  std::cout.flush();
	cell_arrays1.resize(N);
	cell_arrays2.resize(N);

	std::cout<<"fields info "<<std::endl;  std::cout.flush();
	print_field_info( cell_arrays1, rx );
	print_field_info( cell_arrays1, ry );
	print_field_info( cell_arrays1, rz );
	print_field_info( cell_arrays1, e );
	print_field_info( cell_arrays1, dist );
	print_field_info( cell_arrays2, atype );
	print_field_info( cell_arrays2, mid );

	std::cout<<"get pointers"<<std::endl;  std::cout.flush();
	auto rx_ptr = cell_arrays1[rx] ;
	auto ry_ptr = cell_arrays1[ry] ;
	auto rz_ptr = cell_arrays1[rz] ;
	auto rx2_ptr = cell_arrays2[rx] ;
	auto ry2_ptr = cell_arrays2[ry] ;
	auto rz2_ptr = cell_arrays2[rz] ;
	auto e_ptr = cell_arrays1[e] ;
	auto at_ptr = cell_arrays2[atype] ;
	auto m_ptr = cell_arrays2[mid] ;
	auto dist_ptr = cell_arrays1[dist] ;

	std::cout<<"initialize values"<<std::endl;  std::cout.flush();
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
		dist_ptr[i] = 0.0;
	}
	
	std::cout<<"apply compute_distance to arrays (float)"<<std::endl;  std::cout.flush();
	apply_to_arrays( compute_distance , 0, N, dist_ptr, rx_ptr, ry_ptr, rz_ptr, rx2_ptr, ry2_ptr, rz2_ptr );
	for(size_t j=0;j<N;j++)
	{
		std::cout<<"dist["<<j<<"]="<<dist_ptr[j]<<std::endl; std::cout.flush();
	}

	std::cout<<"apply compute_distance to arrays (double)"<<std::endl;  std::cout.flush();
	apply_to_arrays( compute_distance_d , 0, N, e_ptr, rx_ptr, ry_ptr, rz_ptr, rx2_ptr, ry2_ptr, rz2_ptr );
	for(size_t j=0;j<N;j++)
	{
		std::cout<<"e["<<j<<"]="<<e_ptr[j]<<std::endl; std::cout.flush();
	}


	std::cout<<"borrow array pointers"<<std::endl; std::cout.flush();
	auto zip = soatl::make_field_pointers( N, soatl::cst::align<64>(), soatl::cst::chunk<8>(), rx, ry, rz, e, dist );

	zip[rx] = cell_arrays2[rx];
	zip[ry] = cell_arrays2[ry];
	zip[rz] = cell_arrays2[rz];
	zip[dist] = cell_arrays1[dist];
	zip[e] = cell_arrays1[e];

	std::cout<<"copy arrays"<<std::endl;
	soatl::copy( zip , cell_arrays1 , 0, N, rx, ry, rz );

	double ax = rx2_ptr[0];
	double ay = ry2_ptr[0];
	double az = rz2_ptr[0];

	std::cout<<"apply lambda (float)"<<std::endl;
	apply_to_arrays( [ax,ay,az](float& d, double x, double y, double z)
			 {
			     x = ax - x;
			     y = ay - y;
			     z = az - z;
			     d = /*std::sqrt*/ (x*x+y*y+z*z);
			 }
			 , 0, N,
			 zip[dist], zip[rx], zip[ry], zip[rz]
			 );
	for(size_t j=0;j<N;j++)
	{
		std::cout<<"dist["<<j<<"]="<<dist_ptr[j]<<std::endl; std::cout.flush();
	}

	std::cout<<"apply lambda (double)"<<std::endl;
	apply_to_arrays( [ax,ay,az](double& d, double x, double y, double z)
			 {
			     x = ax - x;
			     y = ay - y;
			     z = az - z;
			     d = /*std::sqrt*/ (x*x+y*y+z*z);
			 }
			 , 0, N,
			 zip[e], zip[rx], zip[ry], zip[rz]
			 );
	for(size_t j=0;j<N;j++)
	{
		std::cout<<"e["<<j<<"]="<<e_ptr[j]<<std::endl; std::cout.flush();
	}

}


