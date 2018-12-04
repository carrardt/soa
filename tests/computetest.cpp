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
#include "soatl/static_packed_field_arrays.h"
#include "soatl/compute.h"

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
	dist = /*std::exp(x2+y2+z2) / std::sqrt*/ ( x*x + y*y + z*z );
}

template<typename T, size_t id>
static inline void print_field_info(const T& arrays, soatl::FieldId<id> f)
{
	auto ptr = arrays[f];
	std::cout<<soatl::FieldDescriptor<id>::name()<<" : array type is "<<typeid(ptr).name()<<std::endl;
}


template<size_t A, size_t C, size_t N, size_t... ids>
static inline void check_static_fields( soatl::StaticPackedFieldArrays<A,C,N,ids...>& field_arrays)
{
	size_t alignment = field_arrays.alignment();
	const void* ptr;
	size_t addr;
	TEMPLATE_LIST_BEGIN
		ptr = field_arrays[soatl::FieldId<ids>()] ,
		addr = reinterpret_cast<size_t>(ptr) ,
		assert( ( addr % alignment ) == 0 ) 
	TEMPLATE_LIST_END
	
	size_t count = field_arrays.size();
	size_t k = 0;
	for(size_t i=0;i<count;i++)
        {
		TEMPLATE_LIST_BEGIN
			field_arrays[soatl::FieldId<ids>()][i] = static_cast<typename soatl::FieldDescriptor<ids>::value_type>( k ) ,
			++k
		TEMPLATE_LIST_END
	}
	// read from and write to the area beyond size() and up to next chunk boundary
	for(size_t i=count;i<field_arrays.chunk_ceil();i++)
        {
		TEMPLATE_LIST_BEGIN
			k += static_cast<size_t>( field_arrays[soatl::FieldId<ids>()][i] ) ,
			field_arrays[soatl::FieldId<ids>()][i] = static_cast<typename soatl::FieldDescriptor<ids>::value_type>( k ) 
		TEMPLATE_LIST_END
	}

	// read back values in [0;size()[ to check values are still correct
	k=0;
	for(size_t i=0;i<count;i++)
        {
		bool value_ok = false;
		size_t findex = 0;
		TEMPLATE_LIST_BEGIN
			value_ok = field_arrays[soatl::FieldId<ids>()][i] == static_cast<typename soatl::FieldDescriptor<ids>::value_type>( k ) ,
			//std::cout<<"value["<<findex<<"]="<< (size_t)(field_arrays[fids][i] )<<std::endl ,
			assert(value_ok) ,
			++ findex ,
			++k
		TEMPLATE_LIST_END
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
	auto cell_arrays2 = soatl::make_packed_field_arrays( atype,rx,mid,ry,rz );

	std::cout<<"arrays1: alignment="<<cell_arrays1.alignment()<<", chunksize="<<cell_arrays1.chunksize() <<std::endl;  std::cout.flush();
	std::cout<<"arrays2: alignment="<<cell_arrays2.alignment()<<", chunksize="<<cell_arrays2.chunksize() <<std::endl;  std::cout.flush();

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
	soatl::apply_simd( compute_distance , 0, N, dist_ptr, rx_ptr, ry_ptr, rz_ptr, rx2_ptr, ry2_ptr, rz2_ptr );
	for(size_t j=0;j<N;j++)
	{
		if(rdist(rng)>0.999) std::cout<<"dist["<<j<<"]="<<dist_ptr[j]<<std::endl; std::cout.flush();
	}

	std::cout<<"apply compute_distance to arrays (double)"<<std::endl;  std::cout.flush();
	soatl::apply_simd( compute_distance_d , N, e_ptr, rx_ptr, ry_ptr, rz_ptr, rx2_ptr, ry2_ptr, rz2_ptr );
	for(size_t j=0;j<N;j++)
	{
		if(rdist(rng)>0.999) std::cout<<"e["<<j<<"]="<<e_ptr[j]<<std::endl; std::cout.flush();
	}


	std::cout<<"borrow array pointers"<<std::endl; std::cout.flush();
	// make a bunch of unmanaged pointers to fields rx, ry, rz, e and dist using size, alignment and chunksize from cell_arrays1
	auto zip = soatl::make_field_pointers( cell_arrays1, rx, ry, rz, e, dist );

	zip.set_pointer(rx , cell_arrays2[rx] );
	zip.set_pointer(ry , cell_arrays2[ry] );
	zip.set_pointer(rz , cell_arrays2[rz] );
	zip.set_pointer(dist , cell_arrays1[dist] );
	zip.set_pointer(e , cell_arrays1[e] );

	std::cout<<"copy arrays"<<std::endl;
	soatl::copy( zip , cell_arrays1 , 0, N, rx, ry, rz );

	const double ax = rx2_ptr[0];
	const double ay = ry2_ptr[0];
	const double az = rz2_ptr[0];

	std::cout<<"apply lambda (float)"<<std::endl;
	soatl::apply_simd( [ax,ay,az](float& d, double x, double y, double z)
		{
			x = ax - x;
			y = ay - y;
			z = az - z;
			d = /*std::sqrt*/ (x*x+y*y+z*z);
		}
		, zip, dist, rx, ry, rz );

	for(size_t j=0;j<N;j++)
	{
		if(rdist(rng)>0.999) std::cout<<"dist["<<j<<"]="<<zip[dist][j]<<std::endl; std::cout.flush();
	}

	std::cout<<"apply lambda (double)"<<std::endl;
	soatl::apply_simd( [ax,ay,az](double& d, double x, double y, double z)
		{
			x = ax - x;
			y = ay - y;
			z = az - z;
			d = std::sqrt(x*x+y*y+z*z);
		}
		, 16, N-16, zip, e,rx,ry,rz );
	for(size_t j=0;j<N;j++)
	{
		if(rdist(rng)>0.999) std::cout<<"e["<<j<<"]="<<zip[e][j]<<std::endl; std::cout.flush();
	}

	auto unalignedzip = soatl::make_field_pointers( cell_arrays1.size(), soatl::cst::unaligned, soatl::cst::no_chunk, e,rx,ry,rz );
	unalignedzip.set_pointers( cell_arrays1, e,rx,ry,rz );
	std::cout<<"apply_simd lambda (double) to minmally aligned pointers with no chunk"<<std::endl;
	soatl::apply_simd( [ax,ay,az](double& d, double x, double y, double z)
		{
			x = ax - x;
			y = ay - y;
			z = az - z;
			d = std::sqrt(x*x+y*y+z*z);
		}
		, soatl::SimdRequirements<double>::chunksize
		, N - soatl::SimdRequirements<double>::chunksize
		, unalignedzip, e,rx,ry,rz );
	for(size_t j=0;j<N;j++)
	{
		if(rdist(rng)>0.999) std::cout<<"e["<<j<<"]="<<zip[e][j]<<std::endl; std::cout.flush();
	}

	// ok, non-simd version does not require any alignment
	std::cout<<"apply lambda (double) to unaligned pointers"<<std::endl;
	soatl::apply( [ax,ay,az](double& d, double x, double y, double z)
		{
			x = ax - x;
			y = ay - y;
			z = az - z;
			d = std::sqrt(x*x+y*y+z*z);
		}
		, 1, N-1, unalignedzip, e,rx,ry,rz );
	for(size_t j=0;j<N;j++)
	{
		if(rdist(rng)>0.999) std::cout<<"e["<<j<<"]="<<zip[e][j]<<std::endl; std::cout.flush();
	}

	auto vect = soatl::make_static_packed_field_arrays( soatl::cst::count<128>(), rx, ry, rz, dist, e );
	std::cout<<"sizeof(vect)="<<sizeof(vect)<<", data_size="<<vect.data_size()<<", size="<<vect.size()<<std::endl; std::cout.flush();

	std::cout<<"check vect alignment and aliasing"<<std::endl; std::cout.flush();
	check_static_fields( vect );
	if( N >= vect.size() )
	{
		std::cout<<"copy to vect"<<std::endl; std::cout.flush();
		soatl::copy( vect, cell_arrays1, rx,ry,rz,dist,e );
		std::cout<<"compute on vect"<<std::endl; std::cout.flush();
		soatl::apply_simd( [ax,ay,az](double& d, double x, double y, double z)
			{
				x = ax - x;
				y = ay - y;
				z = az - z;
				d = /*std::sqrt*/ (x*x+y*y+z*z);
			}
			, vect, e, rx, ry, rz );
		
		std::cout<<"print vect"<<std::endl; std::cout.flush();
		for(size_t j=0;j<vect.size();j++)
		{
			if(rdist(rng)>0.99) std::cout<<"e["<<j<<"]="<<vect[e][j]<<std::endl; std::cout.flush();
		}
	}

	return 0;
}


