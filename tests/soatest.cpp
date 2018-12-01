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

// WARNING: assumes that elements in arrays support 'operator = (const size_t&)' and 'operator == (const size_t&)'
template<typename FieldsT, size_t... _ids>
static inline void check_field_arrays_aliasing( size_t N, FieldsT& field_arrays, const soatl::FieldId<_ids>& ... fids )
{
	std::uniform_int_distribution<> rndist(0,N*2);

	/*
	std::cout<<"check_field_arrays_aliasing: align="<<field_arrays.alignment()<<", chunksize="<<field_arrays.chunksize()<<", using fields :"<<std::endl;
	TEMPLATE_LIST_BEGIN
		std::cout<<"\t"<< soatl::FieldDescriptor<_ids>::name() <<std::endl
	TEMPLATE_LIST_END
	*/

        for(size_t j=1;j<=N;j++)
        {
                field_arrays.resize(j);

		// test pointer alignment
		{
			size_t alignment = field_arrays.alignment();
			void* ptr;
			size_t addr;
			TEMPLATE_LIST_BEGIN
				ptr = field_arrays.get( fids ) ,
				addr = reinterpret_cast<size_t>(ptr) ,
				assert( ( addr % alignment ) == 0 ) 
			TEMPLATE_LIST_END
		}

		// container 'c' must guarantee that acces beyond c.size() is valid up to the next chunk boundary
		// container only guarantees that read or write access beyond c.size() (up to the next chunk) but values in this area are undefined
		//std::cout<<"check_field_arrays_aliasing: size="<<field_arrays.size()<<", capacity="<<field_arrays.capacity()<<", chunk boundary="<<field_arrays.chunk_ceil()<<std::endl;
		size_t k = 0;
		for(size_t i=0;i<j;i++)
                {
			TEMPLATE_LIST_BEGIN
				field_arrays.get( fids )[i] = static_cast<typename soatl::FieldDescriptor<_ids>::value_type>( k ) ,
				++k
			TEMPLATE_LIST_END
		}
		// read from and write to the area beyond size() and up to next chunk boundary
		for(size_t i=j;i<field_arrays.chunk_ceil();i++)
                {
			TEMPLATE_LIST_BEGIN
				k += static_cast<size_t>( field_arrays.get( fids )[i] ) ,
				field_arrays.get( fids )[i] = static_cast<typename soatl::FieldDescriptor<_ids>::value_type>( k ) 
			TEMPLATE_LIST_END
		}

		// read back values in [0;size()[ to check values are still correct
		k=0;
		for(size_t i=0;i<j;i++)
                {
			bool value_ok = false;
			size_t findex = 0;
			TEMPLATE_LIST_BEGIN
				value_ok = field_arrays.get( fids )[i] == static_cast<typename soatl::FieldDescriptor<_ids>::value_type>( k ) ,
				//std::cout<<"value["<<findex<<"]="<< (size_t)(field_arrays.get( FieldDescriptors() )[i] )<<std::endl ,
				assert(value_ok) ,
				++ findex ,
				++k
			TEMPLATE_LIST_END
		}
		
		// test that a resize keeps values
		size_t ns = rndist(rng);
		field_arrays.resize( ns );
		size_t cs = std::min( ns , j );
		//std::cout<<"resize from "<<j<<" to "<<ns<<", cs="<<cs<< std::endl;
		k=0;
		for(size_t i=0;i<cs;i++)
                {
			bool value_ok = false;
			size_t findex = 0;
			TEMPLATE_LIST_BEGIN
				value_ok = field_arrays.get( fids )[i] == static_cast<typename soatl::FieldDescriptor<_ids>::value_type>( k ) ,
				//std::cout<<"value["<<findex<<"]="<< (size_t)( field_arrays.get( FieldDescriptors() )[i] )<<std::endl ,
				assert(value_ok) ,
				++ findex ,
				++k
			TEMPLATE_LIST_END
		}

	}

	// test that all pointers are 0 (nullptr) when capacity is 0
	{
		field_arrays.resize(0); // assumes that capacity is adjusted to 0 when container is resized to 0
		void* ptr;
		size_t findex = 0;
		TEMPLATE_LIST_BEGIN
			ptr = field_arrays.get( fids ) ,
			//std::cout<<"ptr["<<findex<<"]="<<ptr<<std::endl ,
			++findex , 
			assert(ptr==nullptr) 
		TEMPLATE_LIST_END
	}
}

template<size_t A, size_t C>
static inline void test_packed_field_arrays_aliasing()
{
	std::cout<<"test_pack_field_arrays_aliasing<"<<A<<","<<C<<">"<<std::endl;

	auto rx = particle_rx;
	auto ry = particle_ry;
	auto rz = particle_rz;
	auto e = particle_e;
	auto atype = particle_atype;
	auto mid = particle_mid;
	auto tmp1 = particle_tmp1;
	auto tmp2 = particle_tmp2;
	auto dist = particle_dist;

	{
		auto field_arrays = soatl::make_packed_field_arrays( soatl::cst::align<A>(), soatl::cst::chunk<C>(), atype,rx,mid,ry,dist,rz,tmp1 );
		assert( field_arrays.alignment()==A && field_arrays.chunksize()==C );
		check_field_arrays_aliasing(1063, field_arrays , atype,rx,mid,ry,dist,rz,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, rx,ry,rz,atype,mid,dist,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, atype,mid,dist,tmp1,rx,ry,rz );
	}

	{
		auto field_arrays = soatl::make_packed_field_arrays( soatl::cst::align<A>(), soatl::cst::chunk<C>(), atype,mid,dist,tmp2,tmp1,rx,ry,rz );
		assert( field_arrays.alignment()==A && field_arrays.chunksize()==C );
		check_field_arrays_aliasing(1063,field_arrays, atype,rx,mid,ry,dist,rz,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, rx,ry,rz,atype,mid,tmp2,dist,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, atype,mid,dist,tmp1,rx,ry,tmp2,rz );
	}

	{
		auto field_arrays = soatl::make_packed_field_arrays( soatl::cst::align<A>(), soatl::cst::chunk<C>(), rx,ry,rz,atype,mid,dist,tmp2,tmp1 );
		assert( field_arrays.alignment()==A && field_arrays.chunksize()==C );
		check_field_arrays_aliasing(1063,field_arrays, atype,tmp2,rx,mid,ry,dist,rz,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, rx,ry,tmp2,rz,atype,mid,dist,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, atype,mid,dist,tmp1,rx,ry,tmp2,rz );
	}

}

template<size_t A, size_t C>
static inline void test_field_arrays_aliasing()
{
	std::cout<<"test_field_arrays_aliasing<"<<C<<">"<<std::endl;	

	auto rx = particle_rx;
	auto ry = particle_ry;
	auto rz = particle_rz;
	auto e = particle_e;
	auto atype = particle_atype;
	auto mid = particle_mid;
	auto tmp1 = particle_tmp1;
	auto tmp2 = particle_tmp2;
	auto dist = particle_dist;

	{
		auto field_arrays = soatl::make_field_arrays( soatl::cst::align<A>(), soatl::cst::chunk<C>(), atype,rx,mid,ry,dist,rz,tmp1 );
		assert( field_arrays.alignment()==A && field_arrays.chunksize()==C );
		check_field_arrays_aliasing(1063,field_arrays, atype,rx,mid,ry,dist,rz,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, rx,ry,rz,atype,mid,dist,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, atype,mid,dist,tmp1,rx,ry,rz );
	}

	{
		auto field_arrays = soatl::make_field_arrays( soatl::cst::align<A>(), soatl::cst::chunk<C>(), atype,mid,dist,tmp2,tmp1,rx,ry,rz );
		assert( field_arrays.alignment()==A && field_arrays.chunksize()==C );
		check_field_arrays_aliasing(1063,field_arrays, atype,rx,mid,ry,dist,rz,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, rx,ry,rz,atype,mid,tmp2,dist,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, atype,mid,dist,tmp1,rx,ry,tmp2,rz );
	}

	{
		auto field_arrays = soatl::make_field_arrays( soatl::cst::align<A>(), soatl::cst::chunk<C>(), rx,ry,rz,atype,mid,dist,tmp2,tmp1 );
		assert( field_arrays.alignment()==A && field_arrays.chunksize()==C );
		check_field_arrays_aliasing(1063,field_arrays, atype,tmp2,rx,mid,ry,dist,rz,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, rx,ry,tmp2,rz,atype,mid,dist,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, atype,mid,dist,tmp1,rx,ry,tmp2,rz );
	}
}


int main(int argc, char* argv[])
{

	int seed = 0;
	size_t N = 10000;

	if(argc>=2) { N=atoi(argv[1]); }
	if(argc>=3) { seed=atoi(argv[2]); }

	rng.seed( seed );

	test_packed_field_arrays_aliasing<1,1>();
	test_packed_field_arrays_aliasing<1,2>();
	test_packed_field_arrays_aliasing<1,3>();
	test_packed_field_arrays_aliasing<1,6>();
	test_packed_field_arrays_aliasing<1,16>();

	test_packed_field_arrays_aliasing<16,1>();
	test_packed_field_arrays_aliasing<16,2>();
	test_packed_field_arrays_aliasing<16,3>();
	test_packed_field_arrays_aliasing<16,6>();
	test_packed_field_arrays_aliasing<16,16>();

	test_packed_field_arrays_aliasing<64,1>();
	test_packed_field_arrays_aliasing<64,2>();
	test_packed_field_arrays_aliasing<64,3>();
	test_packed_field_arrays_aliasing<64,6>();
	test_packed_field_arrays_aliasing<64,16>();

	test_field_arrays_aliasing<1,1>();
	test_field_arrays_aliasing<1,2>();
	test_field_arrays_aliasing<1,3>();
	test_field_arrays_aliasing<1,6>();
	test_field_arrays_aliasing<1,16>();

	test_field_arrays_aliasing<8,1>();
	test_field_arrays_aliasing<8,2>();
	test_field_arrays_aliasing<8,3>();
	test_field_arrays_aliasing<8,6>();
	test_field_arrays_aliasing<8,16>();

	test_field_arrays_aliasing<64,1>();
	test_field_arrays_aliasing<64,2>();
	test_field_arrays_aliasing<64,3>();
	test_field_arrays_aliasing<64,6>();
	test_field_arrays_aliasing<64,16>();

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

	//check_field_arrays_aliasing(N,cell_arrays2 , atype,rx,mid,ry,rz );

	// rebind operator ?
	// zip operator ?
	// embed several field_arrays ?

	// zip arrays
	cell_arrays1.resize(N);
	cell_arrays2.resize(N);

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


	auto zip = soatl::make_field_pointers( N, soatl::cst::align<64>(), soatl::cst::chunk<8>(), rx, ry, rz, dist );

	zip.get(rx) = cell_arrays2.get(rx);
	zip.get(ry) = cell_arrays2.get(ry);
	zip.get(rz) = cell_arrays2.get(rz);
	zip.get(dist) = cell_arrays1.get(dist);

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
			 zip.get(dist), zip.get(rx), zip.get(ry), zip.get(rz)
			 );
}


