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

// WARNING: assumes that elements in arrays support 'operator = (const size_t&)' and 'operator == (const size_t&)'
template<typename FieldsT, typename... _ids>
static inline void check_field_arrays_aliasing( size_t N, FieldsT& field_arrays, const soatl::FieldId<_ids>& ... fids )
{
	std::uniform_int_distribution<> rndist(0,N*2);

#ifdef DEBUG
	std::cout<<"check_field_arrays_aliasing: align="<<field_arrays.alignment()<<", chunksize="<<field_arrays.chunksize()<<", using fields :"<<std::endl;
	TEMPLATE_LIST_BEGIN
		std::cout<<"\t"<< soatl::FieldDescriptor<_ids>::name() <<std::endl
	TEMPLATE_LIST_END
#endif
        for(size_t j=1;j<=N;j++)
        {
                field_arrays.resize(j);

		// test pointer alignment
		{
			size_t alignment = field_arrays.alignment();
			void* ptr;
			size_t addr;
			TEMPLATE_LIST_BEGIN
				ptr = field_arrays[fids] ,
				addr = reinterpret_cast<size_t>(ptr) ,
				assert( ( addr % alignment ) == 0 ) 
			TEMPLATE_LIST_END
		}

		// container 'c' must guarantee that acces beyond c.size() is valid up to the next chunk boundary
		// container only guarantees that read or write access beyond c.size() (up to the next chunk) but values in this area are undefined
#ifdef DEBUG
		std::cout<<"check_field_arrays_aliasing: a="<<field_arrays.alignment()
			 <<", c="<<field_arrays.chunksize()
			 <<", size="<<field_arrays.size()
			 <<", capacity="<<field_arrays.capacity()
			 <<", chunk boundary="<<field_arrays.chunk_ceil()<<std::endl;
#endif
		size_t k = 0;
		for(size_t i=0;i<j;i++)
                {
			TEMPLATE_LIST_BEGIN
				field_arrays[fids][i] = static_cast<typename soatl::FieldDescriptor<_ids>::value_type>( k ) ,
				++k
			TEMPLATE_LIST_END
		}
		// read from and write to the area beyond size() and up to next chunk boundary
		for(size_t i=j;i<field_arrays.chunk_ceil();i++)
                {
			TEMPLATE_LIST_BEGIN
				k += static_cast<size_t>( field_arrays[fids][i] ) ,
				field_arrays[fids][i] = static_cast<typename soatl::FieldDescriptor<_ids>::value_type>( k ) 
			TEMPLATE_LIST_END
		}

		// read back values in [0;size()[ to check values are still correct
		k=0;
		for(size_t i=0;i<j;i++)
                {
			bool value_ok = false;
			size_t findex = 0;
			TEMPLATE_LIST_BEGIN
				value_ok = field_arrays[fids][i] == static_cast<typename soatl::FieldDescriptor<_ids>::value_type>( k ) ,
				//std::cout<<"value["<<findex<<"]="<< (size_t)(field_arrays[fids][i] )<<std::endl ,
				assert(value_ok) ,
				++ findex ,
				++k
			TEMPLATE_LIST_END
		}
		
		// test that a resize keeps values
		size_t ns = rndist(rng);
		size_t cs = std::min( ns , j );
#ifdef DEBUG
		std::cout<<"resize from "<<j<<" to "<<ns<<", cs="<<cs<< std::endl;
#endif
		field_arrays.resize( ns );
		k=0;
		for(size_t i=0;i<cs;i++)
                {
			bool value_ok = false;
			size_t findex = 0;
			TEMPLATE_LIST_BEGIN
				value_ok = field_arrays[fids][i] == static_cast<typename soatl::FieldDescriptor<_ids>::value_type>( k ) ,
#ifdef DEBUG
				std::cout<<"value["<<_ids<<"]["<<i<<"] ="<< (size_t)( field_arrays[fids][i] )<<std::endl ,
#endif
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
			ptr = field_arrays[fids] ,
#ifdef DEBUG
			std::cout<<"ptr["<<findex<<"]="<<ptr<<std::endl ,
#endif
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
	std::cout<<"test_field_arrays_aliasing<"<<A<<","<<C<<">"<<std::endl;	

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

	test_packed_field_arrays_aliasing<2,1>();
	test_packed_field_arrays_aliasing<2,2>();
	test_packed_field_arrays_aliasing<2,3>();
	test_packed_field_arrays_aliasing<2,6>();
	test_packed_field_arrays_aliasing<2,16>();

	test_packed_field_arrays_aliasing<4,1>();
	test_packed_field_arrays_aliasing<4,2>();
	test_packed_field_arrays_aliasing<4,3>();
	test_packed_field_arrays_aliasing<4,6>();
	test_packed_field_arrays_aliasing<4,16>();

	test_packed_field_arrays_aliasing<8,1>();
	test_packed_field_arrays_aliasing<8,2>();
	test_packed_field_arrays_aliasing<8,3>();
	test_packed_field_arrays_aliasing<8,6>();
	test_packed_field_arrays_aliasing<8,16>();

	test_packed_field_arrays_aliasing<16,1>();
	test_packed_field_arrays_aliasing<16,2>();
	test_packed_field_arrays_aliasing<16,3>();
	test_packed_field_arrays_aliasing<16,6>();
	test_packed_field_arrays_aliasing<16,16>();

	test_packed_field_arrays_aliasing<32,1>();
	test_packed_field_arrays_aliasing<32,2>();
	test_packed_field_arrays_aliasing<32,3>();
	test_packed_field_arrays_aliasing<32,6>();
	test_packed_field_arrays_aliasing<32,16>();

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

	return 0;
}


