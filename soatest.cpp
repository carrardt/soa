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

enum ParticleField
{
	PARTICLE_RX,
	PARTICLE_RY,
	PARTICLE_RZ,
	PARTICLE_E,
	PARTICLE_ATOM_TYPE,
	PARITCLE_MOLECULE_ID,

	PARTICLE_PAIR_DIST,

	PARTICLE_TMP1,
	PARTICLE_TMP2,
	PARTICLE_TMP3,

	PARTICLE_FIELD_COUNT
};

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

template<typename OperatorT, size_t A, size_t C, typename... FieldDescriptors >
static inline void apply_to_field_arrays( OperatorT f, size_t first, size_t N, const soatl::FieldArrays<A,C,FieldDescriptors...>& field_arrays)
{
	size_t last = first+N;
	#pragma omp simd
	for(size_t i=first;i<last;i++)
	{
		f( field_arrays.get(FieldDescriptors())[i] ... );
	}
}

// WARNING: assumes that elements in arrays support 'operator = (const size_t&)' and 'operator == (const size_t&)'
template<typename FieldsT, typename... FieldDescriptors>
static inline void check_field_arrays_aliasing( size_t N, FieldsT& field_arrays, const FieldDescriptors & ... fd )
{
	std::cout<<"check_field_arrays_aliasing: align="<<field_arrays.alignment()<<", chunksize="<<field_arrays.chunksize()<<std::endl;
        for(size_t j=1;j<=N;j++)
        {
                field_arrays.resize(j);
		//std::cout<<"check_field_arrays_aliasing: size="<<field_arrays.size()<<", capacity="<<field_arrays.capacity()<<std::endl;
		size_t k=0;
		for(size_t i=0;i<j;i++)
                {
			TEMPLATE_LIST_BEGIN
				field_arrays.get( FieldDescriptors() )[i] = static_cast<typename FieldDescriptors::value_type>( k ) ,
				++k
			TEMPLATE_LIST_END
		}
		k=0;
		for(size_t i=0;i<j;i++)
                {
			bool value_ok = false;
			TEMPLATE_LIST_BEGIN
				value_ok = field_arrays.get( FieldDescriptors() )[i] == static_cast<typename FieldDescriptors::value_type>( k ) ,
				assert(value_ok) ,
				++k
			TEMPLATE_LIST_END
		}
		//std::cout<<k<<" values tested"<<std::endl;
	}

	field_arrays.resize(0);
	void* ptr;
	TEMPLATE_LIST_BEGIN
		ptr = field_arrays.get( FieldDescriptors() ) ,
		assert(ptr==nullptr) 
	TEMPLATE_LIST_END
}

template<size_t A, size_t C>
static inline void test_packed_field_arrays_aliasing()
{
	std::cout<<"test_pack_field_arrays_aliasing<"<<A<<","<<C<<">"<<std::endl;

	FieldDataDescriptor<double,PARTICLE_RX> rx("rx");
	FieldDataDescriptor<double,PARTICLE_RY> ry("ry");
	FieldDataDescriptor<double,PARTICLE_RZ> rz("rz");
	FieldDataDescriptor<float,PARTICLE_E> e("e");
	FieldDataDescriptor<unsigned char,PARTICLE_ATOM_TYPE> atype("type");
	FieldDataDescriptor<int32_t,PARITCLE_MOLECULE_ID> mid("molecule");
	FieldDataDescriptor<float,PARTICLE_PAIR_DIST> dist("distance");
	FieldDataDescriptor<int16_t,PARTICLE_TMP1> tmp1("temporary 1");
	FieldDataDescriptor<int8_t,PARTICLE_TMP2> tmp2("temporary 2");

	{
		auto field_arrays = soatl::make_packed_field_arrays( soatl::cst::align<A>(), soatl::cst::chunk<C>(), atype,rx,mid,ry,dist,rz,tmp1 );
		check_field_arrays_aliasing(1063, field_arrays , atype,rx,mid,ry,dist,rz,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, rx,ry,rz,atype,mid,dist,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, atype,mid,dist,tmp1,rx,ry,rz );
	}

	{
		auto field_arrays = soatl::make_packed_field_arrays( soatl::cst::align<A>(), soatl::cst::chunk<C>(), atype,mid,dist,tmp2,tmp1,rx,ry,rz );
		check_field_arrays_aliasing(1063,field_arrays, atype,rx,mid,ry,dist,rz,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, rx,ry,rz,atype,mid,tmp2,dist,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, atype,mid,dist,tmp1,rx,ry,tmp2,rz );
	}

	{
		auto field_arrays = soatl::make_packed_field_arrays( soatl::cst::align<A>(), soatl::cst::chunk<C>(), rx,ry,rz,atype,mid,dist,tmp2,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, atype,tmp2,rx,mid,ry,dist,rz,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, rx,ry,tmp2,rz,atype,mid,dist,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, atype,mid,dist,tmp1,rx,ry,tmp2,rz );
	}

}

template<size_t A, size_t C>
static inline void test_field_arrays_aliasing()
{
	std::cout<<"test_field_arrays_aliasing<"<<C<<">"<<std::endl;	

	FieldDataDescriptor<double,PARTICLE_RX> rx("rx");
	FieldDataDescriptor<double,PARTICLE_RY> ry("ry");
	FieldDataDescriptor<double,PARTICLE_RZ> rz("rz");
	FieldDataDescriptor<float,PARTICLE_E> e("e");
	FieldDataDescriptor<unsigned char,PARTICLE_ATOM_TYPE> atype("type");
	FieldDataDescriptor<int32_t,PARITCLE_MOLECULE_ID> mid("molecule");
	FieldDataDescriptor<float,PARTICLE_PAIR_DIST> dist("distance");
	FieldDataDescriptor<int16_t,PARTICLE_TMP1> tmp1("temporary 1");
	FieldDataDescriptor<int8_t,PARTICLE_TMP2> tmp2("temporary 2");

	{
		auto field_arrays = soatl::make_field_arrays( soatl::cst::align<A>(), soatl::cst::chunk<C>(), atype,rx,mid,ry,dist,rz,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, atype,rx,mid,ry,dist,rz,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, rx,ry,rz,atype,mid,dist,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, atype,mid,dist,tmp1,rx,ry,rz );
	}

	{
		auto field_arrays = soatl::make_field_arrays( soatl::cst::align<A>(), soatl::cst::chunk<C>(), atype,mid,dist,tmp2,tmp1,rx,ry,rz );
		check_field_arrays_aliasing(1063,field_arrays, atype,rx,mid,ry,dist,rz,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, rx,ry,rz,atype,mid,tmp2,dist,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, atype,mid,dist,tmp1,rx,ry,tmp2,rz );
	}

	{
		auto field_arrays = soatl::make_field_arrays( soatl::cst::align<A>(), soatl::cst::chunk<C>(), rx,ry,rz,atype,mid,dist,tmp2,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, atype,tmp2,rx,mid,ry,dist,rz,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, rx,ry,tmp2,rz,atype,mid,dist,tmp1 );
		check_field_arrays_aliasing(1063,field_arrays, atype,mid,dist,tmp1,rx,ry,tmp2,rz );
	}
}


int main(int argc, char* argv[])
{
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

	int seed = 0;
	size_t N = 10000;

	if(argc>=2) { N=atoi(argv[1]); }
	if(argc>=3) { seed=atoi(argv[2]); }

	FieldDataDescriptor<double,PARTICLE_RX> rx("rx");
	FieldDataDescriptor<double,PARTICLE_RY> ry("ry");
	FieldDataDescriptor<double,PARTICLE_RZ> rz("rz");
	FieldDataDescriptor<float,PARTICLE_E> e("e");
	FieldDataDescriptor<unsigned char,PARTICLE_ATOM_TYPE> atype("type");
	FieldDataDescriptor<int32_t,PARITCLE_MOLECULE_ID> mid("molecule");
	FieldDataDescriptor<float,PARTICLE_PAIR_DIST> dist("distance");
	FieldDataDescriptor<int16_t,PARTICLE_TMP1> tmp1("temporary 1");
	FieldDataDescriptor<int8_t,PARTICLE_TMP2> tmp2("temporary 2");

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

	auto zip = soatl::make_field_pointers( N, soatl::cst::align<64>(), soatl::cst::chunk<8>(), atype,rx,mid,ry,rz );

}


