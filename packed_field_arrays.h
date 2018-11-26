#pragma once

#include <cstdint>

// helper to iterate over a tuple
template<int N> struct _PackedFieldArraysTupleHelper;
template<int N> struct _PackedFieldArraysTupleHelper
{
};

template<> struct _PackedFieldArraysTupleHelper<0>
{
};

template<typename... FieldDescriptors>
struct PackedFieldArrays
{
	using ArrayTuple = std::tuple< typename FieldDescriptors::pointer_type ... > ;
	using TupleHelper = _FieldArraysTupleHelper< sizeof...(FieldDescriptors) >;

	template<typename _T,int _Id>
	inline typename FieldDataDescriptor<_T,_Id>::pointer_type get(const FieldDataDescriptor<_T,_Id>&)
	{
		static constexpr int index = FindFieldIndex<_Id,FieldDescriptors...>::index;
		//std::cout<<f.name()<<" found at index "<<index<<std::endl;
		return std::get<index>(m_field_arrays);
	}

	inline void allocate(size_t s)
	{
		//TupleHelper::allocate( m_field_arrays , s );
	}

private:
	
	uint8_t* m_storage_ptr = nullptr;
};


template<typename... FieldDescriptors>
inline
PackedFieldArrays<FieldDescriptors...> make_packed_field_arrays(const FieldDescriptors&... fdt)
{
	return PackedFieldArrays<FieldDescriptors...>();
}

