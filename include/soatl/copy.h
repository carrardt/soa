#pragma once

#include "soatl/field_descriptor.h"
#include <cstdlib> // for size_t
#include <cstring>

namespace soatl
{
	template<typename DstArrays, typename SrcArrays, size_t... _ids> struct FieldArraysCopyHelper;
	template<typename DstArrays, typename SrcArrays, size_t id, size_t... _ids>
	struct FieldArraysCopyHelper<DstArrays,SrcArrays, id, _ids...>
	{
		static inline void copy( DstArrays& dst, const SrcArrays& src, size_t start, size_t count )
		{
			using ValueType = typename FieldDescriptor<id>::value_type;
			//std::cout<<"copy : field='"<< FieldDescriptor<id>::name() <<"', range=["<<start<<";"<<start+count<<"[, d="<<(void*)d<<", s="<<(void*)s<< std::endl;

			// copy, version 1 : always work
			std::memcpy( dst[FieldId<id>()]+start, src[ FieldId<id>() ]+start, sizeof(ValueType)*count );			

			// copy, version 2 : crashes if compiler generates aligned move instructions and arrays alignment is not sufficient. It happens with gcc 5.4 (-O3)
			//auto d = dst[ FieldId<id>() ];
			//auto s = src[ FieldId<id>() ];
			//for(size_t i=start; i<(start+count); i++) { d[i] = s[i]; }
			
			// copy, version 3 : always work, might be less efficient than memcpy
			//uint8_t* d = reinterpret_cast<uint8_t*>( dst[ FieldId<id>() ] );
			//const uint8_t* s = reinterpret_cast<uint8_t*>( src[ FieldId<id>() ] );
			//for(size_t i=start*sizeof(ValueType); i<(start+count)*sizeof(ValueType); i++) { d[i] = s[i]; }

			FieldArraysCopyHelper<DstArrays,SrcArrays,_ids...>::copy(dst,src,start,count);
		}
	};

	template<typename DstArrays, typename SrcArrays>
	struct FieldArraysCopyHelper<DstArrays,SrcArrays>
	{
		static inline void copy(DstArrays&,const SrcArrays&,size_t,size_t) {}
	};

	template<typename DstArrays, typename SrcArrays, size_t... _ids>
	static inline void copy( DstArrays& dst, const SrcArrays& src, size_t start, size_t count, const FieldId<_ids>&... )
	{
		assert( (start+count) <= dst.size() );
		assert( (start+count) <= src.size() );
		FieldArraysCopyHelper<DstArrays,SrcArrays,_ids...>::copy(dst,src,start,count);
	}

}

