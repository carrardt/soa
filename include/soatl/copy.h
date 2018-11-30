#pragma once

#include <cstdlib> // for size_t

namespace soatl
{
	template<typename DstArrays, typename SrcArrays, typename... FDS> struct FieldArraysCopyHelper;
	template<typename DstArrays, typename SrcArrays, typename FD, typename... FDS>
	struct FieldArraysCopyHelper<DstArrays,SrcArrays, FD, FDS...>
	{
		static inline void copy( DstArrays& dst, const SrcArrays& src, size_t start, size_t count )
		{
			auto d = dst.get( FD() );
			auto s = src.get( FD() );
			//std::cout<<"copy field #"<<FD::FieldId<<" from "<<start<<" to "<<start+count<<", d="<<(void*)d<<", s="<<(void*)s<< std::endl;
			for(size_t i=start; i<(start+count); i++)
			{
				d[i] = s[i];
				//std::cout<<"d["<<i<<"]="<<(size_t)(d[i])<<", s["<<i<<"]="<<(size_t)(s[i]) <<std::endl;
			}
			FieldArraysCopyHelper<DstArrays,SrcArrays,FDS...>::copy(dst,src,start,count);
		}
	};

	template<typename DstArrays, typename SrcArrays>
	struct FieldArraysCopyHelper<DstArrays,SrcArrays>
	{
		static inline void copy(DstArrays&,const SrcArrays&,size_t,size_t) {}
	};

	template<typename DstArrays, typename SrcArrays, typename... FDS>
	static inline void copy( DstArrays& dst, const SrcArrays& src, size_t start, size_t count, const FDS&... )
	{
		assert( (start+count) <= dst.size() );
		assert( (start+count) <= src.size() );
		FieldArraysCopyHelper<DstArrays,SrcArrays,FDS...>::copy(dst,src,start,count);
	}

}

