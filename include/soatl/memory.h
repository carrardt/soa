#pragma once

namespace soatl
{

	struct DefaultAllocationStrategy
	{
		static inline size_t update_capacity(size_t s, size_t capacity, size_t chunksize)
		{
			if( s>capacity || ( (s+2*chunksize) <= capacity ) || s==0 )
			{
				return ( (s+chunksize-1)/chunksize ) * chunksize ;
			}
			else
			{
				return capacity;
			}
		}
	};

}


