#pragma once

namespace soatl
{

	struct ChunkIncrementalAllocationStrategy
	{
		static inline size_t update_capacity(size_t s, size_t capacity, size_t chunksize)
		{
			size_t newCapacity = 0;
			if( s>capacity || ( (s+2*chunksize) <= capacity ) || s==0 )
			{
				newCapacity = ( (s+chunksize-1)/chunksize ) * chunksize ;
			}
			else
			{
				newCapacity = capacity;
			}
			assert( newCapacity >= s );
			assert( (newCapacity % chunksize) == 0 );
			return newCapacity;
		}
	};

	using DefaultAllocationStrategy = ChunkIncrementalAllocationStrategy;

}


