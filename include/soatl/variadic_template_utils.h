#pragma once


namespace soatl {

struct __pass
{
    template<typename ...T> __pass(T...) {}
};

}

#define TEMPLATE_LIST_BEGIN 	soatl::__pass{(
#define TEMPLATE_LIST_END 	,1)...};

