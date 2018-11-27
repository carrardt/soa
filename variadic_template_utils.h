#pragma once

struct __pass
{
    template<typename ...T> __pass(T...) {}
};

#define TEMPLATE_LIST_BEGIN 	__pass{(
#define TEMPLATE_LIST_END 	,1)...};

