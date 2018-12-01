#pragma once

#include <soatl/field_descriptor.h>
#include <cstdint>

namespace private_declarations {
	enum FieldIdEnum {

#define SOATL_DECLARE_FIELD(type,name,desc) name,

#include "data_fields_db.h"

#undef SOATL_DECLARE_FIELD

	FIELD_COUNT
	};
}

#define SOATL_DECLARE_FIELD(__type,__name,__desc) \
static soatl::FieldId<private_declarations::__name> __name; \
namespace soatl { \
template<> struct FieldDescriptor<private_declarations::__name> { \
	using value_type = __type; \
	static constexpr size_t Id = private_declarations::__name; \
	static const char* name() { return __desc ; } \
}; }

#include "data_fields_db.h"

#undef SOATL_DECLARE_FIELD

