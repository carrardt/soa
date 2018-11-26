#pragma once

#include <cstdlib> // for size_t

template<typename _T,int _Id>
struct FieldDataDescriptor
{
	using value_type = _T;
	static constexpr int FieldId = _Id;

	inline FieldDataDescriptor(const std::string& n) : m_name(n) {}
	inline const std::string& name() const { return m_name; }

	// permet de decrire la facon dont on veut acceder a un element
	inline const FieldDataDescriptor& read_only() const { return *this; }

private:
	const std::string m_name;
};


// helper to get tuple index from field Id
template<bool found, int _Index, int _Id, typename... FDS >
struct FieldIdx {};
template<int _Index, int _Id, typename... FDS> 
struct FieldIdx<true,_Index,_Id,FDS...> { static constexpr int index = _Index; };
template<int _Index, int _Id, typename FD, typename... FDS> 
struct FieldIdx<false,_Index,_Id,FD,FDS...> { static constexpr int index = FieldIdx<_Id==FD::FieldId,_Index+1,_Id,FDS...>::index; };
template<int _Index, int _Id>
struct FieldIdx<false,_Index,_Id> { static constexpr int index = -1; };
template<int _Id, typename FD, typename... FDS>
struct FindFieldIndex { static constexpr int index = FieldIdx<false,-1,_Id,FD,FDS...>::index; };

