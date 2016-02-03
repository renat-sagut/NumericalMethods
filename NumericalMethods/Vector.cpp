#include "Warnings.h"
#include "Vector.h"

#pragma warning(push)

#include <type_traits>
#include <assert.h>

#pragma warning(pop)

namespace numerical { namespace basic {

	template class Vector<double, 2>;

	template<class TDataType, unsigned TSize>
	auto Vector<TDataType, TSize>::unitVector(unsigned const i) -> VectorType
	{
		assert(i < TSize);
		VectorType result;
		result[i] = static_cast<TDataType>(1.0);

		return result;
	}

	template <class TDataType, unsigned TSize>
	Vector<TDataType, TSize>::Vector()
	{
		static_assert(TSize > 0, "Invalid TSize template parameter used");

		for (unsigned i = 0; i < TSize; ++i)
		{
			this->data[i] = static_cast<TDataType>(0.0);
		}
	}

	template<class TDataType, unsigned TSize>
	Vector<TDataType, TSize>::Vector(TDataType initialValue)
	{
		static_assert(TSize > 0, "Invalid TSize template parameter used");

		for (unsigned i = 0; i < TSize; ++i)
		{
			this->data[i] = initialValue;
		}
	}

	template<class TDataType, unsigned TSize>
	Vector<TDataType, TSize>::Vector(VectorType const& vec)
	{
		for (unsigned i = 0; i < TSize; ++i)
		{
			this->data[i] = vec.data[i];
		}
	}

	template<class TDataType, unsigned TSize>
	Vector<TDataType, TSize>::Vector(TDataType const& val0, TDataType const& val1)
	{
		static_assert(TSize == 2, "Invalid constructor of size 2 used");

		this->data[0] = val0;
		this->data[1] = val1;
	}

	template<class TDataType, unsigned TSize>
	auto Vector<TDataType, TSize>::operator[](unsigned const i) -> TDataType&
	{
		assert(i < TSize);
		return this->data[i];
	}

	template<class TDataType, unsigned TSize>
	auto Vector<TDataType, TSize>::operator[](unsigned const i) const -> TDataType const&
	{
		//return this->operator[i];
		assert(i < TSize);
		return this->data[i];

	}

	template<class TDataType, unsigned TSize>
	auto Vector<TDataType, TSize>::operator=(VectorType const& rhs) -> VectorType&
	{
		for (unsigned i = 0; i < TSize; ++i)
		{
			this->data[i] = rhs.data[i];
		}

		return *this;
	}

	template<class TDataType, unsigned TSize>
	auto Vector<TDataType, TSize>::operator+(VectorType const& rhs) const -> VectorType
	{
		VectorType result(*this);
		for (unsigned i = 0; i < TSize; ++i)
		{
			result.data[i] += rhs.data[i];
		}

		return result;
	}

	template<class TDataType, unsigned TSize>
	auto Vector<TDataType, TSize>::operator+=(VectorType const& rhs) -> VectorType&
	{
		for (unsigned i = 0; i < TSize; ++i)
		{
			this->data[i] += rhs.data[i];
		}

		return *this;
	}

	template<class TDataType, unsigned TSize>
	auto Vector<TDataType, TSize>::operator-(VectorType const& rhs) const -> VectorType
	{
		VectorType result(*this);
		for (unsigned i = 0; i < TSize; ++i)
		{
			result.data[i] -= rhs.data[i];
		}

		return result;
	}

	template<class TDataType, unsigned TSize>
	auto Vector<TDataType, TSize>::operator*(TDataType const value) const -> VectorType
	{
		VectorType result(*this);
		for (unsigned i = 0; i < TSize; ++i)
		{
			result.data[i] *= value;
		}

		return result;
	}

	template<class TDataType, unsigned TSize>
	auto Vector<TDataType, TSize>::operator*=(TDataType const value) -> VectorType&
	{
		for (unsigned i = 0; i < TSize; ++i)
		{
			this->data[i] *= value;
		}

		return *this;
	}
}}

#pragma warning(disable : WARNING_UNREFERENCED_INLINE_FUNCTION_HAS_BEEN_REMOVED)