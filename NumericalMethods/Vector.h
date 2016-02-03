#pragma once

namespace numerical { namespace basic {
	template<class TDataType, unsigned TSize>
	class Vector final
	{
	public:
		using VectorType = Vector<TDataType, TSize>;

		static auto unitVector(unsigned const i) -> VectorType;

		Vector();

		explicit Vector(TDataType initialValue);

		Vector(VectorType const& vec);

		Vector(TDataType const& val0, TDataType const& val1);

		auto operator[](unsigned const i) -> TDataType&;

		auto operator[](unsigned const i) const -> TDataType const&;

		auto operator=(VectorType const& rhs) -> VectorType&;

		auto operator+(VectorType const& rhs) const -> VectorType;

		auto operator+=(VectorType const& rhs) -> VectorType&;

		auto operator-(VectorType const& rhs) const -> VectorType;

		auto operator*(TDataType const value) const -> VectorType;

		auto operator*=(TDataType const value) -> VectorType&;

	private:
		TDataType data[TSize];
	};

	using Vector2d = Vector<double, 2>;
}}