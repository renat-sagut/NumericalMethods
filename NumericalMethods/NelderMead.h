#pragma once

#include "Vector.h"
using numerical::basic::Vector;
#include "Warnings.h"

#pragma warning(push)
#pragma warning(disable : WARNING_BYTES_PADDING_ADDED)

namespace numerical { namespace optimization {
	template<unsigned TSize>
	class NelderMead final
	{
	public:
		using PointType = Vector<double, TSize>;

		using FunctionType = auto (*)(PointType const&) -> double;

		NelderMead() = delete;

		NelderMead(double const functionTolerance);

		NelderMead(NelderMead const&) = delete;

		auto minimize(PointType const& initialPoint, double const initialDelta, FunctionType const targetFunction) -> PointType;

		auto operator=(NelderMead const&) -> NelderMead& = delete;

	private:
		class Simplex final
		{
		public:
			class PointValuePair final
			{
			public:
				PointValuePair();

				PointValuePair(PointType const& point, double const value);

				auto operator<(PointValuePair const& rhs) const -> bool;

				auto operator<=(PointValuePair const& rhs) const -> bool;

				auto operator>(PointValuePair const& rhs) const -> bool;

				auto operator>=(PointValuePair const& rhs) const -> bool;

				auto getValue() const -> double;

				auto getPoint() const -> PointType const&;

			private:
				PointType point;
				double value;
			};

			enum class TransformType : unsigned
			{
				REFLECTION,
				EXPANSION,
				CONTRACTION,
				MULTIPLE_CONTRACTION
			};

			enum class TransformResult : unsigned
			{
				EXCELLENT,
				GOOD,
				BAD,
				ANY
			};

			Simplex() = delete;

			Simplex(PointType const& initialPoint, double const initialDelta, FunctionType const targetFunction);

			Simplex(Simplex const&) = delete;

			auto operator=(Simplex const&) -> Simplex& = delete;

			//auto operator[](unsigned const i) -> PointType&;

			auto getRelativeTolerance() -> double;

			auto getPointWithLowestValue() const -> PointType const&;

			auto transform(TransformType transformType) -> TransformResult;

			auto extrapolateByFactor(unsigned const targetIndex, double const factor, PointType& calcPoint) const -> void;

			auto contractAroundLowestPoint(double const factor) -> void;

			auto updateIndices() -> void;

		private:
			auto evaluateTargetFunction(PointType const& point) -> double;

			auto setPointValuePair(unsigned const i, PointType const& point) -> void;

			auto setPointValuePair(unsigned const i, PointType const& point, double const value) -> void;

			auto getPointValuePair(unsigned const i) const -> PointValuePair const&;

			unsigned const size = TSize + 1;
			PointValuePair pointValuePairs[TSize + 1];
			FunctionType const targetFunction;
			unsigned indexHighest;
			unsigned indexNextHighest;
			unsigned indexLowest;
			unsigned funtionCallCount = 0;
			unsigned const maxFunctionCallCount = 1000;
		};

		double const functionTolerance;
	};

	using NelderMead2d = NelderMead<2>;

}}

#pragma warning(pop)