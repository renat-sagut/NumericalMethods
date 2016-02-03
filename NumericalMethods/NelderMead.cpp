#include "NelderMead.h"

#pragma warning(push)

#include <type_traits>
#include <assert.h>
#include <cmath>
#include <exception>

#pragma warning(pop)

#pragma warning(disable : WARNING_FUNCTION_NOT_INLINED)

namespace numerical { namespace optimization {

	template class NelderMead<2>;

	template<unsigned TSize>
	NelderMead<TSize>::Simplex::PointValuePair::PointValuePair() = default;

	template<unsigned TSize>
	NelderMead<TSize>::Simplex::PointValuePair::PointValuePair(PointType const& point, double const value)
		: PointValuePair::point{point}, PointValuePair::value{value}
	{

	}

	template<unsigned TSize>
	auto NelderMead<TSize>::Simplex::PointValuePair::operator<(PointValuePair const& rhs) const -> bool
	{
		return this->value < rhs.value;
	}

	template<unsigned TSize>
	auto NelderMead<TSize>::Simplex::PointValuePair::operator<=(PointValuePair const& rhs) const -> bool
	{
		return this->value <= rhs.value;
	}

	template<unsigned TSize>
	auto NelderMead<TSize>::Simplex::PointValuePair::operator>(PointValuePair const& rhs) const -> bool
	{
		return this->value > rhs.value;
	}

	template<unsigned TSize>
	auto NelderMead<TSize>::Simplex::PointValuePair::operator>=(PointValuePair const& rhs) const -> bool
	{
		return this->value >= rhs.value;
	}

	template<unsigned TSize>
	auto NelderMead<TSize>::Simplex::PointValuePair::getValue() const -> double
	{
		return this->value;
	}

	template<unsigned TSize>
	auto NelderMead<TSize>::Simplex::PointValuePair::getPoint() const -> PointType const&
	{
		return this->point;
	}

	template<unsigned TSize>
	NelderMead<TSize>::Simplex::Simplex(PointType const& initialPoint, double const initialDelta, FunctionType const targetFunction)
		: Simplex::targetFunction{targetFunction}
	{
		setPointValuePair(0, initialPoint);

		for (unsigned i = 1; i < this->size; ++i)
		{
			setPointValuePair(i, initialPoint + PointType::unitVector(i - 1) * initialDelta);
		}
	}

	/*template<unsigned TSize>
	auto NelderMead<TSize>::Simplex::operator[](unsigned const i) -> PointType&
	{
		assert(i < this->size);
		return this->points[i];
	}*/

	template<unsigned TSize>
	auto NelderMead<TSize>::Simplex::contractAroundLowestPoint(double const factor) -> void
	{
		auto const& lowestPoint = getPointValuePair(this->indexLowest).getPoint();

		for (unsigned i = 0; i < this->size; ++i)
		{
			if (i == this->indexLowest)
			{
				continue;
			}

			auto contractedPoint = (getPointValuePair(i).getPoint() + lowestPoint) * factor;
			
			setPointValuePair(i, contractedPoint);
		}
	}

	template<unsigned TSize>
	auto NelderMead<TSize>::Simplex::extrapolateByFactor(unsigned const targetIndex, double const factor, PointType& calcPoint) const -> void
	{
		PointType middlePoint(0.0);
		for (unsigned i = 0; i < this->size; ++i)
		{
			if (i == targetIndex)
			{
				continue;
			}

			middlePoint += getPointValuePair(i).getPoint();
		}

		middlePoint *= 1.0 / TSize;

		auto targetPoint = getPointValuePair(targetIndex).getPoint();
		calcPoint = middlePoint + (targetPoint - middlePoint) * factor;
	}

	template<unsigned TSize>
	auto NelderMead<TSize>::Simplex::transform(TransformType transformType) -> TransformResult
	{
		double const reflectionFactor = -1.0;
		double const expansionFactor = 2.0;
		double const contractionFactor = 0.5;

		auto result = TransformResult::BAD;

		switch (transformType)
		{
		case TransformType::REFLECTION:
		{
			PointType reflectedPoint;
			extrapolateByFactor(this->indexHighest, reflectionFactor, reflectedPoint);
			auto valueInReflectedPoint = evaluateTargetFunction(reflectedPoint);

			if (valueInReflectedPoint >= getPointValuePair(this->indexHighest).getValue())
			{
				break;
			}

			if (valueInReflectedPoint < getPointValuePair(this->indexLowest).getValue())
			{
				result = TransformResult::EXCELLENT;
			}
			else if (valueInReflectedPoint < getPointValuePair(this->indexNextHighest).getValue())
			{
				result = TransformResult::GOOD;
			}

			setPointValuePair(this->indexHighest, reflectedPoint, valueInReflectedPoint);
		}
		break;
		case TransformType::EXPANSION:
		{
			PointType expandedPoint;
			extrapolateByFactor(this->indexHighest, expansionFactor, expandedPoint);
			auto valueInExpandedPoint = evaluateTargetFunction(expandedPoint);
			if (valueInExpandedPoint < getPointValuePair(this->indexHighest).getValue())
			{
				setPointValuePair(this->indexHighest, expandedPoint, valueInExpandedPoint);
				result = TransformResult::EXCELLENT;
			}
		}
		break;
		case TransformType::CONTRACTION:
		{
			PointType contractedPoint;
			extrapolateByFactor(this->indexHighest, contractionFactor, contractedPoint);
			auto valueInContractedPoint = evaluateTargetFunction(contractedPoint);
			if (valueInContractedPoint < getPointValuePair(this->indexNextHighest).getValue())
			{
				setPointValuePair(this->indexHighest, contractedPoint, valueInContractedPoint);
				result = TransformResult::GOOD;
			}
		}
		break;
		case TransformType::MULTIPLE_CONTRACTION:
		{
			contractAroundLowestPoint(contractionFactor);
			result = TransformResult::ANY;
		}
		break;
		default:
			assert(false);
		}

		return result;
	}

	template<unsigned TSize>
	auto NelderMead<TSize>::Simplex::getRelativeTolerance() -> double
	{
		double const tiny = 1.0e-10;

		double highestValue = getPointValuePair(this->indexHighest).getValue();
		double lowestValue = getPointValuePair(this->indexLowest).getValue();

		return 2.0 * abs(highestValue - lowestValue) / (abs(highestValue) + abs(lowestValue) + tiny);
	}

	template<unsigned TSize>
	auto NelderMead<TSize>::Simplex::getPointWithLowestValue() const -> PointType const&
	{
		return getPointValuePair(this->indexLowest).getPoint();
	}

	template<unsigned TSize>
	auto NelderMead<TSize>::Simplex::updateIndices() -> void
	{
		this->indexHighest = 0;
		this->indexNextHighest = 0;
		this->indexLowest = 0;

		for (unsigned i = 0; i < this->size; ++i)
		{
			if (getPointValuePair(i) >= getPointValuePair(this->indexHighest))
			{
				this->indexNextHighest = this->indexHighest;
				this->indexHighest = i;
			}
			else if (getPointValuePair(i) > getPointValuePair(this->indexNextHighest))
			{
				if (i != this->indexHighest)
				{
					this->indexNextHighest = i;
				}
			}

			if (getPointValuePair(i) < getPointValuePair(this->indexLowest))
			{
				this->indexLowest = i;
			}
		}
	}

	template<unsigned TSize>
	auto NelderMead<TSize>::Simplex::evaluateTargetFunction(PointType const& point) -> double
	{
		this->funtionCallCount += 1;
		if (this->funtionCallCount > this->maxFunctionCallCount)
			throw std::exception("Reached max function call limit.");

		return this->targetFunction(point);
	}

	template<unsigned TSize>
	auto NelderMead<TSize>::Simplex::setPointValuePair(unsigned const i, PointType const& point) -> void
	{
		assert(i < this->size);
		this->pointValuePairs[i] = {point, evaluateTargetFunction(point)};
	}

	template<unsigned TSize>
	auto NelderMead<TSize>::Simplex::setPointValuePair(unsigned const i, PointType const& point, double const value) -> void
	{
		assert(i < this->size);
		this->pointValuePairs[i] = {point, value};
	}

	template<unsigned TSize>
	auto NelderMead<TSize>::Simplex::getPointValuePair(unsigned const i) const -> PointValuePair const&
	{
		assert(i < this->size);
		return this->pointValuePairs[i];
	}

	template<unsigned TSize>
	NelderMead<TSize>::NelderMead(double const functionTolerance)
		: NelderMead::functionTolerance{functionTolerance}
	{
		static_assert(TSize >= 2, "Invalid TSize template parameter.");
		assert(this->functionTolerance > 0);
	}

	template<unsigned TSize>
	auto NelderMead<TSize>::minimize(PointType const& initialPoint, double const initialDelta, FunctionType const targetFunction) -> PointType
	{
		try
		{
			Simplex simplex(initialPoint, initialDelta, targetFunction);

			for (;;)
			{
				simplex.updateIndices();

				if (simplex.getRelativeTolerance() < this->functionTolerance)
				{
					break;
				}

				auto result = simplex.transform(Simplex::TransformType::REFLECTION);
				if (result == Simplex::TransformResult::EXCELLENT)
				{
					simplex.transform(Simplex::TransformType::EXPANSION);
				}
				else if (result == Simplex::TransformResult::BAD)
				{
					auto result = simplex.transform(Simplex::TransformType::CONTRACTION);
					if (result == Simplex::TransformResult::BAD)
					{
						simplex.transform(Simplex::TransformType::MULTIPLE_CONTRACTION);
					}
				}
			}

			return simplex.getPointWithLowestValue();
		}
		catch (std::exception& /*e*/)
		{
			return initialPoint;
		}
	}

}}

#pragma warning(disable : WARNING_UNREFERENCED_INLINE_FUNCTION_HAS_BEEN_REMOVED)