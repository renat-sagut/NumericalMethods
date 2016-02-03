#include "Vector.h"
using numerical::basic::Vector2d;
#include "NelderMead.h"
using numerical::optimization::NelderMead2d;
#include "Warnings.h"

#pragma warning(push)
#pragma warning(disable : WARNING_BYTES_PADDING_ADDED)
#pragma warning(disable : WARNING_BEHAVIOR_CHANGE)

#include <cmath>
#include <iostream>

#pragma warning(pop)

int main()
{
	NelderMead2d nelder(1e-6);

	/*// f(x,y) = x^3 + y^3 - 3x - 3y + 5
	NelderMead2d::FunctionType targetFunction = [](Vector2d const& point) -> double {
		return pow(point[0], 3.0) + pow(point[1], 3.0) - 3.0 * point[0] - 3.0 * point[1] + 5.0;
	};*/

	/*// f(x,y) = x^2 + y^2 - 4x - y - xy
	// Точка минимума: (3, 2)
	// Минимальное значение функции: -7
	NelderMead2d::FunctionType targetFunction = [](Vector2d const& point) -> double {
		return pow(point[0], 2.0) + pow(point[1], 2.0) - 4.0 * point[0] - point[1] - point[0] * point[1];
	};*/

	// Himmelblau's function
	// f(x,y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2
	// It has four identical local minima:
	// f(3.0, 2.0) = 0.0,
	// f(-2.805118, 3.131312) = 0.0,
	// f(-3.779310, -3.283186) = 0.0,
	// f(3.584428, -1.848126) = 0.0.
	NelderMead2d::FunctionType targetFunction = [](Vector2d const& point) -> double {
		return pow(pow(point[0], 2.0) + point[1] - 11.0, 2.0) + pow(point[0] + pow(point[1], 2.0) - 7.0, 2.0);
	};

	auto solutionPoint = nelder.minimize({0, 0}, 0.5, targetFunction);
	//auto solutionPoint = nelder.minimize({-5, 5}, 0.5, targetFunction);
	//auto solutionPoint = nelder.minimize({-5, 0}, 0.5, targetFunction);
	//auto solutionPoint = nelder.minimize({3.0, -2.0}, 0.5, targetFunction);

	std::cout << "Solution point is (" << solutionPoint[0] << ", " << solutionPoint[1] << ")" << std::endl;
	std::cout << "Function value is " << targetFunction(solutionPoint) << std::endl;

	return 0;
}

#pragma warning(disable : WARNING_UNREFERENCED_INLINE_FUNCTION_HAS_BEEN_REMOVED)
#pragma warning(disable : WARNING_FUNCTION_NOT_INLINED)