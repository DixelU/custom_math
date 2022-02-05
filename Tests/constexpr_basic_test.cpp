
#include <iostream>
#include <chrono>

#include "../sq_matrix.h"

constexpr std::pair<float, float> some_function()
{
	dixelu::sq_matrix<float, 3> sq;
	sq[2][2] = 2;
	sq[1][1] = 5;
	sq[1][0] = 1;
	sq[0][0] = 3;
	sq ^= -3;
	auto n = sq[0].get_norm(1.f / 3);
	return { sq.determinant(), n };
}

int main()
{
	double x, p;
	std::cin >> x >> p;

	auto b = std::chrono::high_resolution_clock::now();
	auto t = dixelu::utils::constexpr_pow<double>(x, p);
	auto e = std::chrono::high_resolution_clock::now();

	std::cout << " Got " << t << " in " << std::chrono::duration_cast<std::chrono::microseconds>(e - b).count() << " microsec (constexpr)" << std::endl;

	b = std::chrono::high_resolution_clock::now();
	auto t2 = std::pow(x, p);
	e = std::chrono::high_resolution_clock::now();

	std::cout << " Got " << t << " in " << std::chrono::duration_cast<std::chrono::microseconds>(e - b).count() << " microsec (nonconstexpr)" << std::endl;

	constexpr auto pair = some_function();
	constexpr auto asdf_const = dixelu::utils::constexpr_pow<double>(1.05498, 3215.2165);

	std::cout << pair.first << " " << pair.second << " ~ " << asdf_const << std::endl;
	return 0; // :)
}