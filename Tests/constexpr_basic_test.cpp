#include "../math_utils.h"
#include "../integers.h"
#include "../sq_matrix.h"

#include <iostream>
#include <sstream>
#include <chrono>

__DIXELU_RELAXED_CONSTEXPR std::pair<float, float> constexpr_matrix_proof()
{
	dixelu::sq_matrix<float, 3> sq;
	sq[2][2] = 2;
	sq[1][1] = 5;
	sq[1][0] = 1;
	sq[0][0] = 3;
	sq ^= -3;
	sq.transpose();
	auto sq5 = sq.to<5>();
	auto n = sq5[0].get_norm(1.f / 3);
	return { sq.determinant(), n };
}

int constexpr_pow_test()
{
	std::cout << "Test constexpr pow in runtime. Enter any two FP values: " << std::flush;
	double x, p;
	std::cin >> x >> p;

	auto b = std::chrono::high_resolution_clock::now();
	auto t = dixelu::utils::constexpr_pow<double>(x, p);
	auto e = std::chrono::high_resolution_clock::now();

	std::cout << " Got " << t << " in " <<
		std::chrono::duration_cast<std::chrono::microseconds>(e - b).count() <<
		" microsec (constexpr)" << std::endl;

	b = std::chrono::high_resolution_clock::now();
	auto t2 = std::pow(x, p);
	e = std::chrono::high_resolution_clock::now();

	std::cout << " Got " << t << " in " <<
		std::chrono::duration_cast<std::chrono::microseconds>(e - b).count() <<
		" microsec (nonconstexpr)" << std::endl;

	b = std::chrono::high_resolution_clock::now();
	__DIXELU_RELAXED_CONSTEXPR auto pair = constexpr_matrix_proof();
	__DIXELU_RELAXED_CONSTEXPR auto asdf_const = dixelu::utils::constexpr_pow<double>(1.05498, 3215.2165);
	e = std::chrono::high_resolution_clock::now();

	std::cout << "Constexpr matrix test: ";
	std::cout << pair.first << " " << pair.second << " ~ " << asdf_const << std::endl;

	std::cout << "Proof of constexpr evauation: " <<
		std::chrono::duration_cast<std::chrono::microseconds>(e - b).count() <<
		" microseconds between chrono::now calls" << std::endl;

	return 0; // :)
}

using longuint = dixelu::long_uint<10>;

void longuint_test()
{
	//std::string damn(100000, '1');
	std::cout << "Test static long arithmetic with " << std::numeric_limits<longuint>::digits10 << " digits" << std::endl;

	std::cout << "Generating number..." << std::endl;
	auto b = longuint(~0ull, longuint::__fill_fields_tag{});//longuint::__from_decimal_std_string(damn);
	auto a = b;
	std::cout << "Begin max*max multiplication..." << std::endl;

	auto begin = std::chrono::steady_clock::now();
	auto d = a * b;
	std::cout << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - begin).count() << " mcsec" << std::endl;

	std::cout << d << std::endl;

	std::cout << "Testing longuint with constexpr sqrt. Enter any integer value: " << std::endl;
	std::cin >> a;
	auto ans = dixelu::utils::constexpr_sqrt(a);
	std::cout << ans << " ~ when squared equals " << ans * ans << std::endl;

	std::cout << "END longuint test" << std::endl;
	std::cout << std::flush;
}

int main(void)
{
	constexpr_pow_test();
	longuint_test();
}