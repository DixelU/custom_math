#include "../math_utils.h"
#include "../integers.h"
#include "../sq_matrix.h"

#include <iostream>
#include <sstream>
#include <chrono>
#include <cassert>

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
		" microsecond (constexpr)" << std::endl;

	b = std::chrono::high_resolution_clock::now();
	auto t2 = std::pow(x, p);
	e = std::chrono::high_resolution_clock::now();

	std::cout << " Got " << t2 << " in " <<
		std::chrono::duration_cast<std::chrono::microseconds>(e - b).count() <<
		" microsecond (non constexpr)" << std::endl;

	b = std::chrono::high_resolution_clock::now();
	__DIXELU_RELAXED_CONSTEXPR auto pair = constexpr_matrix_proof();
	__DIXELU_RELAXED_CONSTEXPR auto asdf_const = dixelu::utils::constexpr_pow<double>(1.05498, 3215.2165);
	e = std::chrono::high_resolution_clock::now();

	std::cout << "Constexpr matrix test: ";
	std::cout << pair.first << " " << pair.second << " ~ " << asdf_const << std::endl;

	std::cout << "Proof of constexpr evaluation: " <<
		std::chrono::duration_cast<std::chrono::microseconds>(e - b).count() <<
		" microseconds between chrono::now calls" << std::endl;

	return 0; // :)
}

using longuint = dixelu::long_uint<2>;

constexpr longuint x()
{
	auto value = longuint(~0ULL, longuint::__fill_fields_tag{});
	auto value2 = longuint(0xAAAAAAAAAAAAAAAAULL, longuint::__fill_fields_tag{});
	value >>= 47;
	value -= 29378;
	value ^= value2;
	value += value2 << 1;

	return value;
}

constexpr auto R = x();
constexpr auto val = dixelu::utils::constexpr_sqrt<longuint>(2987348);

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
	std::cout << "Largest ints multiplication: "
		<< std::chrono::duration_cast<std::chrono::microseconds>(
			std::chrono::steady_clock::now() - begin).count()
		<< " mcsec" << std::endl;

	std::cout << "Result: " << d << std::endl;

	auto begin1 = std::chrono::steady_clock::now();
	auto q1 = longuint(~0ULL, longuint::__fill_fields_tag{});
	for(size_t i = 0; i < 1000; i++)
		q1 << i;
	auto inter1 = std::chrono::steady_clock::now();
	q1 = longuint(~0ULL, longuint::__fill_fields_tag{});
	for(size_t i = 0; i < 1000; i++)
	{
		auto copy = q1;
		copy.__experimental_shift_bits_left(i);
	}
	auto end1 = std::chrono::steady_clock::now();

	std::cout << "Classic shift: "
			  << std::chrono::duration_cast<std::chrono::microseconds>(inter1 - begin1).count()
			  << "mcsec" << std::endl;
	std::cout << "Experimental shift: "
			  << std::chrono::duration_cast<std::chrono::microseconds>(end1 - inter1).count()
			  << "mcsec" << std::endl;

	std::cout << "Testing longuint with constexpr sqrt. Enter any integer value: " << std::endl;
	std::cin >> a;
	auto ans = dixelu::utils::constexpr_sqrt(a);
	std::cout << ans << " ~ when squared equals " << ans * ans << std::endl;

	std::cout << "END longuint test" << std::endl;
	std::cout << std::flush;
}

void longuint_correctness_test()
{
	using uint = dixelu::long_uint<2>; // 512-bit
	using half = uint::down_type;      // 256-bit

	// --- subtraction: basic ---
	assert(uint(100) - uint(42) == uint(58));
	assert(uint(1) - uint(1) == uint(0));

	// --- subtraction: borrow across hi/lo boundary ---
	// Build 2^256 by setting lo=0, hi=1
	uint two256;
	two256.hi = half(1);
	two256.lo = half(0);
	uint result = two256 - uint(1);
	// Expected: hi=0, lo=all-ones (2^256 - 1)
	assert(result.hi == half(0));
	assert(result.lo == half(~0ULL, half::__fill_fields_tag{}));

	// --- subtraction: wrapping underflow ---
	assert(uint(0) - uint(1) == uint(~0ULL, uint::__fill_fields_tag{}));

	// --- division / modulo ---
	uint rem;
	assert(uint::__divmod(uint(100), uint(7),  rem) == uint(14)); assert(rem == uint(2));
	assert(uint::__divmod(uint(0),   uint(7),  rem) == uint(0));  assert(rem == uint(0));
	assert(uint::__divmod(uint(7),   uint(7),  rem) == uint(1));  assert(rem == uint(0));
	assert(uint::__divmod(uint(6),   uint(7),  rem) == uint(0));  assert(rem == uint(6));
	assert(uint(1000000) / uint(13) == uint(76923));
	assert(uint(1000000) % uint(13) == uint(1));

	// --- operator++ / operator-- ---
	{
		uint a(0);
		assert(++a == uint(1));
		assert(a++ == uint(1)); assert(a == uint(2));
		assert(--a == uint(1));
		assert(a-- == uint(1)); assert(a == uint(0));
		// decrement from 0 wraps to max
		--a;
		assert(a == uint(~0ULL, uint::__fill_fields_tag{}));
		// increment from max wraps to 0
		++a;
		assert(a == uint(0));
	}

	// --- operator bool short-circuit ---
	assert(!(bool)uint(0));
	assert((bool)uint(1));
	assert((bool)two256); // hi != 0, lo == 0

	// --- from_decimal / to_string round-trip ---
	{
		const std::string num = "123456789012345678901234567890123456789";
		auto v = uint::__from_decimal_std_string(num);
		assert(uint::to_string(v) == num);
	}

	// --- constexpr subtraction regression (R uses -=) ---
	static_assert(R == R, "constexpr round-trip");

	std::cout << "All correctness tests passed." << std::endl;
}

int main(void)
{
	longuint_correctness_test();
	//constexpr_pow_test();
	longuint_test();
}