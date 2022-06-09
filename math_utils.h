#ifndef _DIXELU_MATH_UTILS_H_
#define _DIXELU_MATH_UTILS_H_

#if (defined(__cpp_constexpr) && (__cpp_constexpr >= 201304L))
#ifndef __DIXELU_RELAXED_CONSTEXPR
#define __DIXELU_RELAXED_CONSTEXPR constexpr
#endif
#ifndef __DIXELU_ATMOST_ONE_CONSTRUCTION
#define __DIXELU_ATMOST_ONE_CONSTRUCTION __DIXELU_RELAXED_CONSTEXPR
#endif
#else
#ifndef __DIXELU_RELAXED_CONSTEXPR
#define __DIXELU_RELAXED_CONSTEXPR
#endif
#ifndef __DIXELU_ATMOST_ONE_CONSTRUCTION
#define __DIXELU_ATMOST_ONE_CONSTRUCTION static const
#endif
#endif

#ifndef __DIXELU_CONDITIONAL_CPP14_SPECIFIERS
#define __DIXELU_CONDITIONAL_CPP14_SPECIFIERS inline __DIXELU_RELAXED_CONSTEXPR
#endif

namespace dixelu
{
	namespace utils
	{
		template<typename T>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS T constexpr_abs(const T& value)
		{
			return value < T(0) ? -value : value;
		}

		template<typename T>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS const T& constexpr_min(const T& m1, const T& m2)
		{
			return m1 < m2 ? m1 : m2;
		}

		template<typename T>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS const T& constexpr_max(const T& m1, const T& m2)
		{
			return m2 < m1 ? m1 : m2;
		}

		namespace details
		{

			template<typename T, bool is_int>
			struct __try_unsigned {};

			template<typename T>
			struct __try_unsigned<T, true>
			{
				using type = typename std::make_unsigned<T>::type;
			};

			template<typename T>
			struct __try_unsigned<T, false>
			{
				using type = T;
			};

			template<typename T>
			struct try_unsigned
			{
				using type = typename
					__try_unsigned<T, std::numeric_limits<T>::is_integer>::type;
			};


			template<typename T, bool is_int>
			struct __try_signed {};

			template<typename T>
			struct __try_signed<T, true>
			{
				using type = typename std::make_signed<T>::type;
			};

			template<typename T>
			struct __try_signed<T, false>
			{
				using type = T;
			};

			template<typename T>
			struct try_signed
			{
				using type = typename
					__try_signed<T, std::numeric_limits<T>::is_integer>::type;
			};

			template<typename T, bool is_signed>
			struct __opposite_sign_type {};

			template<typename T>
			struct __opposite_sign_type<T, true>
			{
				using type = typename
					__try_unsigned<T, std::numeric_limits<T>::is_integer>::type;
			};

			template<typename T>
			struct __opposite_sign_type<T, false>
			{
				using type = typename
					__try_signed<T, std::numeric_limits<T>::is_integer>::type;
			};

			template<typename T>
			struct opposite_sign_type
			{
				using type = typename
					__opposite_sign_type<T, std::is_signed<T>::value>::type;
			};

#ifndef WITHOUT_CONSTEXPR_FUNCTIONS
			template<typename T>
			__DIXELU_CONDITIONAL_CPP14_SPECIFIERS T __safe_mul(T x, T y)
			{
				if (std::numeric_limits<T>::is_integer)
				{
					using uT = typename try_unsigned<T>::type;
					bool xlz = x < 0;
					bool ylz = y < 0;
					uT xabs(x);
					uT yabs(y);
					return (xlz == ylz) ? T(xabs * yabs) : -T(xabs * yabs);
				}
				else
				{
					return x * y;
				}
			}

			template<typename T>
			__DIXELU_CONDITIONAL_CPP14_SPECIFIERS T __sqr(T x)
			{
				if (std::numeric_limits<T>::is_integer)
				{
					using uT = typename try_unsigned<T>::type;
					auto xabs = uT(constexpr_abs<T>(x));
					return T(xabs * xabs);
				}
				else
				{
					return x * x;
				}
			}

			template<typename T>
			__DIXELU_CONDITIONAL_CPP14_SPECIFIERS T __uintpow(T x, std::size_t n)
			{
				return n == 0 ? T(1) : __safe_mul(__sqr(__uintpow(x, n >> 1)), (n & 1 ? x : T(1)));
			}

			template<typename T>
			__DIXELU_CONDITIONAL_CPP14_SPECIFIERS T __intpow(T x, std::ptrdiff_t n)
			{
				return (n < 0) ? T(1) / __uintpow(x, -n) : __uintpow(x, n);
			}

			template<typename T, unsigned int n>
			struct lookup_table
			{
				static constexpr unsigned int radix = (std::numeric_limits<T>::radix) ? std::numeric_limits<T>::radix : 2;
				static constexpr unsigned int max_degree_flt = (std::numeric_limits<T>::max_exponent / n);
				static constexpr unsigned int max_degree_int = ((std::numeric_limits<T>::digits / n));
				static constexpr unsigned int max_degree = (!max_degree_flt) ? max_degree_int : max_degree_flt;
				T lookup_table_vals[max_degree]{};
				T lookup_table_roots[max_degree]{};
				__DIXELU_RELAXED_CONSTEXPR lookup_table() :
					lookup_table_vals(), lookup_table_roots()
				{
					__DIXELU_RELAXED_CONSTEXPR T base(radix);
					__DIXELU_RELAXED_CONSTEXPR T zero(0);
					T root = base;
					lookup_table_vals[0] = zero;
					lookup_table_roots[0] = zero;
					for (unsigned int i = 1; i < max_degree; ++i)
					{
						lookup_table_vals[i] = __uintpow(root, n);
						lookup_table_roots[i] = root;
						root *= base;
					}
				}
			};

			template<typename T, unsigned int n>
			__DIXELU_CONDITIONAL_CPP14_SPECIFIERS T root_approx(T x)
			{
				__DIXELU_ATMOST_ONE_CONSTRUCTION lookup_table<T, n> table{};
				unsigned int guessed_begin = 0;
				unsigned int guessed_end = lookup_table<T, n>::max_degree - 1;
				do
				{
					auto center = (guessed_end + guessed_begin) >> 1;
					if (table.lookup_table_vals[center] < x)
						guessed_begin = center;
					else
						guessed_end = center;

				} while (guessed_end - guessed_begin > 1);
				return table.lookup_table_roots[guessed_end];
			}

			template<typename T, std::size_t n>
			__DIXELU_CONDITIONAL_CPP14_SPECIFIERS T __uintroot(T x)
			{
				__DIXELU_RELAXED_CONSTEXPR T epsilon = std::numeric_limits<T>::epsilon();
				if (x == T(0))
					return T(0);
				T n_conv(n);
				T one(1);
				T x_k(root_approx<T, n>(x));
				T x_prev_k(x_k);
				T x_p_prev_k(0);
				T x_k_nm1(0);
				T coef(n_conv - 1);
				coef /= n_conv;
				T diff(0), diff_step(0);
				T last_diff_multiplier = one;
				do
				{
					x_k_nm1 = one / __uintpow(x_k, n - 1);
					x_p_prev_k = x_prev_k;
					x_prev_k = x_k;
					x_k = coef * x_k + (x / n_conv) * x_k_nm1;
					diff = constexpr_abs(x_k - x_prev_k);
					diff_step = constexpr_abs(x_k - x_p_prev_k);
					last_diff_multiplier = (x_k > one) ? x_prev_k : one;
				} while (diff > epsilon * last_diff_multiplier && diff_step > epsilon * last_diff_multiplier);
				return x_k;
			}

			template<typename T>
			__DIXELU_CONDITIONAL_CPP14_SPECIFIERS T __positive_pow(T x, T p);

			template<typename T>
			__DIXELU_CONDITIONAL_CPP14_SPECIFIERS T __frac_positive_pow(T x, T p)
			{
				constexpr unsigned int rolling_ppow_bits = 4u;
				constexpr unsigned int rolling_power = ((1u << rolling_ppow_bits));
				__DIXELU_RELAXED_CONSTEXPR T epsilon = std::numeric_limits<T>::epsilon();
				if (p <= epsilon)
					return T(1);
				if (p >= T(1))
					return __positive_pow(x, p);
				return __uintroot<T, rolling_power>(__positive_pow<T>(x, p * rolling_power));
			}

			template<typename T>
			__DIXELU_CONDITIONAL_CPP14_SPECIFIERS T __arctanh_naive(T x)
			{
				T xsq = x * x;
				T xdeg = x;
				T sum = x;
				T rec = 1;
				T prev_sum = 0;
				while (prev_sum != sum)
				{
					prev_sum = sum;
					xdeg *= xsq;
					rec += T(2);
					sum += xdeg / rec;
				}
				return sum;
			}

			template<typename T>
			__DIXELU_CONDITIONAL_CPP14_SPECIFIERS T __log_e_naive(T x)
			{
				x -= T(1);
				return 2 * __arctanh_naive(x / (2 + x));
			}

			template<typename T>
			__DIXELU_CONDITIONAL_CPP14_SPECIFIERS T __log(T x)
			{
				constexpr unsigned int rolling_ppow_bits = 4u;
				constexpr unsigned int rolling_power = ((1u << rolling_ppow_bits));
				constexpr T log_top_edge = T(2);
				constexpr T log_bottom_edge = T(1) / log_top_edge;
				if (x == T(0))
					return -std::numeric_limits<T>::infinity();

				T multiplier = T(1);
				while (x > log_top_edge || x < log_bottom_edge)
				{
					multiplier *= rolling_power;
					x = __uintroot<T, rolling_power>(x);
				}
				return multiplier * __log_e_naive(x);
			}

			template<typename T>
			__DIXELU_CONDITIONAL_CPP14_SPECIFIERS T __exp_naive(T x = 1)
			{
				T sum(1);
				T prev_sum(0);
				T fact(1);
				T xdeg(x);
				std::size_t iter(1);
				while (prev_sum != sum)
				{
					prev_sum = sum;
					sum += xdeg / fact;
					xdeg *= x;
					++iter;
					fact *= iter;
				}
				return sum;
			}

			template<typename T>
			__DIXELU_CONDITIONAL_CPP14_SPECIFIERS T __exp(T x)
			{
				__DIXELU_RELAXED_CONSTEXPR T __e = __exp_naive<T>();
				const std::ptrdiff_t x_z = x;
				const T dx = x - x_z;
				const auto e_x_z = __intpow(__e, x_z);
				const auto e_dx = __exp_naive(dx);
				return e_x_z * e_dx;
			}

			template<typename T>
			__DIXELU_CONDITIONAL_CPP14_SPECIFIERS T __positive_pow(T x, T p)
			{
				std::size_t p_w(p);
				auto d_p = p - T(p_w);
				auto uintpow = __uintpow(x, p_w);
				auto fracpow = __frac_positive_pow(x, d_p);
				return uintpow * fracpow;
			}

			template<typename T>
			__DIXELU_CONDITIONAL_CPP14_SPECIFIERS T __pow_sqroot(T a, T b)
			{
				if (b < 0)
					return T(1) / __positive_pow(a, -b);
				return __positive_pow(a, b);
			}

			template<typename T>
			__DIXELU_CONDITIONAL_CPP14_SPECIFIERS T __pow_explog(T a, T b)
			{
				return __exp(b * __log(a));
			}
#endif
		} // namespace details

		template<typename T>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS T constexpr_sqrt(T x)
		{
#ifdef WITHOUT_CONSTEXPR_FUNCTIONS
			return std::sqrt(x);
#else
			return details::__uintroot<T, 2>(x);
#endif // WITHOUT_CONSTEXPR_FUNCTIONS
		}

		template<typename T>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS T constexpr_pow(T a, T b)
		{
#ifdef WITHOUT_CONSTEXPR_FUNCTIONS
			return std::pow(a, b);
#else
			return details::__pow_explog(a, b);
#endif // WITHOUT_CONSTEXPR_FUNCTIONS
		}

		template<typename T>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS T constexpr_intpow(T a, std::ptrdiff_t b)
		{
			std::ptrdiff_t sign = 1;
			bool inverse = false;
			if (a < 0)
			{
				sign = -1;
				a *= sign;
			}
			if (b < 0)
			{
				b *= -1;
				inverse = true;
			}

			auto v = utils::details::__uintpow(a, b);
			return T(sign) * ((inverse) ? T(1) / v : v);
		}
	} // namespace utils
}


#endif // _DIXELU_MATH_UTILS_H_