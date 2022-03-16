#ifndef _DIXELU_SQ_MATRIX_H_
#define _DIXELU_SQ_MATRIX_H_

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>
#include <limits>   
#include <ostream>
#include <iomanip>
#include <functional>
#include <type_traits>
#include <initializer_list>

#define __DEFAULT_DIXELU_FUNC_SPECIFIERS inline

#if (__cplusplus >= 201402L) || (defined(_MSVC_LANG) && _MSVC_LANG >= 201402L)
#define __DIXELU_COND_CONSTEXPR constexpr
#else
#define __DIXELU_COND_CONSTEXPR
#endif

#define __DIXELU_CONDITIONAL_SPECIFIERS __DEFAULT_DIXELU_FUNC_SPECIFIERS __DIXELU_COND_CONSTEXPR

namespace dixelu
{

	namespace utils
	{

		template<typename T>
		__DIXELU_CONDITIONAL_SPECIFIERS T constexpr_abs(const T& value)
		{
			return value < 0 ? -value : value;
		}

		template<typename T>
		__DIXELU_CONDITIONAL_SPECIFIERS const T& constexpr_min(const T& m1, const T& m2)
		{
			return m1 < m2 ? m1 : m2;
		}

		template<typename T>
		__DIXELU_CONDITIONAL_SPECIFIERS const T& constexpr_max(const T& m1, const T& m2)
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
			__DIXELU_CONDITIONAL_SPECIFIERS T __safe_mul(T x, T y)
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
			__DIXELU_CONDITIONAL_SPECIFIERS T __sqr(T x)
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
			__DIXELU_CONDITIONAL_SPECIFIERS T __uintpow(T x, std::size_t n)
			{
				return n == 0 ? 1 : __safe_mul(__sqr(__uintpow(x, n >> 1)), (n & 1 ? x : T(1)));
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
				__DIXELU_COND_CONSTEXPR lookup_table() :
					lookup_table_vals(), lookup_table_roots()
				{
					__DIXELU_COND_CONSTEXPR T base(radix);
					__DIXELU_COND_CONSTEXPR T zero(0);
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
			__DIXELU_CONDITIONAL_SPECIFIERS T root_approx(T x)
			{
				__DIXELU_COND_CONSTEXPR lookup_table<T, n> table{};
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
			__DIXELU_CONDITIONAL_SPECIFIERS T __uintroot(T x)
			{
				__DIXELU_COND_CONSTEXPR T epsilon = std::numeric_limits<T>::epsilon();
				if (x <= epsilon)
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
			__DIXELU_CONDITIONAL_SPECIFIERS T __positive_pow(T x, T p);

			template<typename T>
			__DIXELU_CONDITIONAL_SPECIFIERS T __frac_positive_pow(T x, T p)
			{
				constexpr unsigned int rolling_ppow_bits = 4u;
				constexpr unsigned int rolling_power = ((1u << rolling_ppow_bits));
				__DIXELU_COND_CONSTEXPR T epsilon = std::numeric_limits<T>::epsilon();
				if (p <= epsilon)
					return T(1);
				if (p >= T(1))
					return __positive_pow(x, p);
				return __uintroot<T, rolling_power>(__positive_pow<T>(x, p * rolling_power));
			}

			template<typename T>
			__DIXELU_CONDITIONAL_SPECIFIERS T __positive_pow(T x, T p)
			{
				std::size_t p_w(p);
				auto d_p = p - T(p_w);
				auto uintpow = __uintpow(x, p_w);
				auto fracpow = __frac_positive_pow(x, d_p);
				return uintpow * fracpow;
			}

			template<typename T>
			__DIXELU_CONDITIONAL_SPECIFIERS T __pow(T a, T b)
			{
				if (b < 0)
					return T(1) / __positive_pow(a, -b);
				return __positive_pow(a, b);
			}

#endif
		} // namespace details

		template<typename T>
		__DIXELU_CONDITIONAL_SPECIFIERS T constexpr_sqrt(T x)
		{
#ifdef WITHOUT_CONSTEXPR_FUNCTIONS
			return std::sqrt(x);
#else
			return details::__uintroot<T, 2>(x);
#endif // WITHOUT_CONSTEXPR_FUNCTIONS
		}

		template<typename T>
		__DIXELU_CONDITIONAL_SPECIFIERS T constexpr_pow(T a, T b)
		{
#ifdef WITHOUT_CONSTEXPR_FUNCTIONS
			return std::pow(a, b);
#else
			return details::__pow(a, b);
#endif // WITHOUT_CONSTEXPR_FUNCTIONS
		}

		template<typename T>
		__DIXELU_CONDITIONAL_SPECIFIERS T constexpr_intpow(T a, std::ptrdiff_t b)
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

	template<typename general_float_type, std::size_t dims>
	struct point
	{
		general_float_type base_array[dims];
		using general_paired_ftype = typename
			utils::details::opposite_sign_type<general_float_type>::type;
		using self_type = point<general_float_type, dims>;
		__DIXELU_COND_CONSTEXPR point() :
			base_array()
		{
			for (std::size_t i = 0; i < dims; ++i)
				base_array[i] = general_float_type();
		}
		__DIXELU_COND_CONSTEXPR explicit point(general_float_type v) :
			base_array()
		{
			for (std::size_t i = 0; i < dims; ++i)
				base_array[i] = v;
		}
		__DIXELU_COND_CONSTEXPR point(const std::initializer_list<general_paired_ftype>& il_d) :
			base_array()
		{
			auto y = il_d.begin();
			for (std::size_t i = 0; i < dims && y != il_d.end(); ++i, ++y)
				base_array[i] = *y;
		}
		__DIXELU_COND_CONSTEXPR point(const std::initializer_list<int>& il) :
			base_array()
		{
			auto y = il.begin();
			for (std::size_t i = 0; i < dims && y != il.end(); ++i, ++y)
				base_array[i] = *y;
		}
		__DIXELU_COND_CONSTEXPR point(const std::vector<general_paired_ftype>& il_d) :
			base_array()
		{
			auto y = il_d.cbegin();
			for (std::size_t i = 0; i < dims && y != il_d.end(); ++i, ++y)
				base_array[i] = *y;
		}
		__DIXELU_COND_CONSTEXPR point(const std::vector<int>& il) :
			base_array()
		{
			auto y = il.cbegin();
			for (std::size_t i = 0; i < dims && y != il.end(); ++i, ++y)
				base_array[i] = *y;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS general_float_type get_norm2() const
		{
			general_float_type sum = general_float_type();
			for (std::size_t i = 0; i < dims; ++i)
				sum += base_array[i] * base_array[i];
			return sum;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS general_float_type get_norm(general_float_type deg) const
		{
			general_float_type sum = general_float_type();
			for (std::size_t i = 0; i < dims; ++i)
				sum += utils::constexpr_pow<general_float_type>(base_array[i], deg);
			return utils::constexpr_pow<general_float_type>(sum, 1. / deg);
		}
		__DIXELU_CONDITIONAL_SPECIFIERS general_float_type get_norm() const
		{
			return utils::constexpr_sqrt(get_norm2());
		}
		__DIXELU_CONDITIONAL_SPECIFIERS std::size_t get_dims() const
		{
			return dims;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS general_float_type& operator[](std::size_t D)
		{
			return base_array[D];
		}
		__DIXELU_CONDITIONAL_SPECIFIERS const general_float_type& operator[](std::size_t D) const
		{
			return base_array[D];
		}
		__DIXELU_CONDITIONAL_SPECIFIERS void swap(self_type& P)
		{
			for (std::size_t i = 0; i < dims; ++i)
			{
				auto value = base_array[i];
				base_array[i] = P.base_array[i];
				P.base_array[i] = value;
			}
		}
		__DIXELU_CONDITIONAL_SPECIFIERS self_type reverse() const
		{
			self_type N;
			for (std::size_t i = 0; i < dims; ++i)
				N[i] = base_array[dims - 1 - i];
			return N;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS self_type operator+(const self_type& P) const
		{
			self_type N;
			for (std::size_t i = 0; i < dims; ++i)
				N[i] = base_array[i] + P[i];
			return N;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS self_type& operator+=(const self_type& P)
		{
			for (std::size_t i = 0; i < dims; ++i)
				base_array[i] += P[i];
			return *this;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS self_type operator-(const self_type& P) const
		{
			self_type N;
			for (std::size_t i = 0; i < dims; ++i)
				N[i] = base_array[i] - P[i];
			return N;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS self_type& operator-=(const self_type& P)
		{
			for (std::size_t i = 0; i < dims; ++i)
				base_array[i] -= P[i];
			return *this;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS self_type operator*(general_float_type M) const
		{
			self_type N;
			for (std::size_t i = 0; i < dims; ++i)
				N[i] = base_array[i] * M;
			return N;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS self_type& operator*=(general_float_type M)
		{
			for (std::size_t i = 0; i < dims; ++i)
				base_array[i] *= M;
			return *this;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS self_type operator/(general_float_type M) const
		{
			self_type N;
			for (std::size_t i = 0; i < dims; ++i)
				N[i] = base_array[i] / M;
			return N;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS self_type& operator/=(general_float_type M)
		{
			return ((*this) *= (general_float_type(1) / M));
		}
		__DIXELU_CONDITIONAL_SPECIFIERS general_float_type operator*(const self_type& P) const
		{
			general_float_type sum = general_float_type();
			for (std::size_t i = 0; i < dims; ++i)
				sum += base_array[i] * P[i];
			return sum;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS self_type& operator-()
		{
			for (std::size_t i = 0; i < dims; ++i)
				base_array[i] = 0 - base_array[i];
			return *this;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS bool operator<(const self_type& P) const
		{
			for (std::size_t i = 0; i < dims; ++i)
				if (base_array[i] >= P.base_array[i])
					return false;
			return true;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS bool operator>=(const self_type& P) const
		{
			return !(*this < P);
		}
		__DIXELU_CONDITIONAL_SPECIFIERS bool operator>(const self_type& P) const
		{
			return P < *this;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS bool operator==(const self_type& P) const
		{
			for (std::size_t i = 0; i < dims; ++i)
				if (base_array[i] != P[i]) return false;
			return true;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS bool operator!=(const self_type& P) const
		{
			return !(*this == P);
		}
		__DIXELU_CONDITIONAL_SPECIFIERS bool operator<=(const self_type& P) const
		{
			return !(P < *this);
		}
		__DIXELU_CONDITIONAL_SPECIFIERS self_type normalize() const
		{
			return (*this) / get_norm();
		}
	};

	template<typename general_float_type, std::size_t dims>
	std::ostream& operator<<(std::ostream& os, const point<general_float_type, dims>& P)
	{
		os << "(";
		for (std::size_t i = 0; i < dims; ++i)
		{
			os << P[i];
			if (i != dims - 1)
				os << ",";
		}
		os << ")";
		return os;
	}

	template<typename general_float_type, std::size_t dims>
	__DIXELU_CONDITIONAL_SPECIFIERS point<general_float_type, dims> operator*(general_float_type M, point<general_float_type, dims> P)
	{
		point<general_float_type, dims> N;
		for (std::size_t i = 0; i < dims; ++i)
			N[i] = P[i] * M;
		return N;
	}

	template<typename general_float_type, std::size_t dims>
	struct sq_matrix
	{
		static constexpr general_float_type GFLOAT_EPSILON = std::numeric_limits<general_float_type>::epsilon();
		general_float_type utilisation = general_float_type();

		using point_type = point<general_float_type, dims>;
		using self_type = sq_matrix<general_float_type, dims>;

		point_type base_array[dims];

		__DIXELU_COND_CONSTEXPR sq_matrix() :
			base_array()
		{
			for (std::size_t i = 0; i < dims; ++i)
				base_array[i] = point_type();
		}
		__DIXELU_COND_CONSTEXPR explicit sq_matrix(general_float_type E_num) :
			base_array()
		{
			for (std::size_t i = 0; i < dims; ++i)
			{
				base_array[i] = point_type();
				base_array[i][i] = E_num;
			}
		}
		__DIXELU_COND_CONSTEXPR sq_matrix(std::initializer_list<point_type> IL) :
			base_array()
		{
			std::size_t id = 0;
			for (auto& p : IL)
			{
				if (id == dims)
					break;
				base_array[id] = p;
				++id;
			}
		}
		void swap(self_type& p)
		{
			for (std::size_t i = 0; i < dims; ++i)
				p[i].swap(base_array[i]);
		}
		__DIXELU_CONDITIONAL_SPECIFIERS point_type operator*(const point_type& p) const
		{
			point_type T;
			for (std::size_t i = 0; i < dims; ++i)
				T[i] = base_array[i] * p;
			return T;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS self_type operator*(general_float_type num) const
		{
			self_type T;
			for (std::size_t i = 0; i < dims; ++i)
				T[i] = base_array[i] * num;
			return T;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS self_type operator*=(general_float_type num)
		{
			for (std::size_t i = 0; i < dims; ++i)
				base_array[i] *= num;
			return *this;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS self_type operator/(general_float_type num) const
		{
			return *this * (1. / num);
		}
		__DIXELU_CONDITIONAL_SPECIFIERS self_type& operator/=(general_float_type num)
		{
			return ((*this) *= (1. / num));
		}
		__DIXELU_CONDITIONAL_SPECIFIERS self_type operator+(const self_type& p) const {
			self_type T;
			for (std::size_t i = 0; i < dims; ++i)
				T[i] = base_array[i] + p[i];
			return T;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS self_type operator-(const self_type& p) const
		{
			self_type T;
			for (std::size_t i = 0; i < dims; ++i)
				T[i] = base_array[i] - p[i];
			return T;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS self_type& operator+=(const self_type& p)
		{
			for (std::size_t i = 0; i < dims; ++i)
				base_array[i] += p[i];
			return *this;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS self_type& operator-=(const self_type& p)
		{
			for (std::size_t i = 0; i < dims; ++i)
				base_array[i] -= p[i];
			return *this;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS const point_type& operator[](std::size_t i) const
		{
			return base_array[i];
		}
		__DIXELU_CONDITIONAL_SPECIFIERS point_type& operator[](std::size_t i)
		{
			return base_array[i];
		}
		__DIXELU_CONDITIONAL_SPECIFIERS general_float_type& at(std::size_t point_id, std::size_t coordinate)
		{
			if (point_id < dims && coordinate < dims)
			{
				return base_array[point_id][coordinate];
			}
			else return utilisation;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS const general_float_type& at(std::size_t point_id, std::size_t coordinate) const
		{
			if (point_id < dims && coordinate < dims)
			{
				return base_array[point_id][coordinate];
			}
			else return utilisation;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS self_type operator*(const self_type& M) const
		{
			self_type P;
			for (std::size_t y = 0; y < dims; ++y) {
				for (std::size_t x = 0; x < dims; ++x) {
					for (std::size_t i = 0; i < dims; ++i) {
						P[y][x] += base_array[y][i] * M[i][x];
					}
				}
			}
			return P;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS self_type inverse() const
		{
			general_float_type max_value = general_float_type();
			general_float_type mul = general_float_type();
			std::size_t id = 0;
			self_type E(1), A(*this);
			for (std::size_t step = 0; step < dims; ++step)
			{
				id = step;
				max_value = 0;
				for (std::size_t coid = step; coid < dims; ++coid)
				{
					if (utils::constexpr_abs(max_value) < utils::constexpr_abs(A.at(coid, step)))
					{
						max_value = A.at(coid, step);
						id = coid;
					}
				}
				if (utils::constexpr_abs(max_value) <= GFLOAT_EPSILON)
					return self_type();
				if (id != step)
				{
					A[id].swap(A[step]);
					E[id].swap(E[step]);
				}
				for (std::size_t sum_id = 0; sum_id < dims; ++sum_id)
				{
					if (sum_id == step)
						continue;
					mul = A.at(sum_id, step) / A.at(step, step);
					A[sum_id] -= (A[step] * mul);
					E[sum_id] -= (E[step] * mul);
				}
				mul = A.at(step, step);
				A[step] /= mul;
				E[step] /= mul;
			}
			return E;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS general_float_type determinant() const
		{
			general_float_type determ = 1;
			general_float_type max_value = 0, mul = 0;
			std::size_t id = 0;
			self_type A(*this);
			for (std::size_t step = 0; step < dims; ++step)
			{
				id = step;
				max_value = 0;
				for (std::size_t coid = step; coid < dims; ++coid)
				{
					if (utils::constexpr_abs(max_value) < utils::constexpr_abs(A.at(coid, step)))
					{
						max_value = A.at(coid, step);
						id = coid;
					}
				}
				if (utils::constexpr_abs(max_value) <= GFLOAT_EPSILON)
					return 0;
				if (id != step)
					A[id].swap(A[step]);
				for (std::size_t sum_id = 0; sum_id < dims; ++sum_id)
				{
					if (sum_id == step)
						continue;
					mul = A.at(sum_id, step) / A.at(step, step);
					A[sum_id] -= (A[step] * mul);
				}
				mul = A.at(step, step);
				A[step] /= mul;
				determ *= mul;
			}
			return determ;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS static point_type solve_using_eulers_method(self_type A, point_type P)
		{
			general_float_type max_value = 0, mul = 0;
			std::size_t id = 0;
			for (std::size_t step = 0; step < dims; ++step)
			{
				id = step;
				max_value = 0;
				for (std::size_t coid = step; coid < dims; ++coid)
				{
					if (utils::constexpr_abs(max_value) < utils::constexpr_abs(A.at(coid, step)))
					{
						max_value = A.at(coid, step);
						id = coid;
					}
				}
				if (utils::constexpr_abs(max_value) <= GFLOAT_EPSILON)
					return point_type();
				if (id != step)
				{
					A[id].swap(A[step]);
					P[id].swap(P[step]);
				}
				for (std::size_t sum_id = 0; sum_id < dims; ++sum_id)
				{
					if (sum_id == step)
						continue;
					mul = A.at(sum_id, step) / A.at(step, step);
					A[sum_id] -= A[step] * mul;
					P[sum_id] -= P[step] * mul;
				}
				mul = A.at(step, step);
				A[step] /= mul;
				P[step] /= mul;
			}
			return P;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS self_type operator^(int degree)
		{
			bool inverse = false;
			if (degree < 0)
			{
				inverse = true;
				degree = -degree;
			}
			else if (!degree)
				return self_type(1.);
			self_type cur_matrix(1.), deg_co_matrix(*this);
			while (degree)
			{
				switch (degree & 1)
				{
				case 1:
					cur_matrix = cur_matrix * deg_co_matrix;
					deg_co_matrix = deg_co_matrix * deg_co_matrix;
					break;
				case 0:
					deg_co_matrix = deg_co_matrix * deg_co_matrix;
					break;
				}
				degree >>= 1;
			}
			return inverse ? cur_matrix.inverse() : cur_matrix;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS self_type& operator^=(int degree)
		{
			return ((*this) = (*this) ^ degree), * this;
		}
		__DIXELU_CONDITIONAL_SPECIFIERS sq_matrix<general_float_type, dims - 1>
			minor_matrix(const std::size_t& x_minor, const std::size_t& y_minor) const
		{
			auto minor_index = [](std::size_t x, std::size_t y, std::size_t minor_x, std::size_t minor_y) ->
				std::pair<int64_t, int64_t>
			{
				if (x == minor_x)
					return { -1, -1 };
				if (y == minor_y)
					return { -1, -1 };
				if (x > minor_x)
					x -= 1;
				if (y > minor_y)
					y -= 1;
				return { x, y };
			};
			sq_matrix<general_float_type, dims - 1> M;
			for (std::size_t y = 0; y < dims; ++y)
			{
				for (std::size_t x = 0; x < dims; ++x)
				{
					auto mxy = minor_index(x, y, x_minor, y_minor);
					int64_t mx = mxy.first;
					int64_t my = mxy.second;
					if (mx < 0 || my < 0)
						continue;
					M.at(my, mx) = base_array[y][x];
				}
			}
			return M;
		}

		template<std::size_t new_dims>
		__DIXELU_CONDITIONAL_SPECIFIERS sq_matrix<general_float_type, new_dims>
			to(std::size_t start_index = 0)
		{
			sq_matrix<general_float_type, new_dims> new_matrix;

			const auto max_iter_pos = utils::constexpr_min(dims, new_dims - start_index);

			for (std::size_t y = start_index; y < max_iter_pos; ++y)
				for (std::size_t x = start_index; x < max_iter_pos; ++x)
				{
					new_matrix[y - start_index][x - start_index] =
						base_array[y][x];
				}

			return new_matrix;
		}

		__DIXELU_CONDITIONAL_SPECIFIERS self_type& transpose()
		{
			for (std::size_t y = 0; y < dims; ++y)
				for (std::size_t x = 0; x < dims; ++x)
				{
					auto t = base_array[y][x];
					base_array[y][x] = base_array[x][y];
					base_array[x][y] = t;
				}
			return *this;
		}

		__DIXELU_CONDITIONAL_SPECIFIERS general_float_type trace() const
		{
			general_float_type sum(0);
			for (std::size_t i = 0; i < dims; ++i)
				sum += base_array[i][i];
			return sum;
		}

		__DIXELU_CONDITIONAL_SPECIFIERS self_type ppow(double p) const
		{
			self_type mx(*this);
			mx.selfppow(p);
			return mx;
		}

		__DIXELU_CONDITIONAL_SPECIFIERS self_type& selfppow(double p)
		{
			for (std::size_t y = 0; y < dims; ++y)
				for (std::size_t x = 0; x < dims; ++x)
					base_array[y][x] = utils::constexpr_pow(base_array[y][x]);
			return *this;
		}

		__DIXELU_CONDITIONAL_SPECIFIERS self_type pabs() const
		{
			self_type mx(*this);
			mx.selfpabs();
			return mx;
		}

		__DIXELU_CONDITIONAL_SPECIFIERS self_type& selfpabs()
		{
			for (std::size_t y = 0; y < dims; ++y)
				for (std::size_t x = 0; x < dims; ++x)
					base_array[y][x] = utils::constexpr_abs(base_array[y][x]);
			return *this;
		}

		__DIXELU_CONDITIONAL_SPECIFIERS general_float_type psum() const
		{
			general_float_type s = general_float_type(0);
			for (std::size_t y = 0; y < dims; ++y)
				for (std::size_t x = 0; x < dims; ++x)
					s += base_array[y][x];
			return s;
		}

		__DIXELU_CONDITIONAL_SPECIFIERS self_type& selfapply(const std::function<void(general_float_type&)>& func)
		{
			for (std::size_t y = 0; y < dims; ++y)
				for (std::size_t x = 0; x < dims; ++x)
					func(base_array[y][x]);
			return *this;
		}

		__DIXELU_CONDITIONAL_SPECIFIERS self_type apply(const std::function<void(general_float_type&)>& func) const
		{
			self_type mx(*this);
			mx.selfapply(func);
			return mx;
		}

		__DIXELU_CONDITIONAL_SPECIFIERS self_type& selfapply_indexed(const std::function<void(general_float_type&, size_t, size_t)>& func)
		{
			for (std::size_t y = 0; y < dims; ++y)
				for (std::size_t x = 0; x < dims; ++x)
					func(base_array[y][x], y, x);
			return *this;
		}

		__DIXELU_CONDITIONAL_SPECIFIERS self_type apply_indexed(const std::function<void(general_float_type&, size_t, size_t)>& func) const
		{
			self_type mx(*this);
			mx.selfapply_indexed(func);
			return mx;
		}
	};

	template<typename general_float_type, std::size_t dims>
	__DIXELU_CONDITIONAL_SPECIFIERS point<general_float_type, dims>
		cross_prod(const std::array<point<general_float_type, dims>, dims - 1>& points)
	{
		if (dims > 1)
		{
			point<general_float_type, dims> answer;
			std::array<point<general_float_type, dims>, dims> mx;
			mx[0] = point<general_float_type, dims>(std::vector<general_float_type>(dims, 0));
			std::copy(points.begin(), points.end(), mx.begin() + 1);
			sq_matrix<general_float_type, dims> M(mx);
			for (int i = 0; i < dims; ++i)
				answer[i] = M.minor_matrix(i, 0).determinant() * ((i & 1) ? (1.) : (-1.));
			return answer;
		}
		else if (dims == 1)
			return { 0 };
		else
			return {};
	}

	template<size_t dims, typename general_float_type>
	__DIXELU_CONDITIONAL_SPECIFIERS std::ostream& operator<<(std::ostream& in,
		const sq_matrix<general_float_type, dims>& M)
	{
		for (size_t y = 0; y < dims; ++y)
		{
			for (size_t x = 0; x < dims; ++x)
			{
				in << std::setfill(' ') << std::setw(15) << M.at(y, x) << " ";
			}
			in << "\n";
		}
		return in;
	}
} // namespace dixelu

#endif // _DIXELU_SQ_MATRIX_H_
