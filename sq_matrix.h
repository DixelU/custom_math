#ifndef _MULTIDIMENTIONAL_POINT_H_
#define _MULTIDIMENTIONAL_POINT_H_

#include <cstdio>
#include <array>
#include <vector>
#include <cmath>
#include <ostream>
#include <iomanip>
#include <limits>   
#include <initializer_list>

#define __DEFAULT_MDP_SPECIFIERS inline

#if __cplusplus >= 201402L
#define __MDP_COND_CONSTEXPR constexpr
#else
#define __MDP_COND_CONSTEXPR
#endif

#define __MDP_CONDITIONAL_SPECIFIERS __DEFAULT_MDP_SPECIFIERS __MDP_COND_CONSTEXPR

namespace dixelu
{

	namespace utils
	{

		template<typename T>
		__MDP_CONDITIONAL_SPECIFIERS T constexpr_abs(const T& value)
		{
			return value < 0 ? -value : value;
		}

		namespace details
		{

#ifndef WITHOUT_CONSTEXPR_FUNCTIONS
			template<typename T>
			__MDP_CONDITIONAL_SPECIFIERS T __sqr(T x) {
				return x * x;
			}

			template<typename T>
			__MDP_CONDITIONAL_SPECIFIERS T __uintpow(T x, std::size_t n)
			{
				return n == 0 ? 1 : __sqr(__uintpow(x, n >> 1)) * (n & 1 ? x : T(1));
			}

			template<typename T, unsigned int n>
			struct lookup_table
			{
				static const __MDP_COND_CONSTEXPR unsigned int max_degree = T(std::numeric_limits<T>::max_exponent / n);
				T lookup_table_vals[max_degree]{};
				T lookup_table_roots[max_degree]{};
				__MDP_COND_CONSTEXPR lookup_table() :
					lookup_table_vals(), lookup_table_roots()
				{
					__MDP_COND_CONSTEXPR T base(2);
					__MDP_COND_CONSTEXPR T n_T(n);
					__MDP_COND_CONSTEXPR T zero(0);
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
			__MDP_CONDITIONAL_SPECIFIERS T root_approx(T x)
			{
				__MDP_COND_CONSTEXPR lookup_table<T, n> table{};
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
			__MDP_CONDITIONAL_SPECIFIERS T __uintroot(T x)
			{
				__MDP_COND_CONSTEXPR T epsilon = std::numeric_limits<T>::epsilon();
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
				do
				{
					x_k_nm1 = one / __uintpow(x_k, n - 1);
					x_p_prev_k = x_prev_k;
					x_prev_k = x_k;
					x_k = coef * x_k + (x / n_conv) * x_k_nm1;
					diff = constexpr_abs(x_k - x_prev_k);
					diff_step = constexpr_abs(x_k - x_p_prev_k);
				} while (diff > epsilon && diff_step > epsilon);
				return x_k;
			}

			template<typename T>
			__MDP_CONDITIONAL_SPECIFIERS T __positive_pow(T x, T p);

			template<typename T>
			__MDP_CONDITIONAL_SPECIFIERS T __frac_positive_pow(T x, T p)
			{
				constexpr unsigned int rolling_ppow_bits = 4u;
				constexpr unsigned int rolling_power = ((1u << rolling_ppow_bits));
				__MDP_COND_CONSTEXPR T epsilon = std::numeric_limits<T>::epsilon();
				if (p <= epsilon)
					return T(1);
				if (p >= T(1))
					return __positive_pow(x, p);
				return __uintroot<T, rolling_power>(__positive_pow(x, p * rolling_power));
			}

			template<typename T>
			__MDP_CONDITIONAL_SPECIFIERS T __positive_pow(T x, T p)
			{
				std::size_t p_w(p);
				auto d_p = p - T(p_w);
				auto uintpow = __uintpow(x, p_w);
				auto fracpow = __frac_positive_pow(x, d_p);
				return uintpow * fracpow;
			}

			template<typename T>
			__MDP_CONDITIONAL_SPECIFIERS T __pow(T a, T b)
			{
				if (b < 0)
					return T(1) / __positive_pow(a, -b);
				return __positive_pow(a, b);
			}

#endif
		} // namespace details

		template<typename T>
		__MDP_CONDITIONAL_SPECIFIERS T constexpr_sqrt(T x)
		{
#ifdef WITHOUT_CONSTEXPR_FUNCTIONS
			return std::sqrt(x);
#else
			return details::__uintroot<T, 2>(x);
#endif // WITHOUT_CONSTEXPR_FUNCTIONS
		}

		template<typename T>
		__MDP_CONDITIONAL_SPECIFIERS T constexpr_pow(T a, T b)
		{
#ifdef WITHOUT_CONSTEXPR_FUNCTIONS
			return std::pow(a, b);
#else
			return details::__pow(a, b);
#endif // WITHOUT_CONSTEXPR_FUNCTIONS
		}

	} // namespace utils

	template<typename general_float_type, size_t dims>
	struct point
	{
		general_float_type base_array[dims];
		using self_type = point<general_float_type, dims>;
		__MDP_COND_CONSTEXPR point() :
			base_array()
		{
			for (int i = 0; i < dims; i++)
				base_array[i] = general_float_type();
		}
		__MDP_COND_CONSTEXPR explicit point(general_float_type v) :
			base_array()
		{
			for (int i = 0; i < dims; i++)
				base_array[i] = v;
		}
		__MDP_COND_CONSTEXPR point(const std::initializer_list<general_float_type>& il_d) :
			base_array()
		{
			auto y = il_d.begin();
			for (size_t i = 0; i < dims && y != il_d.end(); i++, y++)
				base_array[i] = *y;
		}
		__MDP_COND_CONSTEXPR point(const std::initializer_list<int>& il) :
			base_array()
		{
			auto y = il.begin();
			for (size_t i = 0; i < dims && y != il.end(); i++, y++)
				base_array[i] = *y;
		}
		__MDP_COND_CONSTEXPR point(const std::vector<general_float_type>& il_d) :
			base_array()
		{
			auto y = il_d.cbegin();
			for (size_t i = 0; i < dims && y != il_d.end(); i++, y++)
				base_array[i] = *y;
		}
		__MDP_COND_CONSTEXPR point(const std::vector<int>& il) :
			base_array()
		{
			auto y = il.cbegin();
			for (size_t i = 0; i < dims && y != il.end(); i++, y++)
				base_array[i] = *y;
		}
		__MDP_CONDITIONAL_SPECIFIERS general_float_type get_norm2() const
		{
			general_float_type sum = general_float_type();
			for (int i = 0; i < dims; i++)
				sum += base_array[i] * base_array[i];
			return sum;
		}
		__MDP_CONDITIONAL_SPECIFIERS general_float_type get_norm(general_float_type deg) const
		{
			general_float_type sum = general_float_type();
			for (int i = 0; i < dims; i++)
				sum += utils::constexpr_pow<general_float_type>(base_array[i], deg);
			return utils::constexpr_pow<general_float_type>(sum, 1. / deg);
		}
		__MDP_CONDITIONAL_SPECIFIERS general_float_type get_norm() const
		{
			return utils::constexpr_sqrt(get_norm2());
		}
		__MDP_CONDITIONAL_SPECIFIERS size_t get_dims() const
		{
			return dims;
		}
		__MDP_CONDITIONAL_SPECIFIERS general_float_type& operator[](size_t D)
		{
			return base_array[D];
		}
		__MDP_CONDITIONAL_SPECIFIERS const general_float_type& operator[](size_t D) const
		{
			return base_array[D];
		}
		__MDP_CONDITIONAL_SPECIFIERS void swap(self_type& P)
		{
			for (size_t i = 0; i < dims; i++)
			{
				auto value = base_array[i];
				base_array[i] = P.base_array[i];
				P.base_array[i] = value;
			}
		}
		__MDP_CONDITIONAL_SPECIFIERS self_type operator+(const self_type& P) const
		{
			self_type N;
			for (size_t i = 0; i < dims; i++)
				N[i] = base_array[i] + P[i];
			return N;
		}
		__MDP_CONDITIONAL_SPECIFIERS self_type operator+=(const self_type& P)
		{
			for (size_t i = 0; i < dims; i++)
				base_array[i] += P[i];
			return *this;
		}
		__MDP_CONDITIONAL_SPECIFIERS self_type operator-(const self_type& P) const
		{
			self_type N;
			for (size_t i = 0; i < dims; i++)
				N[i] = base_array[i] - P[i];
			return N;
		}
		__MDP_CONDITIONAL_SPECIFIERS self_type operator-=(const self_type& P)
		{
			for (size_t i = 0; i < dims; i++)
				base_array[i] -= P[i];
			return *this;
		}
		__MDP_CONDITIONAL_SPECIFIERS self_type operator*(general_float_type M) const
		{
			self_type N;
			for (size_t i = 0; i < dims; i++)
				N[i] = base_array[i] * M;
			return N;
		}
		__MDP_CONDITIONAL_SPECIFIERS self_type operator*=(general_float_type M)
		{
			for (size_t i = 0; i < dims; i++)
				base_array[i] *= M;
			return *this;
		}
		__MDP_CONDITIONAL_SPECIFIERS self_type operator/(general_float_type M) const
		{
			self_type N;
			for (size_t i = 0; i < dims; i++)
				N[i] = base_array[i] / M;
			return N;
		}
		__MDP_CONDITIONAL_SPECIFIERS self_type operator/=(general_float_type M)
		{
			return ((*this) *= (1. / M));
		}
		__MDP_CONDITIONAL_SPECIFIERS general_float_type operator*(const self_type& P) const
		{
			general_float_type sum = general_float_type();
			for (size_t i = 0; i < dims; i++)
				sum += base_array[i] * P[i];
			return sum;
		}
		__MDP_CONDITIONAL_SPECIFIERS self_type operator-()
		{
			self_type N;
			for (size_t i = 0; i < dims; i++)
				N[i] = 0 - base_array[i];
			return *this;
		}
		__MDP_CONDITIONAL_SPECIFIERS bool operator<(const self_type& P) const
		{
			for (size_t i = 0; i < dims; i++)
				if (base_array[i] >= P.base_array[i])
					return false;
			return true;
		}
		__MDP_CONDITIONAL_SPECIFIERS bool operator>=(const self_type& P) const
		{
			return !(*this < P);
		}
		__MDP_CONDITIONAL_SPECIFIERS bool operator>(const self_type& P) const
		{
			return P < *this;
		}
		__MDP_CONDITIONAL_SPECIFIERS bool operator==(const self_type& P) const
		{
			for (size_t i = 0; i < dims; i++)
				if (base_array[i] != P[i]) return false;
			return true;
		}
		__MDP_CONDITIONAL_SPECIFIERS bool operator!=(const self_type& P) const
		{
			return !(*this == P);
		}
		__MDP_CONDITIONAL_SPECIFIERS bool operator<=(const self_type& P) const
		{
			return !(P < *this);
		}
		__MDP_CONDITIONAL_SPECIFIERS self_type normalize() const
		{
			return (*this) / get_norm();
		}
	};

	template<typename general_float_type, size_t dims>
	std::ostream& operator<<(std::ostream& os, const point<general_float_type, dims>& P)
	{
		os << "(";
		for (size_t i = 0; i < dims; i++)
		{
			os << P[i];
			if (i != dims - 1)
				os << ",";
		}
		os << ")";
		return os;
	}

	template<typename general_float_type, size_t dims>
	__MDP_CONDITIONAL_SPECIFIERS point<general_float_type, dims> operator*(general_float_type M, point<general_float_type, dims> P)
	{
		point<general_float_type, dims> N;
		for (size_t i = 0; i < dims; i++)
			N[i] = P[i] * M;
		return N;
	}

	template<typename general_float_type, size_t dims>
	struct sq_matrix
	{
		static constexpr general_float_type GFLOAT_EPSILON = std::numeric_limits<general_float_type>::epsilon();
		general_float_type utilisation = general_float_type();

		using point_type = point<general_float_type, dims>;
		using self_type = sq_matrix<general_float_type, dims>;

		point_type base_array[dims];

		__MDP_COND_CONSTEXPR sq_matrix() :
			base_array()
		{
			for (size_t i = 0; i < dims; i++)
				base_array[i] = point_type();
		}
		__MDP_COND_CONSTEXPR explicit sq_matrix(general_float_type E_num) :
			base_array()
		{
			for (size_t i = 0; i < dims; i++)
			{
				base_array[i] = point_type();
				base_array[i][i] = E_num;
			}
		}
		__MDP_COND_CONSTEXPR sq_matrix(std::initializer_list<point_type> IL) :
			base_array()
		{
			size_t id = 0;
			for (auto& p : IL)
			{
				if (id == dims)
					break;
				base_array[id] = p;
				id++;
			}
		}
		void swap(self_type& p)
		{
			for (int i = 0; i < dims; i++)
				p[i].swap(base_array[i]);
		}
		__MDP_CONDITIONAL_SPECIFIERS point_type operator*(const point_type& p) const
		{
			point_type T;
			for (int i = 0; i < dims; i++)
				T[i] = base_array[i] * p;
			return T;
		}
		__MDP_CONDITIONAL_SPECIFIERS self_type operator*(general_float_type num) const
		{
			self_type T;
			for (int i = 0; i < dims; i++)
				T[i] = base_array[i] * num;
			return T;
		}
		__MDP_CONDITIONAL_SPECIFIERS self_type operator*=(general_float_type num)
		{
			for (int i = 0; i < dims; i++)
				base_array[i] *= num;
			return *this;
		}
		__MDP_CONDITIONAL_SPECIFIERS self_type operator/(general_float_type num) const
		{
			return *this * (1. / num);
		}
		__MDP_CONDITIONAL_SPECIFIERS self_type operator/=(general_float_type num)
		{
			return ((*this) *= (1. / num));
		}
		__MDP_CONDITIONAL_SPECIFIERS self_type operator+(const self_type& p) const {
			self_type T;
			for (int i = 0; i < dims; i++)
				T[i] = base_array[i] + p[i];
			return T;
		}
		__MDP_CONDITIONAL_SPECIFIERS self_type operator-(const self_type& p) const
		{
			self_type T;
			for (int i = 0; i < dims; i++)
				T[i] = base_array[i] - p[i];
			return T;
		}
		__MDP_CONDITIONAL_SPECIFIERS self_type operator+=(const self_type& p)
		{
			for (int i = 0; i < dims; i++)
				base_array[i] += p[i];
			return *this;
		}
		__MDP_CONDITIONAL_SPECIFIERS self_type operator-=(const self_type& p)
		{
			for (int i = 0; i < dims; i++) {
				base_array[i] -= p[i];
			}
			return *this;
		}
		__MDP_CONDITIONAL_SPECIFIERS const point_type& operator[](size_t i) const
		{
			return base_array[i];
		}
		__MDP_CONDITIONAL_SPECIFIERS point_type& operator[](size_t i)
		{
			return base_array[i];
		}
		__MDP_CONDITIONAL_SPECIFIERS general_float_type& at(size_t point_id, size_t coordinate)
		{
			if (point_id < dims && coordinate < dims) {
				return base_array[point_id][coordinate];
			}
			else return utilisation;
		}
		__MDP_CONDITIONAL_SPECIFIERS const general_float_type& at(size_t point_id, size_t coordinate) const
		{
			if (point_id < dims && coordinate < dims) {
				return base_array[point_id][coordinate];
			}
			else return utilisation;
		}
		__MDP_CONDITIONAL_SPECIFIERS self_type operator*(const self_type& M) const
		{
			self_type P;
			for (size_t y = 0; y < dims; y++) {
				for (size_t x = 0; x < dims; x++) {
					for (size_t i = 0; i < dims; i++) {
						P[y][x] += base_array[y][i] * M[i][x];
					}
				}
			}
			return P;
		}
		__MDP_CONDITIONAL_SPECIFIERS self_type inverse() const
		{
			general_float_type max_value = general_float_type();
			general_float_type mul = general_float_type();
			size_t id = 0;
			self_type E(1), A(*this);
			for (size_t step = 0; step < dims; step++)
			{
				id = step;
				max_value = 0;
				for (size_t coid = step; coid < dims; coid++)
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
				for (size_t sum_id = 0; sum_id < dims; sum_id++)
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
		__MDP_CONDITIONAL_SPECIFIERS general_float_type determinant() const
		{
			general_float_type determ = 1;
			general_float_type temp = 0, max_value = 0, mul = 0;
			size_t id = 0;
			self_type A(*this);
			for (size_t step = 0; step < dims; step++)
			{
				id = step;
				max_value = 0;
				for (size_t coid = step; coid < dims; coid++)
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
				for (size_t sum_id = 0; sum_id < dims; sum_id++)
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
		__MDP_CONDITIONAL_SPECIFIERS static point_type solve_using_eulers_method(self_type A, point_type P)
		{
			general_float_type max_value = 0, mul = 0;
			size_t id = 0;
			for (size_t step = 0; step < dims; step++)
			{
				id = step;
				max_value = 0;
				for (size_t coid = step; coid < dims; coid++)
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
				for (size_t sum_id = 0; sum_id < dims; sum_id++)
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
		__MDP_CONDITIONAL_SPECIFIERS self_type operator^(int degree)
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
		__MDP_CONDITIONAL_SPECIFIERS self_type& operator^=(int degree)
		{
			return ((*this) = (*this) ^ degree), * this;
		}
		__MDP_CONDITIONAL_SPECIFIERS sq_matrix<general_float_type, dims - 1>
			minor_matrix(const size_t& x_minor, const size_t& y_minor) const
		{
			auto minor_index = [](size_t x, size_t y, size_t minor_x, size_t minor_y) -> std::pair<int64_t, int64_t>
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
			for (size_t y = 0; y < dims; y++)
			{
				for (size_t x = 0; x < dims; x++)
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
	};

	template<typename general_float_type, size_t dims>
	__MDP_CONDITIONAL_SPECIFIERS point<general_float_type, dims>
		cross_prod(const std::array<point<general_float_type, dims>, dims - 1>& points)
	{
		if (dims > 1)
		{
			point<general_float_type, dims> answer;
			std::array<point<general_float_type, dims>, dims> mx;
			mx[0] = point<general_float_type, dims>(std::vector<general_float_type>(dims, 0));
			std::copy(points.begin(), points.end(), mx.begin() + 1);
			sq_matrix<general_float_type, dims> M(mx);
			for (int i = 0; i < dims; i++)
				answer[i] = M.minor_matrix(i, 0).determinant() * ((i & 1) ? (1.) : (-1.));
			return answer;
		}
		else if (dims == 1)
			return { 0 };
		else
			return {};
	}

	template<size_t dims, typename general_float_type>
	__MDP_CONDITIONAL_SPECIFIERS std::ostream& operator<<(std::ostream& in,
		const sq_matrix<general_float_type, dims>& M)
	{
		for (size_t y = 0; y < dims; y++)
		{
			for (size_t x = 0; x < dims; x++)
			{
				in << std::setfill(' ') << std::setw(15) << M.at(y, x) << " ";
			}
			in << "\n";
		}
		return in;
	}

} // namespace dixelu

#endif