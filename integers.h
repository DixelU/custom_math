#ifndef _DIXELU_INTEGERS_H_
#define _DIXELU_INTEGERS_H_

#if (defined(__cpp_constexpr) && (__cpp_constexpr >= 201304L))

#ifndef __DIXELU_RELAXED_CONSTEXPR
#define __DIXELU_RELAXED_CONSTEXPR constexpr
#endif
#else
#ifndef __DIXELU_RELAXED_CONSTEXPR
#define __DIXELU_RELAXED_CONSTEXPR
#endif
#endif

#if (defined(__cpp_constexpr) && (__cpp_constexpr >= 200704L))
#ifndef __DIXELU_STRICT_CONSTEXPR
#define __DIXELU_STRICT_CONSTEXPR constexpr
#endif
#else 
#ifndef __DIXELU_STRICT_CONSTEXPR
#define __DIXELU_STRICT_CONSTEXPR
#endif
#endif

#ifndef __DIXELU_CONDITIONAL_CPP14_SPECIFIERS
#define __DIXELU_CONDITIONAL_CPP14_SPECIFIERS inline __DIXELU_RELAXED_CONSTEXPR
#endif

#include <limits.h>
#include <limits>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <type_traits>
#include <algorithm>
#include <cinttypes>
#include <string>

namespace dixelu
{
	namespace details
	{
		template<std::uint64_t deg, typename base_type>
		struct tree_bits
		{
			static constexpr base_type base_bits = sizeof(base_type) * CHAR_BIT;
			static constexpr base_type value = (tree_bits<deg - 1, base_type>::value << 1);
		};

		template<typename base_type>
		struct tree_bits<0, base_type>
		{
			static constexpr base_type base_bits = sizeof(base_type) * CHAR_BIT;
			static constexpr base_type value = base_bits * 2;
		};
	}

	/* long_uint<0> ~ 128 bit unsigned integer */
	template <std::uint64_t deg>
	struct long_uint
	{
		using base_type = std::uint64_t; // largest of the integer types
		using size_type = std::uint64_t;
		using self_type = long_uint<deg>;
		using down_type = typename std::conditional<deg == 0, base_type, long_uint<deg - 1>>::type;

		static constexpr size_type base_bits = details::tree_bits<0, base_type>::base_bits;
		static constexpr size_type bits = details::tree_bits<deg, base_type>::value;
		static constexpr size_type down_type_bits = bits >> 1;
		static constexpr size_type size = down_type_bits / base_bits;

		struct __fill_fields_tag {};

		down_type hi;
		down_type lo;

		__DIXELU_STRICT_CONSTEXPR
			long_uint() :
			hi(), lo()
		{ }

		__DIXELU_STRICT_CONSTEXPR
			long_uint(base_type value) :
			hi(), lo(value)
		{ }

		template<uint64_t __deg = deg>
		__DIXELU_STRICT_CONSTEXPR
			explicit long_uint(base_type value,
				const __fill_fields_tag& fill_fields_tag,
				typename std::enable_if<(__deg > 0), void>::type* = 0) :
			hi(value, typename long_uint<__deg - 1>::__fill_fields_tag{}),
			lo(value, typename long_uint<__deg - 1>::__fill_fields_tag{})
		{ }

		template<uint64_t __deg = deg>
		__DIXELU_STRICT_CONSTEXPR
			explicit long_uint(base_type value,
				const __fill_fields_tag& fill_fields_tag,
				typename std::enable_if<(__deg == 0), void>::type* = 0) :
			hi(value), lo(value)
		{ }

		template<uint64_t __deg>
		__DIXELU_RELAXED_CONSTEXPR
			explicit long_uint(const long_uint<__deg>& value,
				typename std::enable_if<(__deg < deg), void>::type* = 0) : hi(), lo((down_type)value) { }

		template<uint64_t __deg>
		__DIXELU_RELAXED_CONSTEXPR
			explicit long_uint(const long_uint<__deg>& value,
				typename std::enable_if<(__deg > deg), void>::type* = 0) :
			hi(value.lo.hi),
			lo(value.lo.lo)
		{ }

		template<uint64_t __deg>
		__DIXELU_RELAXED_CONSTEXPR
			explicit long_uint(const long_uint<__deg>& value,
				typename std::enable_if<(__deg == deg), void>::type* = 0) : hi(value.hi), lo(value.lo) { }

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			self_type& operator|=(const self_type& rhs)
		{
			hi |= rhs.hi;
			lo |= rhs.lo;
			return *this;
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			self_type& operator&=(const self_type& rhs)
		{
			hi &= rhs.hi;
			lo &= rhs.lo;
			return *this;
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			self_type& operator^=(const self_type& rhs)
		{
			hi ^= rhs.hi;
			lo ^= rhs.lo;
			return *this;
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			self_type operator~() const
		{
			long_uint res;
			res.hi = ~hi;
			res.lo = ~lo;
			return res;
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			self_type operator|(const self_type& rhs) const
		{
			long_uint res(*this);
			res |= rhs;
			return res;
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			self_type operator&(const self_type& rhs) const
		{
			long_uint res(*this);
			res &= rhs;
			return res;
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			self_type operator^(const self_type& rhs) const
		{
			long_uint res(*this);
			res ^= rhs;
			return res;
		}

		template< bool cond, typename U >
		using resolved_return_type = typename std::enable_if<cond, U>::type;

		template<std::uint64_t __deg = deg>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			resolved_return_type<(__deg > 0), bool> get_bit(size_type bit) const
		{
			if (bit < down_type_bits)
				return lo.get_bit(bit);
			bit -= down_type_bits;
			return hi.get_bit(bit);
		}

		template<std::uint64_t __deg = deg>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			resolved_return_type<(__deg == 0), bool> get_bit(size_type bit) const
		{
			if (bit < down_type_bits)
				return (lo >> bit) & 1;
			bit -= down_type_bits;
			return (hi >> bit) & 1;
		}

		template<std::uint64_t __deg = deg>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			resolved_return_type<(__deg > 0), base_type&> operator[](size_type idx)
		{
			if (idx < down_type_bits / base_bits)
				return lo[idx];
			idx -= down_type_bits / base_bits;
			return hi[idx];
		}

		template<std::uint64_t __deg = deg>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			resolved_return_type<(__deg == 0), base_type&> operator[](size_type idx)
		{
			if (idx < size)
				return lo;
			return hi;
		}

		template<std::uint64_t __deg = deg>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			resolved_return_type<(__deg > 0), const base_type&> operator[](size_type idx) const
		{
			if (idx < size)
				return lo[idx];
			idx -= size;
			return hi[idx];
		}

		template<std::uint64_t __deg = deg>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			resolved_return_type<(__deg == 0), const base_type&> operator[](size_type idx) const
		{
			if (idx < size)
				return lo;
			return hi;
		}

		template<std::uint64_t __deg = deg>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			resolved_return_type<(__deg > 0), void>
			set_bit(size_type bit, bool bit_value)
		{
			if (bit < down_type_bits)
				return lo.set_bit(bit, bit_value);
			bit -= down_type_bits;
			return hi.set_bit(bit, bit_value);
		}

		template<std::uint64_t __deg = deg>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			resolved_return_type<(__deg == 0), void>
			set_bit(size_type bit, bool bit_value)
		{
			if (bit < down_type_bits)
			{
				const auto mask = ~(base_type(1) << bit);
				const auto value = (base_type(bit_value) << bit);
				lo &= mask;
				lo |= value;
				return;
			}
			bit -= down_type_bits;
			const auto mask = ~(base_type(1) << bit);
			const auto value = (base_type(bit_value) << bit);
			hi &= mask;
			hi |= value;
			return;
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			self_type& operator+=(const self_type& rhs)
		{
			const bool lhs_lo_carried_bit = get_bit(down_type_bits - 1);
			const bool rhs_lo_carried_bit = rhs.get_bit(down_type_bits - 1);
			const bool carry_possibility = lhs_lo_carried_bit ^ rhs_lo_carried_bit;
			const bool definite_carry = lhs_lo_carried_bit & rhs_lo_carried_bit;

			hi += rhs.hi;
			lo += rhs.lo;

			if (definite_carry)
				hi += down_type(1);
			else if (carry_possibility)
			{
				const auto carried_bit = get_bit(down_type_bits - 1);
				if (!carried_bit)
					hi += down_type(1);
			}

			return *this;
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			self_type operator+(const self_type& rhs) const
		{
			self_type res(*this);
			res += rhs;
			return res;
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			self_type& operator-=(const self_type& rhs)
		{
			return ((*this) += ((~rhs) + self_type(1)));
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			self_type operator-(const self_type& rhs) const
		{
			self_type res(*this);
			res -= rhs;
			return res;
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			explicit operator size_type() const
		{
			return operator[](0);
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			explicit operator bool() const
		{
			return ((bool)lo | (bool)hi);
		}

		///* https://github.com/glitchub/arith64/blob/master/arith64.c *///
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			self_type& operator<<=(size_type rhs)
		{
			rhs &= (bits - 1);

			if (rhs >= down_type_bits)
			{
				hi = (lo <<= (rhs - down_type_bits));
				lo = 0;
			}
			else if (rhs)
			{
				auto lo_copy = lo;
				auto hi_copy = hi;
				hi = (lo_copy >>= (down_type_bits - rhs)) | (hi_copy <<= rhs);
				lo <<= rhs;
			}

			return *this;
		}

		///* https://github.com/glitchub/arith64/blob/master/arith64.c *///
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			self_type& operator>>=(size_type rhs)
		{
			rhs &= (bits - 1);

			if (rhs >= down_type_bits)
			{
				lo = (hi >>= (rhs - down_type_bits));
				hi = 0;
			}
			else if (rhs)
			{
				auto lo_copy = lo;
				auto hi_copy = hi;
				lo = (hi_copy <<= (down_type_bits - rhs)) | (lo_copy >>= rhs);
				hi >>= rhs;
			}

			return *this;
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			self_type operator<<(size_type rhs) const
		{
			self_type res(*this);
			res <<= rhs;
			return res;
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			self_type operator>>(size_type rhs) const
		{
			self_type res(*this);
			res >>= rhs;
			return res;
		}

		template<std::uint64_t __deg>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			explicit operator resolved_return_type<(__deg == deg || !__deg), long_uint<__deg>>() const
		{
			return *this;
		}

		template<std::uint64_t __deg>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			explicit operator resolved_return_type<(__deg < deg && __deg), long_uint<__deg>>() const
		{
			return (long_uint<__deg>)(lo);
		}

		template<std::uint64_t __deg>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			explicit operator resolved_return_type<(__deg > deg&& __deg), long_uint<__deg>>() const
		{
			long_uint<deg + 1> res;
			res.lo = *this;
			return (long_uint<__deg>)(res);
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			self_type operator*(const self_type& rhs) const
		{
			auto up_result = long_uint<deg + 1>::__downtype_mul_long(
				(long_uint<deg + 1>)(*this),
				(long_uint<deg + 1>)rhs
			);
			return up_result.lo;
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			self_type& operator*=(const self_type& rhs)
		{
			return (*this = (*this * rhs));
		}

		template<typename T>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			static T __downtype_mul_call(const T& lhs, const T& rhs,
				typename std::enable_if<(std::is_same<T, typename long_uint<deg>::base_type>::value), void>::type* = 0)
		{
			return lhs * rhs;
		}

		template<std::uint64_t __deg>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			static long_uint<__deg> __downtype_mul_call(const long_uint<__deg>& lhs, const long_uint<__deg>& rhs)
		{
			return __downtype_mul<__deg>(lhs, rhs);
		}



		template<std::uint64_t __deg>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			static typename long_uint<__deg>::down_type
			__downtype_get_lo(const typename long_uint<__deg>::down_type& value,
				typename std::enable_if<(__deg == 0), void>::type* = 0)
		{
			return (value & 0x00000000ffffffffULL);
		}

		template<std::uint64_t __deg>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			static typename long_uint<__deg>::down_type
			__downtype_get_lo(const typename long_uint<__deg>::down_type& value,
				typename std::enable_if<(__deg > 0), void>::type* = 0)
		{
			typename long_uint<__deg>::down_type res;
			res.lo = value.lo;
			return res;
		}

		template<std::uint64_t __deg>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			static typename long_uint<__deg>::down_type
			__downtype_get_hi(const typename long_uint<__deg>::down_type& value,
				typename std::enable_if<(__deg == 0), void>::type* = 0)
		{
			return (value & 0xffffffff00000000ULL) >> 32;
		}

		template<std::uint64_t __deg>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			static typename long_uint<__deg>::down_type
			__downtype_get_hi(const typename long_uint<__deg>::down_type& value,
				typename std::enable_if<(__deg > 0), void>::type* = 0)
		{
			typename long_uint<__deg>::down_type res;
			res.lo = value.hi;
			return res;
		}

		template<std::uint64_t __deg>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			static long_uint<__deg> __downtype_mul_long(const long_uint<__deg>& lhs, const long_uint<__deg>& rhs)
		{
			long_uint<__deg> carry(1);
			long_uint<__deg> result;
			long_uint<__deg> rolling_rhs = rhs;
			long_uint<__deg> diminishing_lhs = lhs;
			while(diminishing_lhs)
			{
				if(diminishing_lhs & carry)
				{
					result += rolling_rhs;
					diminishing_lhs &= ~(carry);
				}
				rolling_rhs <<= 1;
				carry <<= 1;
			}
			return result;
		}

		template<std::uint64_t __deg>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
		// todo: fix
			static long_uint<__deg> __downtype_mul(const long_uint<__deg>& lhs, const long_uint<__deg>& rhs)
		{
			using local_down_type = typename long_uint<__deg>::down_type;
			constexpr size_type half_down_type_mask_bits = (long_uint<__deg>::down_type_bits >> 1);
			const local_down_type himask = (local_down_type(1) << half_down_type_mask_bits - 1);

			const local_down_type lhs_lo = long_uint<__deg>::template __downtype_get_lo<__deg>(lhs.lo);
			const local_down_type rhs_lo = long_uint<__deg>::template __downtype_get_lo<__deg>(rhs.lo);
			const local_down_type lhs_hi = long_uint<__deg>::template __downtype_get_hi<__deg>(rhs.lo);
			const local_down_type rhs_hi = long_uint<__deg>::template __downtype_get_hi<__deg>(rhs.lo);

			const local_down_type hihi = __downtype_mul_call(lhs_hi, rhs_hi);
			const local_down_type lolo = __downtype_mul_call(lhs_lo, rhs_lo);

			bool possible_overflow =
				((rhs_lo |
					lhs_lo |
					rhs_hi |
					lhs_hi) & himask) != 0;

			long_uint<__deg> hilo{};

			if (possible_overflow)
			{
				/*std::cout << "PO " << to_hex_string(lhs) << " " << to_hex_string(rhs) << std::endl;

				bool lhs_diff_positive = true;
				bool rhs_diff_positive = true;
				const local_down_type lhs_abs_difference = (lhs_diff_positive = (lhs_hi > lhs_lo)) ? lhs_hi - lhs_lo : lhs_lo - lhs_hi;
				const local_down_type rhs_abs_difference = (rhs_diff_positive = (rhs_hi > rhs_lo)) ? rhs_hi - rhs_lo : rhs_lo - rhs_hi;
				const local_down_type hihilolo = __downtype_mul_call(lhs_abs_difference, rhs_abs_difference);

				if (lhs_diff_positive != rhs_diff_positive)
				{
					hilo = ((long_uint<__deg>)hihi - (long_uint<__deg>)hihilolo + (long_uint<__deg>)lolo);
					std::cout << "eq " << 
						to_hex_string(hihi) << " " << 
						to_hex_string(lolo) << " " << 
						to_hex_string(hihilolo) << " " << 
						to_hex_string(hilo) << std::endl;
				}
				else
				{
					hilo = ((long_uint<__deg>)hihilolo + (long_uint<__deg>)hihi + (long_uint<__deg>)lolo);
					std::cout << "nq " << 
						to_hex_string(hihi) << " " << 
						to_hex_string(lolo) << " " << 
						to_hex_string(hihilolo) << " " << 
						to_hex_string(hilo) << std::endl;
				}*/
			}
			else
			{
				const local_down_type rhs_sum = rhs_lo + rhs_lo;
				const local_down_type lhs_sum = lhs_lo + lhs_lo;
				const local_down_type hihilolo = __downtype_mul_call(rhs_sum, lhs_sum);
				hilo = (long_uint<__deg>)(hihilolo - hihi - lolo);
			}

			auto res_hi = (hilo >> half_down_type_mask_bits) << (2 * half_down_type_mask_bits);
			auto res_lo = (hilo << (3 * half_down_type_mask_bits)) >> (3 * half_down_type_mask_bits);
/*
			if (possible_overflow)
				std::cout << "gg " << 
							to_hex_string(res_hi) << " " << 
							to_hex_string(res_lo) << " " << 
							std::endl;
*/
			long_uint<__deg> res;

			res.hi = hihi;
			res.lo = lolo;
/*
			if (possible_overflow)
				std::cout << "pr " << 
							to_hex_string(res) << " " << 
							std::endl;
*/
			res += res_hi;
			res += res_lo;
/*
			if (possible_overflow)
				std::cout << "as " << 
							to_hex_string(res) << " " << 
							std::endl;
*/
			return res;
		}


		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			bool operator<(const self_type& rhs) const
		{
			return hi < rhs.hi || (hi == rhs.hi && lo < rhs.lo);
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			bool operator==(const self_type& rhs) const
		{
			return hi == rhs.hi && lo == rhs.lo;
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			bool operator>(const self_type& rhs) const
		{
			return hi > rhs.hi || (hi == rhs.hi && lo > rhs.lo);
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			bool operator>=(const self_type& rhs) const
		{
			return !((*this) < rhs);
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			bool operator<=(const self_type& rhs) const
		{
			return !((*this) > rhs);
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			bool operator!=(const self_type& rhs) const
		{
			return !((*this) == rhs);
		}

		///* https://github.com/glitchub/arith64/blob/master/arith64.c *///
		template<std::uint64_t __deg = deg>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			resolved_return_type<(__deg > 0), size_type> __leading_zeros() const
		{
			if (hi == 0)
				return lo.__leading_zeros() + down_type_bits;
			return hi.__leading_zeros();
		}

		template<std::uint64_t __deg = deg>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			resolved_return_type<(__deg == 0), size_type> __leading_zeros() const
		{
			size_type b = 0, n = 0;
			base_type a = (hi) ? hi : ((n += base_bits), lo);
			b = !(a & 0xffffffff00000000ULL) << 5; n += b; a <<= b;
			b = !(a & 0xffff000000000000ULL) << 4; n += b; a <<= b;
			b = !(a & 0xff00000000000000ULL) << 3; n += b; a <<= b;
			b = !(a & 0xf000000000000000ULL) << 2; n += b; a <<= b;
			b = !(a & 0xc000000000000000ULL) << 1; n += b; a <<= b;
			return n + !(a & 0x8000000000000000ULL);
		}

		///* https://github.com/glitchub/arith64/blob/master/arith64.c *///
		template<std::uint64_t __deg = deg>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			resolved_return_type<(__deg > 0), size_type> __trailing_zeros() const
		{
			if (lo == 0)
				return hi.__trailing_zeros() + down_type_bits;
			return lo.__trailing_zeros();
		}

		template<std::uint64_t __deg = deg>
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			resolved_return_type<(__deg == 0), size_type> __trailing_zeros() const
		{
			size_type b = 0, n = 0;
			base_type a = (lo) ? lo : ((n += base_bits), hi);
			b = !(a & 0x00000000ffffffffULL) << 5; n += b; a >>= b;
			b = !(a & 0x000000000000ffffULL) << 4; n += b; a >>= b;
			b = !(a & 0x00000000000000ffULL) << 3; n += b; a >>= b;
			b = !(a & 0x000000000000000fULL) << 2; n += b; a >>= b;
			b = !(a & 0x0000000000000003ULL) << 1; n += b; a >>= b;
			return n + !(a & 0x0000000000000001ULL);
		}

		///* https://github.com/glitchub/arith64/blob/master/arith64.c *///
		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			static self_type __divmod(self_type lhs,
				const self_type& rhs,
				self_type& rem_out, bool assign_rem = true)
		{
			auto& a = lhs;
			const auto& b = rhs;

			if (b > a)
			{
				if (assign_rem)
					rem_out = a;
				return self_type();
			}
			if (b == 0)
				return self_type(); // shrug

			const self_type one(1);
			size_type iter_count = b.__leading_zeros() - a.__leading_zeros() + 1;
			self_type rem = a >> iter_count;
			a <<= bits - iter_count;
			self_type wrap = 0;
			while (iter_count-- > 0)
			{
				rem = ((rem << 1) | (a >> (bits - 1)));
				a = ((a << 1) | (wrap & one));
				wrap = (b > rem) ? self_type() : (~self_type()); // warning! hot spot?
				rem -= (b & wrap);
			}

			if (assign_rem)
				rem_out = rem;
			return (a << 1) | (wrap & one);
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			self_type operator/(const self_type& rhs) const
		{
			self_type res;
			res = __divmod(*this, rhs, res, false);
			return res;
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			self_type operator%(const self_type& rhs) const
		{
			self_type res;
			__divmod(*this, rhs, res, true);
			return res;
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			self_type& operator/=(const self_type& rhs)
		{
			return (*this = (*this / rhs));
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			self_type& operator%=(const self_type& rhs)
		{
			return (*this = (*this % rhs));
		}

		__DIXELU_CONDITIONAL_CPP14_SPECIFIERS
			self_type operator-() const
		{
			self_type res(*this);
			return (self_type() - res);
		}

		// dumb
		inline static std::string to_string(self_type value)
		{
			constexpr size_type radix_size = 19;
			std::string res_array[bits / (radix_size * 3 /*log2 of 10 ~=~ 3*/) + 1]; // SSO?
			std::string res;
			const self_type conversion_radix(10000000000000000000ull); // max
			//self_type conversion_radix(10000000000ull); // sso optimal
			const self_type zero;
			self_type rem;
			size_type idx = 0;
			while (value != zero)
			{
				value = __divmod(value, conversion_radix, rem);
				res_array[idx] = std::to_string(rem[0]);
				idx++;
			}
			if (!idx)
				return "0";
			bool first_run = true;
			while (idx-- > 0)
			{
				res += ((!first_run && res_array[idx].size() < radix_size) ?
					(std::string(radix_size - res_array[idx].size(), '0')) : "") +
					res_array[idx];
				first_run = false;
			}
			return res;
		}

		template<std::uint64_t __deg, bool frontal_0x = true>
		inline static resolved_return_type<(__deg > 0), std::string>
			to_hex_string(const long_uint<__deg>& value, const bool leading_zeros_flag = false)
		{
			std::string result;
			if (value.hi != 0)
				result += to_hex_string<__deg - 1, false>(value.hi, leading_zeros_flag);
			result += to_hex_string<__deg - 1, false>(value.lo, true);

			if (!leading_zeros_flag && result.size() > 1 && result.front() == '0')
			{
				auto it = std::find_if(result.begin(), result.end(), [](const char& c) {return c > '0'; });
				if(it == result.end())
					--it;
				if (it != result.begin())
					result = std::string(it, result.end());
			}
			if (frontal_0x)
				result = "0x" + result;
			return result;
		}

		template<std::uint64_t __deg, bool frontal_0x = true>
		inline static resolved_return_type<(__deg == 0), std::string>
			to_hex_string(const long_uint<__deg>& value, const bool leading_zeros_flag = false)
		{
			std::stringstream ss;
			if (value.hi)
				ss << std::setfill('0') << std::setw(16) << std::hex << value.hi;
			ss << std::setfill('0') << std::setw(16) << std::hex << value.lo;
			auto result = ss.str();
			
			if (!leading_zeros_flag && result.size() > 1 && result.front() == '0')
			{
				auto it = std::find_if(result.begin(), result.end(), [](const char& c) {return c > '0'; });
				if(it == result.end())
					--it;
				if (it != result.begin())
					result = std::string(it, result.end());
			}
			if (frontal_0x)
				result = "0x" + result;
			return result;
		}

		inline static std::string to_hex_string(const base_type& value)
		{
			std::stringstream ss;
			ss << std::setfill('0') << std::setw(16) << std::hex << value;
			auto result = ss.str();

			if (result.size() > 1 && result.front() == '0')
			{
				auto it = std::find_if(result.begin(), result.end(), [](const char& c) {return c > '0'; });
				if(it == result.end())
					--it;
				if (it != result.begin())
					result = std::string(it, result.end());
			}
			return "0x" + result;
		}

		inline static self_type __from_decimal_char_string(const char* str, size_type str_size)
		{
			self_type value;
			self_type radix = 10;

			while (str_size)
			{
				auto ull_value = *str - '0';
				value *= radix;
				value += self_type(ull_value);

				++str;
				--str_size;
			}

			return value;
		}

		inline static self_type __from_decimal_std_string(const std::string& str)
		{
			return __from_decimal_char_string(str.data(), str.size());
		}

		/*	static self_type __from_hex_string(const char* str, size_type str_size)
			{
				constexpr size_type radix_size = 8;
				if (str_size > 2 && (str[1] == 'x' || str[1] == 'X'))
				{
					str_size -= 2;
					str += 2;
				}
				self_type value;
				size_type leftout_size = std::min(radix_size, str_size);
				char* end = (char*)str + leftout_size;
				while (str_size)
				{
					auto ull_value = std::strtoull(str, &end, 16);
					value <<= 64;
					value |= self_type(ull_value);
					str += leftout_size;
					str_size -= leftout_size;
					leftout_size = std::min(radix_size, str_size);
					end = (char*)str + leftout_size;
				}
				return value;
			}*/
	};

	template <std::uint64_t deg>
	inline std::ostream& operator<<(std::ostream& out, const dixelu::long_uint<deg>& a)
	{
		return (out << a.to_string(a));
	}
}

namespace std
{
	template<dixelu::long_uint<0>::size_type deg>
	class numeric_limits<dixelu::long_uint<deg>>
	{
		using base = dixelu::long_uint<deg>;
	public:
		static constexpr bool is_specialized = true;

		static constexpr base min() noexcept { return base(); }
		static constexpr base max() noexcept { return base(~0, base::__fill_fields_tag()); }
		static constexpr base lowest() noexcept { return base(); }

		static constexpr int  digits = base::bits;
		static constexpr int  digits10 = base::bits / 3.32192809488736234787031942948939017586483139302458061205475639581593477660862521585013974335937015509966;
		static constexpr bool is_signed = false;
		static constexpr bool is_integer = true;
		static constexpr bool is_exact = true;
		static constexpr int  radix = 2;

		static constexpr base epsilon() noexcept { return base(); }
		static constexpr base round_error() noexcept { return base(); }

		static constexpr int min_exponent = 0;
		static constexpr int min_exponent10 = 0;
		static constexpr int max_exponent = 0;
		static constexpr int max_exponent10 = 0;
		static constexpr int max_digits10 = 0;

		static constexpr bool has_infinity = false;
		static constexpr bool has_quiet_NaN = false;
		static constexpr bool has_signaling_NaN = false;

		static constexpr float_denorm_style has_denorm = denorm_absent;
		static constexpr bool               has_denorm_loss = false;

		static constexpr base infinity() noexcept { return base(); }
		static constexpr base quiet_NaN() noexcept { return base(); }
		static constexpr base signaling_NaN() noexcept { return base(); }
		static constexpr base denorm_min() noexcept { return base(); }

		static constexpr bool is_iec559 = false;
		static constexpr bool is_bounded = true;
		static constexpr bool is_modulo = true;

		static constexpr bool traps = false;
		static constexpr bool tinyness_before = false;
		static constexpr float_round_style round_style = round_toward_zero;
	};
}
#endif //_DIXELU_INTEGERS_H_