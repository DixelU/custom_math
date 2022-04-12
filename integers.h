
#ifndef _DIXELU_INTEGERS_H_
#define _DIXELU_INTEGERS_H_

#include <limits.h>
#include <type_traits>
#include <cinttypes>

namespace dixelu
{
	template<std::uintmax_t deg, typename base_type>
	struct bits
	{
		static constexpr base_type base_bits = sizeof(base_type) * CHAR_BIT;
		static constexpr base_type value = (bits<deg - 1, base_type>::value << 1);
	};

	template<typename base_type>
	struct bits<0, base_type>
	{
		static constexpr base_type base_bits = sizeof(base_type) * CHAR_BIT;
		static constexpr base_type value = base_bits * 2;
	};

	/* long_uint<0> ~ 128 bit unsigned integer */
	template <std::uintmax_t deg>
	struct long_uint
	{
		using base_type = std::uintmax_t; // largest of the integer types
		using size_type = std::uintmax_t;
		using self_type = long_uint<deg>;
		using down_type = typename std::conditional<deg == 0, base_type, long_uint<deg - 1>>::type;

		static constexpr size_type base_bits = bits<0, base_type>::base_bits;
		static constexpr size_type bits = bits<deg, base_type>::value;
		static constexpr size_type down_type_bits = bits >> 1;
		static constexpr size_type size = down_type_bits / base_bits;

		down_type hi;
		down_type lo;

		long_uint(base_type value = 0) : hi(), lo(value) { }

		template<uintmax_t __deg>
		explicit long_uint(const long_uint<__deg>& value,
			typename std::enable_if<(__deg < deg), void>::type* = 0) : hi(), lo(value) { }

		template<uintmax_t __deg>
		explicit long_uint(const long_uint<__deg>& value,
			typename std::enable_if<(__deg > deg), void>::type* = 0) :
			hi((typename std::conditional<__deg == 0, base_type, long_uint<__deg - 1>>::type)
				(value.lo >> (down_type_bits >> 1))),
			lo((typename std::conditional<__deg == 0, base_type, long_uint<__deg - 1>>::type)
				(value.lo& (long_uint<__deg>(1) << ((down_type_bits >> 1)))))
		{ }

		template<uintmax_t __deg>
		explicit long_uint(const long_uint<__deg>& value,
			typename std::enable_if<(__deg == deg), void>::type* = 0) : hi(value.hi), lo(value.lo) { }


		self_type& operator|=(const self_type& rhs)
		{
			hi |= rhs.hi;
			lo |= rhs.lo;
			return *this;
		}

		self_type& operator&=(const self_type& rhs)
		{
			hi &= rhs.hi;
			lo &= rhs.lo;
			return *this;
		}

		self_type& operator^=(const self_type& rhs)
		{
			hi ^= rhs.hi;
			lo ^= rhs.lo;
			return *this;
		}

		self_type operator~() const
		{
			long_uint res;
			res.hi = ~hi;
			res.lo = ~lo;
			return res;
		}

		self_type operator|(const self_type& rhs) const
		{
			long_uint res(*this);
			res |= rhs;
			return res;
		}

		self_type operator&(const self_type& rhs) const
		{
			long_uint res(*this);
			res &= rhs;
			return res;
		}

		self_type operator^(const self_type& rhs) const
		{
			long_uint res(*this);
			res ^= rhs;
			return res;
		}

		template< bool cond, typename U >
		using resolved_return_type = typename std::enable_if<cond, U>::type;

		template<std::uintmax_t __deg = deg>
		resolved_return_type<(__deg > 0), bool> get_bit(size_type bit) const
		{
			if (bit < down_type_bits)
				return lo.get_bit(bit);
			bit -= down_type_bits;
			return hi.get_bit(bit);
		}

		template<std::uintmax_t __deg = deg>
		resolved_return_type<(__deg == 0), bool> get_bit(size_type bit) const
		{
			if (bit < down_type_bits)
				return (lo >> bit) & 1;
			bit -= down_type_bits;
			return (hi >> bit) & 1;
		}

		template<std::uintmax_t __deg = deg>
		resolved_return_type<(__deg > 0), base_type&> operator[](size_type idx)
		{
			if (idx < down_type_bits / base_bits)
				return lo[idx];
			idx -= down_type_bits / base_bits;
			return hi[idx];
		}

		template<std::uintmax_t __deg = deg>
		resolved_return_type<(__deg == 0), base_type&> operator[](size_type idx)
		{
			if (idx < size)
				return lo;
			return hi;
		}

		template<std::uintmax_t __deg = deg>
		resolved_return_type<(__deg > 0), const base_type&> operator[](size_type idx) const
		{
			if (idx < size)
				return lo[idx];
			idx -= size;
			return hi[idx];
		}

		template<std::uintmax_t __deg = deg>
		resolved_return_type<(__deg == 0), const base_type&> operator[](size_type idx) const
		{
			if (idx < size)
				return lo;
			return hi;
		}

		template<std::uintmax_t __deg = deg>
		resolved_return_type<(__deg > 0), void>
			set_bit(size_type bit, bool bit_value)
		{
			if (bit < down_type_bits)
				return lo.set_bit(bit, bit_value);
			bit -= down_type_bits;
			return hi.set_bit(bit, bit_value);
		}

		template<std::uintmax_t __deg = deg>
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

		self_type& operator+=(const self_type& rhs)
		{
			const bool lhs_lo_carried_bit = get_bit(down_type_bits - 1);
			const bool rhs_lo_carried_bit = rhs.get_bit(down_type_bits - 1);
			const bool carry_possibility = lhs_lo_carried_bit | rhs_lo_carried_bit;

			hi += rhs.hi;
			lo += rhs.lo;

			const auto carried_bit = get_bit(down_type_bits - 1);
			if (carry_possibility && !carried_bit)
				hi += down_type(1);

			return *this;
		}

		self_type operator+(const self_type& rhs) const
		{
			self_type res(*this);
			res += rhs;
			return res;
		}

		self_type& operator-=(const self_type& rhs)
		{
			return ((*this) += ((~rhs) + self_type(1)));
		}

		self_type operator-(const self_type& rhs) const
		{
			self_type res(*this);
			res -= rhs;
			return res;
		}

		explicit operator size_type() const
		{
			return operator[](0);
		}

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

		self_type operator<<(size_type rhs) const
		{
			self_type res(*this);
			res <<= rhs;
			return res;
		}

		self_type operator>>(size_type rhs) const
		{
			self_type res(*this);
			res >>= rhs;
			return res;
		}

		template<std::uintmax_t __deg>
		explicit operator resolved_return_type<(__deg == deg || !__deg), long_uint<__deg>>() const
		{
			return *this;
		}

		template<std::uintmax_t __deg>
		explicit operator resolved_return_type<(__deg < deg&& __deg), long_uint<__deg>>() const
		{
			return (long_uint<deg - 1>)(lo);
		}

		template<std::uintmax_t __deg>
		explicit operator resolved_return_type<(__deg > deg&& __deg), long_uint<__deg>>() const
		{
			long_uint<deg + 1> res;
			res.lo = *this;
			return (long_uint<__deg>)(res);
		}

		self_type operator*(const self_type& rhs) const
		{
			auto up_result = long_uint<deg + 1>::__downtype_mul(
				(long_uint<deg + 1>)(*this),
				(long_uint<deg + 1>)rhs
			);
			return up_result.lo;
		}

		self_type& operator*=(const self_type& rhs)
		{
			return (*this = (*this * rhs));
		}

		template<std::uintmax_t __deg>
		static long_uint<__deg> __downtype_mul(const long_uint<__deg>& lhs, const long_uint<__deg>& rhs,
			typename std::enable_if<(__deg > 0), void>::type* = 0)
		{
			using local_down_type = typename long_uint<__deg>::down_type;
			constexpr size_type half_down_type_mask_bits = (long_uint<__deg>::down_type_bits >> 1);

			const local_down_type lolo_mask = (~local_down_type()) >> half_down_type_mask_bits;
			const local_down_type hilo_mask = ~lolo_mask;

			const local_down_type lhs_lo = lhs.lo & lolo_mask;
			const local_down_type rhs_lo = rhs.lo & lolo_mask;
			const local_down_type lhs_hi = (lhs.lo >> half_down_type_mask_bits) & lolo_mask;
			const local_down_type rhs_hi = (rhs.lo >> half_down_type_mask_bits) & lolo_mask;

			const local_down_type hihi = __downtype_mul<__deg - 1>(lhs_hi, rhs_hi);
			const local_down_type lolo = __downtype_mul<__deg - 1>(lhs_lo, rhs_lo);
			const local_down_type hilo = __downtype_mul<__deg - 1>((lhs_hi + rhs_lo), (lhs_lo + rhs_hi)) - hihi - lolo;

			long_uint<__deg> res;
			res.hi = hihi;
			res.lo = lolo;
			res.hi += (hilo & hilo_mask) >> half_down_type_mask_bits;
			res.lo += (hilo & lolo_mask) << half_down_type_mask_bits;

			return res;
		}

		template<std::uintmax_t __deg>
		static long_uint<__deg> __downtype_mul(const long_uint<__deg>& lhs, const long_uint<__deg>& rhs,
			typename std::enable_if<(__deg == 0), void>::type* = 0)
		{
			using local_down_type = typename long_uint<__deg>::down_type;
			constexpr size_type half_down_type_mask_bits = (long_uint<__deg>::down_type_bits >> 1);
			const local_down_type lolo_mask = (~local_down_type()) >> half_down_type_mask_bits;
			const local_down_type hilo_mask = ~lolo_mask;

			const local_down_type lhs_lo = lhs.lo & lolo_mask;
			const local_down_type rhs_lo = rhs.lo & lolo_mask;
			const local_down_type lhs_hi = (lhs.lo >> half_down_type_mask_bits) & lolo_mask;
			const local_down_type rhs_hi = (rhs.lo >> half_down_type_mask_bits) & lolo_mask;

			const local_down_type hihi = lhs_hi * rhs_hi;
			const local_down_type lolo = lhs_lo * rhs_lo;
			const local_down_type hilo = (lhs_hi + rhs_lo) * (lhs_lo + rhs_hi) - hihi - lolo;

			long_uint<__deg> res;
			res.hi = hihi;
			res.lo = lolo;
			res.hi += (hilo & hilo_mask) >> half_down_type_mask_bits;
			res.lo += (hilo & lolo_mask) << half_down_type_mask_bits;

			return res;
		}

		bool operator<(const self_type& rhs) const
		{
			size_type idx = size - 1;
			bool less = true;
			do
			{
				less &= (*this)[idx] < rhs[idx];
			} while (!idx && less);
			return !idx;
		}

		bool operator==(const self_type& rhs) const
		{
			size_type idx = size - 1;
			bool less = true;
			do
			{
				less &= (*this)[idx] == rhs[idx];
			} while (!idx && less);
			return !idx;
		}

		bool operator>(const self_type& rhs) const
		{
			size_type idx = size - 1;
			bool less = true;
			do
			{
				less &= (*this)[idx] > rhs[idx];
			} while (!idx && less);
			return !idx;
		}

		bool operator>=(const self_type& rhs) const
		{
			return !((*this) < rhs);
		}

		bool operator<=(const self_type& rhs) const
		{
			return !((*this) > rhs);
		}

		size_type __leading_zeros() const
		{

		}

		static self_type __divmod(const self_type& lhs,
			const self_type& rhs,
			self_type& modres)
		{
			if (lhs < rhs)
			{
				modres = lhs;
				return 0;
			}
			//....
		}
		/*
				self_type& operator/(const self_type& rhs) const
				{
					if ((*this) < rhs)


					return *this;
				}
		*/
	};
}
#endif //_DIXELU_INTEGERS_H_