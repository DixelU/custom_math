#ifndef _DIXELU_INTEGERS_H_
#define _DIXELU_INTEGERS_H_

#include <type_traits>
#include <cinttypes>

#include "math_utils.h"

namespace dixelu
{
#include <type_traits>
#include <cinttypes>

	template<std::size_t deg, typename base_type>
	struct bits
	{
		static constexpr base_type base_bits = sizeof(base_type) * 8;
		static constexpr base_type value = (deg == 0) ? base_bits : bits<deg - 1, base_type>::value * 2;
	};

	template<typename base_type>
	struct bits<0, base_type>
	{
		static constexpr base_type base_bits = sizeof(base_type) * 8;
		static constexpr base_type value = sizeof(base_type) * 8;
	};

	template <std::size_t deg>
	struct long_uint
	{
		using base_type = std::uint64_t;
		using size_type = std::uint64_t;
		using self_type = long_uint<deg>;
		using down_type = typename std::conditional<deg == 0, base_type, long_uint<deg - 1>>::type;

		static constexpr size_type bits = bits<deg, base_type>::value;

		down_type hi;
		down_type lo;

		long_uint() : hi(), lo() { }

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

		self_type operator~()
		{
			long_uint res;
			res.hi = ~hi;
			res.lo = ~lo;
			return res;
		}

		self_type operator|(const self_type& rhs)
		{
			long_uint res(*this);
			res |= rhs;
			return res;
		}

		self_type operator&(const self_type& rhs)
		{
			long_uint res(*this);
			res &= rhs;
			return res;
		}

		self_type operator^(const self_type& rhs)
		{
			long_uint res(*this);
			res ^= rhs;
			return res;
		}

		template< bool cond, typename U >
		using resolved_return_type = typename std::enable_if< cond, U >::type;

		template<std::uint64_t __deg = deg>
		resolved_return_type<(__deg > 0), bool> get_bit(size_type bit)
		{
			if (bit < (bits >> 1))
				return hi.get_bit(bit);
			bit -= (bits >> 1);
			return lo.get_bit(bit);
		}

		template<std::uint64_t __deg = deg>
		resolved_return_type<(__deg == 0), bool> get_bit(size_type bit)
		{
			if (bit < (bits >> 1))
				return (hi >> bit) & 1;
			bit -= (bits >> 1);
			return (hi >> bit) & 1;
		}

		self_type& operator+=(const self_type& rhs)
		{


		}
	};
}

#endif //_DIXELU_INTEGERS_H_