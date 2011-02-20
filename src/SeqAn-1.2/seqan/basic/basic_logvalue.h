 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id: basic_logvalue.h 953 2007-07-27 11:48:23Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_LOGVALUE_H
#define SEQAN_HEADER_BASIC_LOGVALUE_H


namespace SEQAN_NAMESPACE_MAIN
{

	
//////////////////////////////////////////////////////////////////////////////
// LogProb
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TValue = double, typename TSpec = Default>
class LogProb;

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec>
class LogProb
{
public:	
	TValue data_value;

	//////////////////////////////////////////////////////////////////////////////

	LogProb() : data_value(std::log(0.0)) {	}

	//////////////////////////////////////////////////////////////////////////////

	template<typename TOtherValue>
	LogProb(TOtherValue const& _other) {
		data_value = std::log(_other);
	}

	template<typename TValue2, typename TSpec2>
	LogProb(LogProb<TValue2, TSpec2> const & _other) {
		data_value = _other.data_value;
	}

	//////////////////////////////////////////////////////////////////////////////

	~LogProb() {}

	//////////////////////////////////////////////////////////////////////////////

	operator int() const {
		return (int) std::exp(data_value);
	}

	operator float() const {
		return (float) std::exp(data_value);
	}

	operator double() const {
		return (double) std::exp(data_value);
	}


	//////////////////////////////////////////////////////////////////////////////

	template<typename TOtherValue>
	LogProb& operator=(TOtherValue const& rhs) {
		data_value = std::log(rhs);
		return *this;
	}

	template<typename TValue2, typename TSpec2>
	LogProb& operator=(LogProb<TValue2, TSpec2> const& rhs) {
		data_value = rhs.data_value;
		return *this;
	}

	//////////////////////////////////////////////////////////////////////////////
	
	template<typename TOtherValue>
	LogProb& operator*=(TOtherValue const& rhs) {
		data_value += std::log(rhs);
		return *this;
	}

	template<typename TValue2, typename TSpec2>
	LogProb& operator*=(LogProb<TValue2, TSpec2> const& rhs) {
		data_value += rhs.data_value;
		return *this;
	}

	//////////////////////////////////////////////////////////////////////////////

	template<typename TOtherValue>
	LogProb operator*(TOtherValue const& other) {
		LogProb result = *this;
		result *= LogProb(other);
		return result;
	}

	template<typename TValue2, typename TSpec2>
	LogProb operator*(LogProb<TValue2, TSpec2> const& other) {
		LogProb result = *this;
		result *= other;
		return result;
	}

	//////////////////////////////////////////////////////////////////////////////
	
	template<typename TOtherValue>
	LogProb& operator/=(TOtherValue const& rhs) {
		data_value -= std::log(rhs);
		return *this;
	}

	template<typename TValue2, typename TSpec2>
	LogProb& operator/=(LogProb<TValue2, TSpec2> const& rhs) {
		data_value -= rhs.data_value;
		return *this;
	}

	//////////////////////////////////////////////////////////////////////////////

	template<typename TOtherValue>
	LogProb operator/(TOtherValue const& other) {
		LogProb result = *this;
		result /= LogProb(other);
		return result;
	}

	template<typename TValue2, typename TSpec2>
	LogProb operator/(LogProb<TValue2, TSpec2> const& other) {
		LogProb result = *this;
		result /= other;
		return result;
	}

	//////////////////////////////////////////////////////////////////////////////
	
	template<typename TOtherValue>
	LogProb& operator+=(TOtherValue const& rhs) {
		data_value = std::log(std::exp(data_value) + rhs);
		return *this;
	}

	template<typename TValue2, typename TSpec2>
	LogProb& operator+=(LogProb<TValue2, TSpec2> const& rhs) {
		if (data_value > rhs.data_value) {
			if ((rhs.data_value == std::log(0.0)) || (data_value - rhs.data_value > 100)) return *this;
			data_value = data_value + std::log(1 + std::exp(rhs.data_value - data_value));
		} else {
			if ((data_value == std::log(0.0)) || (rhs.data_value - data_value > 100)) {
				data_value = rhs.data_value;
				return *this;
			}
			data_value = rhs.data_value + std::log(1 + std::exp(data_value - rhs.data_value));
		}
		return *this;
	}

	//////////////////////////////////////////////////////////////////////////////

	template<typename TOtherValue>
	LogProb operator+(TOtherValue const& other) {
		LogProb result = *this;
		result += LogProb(other);
		return result;
	}

	template<typename TValue2, typename TSpec2>
	LogProb operator+(LogProb<TValue2, TSpec2> const& other) {
		LogProb result = *this;
		result += other;
		return result;
	}

	//////////////////////////////////////////////////////////////////////////////
	
	template<typename TOtherValue>
	LogProb& operator-=(TOtherValue const& rhs) {
		data_value = std::log(std::exp(data_value) - rhs);
		return *this;
	}

	template<typename TValue2, typename TSpec2>
	LogProb& operator-=(LogProb<TValue2, TSpec2> const& rhs) {
		if (data_value > rhs.data_value) {
			if ((rhs.data_value == std::log(0.0)) || (data_value - rhs.data_value > 100)) return *this;
			data_value = data_value + std::log(1 - std::exp(rhs.data_value - data_value));
		} else {
			if ((data_value == std::log(0.0)) || (rhs.data_value - data_value > 100)) {
				data_value = rhs.data_value;
				return *this;
			}
			data_value = rhs.data_value + std::log(1 - std::exp(data_value - rhs.data_value));
		}
		return *this;
	}

	//////////////////////////////////////////////////////////////////////////////

	template<typename TOtherValue>
	LogProb operator-(TOtherValue const& other) {
		LogProb result = *this;
		result -= LogProb(other);
		return result;
	}

	template<typename TValue2, typename TSpec2>
	LogProb operator-(LogProb<TValue2, TSpec2> const& other) {
		LogProb result = *this;
		result -= other;
		return result;
	}

	//////////////////////////////////////////////////////////////////////////////

	template<typename TOtherValue>
	bool operator==(TOtherValue const& other) const {
		return data_value == std::log(other);
	}

	template<typename TValue2, typename TSpec2>
	bool operator==(LogProb<TValue2, TSpec2> const& other) const {
		return data_value == other.data_value;
	}

	//////////////////////////////////////////////////////////////////////////////

	template<typename TOtherValue>
	bool operator!=(TOtherValue const& other) const {
		return data_value != std::log(other);
	}

	template<typename TValue2, typename TSpec2>
	bool operator!=(LogProb<TValue2, TSpec2> const& other) const {
		return data_value != other.data_value;
	}

	//////////////////////////////////////////////////////////////////////////////

	template<typename TOtherValue>
	bool operator<(TOtherValue const& other) const {
		return data_value < std::log(other);
	}

	template<typename TValue2, typename TSpec2>
	bool operator<(LogProb<TValue2, TSpec2> const& other) const {
		return data_value < other.data_value;
	}

	//////////////////////////////////////////////////////////////////////////////

	template<typename TOtherValue>
	bool operator>(TOtherValue const& other) const {
		return data_value > std::log(other);
	}

	template<typename TValue2, typename TSpec2>
	bool operator>(LogProb<TValue2, TSpec2> const& other) const {
		return data_value > other.data_value;
	}

	//////////////////////////////////////////////////////////////////////////////

	template<typename TOtherValue>
	bool operator<=(TOtherValue const& other) const {
		return data_value <= std::log(other);
	}

	template<typename TValue2, typename TSpec2>
	bool operator<=(LogProb<TValue2, TSpec2> const& other) const {
		return data_value <= other.data_value;
	}

	//////////////////////////////////////////////////////////////////////////////

	template<typename TOtherValue>
	bool operator>=(TOtherValue const& other) const {
		return data_value >= std::log(other);
	}

	template<typename TValue2, typename TSpec2>
	bool operator>=(LogProb<TValue2, TSpec2> const& other) const {
		return data_value >= other.data_value;
	}


};

//////////////////////////////////////////////////////////////////////////////

template<typename TStream, typename TValue, typename TSpec>
TStream& operator<<(TStream& os, LogProb<TValue, TSpec> const& rhs) {
	return os << std::exp(rhs.data_value);
}
	

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
