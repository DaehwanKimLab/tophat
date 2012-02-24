// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

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
