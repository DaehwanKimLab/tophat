
/**
 * \file misc.h
 * This file is part of FSA.
 * \author Source code in this file was written by Robert Bradley.
 */

#ifndef UTIL_MISC_INCLUDED
#define UTIL_MISC_INCLUDED

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <map>
#include <cassert>

#include "regexp.h"

using std::cerr;
using std::cout;
using std::endl;

namespace fsa {

  typedef std::pair<unsigned, unsigned> Interval;  ///< genomic interval

  struct Util {

  public:
      
      inline static std::string trim_right(const std::string &source , const std::string& t = " ");
      
      inline static std::string trim_left( const std::string& source, const std::string& t = " ");
      
      inline static std::string trim(const std::string& source, const std::string& t = " ");

    /**
     * \brief Convert to a string.
     */
    template<typename T>
    inline static std::string to_string (const T& t);

    /**
     * \brief Remove newline character, if it exists, from the end of a string.
     */
    inline static void chomp (std::string& str);

    /**
     * \brief Convert a string to upper-case.
     */
    inline static void toupper (std::string& str);

    /**
     * \brief Strip leading 'chr' if present.
     */
    static void strip_leading_chr (std::string& chromosome);

    /**
     * \brief Tokenize a string according to the specified delimiters.
     *
     * Code from:
     * http://oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html
     */
    inline static std::vector<std::string> tokenize_string (const std::string& str,
							    const std::string& delimiters = " ");

    /**
     * \brief Return a random number in [0, max].
     */
    inline static unsigned rand (const unsigned max);

    /**
     * \brief Does a file exist?
     *
     * Uses stat to see if we can get the file attributes.
     * From http://www.techbytes.ca/techbyte103.html.
     */
    static bool exists_file (const std::string& filename);

    /**
     * \brief Function object for comparing duples.
     *
     * Orders by first coordinate, then second coordinate.
     */
    template<typename T>
    struct Duple_comparator : std::binary_function<T, T, bool> {
    public:
      bool operator() (const std::pair<T, T> l, const std::pair<T, T> r) const {
	if (l.first == r.first)
	  return (l.second < r.second);
	else
	  return (l.first < r.first);
      }
    };

    /**
     * \brief Function object for comparing map entries based on their values.
     *
     * T might be, e.g., std::map<char, size_t>
     */
    template<typename T>
    struct Map_value_less : std::binary_function<T, T, bool> {
    public:
      bool operator() (const typename T::value_type& l, const typename T::value_type& r) const {
	return l.second < r.second;
      }
    };


  private:

    static Regexp re_chr_stripper;                  ///< \see strip_leading_chr



  };


  template<typename T>
    inline std::string Util::to_string (const T& t) {

    std::stringstream ss;
    ss << t;

    return ss.str();

  }

  inline void Util::chomp (std::string& str) {

    const int end = str.length() - 1;
    if (end >= 0 && str[end] == '\n')
      str.erase (end);

  }

  inline void Util::toupper (std::string& str) {

    for (std::string::iterator s = str.begin(); s != str.end(); ++s)
      *s = std::toupper (*s);

  }

  inline std::vector<std::string> Util::tokenize_string (const std::string& str,
							 const std::string& delimiters /* = " " */) {

    std::vector<std::string> tokens;

    // Skip delimiters at beginning.
    std::string::size_type lastpos = str.find_first_not_of (delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos = str.find_first_of (delimiters, lastpos);

    while (std::string::npos != pos || std::string::npos != lastpos) {
      // Found a token, add it to the vector.
      tokens.push_back (str.substr (lastpos, pos - lastpos));
      // Skip delimiters.  Note the "not_of"
      lastpos = str.find_first_not_of (delimiters, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of (delimiters, lastpos);
    }

    return tokens;

  }

  inline unsigned Util::rand (const unsigned max) {
    assert (max < RAND_MAX);
    return std::rand() % (max + 1);
  }

    inline std::string Util::trim_right(const std::string &source , const std::string& t)
    {
        std::string str = source;
        return str.erase( str.find_last_not_of(t) + 1);
    }
    
    inline std::string Util::trim_left( const std::string& source, const std::string& t)
    {
        std::string str = source;
        return str.erase(0 , source.find_first_not_of(t) );
    }
    
    inline std::string Util::trim(const std::string& source, const std::string& t)
    {
        std::string str = source;
        return trim_left( trim_right( str , t) , t );
    }

}

#endif /* UTIL_MISC_INCLUDED */
