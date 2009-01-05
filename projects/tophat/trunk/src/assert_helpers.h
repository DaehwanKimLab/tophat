#ifndef ASSERT_HELPERS_H_
#define ASSERT_HELPERS_H_

#include <stdexcept>
#include <string>
#include <cassert>
#include <iostream>

/**
 * Assertion for release-enabled assertions
 */
class ReleaseAssertException : public std::runtime_error {
public:
	ReleaseAssertException(const std::string& msg = "") : std::runtime_error(msg) {}
};

/**
 * Macros for release-enabled assertions, and helper macros to make
 * all assertion error messages more helpful.
 */
#ifndef NDEBUG
#define ASSERT_ONLY(x...) x
#else
#define ASSERT_ONLY(x...)
#endif

#define rt_assert(b)  \
	if(!(b)) { \
		std::cout << "rt_assert at " << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(); \
	}
#define rt_assert_msg(b,msg)  \
	if(!(b)) { \
		std::cout << msg <<  " at " << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(msg); \
	}

#define rt_assert_eq(ex,ac)  \
	if(!((ex) == (ac))) { \
		std::cout << "rt_assert_eq: expected (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(); \
	}
#define rt_assert_eq_msg(ex,ac,msg)  \
	if(!((ex) == (ac))) { \
		std::cout << "rt_assert_eq: " << msg <<  ": (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(msg); \
	}

#ifndef NDEBUG
#define assert_eq(ex,ac)  \
	if(!((ex) == (ac))) { \
		std::cout << "assert_eq: expected (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#define assert_eq_msg(ex,ac,msg)  \
	if(!((ex) == (ac))) { \
		std::cout << "assert_eq: " << msg <<  ": (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#else
#define assert_eq(ex,ac)
#define assert_eq_msg(ex,ac,msg)
#endif

#define rt_assert_neq(ex,ac)  \
	if(!((ex) != (ac))) { \
		std::cout << "rt_assert_neq: expected not (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(); \
	}
#define rt_assert_neq_msg(ex,ac,msg)  \
	if(!((ex) != (ac))) { \
		std::cout << "rt_assert_neq: " << msg << ": (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(msg); \
	}

#ifndef NDEBUG
#define assert_neq(ex,ac)  \
	if(!((ex) != (ac))) { \
		std::cout << "assert_neq: expected not (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#define assert_neq_msg(ex,ac,msg)  \
	if(!((ex) != (ac))) { \
		std::cout << "assert_neq: " << msg << ": (" << (ex) << ", 0x" << std::hex << (ex) << std::dec << ") got (" << (ac) << ", 0x" << std::hex << (ac) << std::dec << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#else
#define assert_neq(ex,ac)
#define assert_neq_msg(ex,ac,msg)
#endif

#define rt_assert_gt(a,b) \
	if(!((a) > (b))) { \
		std::cout << "rt_assert_gt: expected (" << (a) << ") > (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(); \
	}
#define rt_assert_gt_msg(a,b,msg) \
	if(!((a) > (b))) { \
		std::cout << "rt_assert_gt: " << msg << ": (" << (a) << ") > (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(msg); \
	}

#ifndef NDEBUG
#define assert_gt(a,b) \
	if(!((a) > (b))) { \
		std::cout << "assert_gt: expected (" << (a) << ") > (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#define assert_gt_msg(a,b,msg) \
	if(!((a) > (b))) { \
		std::cout << "assert_gt: " << msg << ": (" << (a) << ") > (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#else
#define assert_gt(a,b)
#define assert_gt_msg(a,b,msg)
#endif

#define rt_assert_geq(a,b) \
	if(!((a) >= (b))) { \
		std::cout << "rt_assert_geq: expected (" << (a) << ") >= (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(); \
	}
#define rt_assert_geq_msg(a,b,msg) \
	if(!((a) >= (b))) { \
		std::cout << "rt_assert_geq: " << msg << ": (" << (a) << ") >= (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(msg); \
	}

#ifndef NDEBUG
#define assert_geq(a,b) \
	if(!((a) >= (b))) { \
		std::cout << "assert_geq: expected (" << (a) << ") >= (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#define assert_geq_msg(a,b,msg) \
	if(!((a) >= (b))) { \
		std::cout << "assert_geq: " << msg << ": (" << (a) << ") >= (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#else
#define assert_geq(a,b)
#define assert_geq_msg(a,b,msg)
#endif

#define rt_assert_lt(a,b) \
	if(!(a < b)) { \
		std::cout << "rt_assert_lt: expected (" << a << ") < (" << b << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(); \
	}
#define rt_assert_lt_msg(a,b,msg) \
	if(!(a < b)) { \
		std::cout << "rt_assert_lt: " << msg << ": (" << a << ") < (" << b << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(msg); \
	}

#ifndef NDEBUG
#define assert_lt(a,b) \
	if(!(a < b)) { \
		std::cout << "assert_lt: expected (" << a << ") < (" << b << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#define assert_lt_msg(a,b,msg) \
	if(!(a < b)) { \
		std::cout << "assert_lt: " << msg << ": (" << a << ") < (" << b << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#else
#define assert_lt(a,b)
#define assert_lt_msg(a,b,msg)
#endif

#define rt_assert_leq(a,b) \
	if(!((a) <= (b))) { \
		std::cout << "rt_assert_leq: expected (" << (a) << ") <= (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(); \
	}
#define rt_assert_leq_msg(a,b,msg) \
	if(!((a) <= (b))) { \
		std::cout << "rt_assert_leq: " << msg << ": (" << (a) << ") <= (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		throw ReleaseAssertException(msg); \
	}

#ifndef NDEBUG
#define assert_leq(a,b) \
	if(!((a) <= (b))) { \
		std::cout << "assert_leq: expected (" << (a) << ") <= (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#define assert_leq_msg(a,b,msg) \
	if(!((a) <= (b))) { \
		std::cout << "assert_leq: " << msg << ": (" << (a) << ") <= (" << (b) << ")" << std::endl; \
		std::cout << __FILE__ << ":" << __LINE__ << std::endl; \
		assert(0); \
	}
#else
#define assert_leq(a,b)
#define assert_leq_msg(a,b,msg)
#endif

#endif /*ASSERT_HELPERS_H_*/
