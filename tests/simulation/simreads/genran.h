#ifndef LH3_GENRAN_H_
#define LH3_GENRAN_H_

#include <stdlib.h>
#include <time.h>

#ifndef _WIN32 /* POSIX: rand48 family */

#include <sys/types.h>
#include <unistd.h>

#define ran_seed() srand48(time(0) * (long)getpid())
#define ran_uniform() drand48()

#else /* Windows: this will be pretty BAD. */

#define ran_seed() srand(time(0))
#define ran_uniform() ((double)rand() / RAND_MAX)

#endif


extern "C" {
	double ran_normal();
}


#endif
