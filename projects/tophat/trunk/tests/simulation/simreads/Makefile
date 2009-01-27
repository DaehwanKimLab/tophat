DEBUG_FLAGS = -O0 -g3
RELEASE_FLAGS = -O3 
NOASSERT_FLAGS = -DNDEBUG
64BIT_FLAGS = -m64


simreads: simulate.c seq.c const.c genran.c
	$(CXX) $(DEBUG_FLAGS)	-lz -Wall -o $@ $< seq.c const.c genran.c