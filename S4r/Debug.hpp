#ifndef S4R_DEBUG_HPP_INCLUDED
#define S4R_DEBUG_HPP_INCLUDED

#include <cstdio>
#ifdef S4R_ENABLE_TRACE
# define S4R_TRACE(...) fprintf(stderr, __VA_ARGS__)
# define S4R_ASSERT(COND)
#else
# define S4R_TRACE(...)
# define S4R_ASSERT(COND)
#endif

#define S4R_VERB(verb,...) do{          \
	if(S->options.verbosity >= verb){   \
		fprintf(stdout, "[%d] ", verb); \
		fprintf(stdout, __VA_ARGS__);   \
	}                                   \
}while(0)

#endif // S4R_DEBUG_HPP_INCLUDED
