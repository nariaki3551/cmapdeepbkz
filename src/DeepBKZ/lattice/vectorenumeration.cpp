#ifndef _inc_vector_enumeration
#define _inc_vector_enumeration

//Definitions to control the behaviours of enumeration

#define enum_mode_find_shortest 0x01
//update pruning radius during enumeration

#define enum_mode_all_vectors 0x02  
//enumerate all vectors s.t. |v|<bound and under the pruning function

#define enum_mode_find_abort 0x04
#define enum_mode_find_then_abort 0x04
//Abort the enumeration when a vector is found

#define enum_mode_find_short_projection 0x08 

#define enum_mode_except_trivial 0x10

#define enum_mode_count_only 0x20

#define finish_mode_nowaste 0x01
//In multithread-mode, when a non-working thread exists, then halt the subroutine

#define finish_mode_exact 0x02

//#slots (#holding vectors found in ENUM1a)        
#define slotmax 1000

//Activate if precise real cost is needed
#define depthcost       

#ifdef SSE_ROUND
#include <emmintrin.h>
#include <smmintrin.h>
#endif

#ifdef AVX_PRD
#include <immintrin.h>
#endif


#endif

#endif
