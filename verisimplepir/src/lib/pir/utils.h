#include "head.h"
#include <iostream>
#include <functional>
#include "absl/numeric/int128.h"

#define absu128 absl::uint128

double currentDateTime();
void solinas30LFSR();
void solinas32LFSR();
void solinas45LFSR();
void solinas51LFSR();
void solinas54LFSR();
// Identity function to be used when we don't use primes
inline uint64_t fast_mod_solinas_0(absu128 v_big);
inline uint64_t fast_mod_solinas_0_max_64(uint64_t v);
// When primes are used
inline uint64_t fast_mod_solinas_30(absu128 v_big);
inline uint64_t fast_mod_solinas_30_max_64(uint64_t v);
inline uint64_t fast_mod_solinas_32(absu128 v_big);
inline uint64_t fast_mod_solinas_32_max_64(uint64_t v);
inline uint64_t fast_mod_solinas_40(absu128 v);
inline uint64_t fast_mod_solinas_40_max_64(uint64_t v);
inline uint64_t fast_mod_solinas_45(absu128 v);
inline uint64_t fast_mod_solinas_45_max_64(uint64_t v);
inline uint64_t fast_mod_solinas_51(absu128 v);
inline uint64_t fast_mod_solinas_51_max_64(uint64_t v);
inline uint64_t fast_mod_solinas_54(absu128 v);
inline uint64_t fast_mod_solinas_54_max_64(uint64_t v);
void test_mod_solinas_30();
void test_mod_solinas_32();
void test_mod_solinas_40();
void test_mod_solinas_45();
void test_mod_solinas_51();
void test_mod_solinas_54();

// ============ SETTINGS ============

#define RED "\033[31m"
#define BLUE "\033[34m"
#define YELLOW "\033[33m"
#define END "\033[0m"

#ifndef BSGS
#define BSGS
#endif

#ifndef FAKE_RUN
//#define FAKE_RUN
#endif

#ifndef SH_RUN // Run semi-honest (no verifiability)
//#define SH_RUN
#endif

// ====================

#ifdef SH_RUN
#define SOL_NTT_PRIME_BITS 0 
#define SOL_NTT_PRIME 0ULL // 0-bit prime
#define LOG_Q 32
#else
#define SOL_NTT_PRIME_BITS 32 // TODO SEE ABOVE
#define SOL_NTT_PRIME 4278255617ULL // 32-bit prime 
#define LOG_Q 31
#endif

#define VSPIR_LHE_PRIME_BITS 54
#define VSPIR_LHE_PRIME 18014398492704769ULL // 54-bit prime
#define VSPIR_LHE_LOG_Q 53

#define SOL_NTT_PRIME_BITS_SEC 30 // Not used by code currently
#define SOL_NTT_PRIME_SEC 1040220161ULL // 30-bit prime 
#define LOG_Q_SEC 29

#ifndef FAST_MOD_FUNC__
#define FAST_MOD_FUNC__
extern const std::function<uint64_t(uint64_t)> fast_mod_max_64; // fast_mod calls the 32-bit prime reduction function
extern const std::function<uint64_t(absu128)> fast_mod; // fast_mod calls the 32-bit prime reduction function

extern const std::function<uint64_t(uint64_t)> fast_mod_vspir_max_64;
extern const std::function<uint64_t(absu128)> fast_mod_vspir;

extern const std::function<uint64_t(uint64_t)> fast_mod_secondary_max_64; // fast_mod for secondary prime in CRT (not used in the current code)
extern const std::function<uint64_t(absu128)> fast_mod_secondary;
#endif

// --------

// To run native datatypes (power of two rings)
//#define SOL_NTT_PRIME_BITS 0 
//#define SOL_NTT_PRIME 0ULL // 0-bit prime
//#define LOG_Q 32/64 //8*sizeof(Elem)

//#define SOL_NTT_PRIME_BITS 32
#define SOL_NTT_PRIME_32 4278255617ULL // 32-bit prime
//#define LOG_Q 31

//#define SOL_NTT_PRIME_BITS 30
#define SOL_NTT_PRIME_30 1040220161ULL // 30-bit prime
//#define LOG_Q 29

//#define SOL_NTT_PRIME_BITS 35
#define SOL_NTT_PRIME_35 34630287361ULL // 35-bit prime
//#define LOG_Q 34

//#define SOL_NTT_PRIME_BITS 40
#define SOL_NTT_PRIME_40 1098438934529ULL // 40-bit prime
//#define LOG_Q 39

//#define SOL_NTT_PRIME_BITS 45
#define SOL_NTT_PRIME_45 35184372088321ULL // 45-bit prime
//#define LOG_Q 44

//#define SOL_NTT_PRIME_BITS 51
#define SOL_NTT_PRIME_51 2251799813554177ULL // 51-bit prime
//#define LOG_Q 50

//#define SOL_NTT_PRIME_BITS 54
#define SOL_NTT_PRIME_54 18014398492704769ULL // 54-bit prime
//#define LOG_Q 53