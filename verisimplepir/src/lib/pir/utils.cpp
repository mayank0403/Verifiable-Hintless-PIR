#include <iostream>
#include <fstream>
#include "time.h"
#include <chrono>
#include "utils.h"
#include <random>
#include <cmath>

#if(SOL_NTT_PRIME_BITS == 0)
    const std::function<uint64_t(uint64_t)> fast_mod_max_64 = fast_mod_solinas_0_max_64;
    const std::function<uint64_t(absu128)> fast_mod = fast_mod_solinas_0;
#elif(SOL_NTT_PRIME_BITS == 32)
    const std::function<uint64_t(uint64_t)> fast_mod_max_64 = fast_mod_solinas_32_max_64;
    const std::function<uint64_t(absu128)> fast_mod = fast_mod_solinas_32;
#elif(SOL_NTT_PRIME_BITS == 40)
    const std::function<uint64_t(uint64_t)> fast_mod_max_64 = fast_mod_solinas_40_max_64;
    const std::function<uint64_t(absu128)> fast_mod = fast_mod_solinas_40;
#elif(SOL_NTT_PRIME_BITS == 45)
    const std::function<uint64_t(uint64_t)> fast_mod_max_64 = fast_mod_solinas_45_max_64;
    const std::function<uint64_t(absu128)> fast_mod = fast_mod_solinas_45;
#elif(SOL_NTT_PRIME_BITS == 51)
    const std::function<uint64_t(uint64_t)> fast_mod_max_64 = fast_mod_solinas_51_max_64;
    const std::function<uint64_t(absu128)> fast_mod = fast_mod_solinas_51;
#elif(SOL_NTT_PRIME_BITS == 54)
    const std::function<uint64_t(uint64_t)> fast_mod_max_64 = fast_mod_solinas_54_max_64;
    const std::function<uint64_t(absu128)> fast_mod = fast_mod_solinas_54;
#endif

#if(VSPIR_LHE_PRIME_BITS == SOL_NTT_PRIME_BITS)
    const std::function<uint64_t(uint64_t)> fast_mod_vspir_max_64 = fast_mod_max_64;
    const std::function<uint64_t(absu128)> fast_mod_vspir = fast_mod;
#else
    const std::function<uint64_t(uint64_t)> fast_mod_vspir_max_64 = fast_mod_solinas_54_max_64;
    const std::function<uint64_t(absu128)> fast_mod_vspir = fast_mod_solinas_54;
#endif

#if(SOL_NTT_PRIME_BITS_SEC == 0)
    const std::function<uint64_t(uint64_t)> fast_mod_secondary_max_64 = fast_mod_solinas_0_max_64;
    const std::function<uint64_t(absu128)> fast_mod_secondary = fast_mod_solinas_0;
#elif(SOL_NTT_PRIME_BITS_SEC == 30)
    const std::function<uint64_t(uint64_t)> fast_mod_secondary_max_64 = fast_mod_solinas_30_max_64;
    const std::function<uint64_t(absu128)> fast_mod_secondary = fast_mod_solinas_30;
#endif

double currentDateTime()
{

    std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();

    time_t tnow = std::chrono::system_clock::to_time_t(now);
    tm *date = localtime(&tnow); // todo: dperecated use localtime_s
    date->tm_hour = 0;
    date->tm_min = 0;
    date->tm_sec = 0;

    auto midnight = std::chrono::system_clock::from_time_t(mktime(date));

    return std::chrono::duration<double, std::milli>(now - midnight).count();
}

// Solinas LFSR matrix
// This is for the prime p = 2^30 - 2^25 + 2^15 + 1
// Souce: https://crysp.uwaterloo.ca/software/gmnt/gmz_mcs_thesis.pdf [Greg Zaverucha's thesis]
void solinas30LFSR(){
    int d = 6;
    std::vector<std::vector<int32_t>> lfsr(d, std::vector<int32_t>(d, 0));
    std::vector<int32_t> coeffs(d+1, 0); // to get indexing from 1 to avoid mistakes
    coeffs[1] = 1;
    coeffs[3] = -1;
    coeffs[d] = -1;
    for (int i = 0; i < d; i++){
        lfsr[0][i] = coeffs[d-i];
    }

    for (int i = 1; i < d; i++){
        lfsr[i][0] = lfsr[i-1][d-1] * coeffs[d];
        for (int j = 1; j < d; j++){
            lfsr[i][j] = lfsr[i-1][j-1] + lfsr[i-1][d-1] * coeffs[d-j];
        }
    }

    std::cout << "Solinas 30 LFSR matrix:\n";
    for (int i = 0; i < d; i++){
        for (int j = 0; j < d; j++){
            std::cout << lfsr[i][j] << " ";
        }
        std::cout << std::endl;
    }

}

// Solinas LFSR matrix
// This is for the prime p = 2^32 - 2^24 + 2^16 + 1
// Souce: https://crysp.uwaterloo.ca/software/gmnt/gmz_mcs_thesis.pdf [Greg Zaverucha's thesis]
void solinas32LFSR(){
    int d = 4;
    std::vector<std::vector<int32_t>> lfsr(d, std::vector<int32_t>(d, 0));
    std::vector<int32_t> coeffs(d+1, 0); // to get indexing from 1 to avoid mistakes
    coeffs[1] = 1;
    coeffs[2] = -1;
    coeffs[d] = -1;
    for (int i = 0; i < d; i++){
        lfsr[0][i] = coeffs[d-i];
    }

    for (int i = 1; i < d; i++){
        lfsr[i][0] = lfsr[i-1][d-1] * coeffs[d];
        for (int j = 1; j < d; j++){
            lfsr[i][j] = lfsr[i-1][j-1] + lfsr[i-1][d-1] * coeffs[d-j];
        }
    }

    std::cout << "Solinas 32 LFSR matrix:\n";
    for (int i = 0; i < d; i++){
        for (int j = 0; j < d; j++){
            std::cout << lfsr[i][j] << " ";
        }
        std::cout << std::endl;
    }

}


// Solinas LFSR matrix
// This is for the prime p = 2^45 - 2^9 + 1
// This can be written as the polynomial f(t) = t^5 - t + 1
// Souce: https://crysp.uwaterloo.ca/software/gmnt/gmz_mcs_thesis.pdf [Greg Zaverucha's thesis]
void solinas45LFSR(){
    int d = 5;
    std::vector<std::vector<int32_t>> lfsr(d, std::vector<int32_t>(d, 0));
    std::vector<int32_t> coeffs(d+1, 0); // to get indexing from 1 to avoid mistakes
    coeffs[d-1] = 1;
    coeffs[d] = -1;
    for (int i = 0; i < d; i++){
        lfsr[0][i] = coeffs[d-i];
    }

    for (int i = 1; i < d; i++){
        lfsr[i][0] = lfsr[i-1][d-1] * coeffs[d];
        for (int j = 1; j < d; j++){
            lfsr[i][j] = lfsr[i-1][j-1] + lfsr[i-1][d-1] * coeffs[d-j];
        }
    }

    std::cout << "Solinas 45 LFSR matrix:\n";
    for (int i = 0; i < d; i++){
        for (int j = 0; j < d; j++){
            std::cout << lfsr[i][j] << " ";
        }
        std::cout << std::endl;
    }

}

void solinas51LFSR(){
    int d = 3;
    std::vector<std::vector<int32_t>> lfsr(d, std::vector<int32_t>(d, 0));
    std::vector<int32_t> coeffs(d+1, 0); // to get indexing from 1 to avoid mistakes
    coeffs[d-1] = 1;
    coeffs[d] = -1;
    for (int i = 0; i < d; i++){
        lfsr[0][i] = coeffs[d-i];
    }

    for (int i = 1; i < d; i++){
        lfsr[i][0] = lfsr[i-1][d-1] * coeffs[d];
        for (int j = 1; j < d; j++){
            lfsr[i][j] = lfsr[i-1][j-1] + lfsr[i-1][d-1] * coeffs[d-j];
        }
    }

    std::cout << "Solinas 51 LFSR matrix:\n";
    for (int i = 0; i < d; i++){
        for (int j = 0; j < d; j++){
            std::cout << lfsr[i][j] << " ";
        }
        std::cout << std::endl;
    }

}

// (1<<54) - (1<<24) + 1 [18014398492704769ULL]
void solinas54LFSR(){
    int d = 9;
    std::vector<std::vector<int32_t>> lfsr(d, std::vector<int32_t>(d, 0));
    std::vector<int32_t> coeffs(d+1, 0); // to get indexing from 1 to avoid mistakes
    coeffs[5] = 1;
    coeffs[d] = -1;
    for (int i = 0; i < d; i++){
        lfsr[0][i] = coeffs[d-i];
    }

    for (int i = 1; i < d; i++){
        lfsr[i][0] = lfsr[i-1][d-1] * coeffs[d];
        for (int j = 1; j < d; j++){
            lfsr[i][j] = lfsr[i-1][j-1] + lfsr[i-1][d-1] * coeffs[d-j];
        }
    }

    std::cout << "Solinas 54 LFSR matrix:\n";
    for (int i = 0; i < d; i++){
        for (int j = 0; j < d; j++){
            std::cout << lfsr[i][j] << " ";
        }
        std::cout << std::endl;
    }

}

// Identity function to be used when we don't use primes
inline uint64_t fast_mod_solinas_0(absu128 v_big){
    // Always fits in 64 bits
    uint64_t v = absl::Uint128Low64(v_big);
    return v;
}

// Identity function to be used when we don't use primes
inline uint64_t fast_mod_solinas_0_max_64(uint64_t v){
    return v;
}

// NOTE: this is not constant time!
// prime modulus is 2^30 - 2^25 + 2^15 + 1
inline uint64_t fast_mod_solinas_30(absu128 v_big) {
    // v is split into 12x 5-bit chunks a_0, a_1, ..., a_11 {MSB}

    // Always fits in 64 bits
    uint64_t v = absl::Uint128Low64(v_big);
    
    if (v < SOL_NTT_PRIME_30) {return v;}
    //else if (v == SOL_NTT_PRIME_30) {return 0;}
    else if (v < (SOL_NTT_PRIME_30 << 1)) {return v - SOL_NTT_PRIME_30;}
    
    uint64_t mask_5 = (1ULL << 5) - 1;
    uint64_t mask_10 = (1ULL << 10) - 1;
    uint64_t mask_15 = (1ULL << 15) - 1;
    uint64_t mask_30 = (1ULL << 30) - 1;

    uint64_t a_lower = v & mask_30;
    uint64_t a_upper = v >> 30;

    uint64_t a_6 = a_upper & mask_5;
    uint64_t a_7 = (a_upper >> 5) & mask_5;
    uint64_t a_67 = a_upper & mask_10;
    uint64_t a_79 = (a_upper >> 5) & mask_15;
    uint64_t a_810 = (a_upper >> 10) & mask_15;
    uint64_t a_10 = (a_upper >> 20) & mask_5;
    uint64_t a_1011 = (a_upper >> 20) & mask_10;
    uint64_t a_11 = a_upper >> 25;

    uint64_t v_reduced = a_lower - a_upper - a_79 - a_810 + a_1011 + (a_11 << 1) - (a_67 << 15) - (a_79 << 15) - (a_810 << 15) + (a_11 << 15) + (a_6 << 25) + (a_7 << 25) - (a_10 << 25) - (a_11 << 26); 

    // These loops run for 2-3 times max
    while (int64_t(v_reduced) < 0) {
        v_reduced += SOL_NTT_PRIME_30;
        //std::cout << "adds ";
    }
    //std::cout<<std::endl;

    while (v_reduced >= SOL_NTT_PRIME_30) {
        //std::cout << "subs ";
        v_reduced -= SOL_NTT_PRIME_30;
    }
    //std::cout<<std::endl;

    return v_reduced;
}

// prime modulus is 2^30 - 2^25 + 2^15 + 1
inline uint64_t fast_mod_solinas_30_max_64(uint64_t v) {
    // v is split into 12x 5-bit chunks a_0, a_1, ..., a_11 {MSB}
    if (v < SOL_NTT_PRIME_30) {return v;}
    //else if (v == SOL_NTT_PRIME_30) {return 0;}
    else if (v < (SOL_NTT_PRIME_30 << 1)) {return v - SOL_NTT_PRIME_30;}
    
    uint64_t mask_5 = (1ULL << 5) - 1;
    uint64_t mask_10 = (1ULL << 10) - 1;
    uint64_t mask_15 = (1ULL << 15) - 1;
    uint64_t mask_30 = (1ULL << 30) - 1;

    uint64_t a_lower = v & mask_30;
    uint64_t a_upper = v >> 30;

    uint64_t a_6 = a_upper & mask_5;
    uint64_t a_7 = (a_upper >> 5) & mask_5;
    uint64_t a_67 = a_upper & mask_10;
    uint64_t a_79 = (a_upper >> 5) & mask_15;
    uint64_t a_810 = (a_upper >> 10) & mask_15;
    uint64_t a_10 = (a_upper >> 20) & mask_5;
    uint64_t a_1011 = (a_upper >> 20) & mask_10;
    uint64_t a_11 = a_upper >> 25;

    uint64_t v_reduced = a_lower - a_upper - a_79 - a_810 + a_1011 + (a_11 << 1) - (a_67 << 15) - (a_79 << 15) - (a_810 << 15) + (a_11 << 15) + (a_6 << 25) + (a_7 << 25) - (a_10 << 25) - (a_11 << 26); 

    // These loops run for 2-3 times max
    while (int64_t(v_reduced) < 0) {
        v_reduced += SOL_NTT_PRIME_30;
        //std::cout << "adds ";
    }
    //std::cout<<std::endl;

    while (v_reduced >= SOL_NTT_PRIME_30) {
        //std::cout << "subs ";
        v_reduced -= SOL_NTT_PRIME_30;
    }
    //std::cout<<std::endl;

    return v_reduced;
}

// prime modulus is 2^32 - 2^24 + 2^16 + 1
inline uint64_t fast_mod_solinas_32(absu128 v_big) {
    // v is split into 8x 8-bit chunks a_0, a_1, ..., a_7 {MSB}

    // Always fits in 64 bits
    uint64_t v = absl::Uint128Low64(v_big);

    if (v < SOL_NTT_PRIME_32) return v;
    if (v == SOL_NTT_PRIME_32) return 0;
    
    uint64_t mask_8 = (1ULL << 8) - 1;
    uint64_t mask_16 = (1ULL << 16) - 1;
    uint64_t mask_24 = (1ULL << 24) - 1;
    uint64_t mask_32 = (1ULL << 32) - 1;

    uint64_t a_lower = v & mask_32;
    uint64_t a_upper = v >> 32;

    uint64_t a_4 = a_upper & mask_8;
    uint64_t a_7 = (a_upper >> 24) & mask_8;
    uint64_t a_56 = (a_upper >> 8) & mask_16;

    uint64_t v_reduced = a_lower - (a_7 << 24) - a_upper - (a_56 << 16) - a_56 + a_7 + (a_4 << 24) - (a_4 << 16); 

    //int64_t t1 = a_7 - a_4 - a_5;
    //int64_t t2 = 0 - a_5 - a_6;
    //int64_t t3 = 0 - a_4 - a_5 - a_6;
    //int64_t t4 = a_4 - a_6 - a_7 - a_7;
    //int64_t t5 = t1 + t2*(1<<8) + t3*(1<<16) + t4*(1<<24);
    //std::cout<<a_lower<<" + "<<t5<<" = "<<a_lower + t5<<std::endl;
    //uint64_t v_reduced = a_lower + t5;

    // These loops run for 2-3 times max
    while (int64_t(v_reduced) < 0) {
        v_reduced += SOL_NTT_PRIME_32;
        //std::cout << "adds ";
    }
    //std::cout<<std::endl;

    while (v_reduced >= SOL_NTT_PRIME_32) {
        //std::cout << "subs ";
        v_reduced -= SOL_NTT_PRIME_32;
    }
    //std::cout<<std::endl;

    return v_reduced;
}

// prime modulus is 2^32 - 2^24 + 2^16 + 1
inline uint64_t fast_mod_solinas_32_max_64(uint64_t v) {
    // v is split into 8x 8-bit chunks a_0, a_1, ..., a_7 {MSB}
    if (v < SOL_NTT_PRIME_32) {return v;}
    //else if (v == SOL_NTT_PRIME_32) {return 0;}
    else if (v < (SOL_NTT_PRIME_32 << 1)) {return v - SOL_NTT_PRIME_32;}
    //else if (v < (3 * SOL_NTT_PRIME_32)) {return v - (SOL_NTT_PRIME_32 << 1);}
    
    uint64_t mask_8 = (1ULL << 8) - 1;
    uint64_t mask_16 = (1ULL << 16) - 1;
    uint64_t mask_24 = (1ULL << 24) - 1;
    uint64_t mask_32 = (1ULL << 32) - 1;

    uint64_t a_lower = v & mask_32;
    uint64_t a_upper = v >> 32;

    uint64_t a_4 = a_upper & mask_8;
    uint64_t a_7 = (a_upper >> 24) & mask_8;
    uint64_t a_56 = (a_upper >> 8) & mask_16;

    uint64_t v_reduced = a_lower - (a_7 << 24) - a_upper - (a_56 << 16) - a_56 + a_7 + (a_4 << 24) - (a_4 << 16); 

    //int64_t t1 = a_7 - a_4 - a_5;
    //int64_t t2 = 0 - a_5 - a_6;
    //int64_t t3 = 0 - a_4 - a_5 - a_6;
    //int64_t t4 = a_4 - a_6 - a_7 - a_7;
    //int64_t t5 = t1 + t2*(1<<8) + t3*(1<<16) + t4*(1<<24);
    //std::cout<<a_lower<<" + "<<t5<<" = "<<a_lower + t5<<std::endl;
    //uint64_t v_reduced = a_lower + t5;

    // These loops run for 2-3 times max
    while (int64_t(v_reduced) < 0) {
        v_reduced += SOL_NTT_PRIME_32;
        //std::cout << "adds ";
    }
    //std::cout<<std::endl;

    while (v_reduced >= SOL_NTT_PRIME_32) {
        //std::cout << "subs ";
        v_reduced -= SOL_NTT_PRIME_32;
    }
    //std::cout<<std::endl;

    return v_reduced;
}

// prime modulus is 2^40 - 2^30 + 2^20 + 1
inline uint64_t fast_mod_solinas_40(absu128 v) {
    // v is split into 8x 10-bit chunks a_0, a_1, ..., a_7 {MSB}
    if (v < SOL_NTT_PRIME_40) return static_cast<uint64_t>(v);
    if (v == SOL_NTT_PRIME_40) return 0;

    uint64_t high64 = absl::Uint128High64(v);
    uint64_t low64 = absl::Uint128Low64(v);
    
    uint64_t mask_10 = (1ULL << 10) - 1;
    uint64_t mask_20 = (1ULL << 20) - 1;
    uint64_t mask_30 = (1ULL << 30) - 1;
    uint64_t mask_40 = (1ULL << 40) - 1;

    uint64_t a_lower = low64 & mask_40;
    uint64_t a_upper = (low64 >> 40) | (high64 << 24);

    uint64_t a_4 = a_upper & mask_10;
    uint64_t a_7 = (a_upper >> 30) & mask_10;
    uint64_t a_56 = (a_upper >> 10) & mask_20;

    uint64_t v_reduced = a_lower - (a_7 << 30) - a_upper - (a_56 << 20) - a_56 + a_7 + (a_4 << 30) - (a_4 << 20); 

    // These loops run for 2-3 times max
    while (int64_t(v_reduced) < 0) {
        v_reduced += SOL_NTT_PRIME_40;
        //std::cout << "adds ";
    }
    //std::cout<<std::endl;

    while (v_reduced >= SOL_NTT_PRIME_40) {
        //std::cout << "subs ";
        v_reduced -= SOL_NTT_PRIME_40;
    }
    //std::cout<<std::endl;

    return v_reduced;
}

// prime modulus is 2^40 - 2^30 + 2^20 + 1
inline uint64_t fast_mod_solinas_40_max_64(uint64_t v) {
    // v is split into 8x 10-bit chunks a_0, a_1, ..., a_7 {MSB}
    // However, since input is at most 64 bits, a_7 = 0.
    if (v < SOL_NTT_PRIME_40) {return v;}
    //else if (v == SOL_NTT_PRIME_40) {return 0;}
    else if (v < (SOL_NTT_PRIME_40 << 1)) {return v - SOL_NTT_PRIME_40;}
    //else if (v < (3 * SOL_NTT_PRIME_40)) {return v - (SOL_NTT_PRIME_40 << 1);}
    
    uint64_t mask_10 = (1ULL << 10) - 1;
    uint64_t mask_20 = (1ULL << 20) - 1;
    uint64_t mask_30 = (1ULL << 30) - 1;
    uint64_t mask_40 = (1ULL << 40) - 1;

    uint64_t a_lower = v & mask_40;
    uint64_t a_upper = v >> 40;

    uint64_t a_4 = a_upper & mask_10;
    //uint64_t a_7 = (a_upper >> 30) & mask_10;
    uint64_t a_56 = (a_upper >> 10) & mask_20;

    uint64_t v_reduced = a_lower - a_upper - (a_56 << 20) - a_56 + (a_4 << 30) - (a_4 << 20);

    //uint64_t v_reduced = a_lower - (a_7 << 30) - a_upper - (a_56 << 20) - a_56 + a_7 + (a_4 << 30) - (a_4 << 20); 

    // These loops run for 2-3 times max
    while (int64_t(v_reduced) < 0) {
        v_reduced += SOL_NTT_PRIME_40;
        //std::cout << "adds ";
    }
    //std::cout<<std::endl;

    while (v_reduced >= SOL_NTT_PRIME_40) {
        //std::cout << "subs ";
        v_reduced -= SOL_NTT_PRIME_40;
    }
    //std::cout<<std::endl;

    return v_reduced;
}

// prime modulus is 2^45 - 2^9 + 1
inline uint64_t fast_mod_solinas_45(absu128 v) {
    // v is split into 10x 9-bit chunks a_0, a_1, ..., a_9 {MSB}
    if (v < SOL_NTT_PRIME_45) return static_cast<uint64_t>(v);
    if (v == SOL_NTT_PRIME_45) return 0;

    uint64_t high64 = absl::Uint128High64(v);
    uint64_t low64 = absl::Uint128Low64(v);
    
    uint64_t mask_9 = (1ULL << 9) - 1;
    uint64_t mask_36 = (1ULL << 36) - 1;
    uint64_t mask_45 = (1ULL << 45) - 1;

    uint64_t a_lower = low64 & mask_45;
    uint64_t a_upper = (low64 >> 45) | (high64 << 19);

    uint64_t a_58 = a_upper & mask_36;
    uint64_t a_9 = (a_upper >> 36) & mask_9;

    uint64_t v_reduced = a_lower + (a_58 << 9) - a_upper + (a_9 << 9) - a_9; 

    // These loops run for 2-3 times max
    while (int64_t(v_reduced) < 0) {
        v_reduced += SOL_NTT_PRIME_45;
        //std::cout << "adds ";
    }
    //std::cout<<std::endl;

    while (v_reduced >= SOL_NTT_PRIME_45) {
        //std::cout << "subs ";
        v_reduced -= SOL_NTT_PRIME_45;
    }
    //std::cout<<std::endl;

    return v_reduced;
}

// prime modulus is 2^45 - 2^9 + 1
inline uint64_t fast_mod_solinas_45_max_64(uint64_t v) {
    // v is split into 8x 10-bit chunks a_0, a_1, ..., a_7 {MSB}
    // However, since input is at most 64 bits, a_7 = 0.
    if (v < SOL_NTT_PRIME_45) {return v;}
    //else if (v == SOL_NTT_PRIME_45) {return 0;}
    else if (v < (SOL_NTT_PRIME_45 << 1)) {return v - SOL_NTT_PRIME_45;}
    //else if (v < (3 * SOL_NTT_PRIME_45)) {return v - (SOL_NTT_PRIME_45 << 1);}
    
    //uint64_t mask_9 = (1ULL << 9) - 1;
    uint64_t mask_36 = (1ULL << 36) - 1;
    uint64_t mask_45 = (1ULL << 45) - 1;

    uint64_t a_lower = v & mask_45;
    uint64_t a_upper = (v >> 45);

    //uint64_t a_58 = a_upper & mask_36;
    //uint64_t a_9 = (a_upper >> 36) & mask_9;

    uint64_t v_reduced = a_lower + (a_upper << 9) - a_upper; 

    // These loops run for 2-3 times max
    while (int64_t(v_reduced) < 0) {
        v_reduced += SOL_NTT_PRIME_45;
        //std::cout << "adds ";
    }
    //std::cout<<std::endl;

    while (v_reduced >= SOL_NTT_PRIME_45) {
        //std::cout << "subs ";
        v_reduced -= SOL_NTT_PRIME_45;
    }
    //std::cout<<std::endl;

    return v_reduced;
}

// prime modulus is 2^51 - 2^17 + 1
inline uint64_t fast_mod_solinas_51(absu128 v) {
    // Input would be expected to be 102 bits long and we would split it into 6x 17-bit chunks a_0, a_1, ..., a_5 {MSB}
    if (v < SOL_NTT_PRIME_51) return static_cast<uint64_t>(v);
    if (v == SOL_NTT_PRIME_51) return 0;
    
    uint64_t high64 = absl::Uint128High64(v);
    uint64_t low64 = absl::Uint128Low64(v);

    uint64_t mask_51 = (1ULL << 51) - 1;
    uint64_t mask_34 = (1ULL << 34) - 1;

    uint64_t a_lower = low64 & mask_51;
    uint64_t a_upper = (low64 >> 51) | (high64 << 13);
    uint64_t a_34 = a_upper & mask_34;
    uint64_t a_5 = (a_upper >> 34);

    uint64_t v_reduced = a_lower + (a_34 << 17) + (a_5 << 17) - a_upper - a_5; 
    //std::cout << "51 bit mod" << std::endl;

    // These loops run for 2-3 times max
    while (int64_t(v_reduced) < 0) {
        v_reduced += SOL_NTT_PRIME_51;
        //std::cout << "adds ";
    }
    //std::cout<<std::endl;

    while (v_reduced >= SOL_NTT_PRIME_51) {
        //std::cout << "subs ";
        v_reduced -= SOL_NTT_PRIME_51;
    }
    //std::cout<<std::endl;

    return v_reduced;
}

// prime modulus is 2^51 - 2^17 + 1
inline uint64_t fast_mod_solinas_51_max_64(uint64_t v) {
    // Typically input would be expected to be 102 bits long and we would split it into 6x 17-bit chunks a_0, a_1, ..., a_5 {MSB}
    // But we will only get inputs at most 64 bits long, a_4, a_5 = 0, and a_3 is only 13 bits instead of 17
    if (v < SOL_NTT_PRIME_51) {return v;}
    //if (v == SOL_NTT_PRIME_51) return 0;
    else if (v < (SOL_NTT_PRIME_51 << 1)) {return v - SOL_NTT_PRIME_51;}
    //else if (v < (3 * SOL_NTT_PRIME_51)) {return v - (SOL_NTT_PRIME_51 << 1);}
    
    uint64_t mask_51 = (1ULL << 51) - 1;

    uint64_t a_lower = v & mask_51;
    uint64_t a_upper = v >> 51;

    uint64_t v_reduced = a_lower + (a_upper << 17) - a_upper; 
    //std::cout << "51 bit mod" << std::endl;

    // These loops run for 2-3 times max
    while (int64_t(v_reduced) < 0) {
        v_reduced += SOL_NTT_PRIME_51;
        //std::cout << "adds ";
    }
    //std::cout<<std::endl;

    while (v_reduced >= SOL_NTT_PRIME_51) {
        //std::cout << "subs ";
        v_reduced -= SOL_NTT_PRIME_51;
    }
    //std::cout<<std::endl;

    return v_reduced;
}

// prime modulus is (1<<54) - (1<<24) + 1 [18014398492704769ULL]
inline uint64_t fast_mod_solinas_54(absu128 v) {
    // Input would be expected to be 108 bits long and we would split it into 18x 6-bit chunks a_0, a_1, ..., a_17 {MSB}
    if (v < SOL_NTT_PRIME_54) return static_cast<uint64_t>(v);
    if (v == SOL_NTT_PRIME_54) return 0;
    
    uint64_t high64 = absl::Uint128High64(v);
    uint64_t low64 = absl::Uint128Low64(v);

    uint64_t mask_54 = (1ULL << 54) - 1;
    uint64_t mask_30 = (1ULL << 30) - 1;

    uint64_t a_lower = low64 & mask_54;
    uint64_t a_upper = (low64 >> 54) | (high64 << 10);
    uint64_t a_9_13 = a_upper & mask_30;
    uint64_t a_14_17 = a_upper >> 30;

    uint64_t v_reduced = a_lower - a_upper - a_14_17 + (a_9_13 << 24) + (a_14_17 << 24); 
    //std::cout << "54 bit mod" << std::endl;

    // These loops run for 2-3 times max
    while (int64_t(v_reduced) < 0) {
        v_reduced += SOL_NTT_PRIME_54;
        //std::cout << "adds ";
    }
    //std::cout<<std::endl;

    while (v_reduced >= SOL_NTT_PRIME_54) {
        //std::cout << "subs ";
        v_reduced -= SOL_NTT_PRIME_54;
    }
    //std::cout<<std::endl;

    return v_reduced;
}

// prime modulus is 2^54 - 2^24 + 1
inline uint64_t fast_mod_solinas_54_max_64(uint64_t v) {
    // Typically input would be expected to be 108 bits long and we would split it into 18x 17-bit chunks a_0, a_1, ..., a_17 {MSB}
    // But we will only get inputs at most 64 bits long, a_11, ..., a_17 = 0, and a_10 is only 4 bits instead of 6
    if (v < SOL_NTT_PRIME_54) {return v;}
    //if (v == SOL_NTT_PRIME_54) return 0;
    else if (v < (SOL_NTT_PRIME_54 << 1)) {return v - SOL_NTT_PRIME_54;}
    //else if (v < (3 * SOL_NTT_PRIME_54)) {return v - (SOL_NTT_PRIME_54 << 1);}
    
    uint64_t mask_54 = (1ULL << 54) - 1;

    uint64_t a_lower = v & mask_54;
    uint64_t a_upper = v >> 54;

    uint64_t v_reduced = a_lower - a_upper + (a_upper << 24); 
    //std::cout << "54 bit mod" << std::endl;

    // These loops run for 2-3 times max
    
    while (int64_t(v_reduced) < 0) {
        v_reduced += SOL_NTT_PRIME_54;
        //std::cout << "adds ";
    }
    //std::cout<<std::endl;

    while (v_reduced >= SOL_NTT_PRIME_54) {
        //std::cout << "subs ";
        v_reduced -= SOL_NTT_PRIME_54;
    }
    //std::cout<<std::endl;
    

    return v_reduced;
}

void test_mod_solinas_30() {
    absu128 modl = static_cast<absu128>((1ULL << 30) - (1ULL << 25) + (1ULL << 15) + 1);
    absu128 modl_double = modl * modl;
    std::random_device rd;
    std::mt19937_64 e2(rd());
    std::uniform_int_distribution<long long int> dist(std::llround(std::pow(2,61)), std::llround(std::pow(2,62)));

    std::cout << "Testing Solinas 30 with max 64" << std::endl;
    // Generate random values of at most 64 bits
    for (int i = 0; i < 20000; i++){
        uint64_t v = static_cast<uint64_t>(static_cast<absu128>(dist(e2)) % modl_double);
        //std::cout << "v: " << v << std::endl;
        uint64_t v_reduced = fast_mod_solinas_30_max_64(v);
        if (v_reduced != static_cast<uint64_t>(static_cast<absu128>(v) % modl)) {
            std::cout << "mod_solinas_30 failed" << std::endl; 
            //std::cout << "mod_solinas_30 failed: " << (static_cast<absu128>(v) % modl) << " " << v_reduced << " " << modl << std::endl;
        }
    }

    std::cout << "Testing Solinas 30 with full 128 size" << std::endl;
    // Generate random values of at most double bits
    for (int i = 0; i < 20000; i++){
        uint64_t v1 = dist(e2);
        uint64_t v2 = dist(e2);
        absu128 v = absl::MakeUint128(v1, v2);
        v %= modl_double;
        //std::cout << "v: " << v << std::endl;
        uint64_t v_reduced = fast_mod_solinas_30(v);
        if (v_reduced != static_cast<uint64_t>(static_cast<absu128>(v) % modl)) {
            std::cout << "mod_solinas_30 failed" << std::endl;
            //std::cout << "mod_solinas_30 failed: " << (static_cast<absu128>(v) % modl) << " " << v_reduced << " " << modl << std::endl;
        }
    }
}

void test_mod_solinas_32() {
    absu128 modl = static_cast<absu128>((1ULL << 32) - (1ULL << 24) + (1ULL << 16) + 1);
    absu128 modl_double = modl * modl;
    std::random_device rd;
    std::mt19937_64 e2(rd());
    std::uniform_int_distribution<long long int> dist(std::llround(std::pow(2,61)), std::llround(std::pow(2,62)));

    std::cout << "Testing Solinas 32 with max 64" << std::endl;
    // Generate random values of at most 64 bits
    for (int i = 0; i < 20000; i++){
        uint64_t v = static_cast<uint64_t>(static_cast<absu128>(dist(e2)) % modl_double);
        //std::cout << "v: " << v << std::endl;
        uint64_t v_reduced = fast_mod_solinas_32_max_64(v);
        if (v_reduced != static_cast<uint64_t>(static_cast<absu128>(v) % modl)) {
            std::cout << "mod_solinas_32 failed" << std::endl; 
            //std::cout << "mod_solinas_32 failed: " << (static_cast<absu128>(v) % modl) << " " << v_reduced << " " << modl << std::endl;
        }
    }

    std::cout << "Testing Solinas 32 with full 128 size" << std::endl;
    // Generate random values of at most double bits
    for (int i = 0; i < 20000; i++){
        uint64_t v1 = dist(e2);
        uint64_t v2 = dist(e2);
        absu128 v = absl::MakeUint128(v1, v2);
        v %= modl_double;
        //std::cout << "v: " << v << std::endl;
        uint64_t v_reduced = fast_mod_solinas_32(v);
        if (v_reduced != static_cast<uint64_t>(static_cast<absu128>(v) % modl)) {
            std::cout << "mod_solinas_32 failed" << std::endl;
            //std::cout << "mod_solinas_32 failed: " << (static_cast<absu128>(v) % modl) << " " << v_reduced << " " << modl << std::endl;
        }
    }
}

void test_mod_solinas_40() {
    absu128 modl = static_cast<absu128>((1ULL << 40) - (1ULL << 30) + (1ULL << 20) + 1);
    absu128 modl_double = modl * modl;
    std::random_device rd;
    std::mt19937_64 e2(rd());
    std::uniform_int_distribution<long long int> dist(std::llround(std::pow(2,61)), std::llround(std::pow(2,62)));

    std::cout << "Testing Solinas 40 with max 64" << std::endl;
    // Generate random values of at most 64 bits
    for (int i = 0; i < 20000; i++){
        uint64_t v = static_cast<uint64_t>(static_cast<absu128>(dist(e2)) % modl_double);
        //std::cout << "v: " << v << std::endl;
        uint64_t v_reduced = fast_mod_solinas_40_max_64(v);
        if (v_reduced != static_cast<uint64_t>(static_cast<absu128>(v) % modl)) {
            std::cout << "mod_solinas_40 failed" << std::endl; 
            //std::cout << "mod_solinas_40 failed: " << (static_cast<absu128>(v) % modl) << " " << v_reduced << " " << modl << std::endl;
        }
    }

    std::cout << "Testing Solinas 40 with full 128 size" << std::endl;
    // Generate random values of at most double bits
    for (int i = 0; i < 20000; i++){
        uint64_t v1 = dist(e2);
        uint64_t v2 = dist(e2);
        absu128 v = absl::MakeUint128(v1, v2);
        v %= modl_double;
        //std::cout << "v: " << v << std::endl;
        uint64_t v_reduced = fast_mod_solinas_40(v);
        if (v_reduced != static_cast<uint64_t>(static_cast<absu128>(v) % modl)) {
            std::cout << "mod_solinas_40 failed" << std::endl;
            //std::cout << "mod_solinas_40 failed: " << (static_cast<absu128>(v) % modl) << " " << v_reduced << " " << modl << std::endl;
        }
    }
}

void test_mod_solinas_45() {
    absu128 modl = static_cast<absu128>((1ULL << 45) - (1ULL << 9) + 1);
    absu128 modl_double = modl * modl;
    std::random_device rd;
    std::mt19937_64 e2(rd());
    std::uniform_int_distribution<long long int> dist(std::llround(std::pow(2,61)), std::llround(std::pow(2,62)));

    std::cout << "Testing Solinas 45 with max 64" << std::endl;
    // Generate random values of at most 64 bits
    for (int i = 0; i < 20000; i++){
        uint64_t v = static_cast<uint64_t>(static_cast<absu128>(dist(e2)) % modl_double);
        //std::cout << "v: " << v << std::endl;
        uint64_t v_reduced = fast_mod_solinas_45_max_64(v);
        if (v_reduced != static_cast<uint64_t>(static_cast<absu128>(v) % modl)) {
            std::cout << "mod_solinas_45 failed" << std::endl; 
            //std::cout << "mod_solinas_45 failed: " << (static_cast<absu128>(v) % modl) << " " << v_reduced << " " << modl << std::endl;
        }
    }

    std::cout << "Testing Solinas 45 with full 128 size" << std::endl;
    // Generate random values of at most double bits
    for (int i = 0; i < 20000; i++){
        uint64_t v1 = dist(e2);
        uint64_t v2 = dist(e2);
        absu128 v = absl::MakeUint128(v1, v2);
        v %= modl_double;
        //std::cout << "v: " << v << std::endl;
        uint64_t v_reduced = fast_mod_solinas_45(v);
        if (v_reduced != static_cast<uint64_t>(static_cast<absu128>(v) % modl)) {
            std::cout << "mod_solinas_45 failed" << std::endl;
            //std::cout << "mod_solinas_45 failed: " << (static_cast<absu128>(v) % modl) << " " << v_reduced << " " << modl << std::endl;
        }
    }
}

void test_mod_solinas_51() {
    absu128 modl = static_cast<absu128>((1ULL << 51) - (1ULL << 17) + 1);
    absu128 modl_double = modl * modl;
    std::random_device rd;
    std::mt19937_64 e2(rd());
    std::uniform_int_distribution<long long int> dist(std::llround(std::pow(2,61)), std::llround(std::pow(2,62)));

    std::cout << "Testing Solinas 51 with max 64" << std::endl;
    // Generate random values of at most 64 bits
    for (int i = 0; i < 20000; i++){
        uint64_t v = static_cast<uint64_t>(static_cast<absu128>(dist(e2)) % modl_double);
        //std::cout << "v: " << v << std::endl;
        uint64_t v_reduced = fast_mod_solinas_51_max_64(v);
        if (v_reduced != static_cast<uint64_t>(static_cast<absu128>(v) % modl)) {
            std::cout << "mod_solinas_51 failed" << std::endl; 
            //std::cout << "mod_solinas_51 failed: " << (static_cast<absu128>(v) % modl) << " " << v_reduced << " " << modl << std::endl;
        }
    }

    std::cout << "Testing Solinas 51 with full 128 size" << std::endl;
    // Generate random values of at most double bits
    for (int i = 0; i < 20000; i++){
        uint64_t v1 = dist(e2);
        uint64_t v2 = dist(e2);
        absu128 v = absl::MakeUint128(v1, v2);
        v %= modl_double;
        //std::cout << "v: " << v << std::endl;
        uint64_t v_reduced = fast_mod_solinas_51(v);
        if (v_reduced != static_cast<uint64_t>(static_cast<absu128>(v) % modl)) {
            std::cout << "mod_solinas_51 failed" << std::endl;
            //std::cout << "mod_solinas_51 failed: " << (static_cast<absu128>(v) % modl) << " " << v_reduced << " " << modl << std::endl;
        }
    }
}

void test_mod_solinas_54() {
    absu128 modl = static_cast<absu128>((1ULL << 54) - (1ULL << 24) + 1);
    absu128 modl_double = modl * modl;
    std::random_device rd;
    std::mt19937_64 e2(rd());
    std::uniform_int_distribution<long long int> dist(std::llround(std::pow(2,61)), std::llround(std::pow(2,62)));

    std::cout << "Testing Solinas 54 with max 64" << std::endl;
    // Generate random values of at most 64 bits
    for (int i = 0; i < 20000; i++){
        uint64_t v = static_cast<uint64_t>(static_cast<absu128>(dist(e2)) % modl_double);
        //std::cout << "v: " << v << std::endl;
        uint64_t v_reduced = fast_mod_solinas_54_max_64(v);
        if (v_reduced != static_cast<uint64_t>(static_cast<absu128>(v) % modl)) {
            std::cout << "mod_solinas_54 failed" << std::endl; 
            //std::cout << "mod_solinas_54 failed: " << (static_cast<absu128>(v) % modl) << " " << v_reduced << " " << modl << std::endl;
        }
    }

    std::cout << "Testing Solinas 54 with full 128 size" << std::endl;
    // Generate random values of at most double bits
    for (int i = 0; i < 20000; i++){
        uint64_t v1 = dist(e2);
        uint64_t v2 = dist(e2);
        absu128 v = absl::MakeUint128(v1, v2);
        v %= modl_double;
        //std::cout << "v: " << v << std::endl;
        uint64_t v_reduced = fast_mod_solinas_54(v);
        if (v_reduced != static_cast<uint64_t>(static_cast<absu128>(v) % modl)) {
            std::cout << "mod_solinas_54 failed" << std::endl;
            //std::cout << "mod_solinas_54 failed: " << (static_cast<absu128>(v) % modl) << " " << v_reduced << " " << modl << std::endl;
        }
    }
}