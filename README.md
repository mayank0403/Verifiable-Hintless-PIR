# Verifiable PIR with small client storage
This library builds upon the code of HintlessPIR and VeriSimplePIR. 

### Dependencies
To build the project, do the following:
1. `sudo apt install make clang libssl-dev` (from VeriSimplePIR)
2. Install Bazel for your platform by following the instructions [here](https://bazel.build/install) (from HintlessPIR)

### Building and running
The main test file you need to run for verifiable PIR is `new_pir_test`. Run it using:
```
bazel run -c opt --cxxopt='-std=c++17' //hintless_simplepir:new_pir_test --cxxopt="-w" --copt="-w"
``` 

The above test might fail on non x86-64 CPUs like MacBook Pro. This issue comes from the use Highway SIMD library. In this case, you can disable this feature by running the following:
```
bazel run -c opt --cxxopt='-std=c++17' --cxxopt='-DHWY_COMPILE_ONLY_SCALAR=1' //hintless_simplepir:new_pir_test --cxxopt="-w" --copt="-w"
``` 

To debug (include `--cxxopt='-DHWY_COMPILE_ONLY_SCALAR=1'` on ARM or non-x86):
```
bazel test -c opt --test_output=all --compilation_mode=dbg --cxxopt='-std=c++17' //hintless_simplepir:new_pir_test --cxxopt="-w" --copt="-w"
```

## Parameters
This code uses many parameters across LWE, RLWE, SIS, arranging DB as a 2D matrix, etc. 
### Runtime
The following runtime parameters are set in `new_pir_test.cc`:

| Parameter | Example value for 512 MiB DB | Meaning |
| ----------- | ----------- | ----------- |
| `rows_db` | 32768 | Number of rows in 2D DB matrix |
| `cols_db` | 16384 | Number of columns in 2D DB matrix |
| `use_static_db` | true | Faster DB initialization and transpose with a static DB rather than sampling a random DB |
| `db_record_bit_size` | 8 | Size of each record in 2D DB matrix * number of shards (see below) |
| `db_stack_cells` | 8 | Size of each stack when doing Stacking. Use 1 if not used |
| `lwe_secret_dim` | 1408 | LWE dimension for online phase |
| `offline_lwe_secret_dim` | 2816 | LWE dimension for offline phase to compute Z<sup>T</sup> = D<sup>T</sup> * C<sup>T</sup> |
| `lwe_modulus_bit_size` | 32 | Bitsize of LWE modulus. If over 32, increase size of `Elem` in `mat.h` |
| `lwe_plaintext_bit_size` | 8 | Fixed to 8 (1-byte). You can also probably use 1, 2, 4. Anything > 8 won't work |
| `lwe_error_variance` | 8 | Variance of LWE error sampled from Centered Binomial distribution. LWE secret key is ternary |
| `log_n` | 12 | RLWE ring dimension |
| `error_variance` | 8 | Variance of RLWE error sampled from Centered Binomial. RLWE secret key also from Centered Binomial |
| `qs` | see code | A set of RLWE ciphertext moduli specifying each RNS ciphertext limb |
| `ts` | see code | A set of RLWE plaintext moduli specifying each RNS plaintext limb. We typically use a single prime here for efficiency |
| `rows_per_block` | see code | An artifact of our RLWE packing scheme for hint matrix H. Typically half of RLWE ring dimension |

For other paramters, refer to our paper and the code. For some common params, you can refer to HintlessPIR and VeriSimplePIR repos for another explanation.


Default value is for 512 MiB database from our paper. The code lays out the DB as a 2D matrix with 1-byte records. To represent databases with larger entries, you have 2 options:
1. Stacking: For each PIR request, our code fetches an entire column of the 2D DB matrix. Therefore, to reduce the number of such requests, one can stack each entry of original DB as a stack of 1-byte entries in each column. Each column can have multiple of these stacks. But be careful to never have a stack overflow a column. This is important of security because one would two PIR requests when an overflow happens which leaks private information about the column's identity to the server. Benefits of Stacking are that the client download and server computation are reduced compared to Sharding (see next bullet). We use Stacking in our experiments in the paper.
2. Starding: Another approach to deal with large entries in the DB is to do Sharding. In Sharding, the database is split into multiple databases with each containing only a 1-byte entries. If original DB has 2-byte entries, then we can create two shards with 1-byte entries each. These act as two separate 2D DBs now. The benefit of sharding is reduced client upload because same query can be used for both databases. We do not use Sharding in our experiments in the paper.

### LWE ciphertext moduli and offline LWE vs online LWE
There are 2 instances of LWE: one is used only in the offline phase to compute Z, while the other dominates rest of the computation. Since the magnitude of values in Z is larger than D because Z = C*D, where C has binary values, we need to use LWE with larger plaintext space to compute Z homomorphically. This difference in plaintext space is reflected in the difference between the variables `Delta` and `Delta_preproc` in `verisimplepir/src/lib/pir/lhe.h`. The ciphertext space also differs across the two instances. This is why we have `offline_lwe_secret_dim` and `lwe_secret_dim` separately in the table above. The code uses NTT-friendly Solinas primes for ciphertext modulus and provides efficient modular reduction routines in `verisimplepir/src/lib/pir/utils.cpp`. These primes are specified in `verisimplepir/src/lib/pir/utils.h` as follows:

Online LWE modulus (the main LWE instance):
- `SOL_NTT_PRIME_BITS`: Bitsize of the prime
- `SOL_NTT_PRIME`: The prime modulus
- `LOG_Q`: Bitsize - 1

Offline LWE modulus (only used for computing Z in offline phase):
- `VSPIR_LHE_PRIME_BITS`: Bitsize of offline LWE prime
- `VSPIR_LHE_PRIME`: The prime modulus that defines ciphertext space for offline LHE
- `VSPIR_LHE_LOG_Q`: Bitsize - 1

There are several options for Solinas primes of differetn sizes included in `utils.h` (e.g. `SOL_NTT_PRIME_40`). Their correspdoning fast modular reductions functions are implemented in `utils.cpp`. If you want to use a different prime, you will need to implement a modular reduction function for it and add that function to the defines at top of `utils.cpp`.

### Other controls
You can use the following flags in `verisimplepir/src/lib/pir/utils.h` as well:
1. `SH_RUN`: If you don't want vPIR, and only want to run standard (semi-honest) PIR. The code will change the RLWE plaintext modulus in this case and also don't use prime modulus anymore - the code sets `SOL_NTT_PRIME` = 0 which reverts to using integer rings like HintlessPIR and VeriSimplePIR.
2. `BSGS`: Toggle between Baby-Step-Giant-Step based homomorphic matrix multiplication for query-level preprocessing (H*s computation).
3. `FAKE_RUN`: Since the preprocessing phase can take a while to run for large databases and has no impact on the online metrics, if you want to get online metrics, use this flag to make preprocessing go faster. Online metrics will stay true to a real run.

## Interpreting captured matrics
The code captures many useful metrics. The table explains each.
| Metric (unit) | Phase | Meaning |
| ----------- | ----------- | ----------- |
| Hints (MiB) | Global + client-specific preprocessing | Size of hints. We have 2 hints in offline phase - hint of decrypting responses of offline LWE (to compute Z), and the hint (SIS hash) of online LWE. Both are deleted at the end of the offline phase |
| Offline Up (KiB) | Global + client-specific preprocessing | Client upload in offline phase |
| Offline Down (KiB) | Global + client-specific preprocessing | Client download in offline phase |
| Global Prepr (s) | Global preprocessing | Server's global preprocessing time |
| Server Per-Client Prepr (s) | Client-specific preprocessing | Server's time for client-specific preprocessing |
| Client Local Prepr (s) | Global + client-specific preprocessing | Client's compute time for offline phase |
| Prep Up (Kib) | Query-level preprocessing | Client upload for query-level preprocessing |
| Prep Down (KiB) | Query-level preprocessing | Client download for query-level preprocessing |
| Client Prepa Pre Req (s) | Query-level preprocessing | Client's compute time for query-level preprocessing before talking to server |
| Server Prepa Comp (s) | Query-level preprocessing | Server's time for query-level preprecessing per client |
| Client Prepa Post Req (s) | Query-level preprocessing | Client's compute time for query-level preprocessing after getting response from server |
| Query: Client Req Gen (ms) | Online-only | Client's compute time to generate online PIR request |
| Query: Server Comp (s) | Online-only | Server's time to generate response to online PIR request |
| Query: Client Decryption (ms) | Online-only | Client's time to decrypt the server's response to online PIR request |
| Query: Client Verification (ms) | Online-only | Client's time to verify the online PIR response |
| Query Up (KiB) | Online-only | Client upload for online PIR request |
| Query Down (KiB) | Online-only | Client download for online PIR request |
| Online State (KiB) | Online-only + query-level preprocessing | Rough estimate of client's state during online and query-level preprocessing phases |
| Long Term State (KiB) | - | Persistent state stored by the client |

We use "Online-only" to denote online phase without query-level preprocessing in the "Phase" column above. Note that typically, query-level preprocessing (aka "prepare phase") is included in the online phase, not in the offline phase. Prior works don't preprocess this part and do it alongside online query. However, this phase can be done before the online phase starts. Please refer to our paper for more details. "Prepa": Prepare; "Prepr": Preprocess.

Note that the code doesn't do any communication. So all the timings reported are without the time it takes to communicate the message. It is easy to estimate time spent in communication using the upload and download figures.

## Communication vs Client Storage
A shorter and wider database gives lower communication (fewer rows means fewer RLWE ciphertexts in server's response), while the opposite leads to lower persistent state (only depends on the number of columns).

## Warnings
The code uses a loose bound for SIS hardness and can deem some parameters not hard by a small margin. To check SIS hardness in this case, use [lattice estimator tool](https://github.com/malb/lattice-estimator). 

This code is a proof-of-concept for an academic work. It should *NOT* be used as it is for production. The code comes with no guarantees about its security or correctness, and has not been vetted for the same.

## Out of scope of this code
1. Optimizations from HintlessPIR (section E.1) to reduce the server's response size are not implemented. This includes modulus switching and preprocessing of "a" component of RLWE response.
2. vPIR without honest hint assumption. You can adjust parameters (like increasing primes to allow for larger norm of `Z`) to support this. 
3. cvPIR.

## Acknowledgement
This code is built over two separate repositories:
1. [HintlessPIR](https://github.com/google/hintless_pir/): Efficient standard PIR (semi-honest) without persistent client storage. Commit `4be2ae8`.
2. [VeriSimplePIR](https://github.com/leodec/VeriSimplePIR): Efficient Verifiable PIR (vPIR) with large client storage. Commit `3643bb7`.
