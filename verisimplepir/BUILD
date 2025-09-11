# Core library

cc_library(
    name = "verisimplepir_lib",
    srcs = glob([
        "src/lib/**/*.cpp",
        "src/lib/**/*.cc",
    ]),
    hdrs = glob([
        "src/lib/**/*.h",
        "src/lib/**/*.hpp",
    ]),
    copts = [
        "-std=c++17",
        "-fPIC",
        "-O3",
        "-maes",
        "-msse4",
        "-Wall",
        "-fno-omit-frame-pointer",
        "-Wno-ignored-attributes",
    ],
    includes = [
        "src/lib",
        "test",
    ],
    linkopts = [
        "-lssl",
        "-lcrypto",
        "-lgsl",
        "-lgslcblas",
        "-lrt",
        "-L/usr/lib/x86_64-linux-gnu",
        "-L/usr/local/lib",
        "-Wl,--export-dynamic,-z,defs",
    ],
    visibility = ["//visibility:public"],
    deps = ["@com_google_absl//absl/numeric:int128"],
)