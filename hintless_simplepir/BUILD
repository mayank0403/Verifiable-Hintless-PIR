# Copyright 2024 Google LLC.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# An implementation of the Hintless SimplePIR scheme that uses RLWE-based
# homomorphic encryption to remove the database-dependent "hint" in SimplePIR.

load("@io_bazel_rules_go//proto:def.bzl", "go_proto_library")
load("@rules_cc//cc:defs.bzl", "cc_library")
load("@rules_proto//proto:defs.bzl", "proto_library")

package(
    default_visibility = ["//visibility:public"],
)

licenses(["notice"])

# SimplePIR implementation.
cc_library(
    name = "simplepir",
    hdrs = [
        "simplepir.h",
    ],
    deps = [
        "//lwe:encode",
        "//lwe:lwe_symmetric_encryption",
        "//lwe:sample_error",
        "//lwe:types",
        "@com_github_google_shell-encryption//shell_encryption:statusor_fork",
        "@com_github_google_shell-encryption//shell_encryption/prng:single_thread_hkdf_prng",
        "@com_gitlab_libeigen-eigen//:eigen3",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
    ],
)

# Hintless SimplePIR request and response types.
proto_library(
    name = "serialization_proto",
    srcs = ["serialization.proto"],
    deps = [
        "//linpir:serialization_proto",
        "@com_github_google_shell-encryption//shell_encryption/rns:serialization_proto",
    ],
)

cc_proto_library(
    name = "serialization_cc_proto",
    deps = [":serialization_proto"],
)

# Hintless SimplePIR parameters.
cc_library(
    name = "parameters",
    hdrs = ["parameters.h"],
    deps = [
        "//linpir:parameters",
        "//lwe:types",
        "@com_github_google_shell-encryption//shell_encryption:serialization_cc_proto",
    ],
)

# Test helpers.
cc_library(
    name = "testing",
    hdrs = ["testing.h"],
    deps = [
        ":parameters",
        ":utils",
        "@com_google_absl//absl/random",
    ],
)

# Utilities
cc_library(
    name = "utils",
    hdrs = ["utils.h"],
    deps = [
        ":parameters",
        ":serialization_cc_proto",
        "//lwe:types",
        "@com_gitlab_libeigen-eigen//:eigen3",
        "@com_google_absl//absl/strings:string_view",
        "@com_google_absl//absl/types:span",
        "//verisimplepir:verisimplepir_lib",
    ],
)

# Hintless SimplePIR database and the preprocessed hint states
cc_library(
    name = "database",
    srcs = ["database.cc"],
    hdrs = ["database.h"],
    deps = [
        ":parameters",
        ":utils",
        "//lwe:types",
        "@com_gitlab_libeigen-eigen//:eigen3",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings:string_view",
        "@com_google_absl//absl/types:span",
    ],
)

# highway-based matrix-vector multiplication.
cc_library(
    name = "inner_product_hwy",
    srcs = ["inner_product_hwy.cc"],
    hdrs = ["inner_product_hwy.h"],
    deps = [
        "//lwe:types",
        "@com_github_google_highway//:hwy",
        "@com_google_absl//absl/base:core_headers",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/types:span",
        "//verisimplepir:verisimplepir_lib",
    ],
)

# Build and test :inner_product_hwy on platforms without vector intrinsics.
cc_library(
    name = "inner_product_hwy_scalar",
    srcs = ["inner_product_hwy.cc"],
    hdrs = ["inner_product_hwy.h"],
    local_defines = ["HWY_COMPILE_ONLY_SCALAR"],
    deps = [
        "//lwe:types",
        "@com_github_google_highway//:hwy",
        "@com_google_absl//absl/base:core_headers",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/types:span",
    ],
)

# highway-based database implementation.
cc_library(
    name = "database_hwy",
    srcs = ["database_hwy.cc"],
    hdrs = ["database_hwy.h"],
    deps = [
        ":inner_product_hwy",
        ":parameters",
        ":utils",
        "//lwe:types",
        "@com_github_google_shell-encryption//shell_encryption:statusor_fork",
        "@com_gitlab_libeigen-eigen//:eigen3",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings:string_view",
        "@com_google_absl//absl/types:span",
    ],
)

# Hintless SimplePIR server.
cc_library(
    name = "server",
    srcs = ["server.cc"],
    hdrs = ["server.h"],
    deps = [
        ":database_hwy",
        ":parameters",
        ":serialization_cc_proto",
        ":utils",
        "//linpir:database",
        "//linpir:server",
        "//lwe:lwe_symmetric_encryption",
        "//lwe:types",
        "@com_github_google_shell-encryption//shell_encryption:montgomery",
        "@com_github_google_shell-encryption//shell_encryption:statusor_fork",
        "@com_github_google_shell-encryption//shell_encryption/prng:single_thread_chacha_prng",
        "@com_github_google_shell-encryption//shell_encryption/prng:single_thread_hkdf_prng",
        "@com_github_google_shell-encryption//shell_encryption/rns:rns_context",
        "@com_gitlab_libeigen-eigen//:eigen3",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/strings:string_view",
    ],
)

# Hintless SimplePIR client.
cc_library(
    name = "client",
    srcs = ["client.cc"],
    hdrs = ["client.h"],
    deps = [
        ":parameters",
        ":serialization_cc_proto",
        ":utils",
        "//linpir:client",
        "//lwe:encode",
        "//lwe:lwe_symmetric_encryption",
        "//lwe:types",
        "//lwe:matrix_conversion",
        "@com_github_google_shell-encryption//shell_encryption:int256",
        "@com_github_google_shell-encryption//shell_encryption:montgomery",
        "@com_github_google_shell-encryption//shell_encryption:statusor_fork",
        "@com_github_google_shell-encryption//shell_encryption/prng",
        "@com_github_google_shell-encryption//shell_encryption/prng:single_thread_chacha_prng",
        "@com_github_google_shell-encryption//shell_encryption/prng:single_thread_hkdf_prng",
        "@com_github_google_shell-encryption//shell_encryption/rns:crt_interpolation",
        "@com_github_google_shell-encryption//shell_encryption/rns:rns_context",
        "@com_github_google_shell-encryption//shell_encryption/rns:rns_modulus",
        "@com_gitlab_libeigen-eigen//:eigen3",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/strings:string_view",
        "//verisimplepir:verisimplepir_lib",
    ],
)

# Bridge connection VeriSimplePIR to HintlessPIR.
cc_test(
    name = "new_pir_test",
    srcs = ["new_pir_test.cc"],
    deps = [
        ":client",
        ":database_hwy",
        ":parameters",
        ":server",
        "//linpir:parameters",
        ":simplepir",
        "//lwe:types",
        "//lwe:matrix_conversion",
        "@com_github_google_googletest//:gtest_main",
        "@com_github_google_shell-encryption//shell_encryption/prng:single_thread_hkdf_prng",
        "@com_github_google_shell-encryption//shell_encryption/testing:matchers",
        "@com_github_google_shell-encryption//shell_encryption/testing:status_testing",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "//verisimplepir:verisimplepir_lib",
    ],
)