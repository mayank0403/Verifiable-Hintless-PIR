# Use --platform to force x86_64
FROM --platform=linux/amd64 ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install essential packages including all development libraries
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    clang \
    llvm \
    libssl-dev \
    pkg-config \
    curl \
    gpg \
    gnupg \
    apt-transport-https \
    ca-certificates \
    software-properties-common \
    wget \
    python3 \
    g++ \
    unzip \
    openssl \
    libgsl-dev \
    libgslcblas0 \
    && rm -rf /var/lib/apt/lists/*


# Install Microsoft's GSL
RUN git clone https://github.com/microsoft/GSL.git && \
    cd GSL && \
    cmake -B build . && \
    cmake --build build --target install
# Add Bazel repository and key
RUN curl -fsSL https://bazel.build/bazel-release.pub.gpg | gpg --dearmor > /usr/share/keyrings/bazel-archive-keyring.gpg && \
    echo "deb [arch=amd64 signed-by=/usr/share/keyrings/bazel-archive-keyring.gpg] https://storage.googleapis.com/bazel-apt stable jdk1.8" | tee /etc/apt/sources.list.d/bazel.list

# Install Bazel
RUN apt-get update && \
    apt-get install -y bazel

WORKDIR /workspace

COPY . .

CMD ["/bin/bash"]