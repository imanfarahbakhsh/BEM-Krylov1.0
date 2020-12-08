# hash:sha256:b0fd68353a368a69bf3c965ff2a30979d263a5795a2a8a6cf88a12fcb3669c89
FROM registry.codeocean.com/codeocean/ubuntu:18.04.5

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        gfortran \
    && rm -rf /var/lib/apt/lists/*
