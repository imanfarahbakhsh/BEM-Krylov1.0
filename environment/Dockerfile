# hash:sha256:372c294565881382c0bf0923c4b4b86fd4bc954b2710692608a8c376ebfb73ae
FROM registry.codeocean.com/codeocean/ubuntu:18.04.5

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        gfortran=4:7.4.0-1ubuntu2.3 \
    && rm -rf /var/lib/apt/lists/*
