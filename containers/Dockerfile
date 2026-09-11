# Reproducible toolchain for the scintillator screening pipeline.
#
# Version policy
# --------------
# Pinned to Geant4 11.4.2 (released 2026-06-17), the latest patch release of the
# 11.4 series. The development machine reports 11.4.0 via geant4-config, but its
# G4Version.hh carries the tag "geant4-11-04-beta-01" -- i.e. it may be a beta
# build. Publishing results from a beta is not defensible, so the container pins
# a released version and the local installation should be brought into line with
# it. The equivalence test (same config + same seed -> same output checksum)
# between host and container is what proves the two agree.
#
# ROOT is deliberately absent. Geant4's own G4AnalysisManager writes .root files
# without a ROOT installation, which keeps this image to one long build instead
# of two.
#
# Build:  docker build -t scint:11.4.2 -f containers/Dockerfile .
# Run:    docker run --rm -v "$PWD":/work -w /work scint:11.4.2 <command>

FROM ubuntu:24.04

ARG GEANT4_VERSION=11.4.2
ARG BUILD_JOBS=8
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        cmake \
        curl \
        expat \
        libexpat1-dev \
        libxerces-c-dev \
        python3 \
        python3-pip \
        python3-venv \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# ---- Geant4 -----------------------------------------------------------------
# GEANT4_INSTALL_DATA downloads the physics datasets at configure time so the
# image is self-contained and needs no network at run time.
RUN set -eux; \
    mkdir -p /tmp/g4; cd /tmp/g4; \
    curl -fsSL -o geant4.tar.gz \
      "https://gitlab.cern.ch/geant4/geant4/-/archive/v${GEANT4_VERSION}/geant4-v${GEANT4_VERSION}.tar.gz"; \
    tar xzf geant4.tar.gz; \
    cmake -S "geant4-v${GEANT4_VERSION}" -B build \
        -DCMAKE_INSTALL_PREFIX=/opt/geant4 \
        -DCMAKE_BUILD_TYPE=Release \
        -DGEANT4_INSTALL_DATA=ON \
        -DGEANT4_BUILD_MULTITHREADED=ON \
        -DGEANT4_USE_GDML=ON \
        -DGEANT4_USE_SYSTEM_EXPAT=ON; \
    cmake --build build -j"${BUILD_JOBS}"; \
    cmake --install build; \
    rm -rf /tmp/g4

ENV PATH=/opt/geant4/bin:$PATH
ENV CMAKE_PREFIX_PATH=/opt/geant4:$CMAKE_PREFIX_PATH

# Geant4 dataset locations. geant4.sh does not survive a non-interactive shell
# reliably, so they are baked in; scint/attenuation.py sets the same variables
# when running against a host installation.
ENV G4LEDATA=/opt/geant4/share/Geant4/data/G4EMLOW8.7 \
    G4LEVELGAMMADATA=/opt/geant4/share/Geant4/data/PhotonEvaporation6.1 \
    G4RADIOACTIVEDATA=/opt/geant4/share/Geant4/data/RadioactiveDecay6.1.2 \
    G4PARTICLEXSDATA=/opt/geant4/share/Geant4/data/G4PARTICLEXS4.1 \
    G4PIIDATA=/opt/geant4/share/Geant4/data/G4PII1.3 \
    G4REALSURFACEDATA=/opt/geant4/share/Geant4/data/RealSurface2.2 \
    G4SAIDXSDATA=/opt/geant4/share/Geant4/data/G4SAIDDATA2.0 \
    G4ABLADATA=/opt/geant4/share/Geant4/data/G4ABLA3.3 \
    G4INCLDATA=/opt/geant4/share/Geant4/data/G4INCL1.2 \
    G4ENSDFSTATEDATA=/opt/geant4/share/Geant4/data/G4ENSDFSTATE3.0 \
    G4NEUTRONHPDATA=/opt/geant4/share/Geant4/data/G4NDL4.7.1

# ---- Python -----------------------------------------------------------------
COPY requirements.txt /tmp/requirements.txt
RUN python3 -m venv /opt/venv \
    && /opt/venv/bin/pip install --no-cache-dir --upgrade pip \
    && /opt/venv/bin/pip install --no-cache-dir -r /tmp/requirements.txt
ENV PATH=/opt/venv/bin:$PATH

WORKDIR /work
CMD ["/bin/bash"]
