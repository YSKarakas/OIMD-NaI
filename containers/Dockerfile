# Reproducible toolchain for the scintillator screening pipeline.
#
# Version policy
# --------------
# Pinned to Geant4 11.4.2 (tag geant4-11-04-patch-02, 12 June 2026).
#
# This pin exists for a concrete reason. The development machine's installation
# reports "11.4.0" from geant4-config, but its G4Version.hh carries the tag
# "geant4-11-04-beta-01" with G4Date 26-June-2025, and its libraries are dated
# September 2025 -- before 11.4.0 was released on 5 December 2025. It is a beta.
# geant4-config cannot reveal this, because G4VERSION_NUMBER is already 1140 in
# the beta. Since the 11.4 cycle refactored G4OpBoundaryProcess and moved
# scintillation to a model-based framework -- exactly the code the light-yield
# geometry study depends on -- results must come from a released version.
#
# ROOT is deliberately absent. Geant4's own G4AnalysisManager writes .root files
# without a ROOT installation, which keeps this to one long build instead of two.
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
        libexpat1-dev \
        libxerces-c-dev \
        python3 \
        python3-pip \
        python3-venv \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# ---- Geant4 -----------------------------------------------------------------
# GEANT4_INSTALL_DATA fetches the physics datasets at configure time so the
# finished image needs no network. Visualisation drivers are left off: nothing
# in this pipeline renders, and they would only lengthen the build.
RUN set -eux; \
    mkdir -p /tmp/g4 && cd /tmp/g4; \
    curl -fsSL -o geant4.tar.gz \
      "https://github.com/Geant4/geant4/archive/refs/tags/v${GEANT4_VERSION}.tar.gz"; \
    tar xzf geant4.tar.gz; \
    cmake -S "geant4-${GEANT4_VERSION}" -B build \
        -DCMAKE_INSTALL_PREFIX=/opt/geant4 \
        -DCMAKE_BUILD_TYPE=Release \
        -DGEANT4_INSTALL_DATA=ON \
        -DGEANT4_BUILD_MULTITHREADED=ON \
        -DGEANT4_USE_GDML=ON \
        -DGEANT4_USE_SYSTEM_EXPAT=ON; \
    cmake --build build -j"${BUILD_JOBS}"; \
    cmake --install build; \
    cd / && rm -rf /tmp/g4

ENV PATH=/opt/geant4/bin:$PATH
ENV CMAKE_PREFIX_PATH=/opt/geant4:$CMAKE_PREFIX_PATH

# Capture the dataset locations from the installation's own geant4.sh rather
# than hard-coding version-numbered paths, which would go stale on an upgrade.
RUN bash -c 'set -a; . /opt/geant4/bin/geant4.sh; set +a; \
             env | grep -E "^G4[A-Z0-9]+DATA=" | sort > /etc/geant4.env' \
    && cat /etc/geant4.env

# ---- Python -----------------------------------------------------------------
COPY requirements.txt /tmp/requirements.txt
RUN python3 -m venv /opt/venv \
    && /opt/venv/bin/pip install --no-cache-dir --upgrade pip \
    && /opt/venv/bin/pip install --no-cache-dir -r /tmp/requirements.txt
ENV PATH=/opt/venv/bin:$PATH

COPY containers/entrypoint.sh /usr/local/bin/entrypoint.sh
WORKDIR /work
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["bash"]
