#!/bin/sh
# Run the simulation campaign inside the pinned container, against the Geant4
# release it carries, with the repository mounted so that every run lands in
# runs/ exactly as a host run would.
#
# Two things the container cannot do for itself are done here first. It has no
# git, so the state of the working tree -- the one fact env.json exists to
# record -- is captured on the host and passed through the environment, and
# scint.registry records it with its provenance. And the binary is the one
# built inside the container (build/container), never the host build: the host
# toolkit is a pre-release whose optical transport differs from the release.
#
#   sh scripts/run_in_container.sh --set all --events 2000 --jobs 2
set -e
cd "$(dirname "$0")/.."

IMAGE=scint:11.4.2
[ -x build/container/scint_optical ] || {
    echo "build/container/scint_optical missing; build it in the container first:" >&2
    echo "  docker run --rm -v \"\$PWD\":/work -w /work $IMAGE sh -c 'cmake -S sim -B build/container && cmake --build build/container -j4'" >&2
    exit 1
}

GIT_COMMIT=$(git rev-parse HEAD)
GIT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
if [ -n "$(git status --porcelain)" ]; then GIT_DIRTY=1; else GIT_DIRTY=0; fi

# The image's entrypoint sources /etc/geant4.env, but that file is empty in the
# image as built, so neither the library path nor the dataset variables are set
# when a command is run directly. geant4.sh from the installation itself is
# sourced here instead; it is the file /etc/geant4.env was meant to be derived
# from. SCINT_G4DATA_BIN points the version record at the container-built tool.
exec docker run --rm -v "$PWD":/work -w /work \
    -e SCINT_GIT_COMMIT="$GIT_COMMIT" -e SCINT_GIT_BRANCH="$GIT_BRANCH" -e SCINT_GIT_DIRTY="$GIT_DIRTY" \
    -e SCINT_G4DATA_BIN=/work/build/container/g4data \
    "$IMAGE" bash -c '. /opt/geant4/bin/geant4.sh && exec /opt/venv/bin/python3 scripts/run_gates.py \
        --binary build/container/scint_optical "$@"' -- "$@"
