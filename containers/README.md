# Reproducible container

Pins Geant4 to a released version so that published results do not depend on the
state of any one developer's machine.

```bash
docker build -t scint:11.4.2 -f containers/Dockerfile .
docker run --rm -v "$PWD":/work -w /work scint:11.4.2 \
    bash -c 'cmake -S sim -B build/container && cmake --build build/container -j8'
```

The build takes on the order of an hour: Geant4 is compiled from source and the
physics datasets (~2 GB) are downloaded at configure time so the image needs no
network afterwards.

## Why no ROOT

Geant4's `G4AnalysisManager` writes `.root` files without a ROOT installation.
Adding ROOT would roughly double an already long build for no benefit to the
pipeline; interactive ROOT analysis happens on the host.

## Version note

`geant4-config --version` on the development machine reports `11.4.0`, but
`G4Version.hh` there carries the tag `geant4-11-04-beta-01`. Whether that
installation is genuinely a beta or a release with a stale tag is unresolved.
The container therefore pins **11.4.2** (released 2026-06-17), and the host
installation should be aligned with it before any results are published.

## Equivalence test

The container's claim is only worth what it is tested against. The check is that
the same configuration and the same RNG seed produce byte-identical outputs on
host and in container; `env.json` in each run records `platform.in_container`,
so the two can always be told apart after the fact.
