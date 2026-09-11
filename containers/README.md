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

## Equivalence: what is and is not reproducible

Measured on 2026-09-12, 20 000 events, 662 keV into 3" x 3" NaI(Tl), identical
configuration and seed.

| Comparison | Result |
|---|---|
| Same build, same seed, run twice | **byte-identical** (host and container both) |
| Host beta vs container release, event by event | 61 % of events differ |
| Host vs container, photopeak efficiency | 0.61090 vs 0.61060 — **0.05 %**, 0.1 sigma |
| Host vs container, total interaction probability | 0.87945 vs 0.87285 — 0.75 %, about 2.9 sigma |

Byte-identity **across** Geant4 versions is impossible by construction: different
model implementations consume random numbers differently, and the two
installations also carry different EM datasets (G4EMLOW 8.7 on the host, 8.8 in
the container). An earlier draft of this file claimed host and container should
produce identical bytes; that was wrong.

The two claims worth making are therefore:

1. **Reproducibility** — a given build reproduces its own results exactly. Verified.
2. **Physics agreement** — the beta and the release agree on the photopeak
   efficiency to 0.05 %. The 0.75 % difference in total interaction probability
   sits at about 2.9 sigma and is plausibly the EM dataset change, since it moves
   "any interaction" without moving full-energy deposition, as a Rayleigh or
   low-energy cross-section update would. More statistics are needed before
   claiming it is real rather than a fluctuation.

`env.json` records `platform.in_container` and the Geant4 tag for every run, so
results from the two builds can always be told apart after the fact.
