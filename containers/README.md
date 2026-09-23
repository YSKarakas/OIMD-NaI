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
configuration and seed; the optical columns re-read on 2026-09-13 after a
referee asked where the runs had been made.

| Comparison | Result |
|---|---|
| Same build, same seed, run twice | **byte-identical** (host and container both) |
| Host beta vs container release, event by event | 61 % of events differ |
| Host vs container, photopeak efficiency | 0.61090 vs 0.61060 -- 0.05 %, 0.1 sigma |
| Host vs container, total interaction probability | 0.87945 vs 0.87285 -- 0.75 %, about 2.9 sigma |
| **Host vs container, detected per generated photon, events above 660 keV** | **release collects 4.3 % more -- 19 sigma (12 217 and 12 209 events, compared unpaired)** |
| **Host vs container, mean detection time** | **3.415 ns vs 2.215 ns -- the release is a third shorter, 70 sigma** |

The light-collection row is the comparison the manuscript quotes, and
`scripts/check_paper_numbers.py` recomputes it: each run on its own events
above 660 keV that generate light, the per-event detected/generated fraction
compared unpaired. Until 23 September 2026 the row gave a different statistic
-- photons detected per event, on the 7477 events above 600 keV in both runs
matched by event number, 4.8 sigma -- which means little when 61 % of the
events differ between the two builds.

An earlier version of this section compared only the two EM quantities, found
them in agreement, and concluded "physics agreement". That conclusion was wrong
for the quantity this project actually measures. The beta and the release
differ in optical transport itself, consistent with the refactoring of
G4OpBoundaryProcess during the 11.4 cycle noted in the Dockerfile. It also said
the host should be aligned with the release "before any results are
published"; that was not done for the first campaign, and a referee found it.

Consequences, in force since 2026-09-13:

1. **Every published run is made in this container** through
   `scripts/run_in_container.sh`, which records `in_container: true`, the
   release tag `geant4-11-04-patch-02`, and the host's git state (the image has
   no git, so the launcher passes it in through the environment).
2. Runs made with the pre-release are config schema 3 and earlier; schema 4 is
   the release campaign that exposed the scintillation cascade (kept as
   evidence); the published campaign is schema 5. The analysis refuses to
   compare across schemas.
3. `/etc/geant4.env`, which the entrypoint sources, is **empty** in the image as
   built (the `grep` in the Dockerfile matched nothing), so `LD_LIBRARY_PATH`
   and the dataset variables are unset for any command run directly. The
   launcher sources `/opt/geant4/bin/geant4.sh` itself. The image was rebuilt
   from this Dockerfile, unchanged, on a new machine on 23 September 2026
   (`scint:11.4.2@sha256:b5e01ff7b473...`); four published runs re-run in it
   reproduced the per-event output checksums they recorded, byte for byte
   (`runs_scratch/diag/rebuilt_image_reproduction.json`). The empty
   `/etc/geant4.env` is left as it is, so that the recipe stays the one the
   published image was built from.

Byte-identity **across** Geant4 versions is impossible by construction:
different model implementations consume random numbers differently, and the two
installations also carry different EM datasets (G4EMLOW 8.7 on the host, 8.8 in
the container).
