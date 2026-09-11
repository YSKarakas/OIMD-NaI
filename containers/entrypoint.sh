#!/bin/bash
# Load the Geant4 dataset locations captured at image build time, then exec.
# The variables are derived from the installation's own geant4.sh rather than
# hard-coded, so a change of Geant4 version cannot leave them silently stale.
set -a
. /etc/geant4.env
set +a
exec "$@"
