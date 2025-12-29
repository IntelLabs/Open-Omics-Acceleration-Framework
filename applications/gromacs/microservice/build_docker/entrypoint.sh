#!/bin/bash
set -eo pipefail

source /opt/conda/bin/activate grms_env
source /grmcs/gmx_mkl_prefix/bin/GMXRC

exec "$@"

