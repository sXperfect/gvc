#!/usr/bin/env bash

# Run 'unittest' unit tests

# Safer Bash script:
#   -e           Exit immediately when a command fails
#   -u           Treat unset variables as an error and exit immediately
#   -x           Print each command before executing it
#   -o pipefail  Set exit code of a pipeline to that of the rightmost command to exit with a non-zero status
set -euo pipefail

# Get git root directory
git rev-parse --git-dir 1>/dev/null # Exit if not inside git repo
readonly git_root_dir="$(git rev-parse --show-toplevel)"

# Run unit tests
(
    cd "${git_root_dir}"
    python3 -m unittest discover --start-directory tests --verbose
)
