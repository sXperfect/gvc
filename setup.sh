#!/usr/bin/env bash
set -euo pipefail

git rev-parse --git-dir >/dev/null
readonly git_root_dir="${1:-$(git rev-parse --show-toplevel)}"

cd "${git_root_dir}"

python3 -m pip install -e ".[all,test]"

cmake -S library/libgvc/src -B library/libgvc/build
cmake --build library/libgvc/build --parallel
