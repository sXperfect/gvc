#!/usr/bin/env bash

# Setup everything for GVC

# Safer Bash script:
#   -e           Exit immediately when a command fails
#   -u           Treat unset variables as an error and exit immediately
#   -x           Print each command before executing it
#   -o pipefail  Set exit code of a pipeline to that of the rightmost command to exit with a non-zero status
set -euxo pipefail

# Get git root directory
if (( $# < 1 )); then
    git rev-parse --git-dir 1>/dev/null # Exit if not inside git repo
    readonly git_root_dir="$(git rev-parse --show-toplevel)"
else
    readonly git_root_dir="${1}"
fi

# Prepare python
readonly venv_dir="${git_root_dir}/tmp/venv"
readonly python_bin="${venv_dir}/bin/python"
{
    cd ${git_root_dir}
    virtualenv -p python3 ${venv_dir}
    ${python_bin} -m pip install -r requirements.txt
    ${python_bin} setup.py build_ext --inplace
}

#? LibGVC
readonly libgvc_dir="${git_root_dir}/library/libgvc"
{
    cd ${libgvc_dir}
    cmake -S src/ -B build
    cd build
    make --jobs
}