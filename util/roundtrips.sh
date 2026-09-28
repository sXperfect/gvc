#!/usr/bin/env bash

# Run and verify GVC encoder/decoder roundtrips

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


# Run 'toy' roundtrip
(
    readonly genotype_matrix_file="${git_root_dir}/data/toy.sample_matrix.txt"
    readonly bitstream_file="${git_root_dir}/data/toy.gvc"
    readonly reconstructed_genotype_matrix_file="${git_root_dir}/data/recon.toy.sample_matrix.txt"

    cd "${git_root_dir}"

    python3 -m gvc \
      -b \
      5 \
      --binarization \
      bit_plane \
      --encoder \
      jbig1 \
      --dist \
      hrl \
      --sort-rows \
      --sort-cols \
      --transpose \
      ${genotype_matrix_file} \
      ${bitstream_file}

    python3 -m gvc \
      --decode \
      --log_level debug \
      ${bitstream_file} \
      ${reconstructed_genotype_matrix_file}

    diff --report-identical-files "${genotype_matrix_file}" "${reconstructed_genotype_matrix_file}"

    rm "${bitstream_file}"
    rm "${reconstructed_genotype_matrix_file}"
)
