from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
import logging as log
from os import makedirs
from os.path import exists, join
from shutil import rmtree
import time
from typing import Optional

import numpy as np

import gvc.common


def _vcf_class():
    try:
        from cyvcf2 import VCF
    except ImportError as exc:
        raise ImportError(
            "VCF support requires the optional dependency group: "
            "python -m pip install 'gvc[vcf]'"
        ) from exc
    return VCF


class FORMAT_ID(IntEnum):
    TXT = 0
    VCF = 1


@dataclass
class MetaHandler:
    vcf_f: object
    metadata_dpath: Optional[str]
    block_size: int

    def init(self):
        self.mkdir_root()
        self.write_header()
        self.min_max_pos_list = []

    @property
    def header_fpath(self):
        return join(self.metadata_dpath, "header.txt")

    @property
    def is_enabled(self):
        return self.metadata_dpath is not None

    def mkdir_root(self):
        if not self.is_enabled:
            return
        if exists(self.metadata_dpath):
            rmtree(self.metadata_dpath)
        makedirs(self.metadata_dpath)

    def init_block(self):
        if self.is_enabled:
            self.pos_arr = np.empty(self.block_size, dtype=np.uint64)

    def write_header(self):
        if not self.is_enabled:
            return
        with open(self.header_fpath, "w") as f:
            f.write(self.vcf_f.raw_header.strip())
        sample_ids = np.array(self.vcf_f.raw_header.strip().split("\n")[-1].split("\t")[9:])
        np.save(join(self.metadata_dpath, "samples"), sample_ids)

    def proc_var(self, i_var, variant):
        if self.is_enabled:
            self.pos_arr[i_var] = variant.POS

    def proc_block(self, block_id, n_vars=None):
        if not self.is_enabled:
            return
        positions = self.pos_arr if n_vars is None else self.pos_arr[:n_vars]
        if positions.size == 0:
            return
        self.min_max_pos_list.append(
            np.array([positions.min(), positions.max()], dtype=np.uint64)
        )
        np.save(join(self.metadata_dpath, str(block_id)), positions)

    def end(self):
        if not self.is_enabled:
            return
        if self.min_max_pos_list:
            root = np.stack(self.min_max_pos_list)
        else:
            root = np.empty((0, 2), dtype=np.uint64)
        np.save(join(self.metadata_dpath, "main"), root)


def reshape_trans_mat(mat, axis):
    block_size, num_samples, p = mat.shape
    if axis == 1:
        return np.concatenate(np.split(mat, num_samples, axis=axis), axis=-1).reshape(
            block_size, -1
        )
    if axis in (2, -1):
        return np.concatenate(np.split(mat, p, axis=axis), axis=1).reshape(
            block_size, -1
        )
    raise ValueError("Invalid axis value!")


def vcf_genotypes_reader(fpath, out_fpath, block_size):
    if block_size <= 0:
        raise ValueError("block_size must be greater than zero")

    VCF = _vcf_class()
    vcf_f = VCF(fpath, strict_gt=True, gts012=True, threads=2)
    num_samples = len(vcf_f.samples)

    metadata_dpath = (
        join("{}{}".format(out_fpath, ".metadata"))
        if out_fpath is not None
        else None
    )
    meta_handler = MetaHandler(vcf_f, metadata_dpath, block_size)
    meta_handler.init()

    def allocate_block(ploidy):
        allele = np.empty(
            (block_size, num_samples, ploidy),
            dtype=gvc.common.SIGNED_ALLELE_DTYPE,
        )
        if ploidy == 1:
            phase = np.empty((block_size, 0), dtype=bool)
        else:
            phase = np.empty(
                (block_size, num_samples, ploidy - 1),
                dtype=bool,
            )
        meta_handler.init_block()
        return allele, phase

    def finalize_block(allele, phase, n_rows, ploidy, block_id):
        allele_block, missing_rep_val, na_rep_val = (
            gvc.binarization.adaptive_max_value(allele[:n_rows])
        )
        allele_block = reshape_trans_mat(allele_block, 1)

        phase_block = phase[:n_rows]
        if ploidy > 1:
            phase_block = reshape_trans_mat(phase_block, 1)

        meta_handler.proc_block(block_id, n_rows)
        return (
            allele_block,
            phase_block,
            ploidy,
            missing_rep_val,
            na_rep_val,
        )

    i_var = 0
    block_id = 0
    current_p = None
    allele_matrix = None
    phase_matrix = None

    try:
        for variant in iter(vcf_f):
            variant_p = int(variant.ploidy)
            if variant_p <= 0:
                raise ValueError("VCF genotype ploidy must be greater than zero")

            # A ParameterSet carries one ploidy value. If ploidy changes before
            # block_size is reached, close the current block so the encoder can
            # emit a new ParameterSet rather than mixing incompatible shapes.
            if i_var and variant_p != current_p:
                yield finalize_block(
                    allele_matrix,
                    phase_matrix,
                    i_var,
                    current_p,
                    block_id,
                )
                block_id += 1
                i_var = 0

            if i_var == 0:
                current_p = variant_p
                allele_matrix, phase_matrix = allocate_block(current_p)

            meta_handler.proc_var(i_var, variant)
            genotypes = variant.genotype.array()
            if genotypes.ndim != 2 or genotypes.shape[0] != num_samples:
                raise ValueError("unexpected cyvcf2 genotype array shape")
            required_columns = current_p + (1 if current_p > 1 else 0)
            if genotypes.shape[1] < required_columns:
                raise ValueError(
                    "cyvcf2 genotype array does not contain the expected "
                    "allele/phasing columns"
                )

            allele_values = genotypes[:, :current_p]
            signed_info = np.iinfo(gvc.common.SIGNED_ALLELE_DTYPE)
            if np.any(allele_values < -2) or np.any(allele_values > signed_info.max):
                raise ValueError(
                    "VCF allele index is outside the supported range -2..{}".format(
                        signed_info.max
                    )
                )
            allele_matrix[i_var, :, :] = allele_values
            if current_p > 1:
                # cyvcf2 uses True for "|" while GVC serializes 0 for "|" and
                # 1 for "/". A single cyvcf2 phase flag is broadcast across
                # separators for polyploid calls, matching the historical GVC
                # convention.
                phase_matrix[i_var, :, :] = np.logical_not(
                    genotypes[:, current_p:]
                )

            i_var += 1
            if i_var == block_size:
                yield finalize_block(
                    allele_matrix,
                    phase_matrix,
                    i_var,
                    current_p,
                    block_id,
                )
                block_id += 1
                i_var = 0

        if i_var:
            yield finalize_block(
                allele_matrix,
                phase_matrix,
                i_var,
                current_p,
                block_id,
            )
    finally:
        vcf_f.close()
        # This must also run when the final block exactly fills block_size.
        meta_handler.end()
