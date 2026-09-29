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
        if not self.is_enabled or not self.min_max_pos_list:
            return
        np.save(join(self.metadata_dpath, "main"), np.stack(self.min_max_pos_list))


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

    metadata_dpath = join("{}{}".format(out_fpath, ".metadata")) if out_fpath is not None else None
    meta_handler = MetaHandler(vcf_f, metadata_dpath, block_size)
    meta_handler.init()

    i_var = 0
    stime = time.time()
    p = 0
    block_id = 0

    for variant in iter(vcf_f):
        p = variant.ploidy
        if i_var == 0:
            allele_matrix = np.empty(
                (block_size, num_samples, p), dtype=gvc.common.SIGNED_ALLELE_DTYPE
            )
            if p == 1:
                phase_matrix = np.empty((block_size, 0), dtype=bool)
            else:
                phase_matrix = np.empty(
                    (block_size, num_samples, p - 1), dtype=bool
                )
            meta_handler.init_block()

        meta_handler.proc_var(i_var, variant)
        genotypes = variant.genotype.array()
        allele_matrix[i_var, :, :] = genotypes[:, :p]
        if p > 1:
            # cyvcf2 uses True for "|" while GVC serializes 0 for "|" and
            # 1 for "/". Convert at the ingestion boundary so every internal
            # path uses the same phase convention as split_genotype_matrix().
            phase_matrix[i_var, :, :] = np.logical_not(genotypes[:, p:])

        i_var += 1
        if i_var == block_size:
            log.debug("Parsing time: %.03f", time.time() - stime)
            allele_matrix, missing_rep_val, na_rep_val = gvc.binarization.adaptive_max_value(
                allele_matrix
            )
            allele_matrix = reshape_trans_mat(allele_matrix, 1)
            if p > 1:
                phase_matrix = reshape_trans_mat(phase_matrix, 1)
            else:
                phase_matrix = phase_matrix[:block_size]
            meta_handler.proc_block(block_id)
            yield allele_matrix, phase_matrix, p, missing_rep_val, na_rep_val
            i_var = 0
            block_id += 1
            stime = time.time()

    vcf_f.close()
    if i_var == 0:
        return

    suballele_matrix, missing_rep_val, na_rep_val = gvc.binarization.adaptive_max_value(
        allele_matrix[:i_var, :]
    )
    subphase_matrix = phase_matrix[:i_var]
    suballele_matrix = reshape_trans_mat(suballele_matrix, 1)
    if p > 1:
        subphase_matrix = reshape_trans_mat(subphase_matrix, 1)
    meta_handler.proc_block(block_id, i_var)
    meta_handler.end()
    yield suballele_matrix, subphase_matrix, p, missing_rep_val, na_rep_val
