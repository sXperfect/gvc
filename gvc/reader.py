from multiprocessing.sharedctypes import Value
import time
import logging as log

import numpy as np
from cyvcf2 import VCF

import gvc.common

from enum import IntEnum

class FORMAT_ID(IntEnum):
    TXT = 0
    VCF = 1
    
def recomp_p(gt):
    return (len(gt)+1)//2

v_recomp_p = np.vectorize(recomp_p)

def get_gt(gt, new_p, i):
    if i<new_p:
        return gt[i]
    else:
        return -2
    
v_get_gt = np.vectorize(get_gt)

def get_phase(gt, new_p, i):
    try:
        return gt[new_p[i]]
    except:
        return False
v_get_phase = np.vectorize(get_phase)

def vcf_genotypes_reader(fpath, block_size):

    vcf_f = VCF(fpath, strict_gt=True, gts012=True)
    num_samples = len(vcf_f.samples)

    assert block_size > 0, "Block size must be >= 1"

    i_var = 0
    stime = time.time()
    for variant in iter(vcf_f):
        p = variant.ploidy

        if i_var == 0:
            if p > 2:
                raise NotImplementedError("Currently support only ploidy <= 2!")

            allele_matrix = np.empty((block_size, p*num_samples), dtype=np.int8)
            phase_matrix = np.empty((block_size, (p-1)*num_samples), dtype=bool)

        # try:
        cyvcf_gt_matrix = np.array(variant.genotypes, dtype=np.int8)
        allele_matrix[i_var, :] = np.concatenate(cyvcf_gt_matrix[:, :p])

        try:
            phase_matrix[i_var, :] = np.concatenate(cyvcf_gt_matrix[:, p:])
        except:
            pass
        
        # except:

        i_var += 1
        
        if i_var == block_size:
            log.info('Parsing time: {:.03f}'.format(time.time() - stime))

             # TODO: handle not_available?
            allele_matrix, any_missing, not_available = gvc.binarization.adaptive_max_value(allele_matrix)

            yield (
                allele_matrix, 
                phase_matrix, 
                p, 
                any_missing, # any_missing
                not_available # not_available
            )

            i_var = 0
            stime = time.time()

    vcf_f.close()

    # TODO: what happer for not_available?
    suballele_matrix, any_missing, not_available = gvc.binarization.adaptive_max_value(allele_matrix[:i_var, :])

    yield (
        suballele_matrix,
        phase_matrix[:i_var, :],
        p, 
        any_missing,
        not_available
    )
    
def raw_vcf_genotypes_reader(fpath, block_size):
    vcf_f = VCF(fpath, strict_gt=True, gts012=True)
    num_samples = len(vcf_f.samples)

    assert block_size > 0, "Block size must be >= 1"

    i_var = 0
    p = 0
    block = [None] * block_size
    for variant in iter(vcf_f):
        p = max(p, variant.ploidy)
        genotypes = variant.genotypes
        
        # assert not (variant.gt_types == 3).any(), '"Missing found!"'
        
        block[i_var] = genotypes
        i_var += 1
        
        if i_var == block_size:
            yield FORMAT_ID.VCF, block[:i_var], num_samples, p
            
            i_var = 0
            p = 0
            block = [None] * block_size

    vcf_f.close()
    yield FORMAT_ID.VCF, block[:i_var], num_samples, p
    
def raw_vcf_to_gt_mat(raw_block):
    block, num_samples, p = raw_block
    block_size = len(block)
    
    try:
        cyvcf_gt_matrix = np.array(block, dtype=np.int8)
        # allele_matrix = cyvcf_gt_matrix[:, :, :p]
        # allele_matrix = np.concatenate(np.split(cyvcf_gt_matrix[:, :, :p], num_samples, axis=1), axis=-1).reshape(block_size, -1)
        allele_mat = np.concatenate(np.split(cyvcf_gt_matrix[:, :, :p], p, axis=-1), axis=1).reshape(block_size, -1)

        try:
            phase_mat = cyvcf_gt_matrix[:, :, p:].reshape(block_size, -1)
        except:
            phase_mat = None
            
        return allele_mat, phase_mat, p, False, False
    
    except ValueError:
        not_available = True
        
        block_np = np.array(block, dtype=object)
        p_per_gt = v_recomp_p(block_np)
        
        assert p_per_gt.max() == p
        
        allele_mat = np.empty((block_size, p*num_samples), dtype=np.int8)
                
        for i in range(p):
            allele_mat[:, i::p] = v_get_gt(block_np, p_per_gt, i)
            
        if p > 1:
            phase_mat = np.empty((block_size, (p-1)*num_samples), dtype=bool)
            for i in range(p-1):
                phase_mat[:, i::(p-1)] = v_get_phase(block_np, p_per_gt, i)
        else:
            phase_mat = None
            
        return allele_mat, phase_mat, p, False, not_available

def txt_genotypes_reader(fpath, block_size):

    lines = []
    is_eof = False

    with open(fpath, 'r') as f:
        for line in f:
            if len(line):
                lines.append(line)

            if len(lines) == block_size:

                # Excecute part 4.1 - splitting genotype matrix
                log.info('Execute part 4.1 - Splitting genotype matrix')
                stime = time.time()
                allele_matrix, phasing_matrix, p = gvc.binarization.split_genotype_matrix(lines)
                log.info('Parsing time: {:.03f}'.format(time.time() - stime))

                assert(np.issubdtype(allele_matrix.dtype, np.signedinteger))

                allele_matrix, any_missing, not_available = gvc.binarization.adaptive_max_value(allele_matrix)

                yield allele_matrix, phasing_matrix, p, any_missing, not_available

                lines.clear()

    if len(lines):

        # Excecute part 4.1 - splitting genotype matrix
        log.info('Execute part 4.1 - Splitting genotype matrix')
        stime = time.time()
        allele_matrix, phasing_matrix, p = gvc.binarization.split_genotype_matrix(lines)
        log.info('Parsing time: {:.03f}'.format(time.time() - stime))

        assert(np.issubdtype(allele_matrix.dtype, np.signedinteger))

        allele_matrix, any_missing, not_available = gvc.binarization.adaptive_max_value(allele_matrix)

        yield allele_matrix, phasing_matrix, p, any_missing, not_available
        
# Transform a sample matrix (extracted as a 2D matrix of GT values from a VCF file) into two matrixes:
#   1) sample value matrix
#   2) phase matrix
#
# Example input:
#   0|1	0|1	0|1	0|0
#   0/0	0/1	0/1	0/0
#   0|0	0|0	0|0	0|0
#
# Example output:
#   gtmat= 01010100    phasemat = 0000
#          00010100               1111
#          00000000               0000

def split_gt_phase(fpath):
    gt_mat = None
    phase_mat = None

    nrows = gvc.utils.line_cnt(fpath)

    with open(fpath) as f:
        for iline, line in enumerate(f):
            line = line.strip()

            variant_data = []
            phase_data = []

            for sample in line.split('\t'):
                if (len(sample) != 3):
                    raise gvc.errors.GvcError('sample matrix seems to contain non-diploid data')

                lcode, sep, rcode = sample
                lcode = int(lcode)
                rcode = int(rcode)

                if sep in ('/', '|'):
                    if sep == '|':
                        is_unphased = False
                    else:
                        is_unphased = True
                    phase_data.append(is_unphased)
                else:
                    raise ValueError('Invalid separator found {}'.format(sep))

                variant_data.append(lcode)
                variant_data.append(rcode)

            if gt_mat is None:
                ncols = len(variant_data)
                gt_mat = np.zeros((nrows, ncols), dtype=np.uint8)

            if phase_mat is None:
                ncols = len(phase_data)
                phase_mat = np.zeros((nrows, ncols), dtype=bool)

            gt_mat[iline, :] = variant_data
            phase_mat[iline, :] = phase_data

    # If all data is phased, do not store
    if not np.any(phase_mat):
        phase_mat = None

    return gt_mat, phase_mat
