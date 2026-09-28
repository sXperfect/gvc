from dataclasses import dataclass
from os import makedirs
from os.path import join, exists
from shutil import rmtree
import time
import logging as log

import numpy as np
from cyvcf2 import VCF

import gvc.common

from enum import IntEnum

class FORMAT_ID(IntEnum):
    TXT = 0
    VCF = 1
    
@dataclass
class MetaHandler(object):
    vcf_f:VCF
    metadata_dpath:str
    block_size:int
    
    def init(self):
        self.mkdir_root()
        self.write_header()
        self.min_max_pos_list = []

    @property
    def header_fpath(self):
        return join(self.metadata_dpath, 'header.txt')
    
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
        if not self.is_enabled:
            return
        
        self.pos_arr = np.empty(self.block_size, np.uint64)
    
    def write_header(self):
        if not self.is_enabled:
            return
        
        with open(self.header_fpath, 'w') as f:
            f.write(self.vcf_f.raw_header.strip())
            
        sample_ids = np.array(self.vcf_f.raw_header.strip().split('\n')[-1].split('\t')[9:])
        np.save(
            join(self.metadata_dpath, "samples"),
            sample_ids
        )
        
    def proc_var(self, i_var, variant):
        if not self.is_enabled:
            return
        
        self.pos_arr[i_var] = variant.POS
        
    def proc_block(self, block_id, n_vars=None):
        if not self.is_enabled:
            return
        
        min_max_pos = np.array([self.pos_arr.min(), self.pos_arr.max()], dtype=np.uint64)
        self.min_max_pos_list.append(min_max_pos)
            
        if n_vars is None:
            np.save(
                join(self.metadata_dpath, f"{block_id}"),
                self.pos_arr
            )
        else:
            np.save(
                join(self.metadata_dpath, f"{block_id}"),
                self.pos_arr[:n_vars]
            )
            
    def end(self):
        if not self.is_enabled:
            return
        
        np.save(
            join(self.metadata_dpath, "main"),
            np.stack(self.min_max_pos_list)
        )

def reshape_trans_mat(mat, axis):
    block_size, num_samples, p = mat.shape
    
    if axis == 1:
        trans_mat = np.concatenate(np.split(mat, num_samples, axis=axis), axis=-1).reshape(block_size, -1)
    elif axis == 2 or axis == -1:
        trans_mat = np.concatenate(np.split(mat, p, axis=axis), axis=1).reshape(block_size, -1)
    else:
        raise ValueError("Invalid axis value!")
        
    return trans_mat

def vcf_genotypes_reader(fpath, out_fpath, block_size):

    vcf_f = VCF(fpath, strict_gt=True, gts012=True, threads=2)
    num_samples = len(vcf_f.samples)
    
    if out_fpath is not None:
        metadata_dpath = join(f'{out_fpath}.metadata')
    else:
        metadata_dpath = None #? Disable metadata handler
    
    meta_handler = MetaHandler(vcf_f, metadata_dpath, block_size)
    meta_handler.init()

    assert block_size > 0, "Block size must be greater than zero"

    i_var = 0
    stime = time.time()
    p = 0
    block_id = 0
    for variant in iter(vcf_f):
        p = variant.ploidy

        if i_var == 0:
            #TODO: is int8 as the data type of allele_matrix sufficient? genotypes is int16
            allele_matrix = np.empty((block_size, num_samples, p), dtype=gvc.common.SIGNED_ALLELE_DTYPE) 
            phase_matrix = np.empty((block_size, num_samples, (p-1)), dtype=bool)
            meta_handler.init_block()

        meta_handler.proc_var(i_var, variant)
        
        #? -2 represent missing allele
        genotypes = variant.genotype.array() #? int16
        allele_matrix[i_var, :, :] = genotypes[:, :p]

        try:
            phase_matrix[i_var, :, :] = genotypes[:, p:]
        except:
            pass

        i_var += 1
        
        if i_var == block_size:
            log.debug('Parsing time: {:.03f}'.format(time.time() - stime))

            allele_matrix, missing_rep_val, na_rep_val = gvc.binarization.adaptive_max_value(allele_matrix)

            allele_matrix = reshape_trans_mat(allele_matrix, 1)
            phase_matrix = reshape_trans_mat(phase_matrix, 1)
            
            meta_handler.proc_block(block_id)
            
            yield (
                allele_matrix, 
                phase_matrix, 
                p, 
                missing_rep_val, # any_missing
                na_rep_val # not_available
            )

            i_var = 0
            block_id += 1
            stime = time.time()

    vcf_f.close()

    if i_var == 0:
        return
    
    suballele_matrix, missing_rep_val, na_rep_val = gvc.binarization.adaptive_max_value(allele_matrix[:i_var, :])
    subphase_matrix = phase_matrix[:i_var, :]
    
    suballele_matrix = reshape_trans_mat(suballele_matrix, 1)
    subphase_matrix = reshape_trans_mat(subphase_matrix, 1)
    
    meta_handler.proc_block(block_id, i_var)
    meta_handler.end()
                
    yield (
        suballele_matrix,
        subphase_matrix,
        p, 
        missing_rep_val,
        na_rep_val
    )