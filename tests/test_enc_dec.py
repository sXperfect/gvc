from os import getcwd
from os.path import join
import itertools as it
import numpy as np
from gvc.common import create_parameter_set
from gvc.sort import sort
from gvc import reader
from gvc import data_structures as ds
from gvc.codec import encode, decode
from gvc.data_structures.consts import BinarizationID, CodecID
from gvc.encoder import run_core, binarize_allele_matrix
from gvc.decoder import decode_encoded_variants
import unittest


class TestEncodeDecode(unittest.TestCase):
    
    def setUp(self) -> None:
        self.vcf01_fpath = 'test_block01.vcf.gz'
        self.block_size = 2048
        self.codec_id = CodecID.JBIG1
        self.transpose = False
        
        self.tsp_params = [
            'ham', 'nn', 1
        ]
        
    def test_roundtrip_manual_vcf01(self):
        fpath = join(getcwd(), 'tests', self.vcf01_fpath)
        for binarization_id in [BinarizationID.ROW_BIN_SPLIT, BinarizationID.BIT_PLANE]:
            for [sort_rows, sort_cols] in it.product([True, True], repeat=2):
                
                ps_params = [
                    binarization_id,
                    self.codec_id,
                    0,
                    sort_rows,
                    sort_cols,
                    self.transpose
                ]
                
                f1_iter = reader.vcf_genotypes_reader(fpath, None, self.block_size)
                
                for block_ID, raw_block in enumerate(f1_iter):

                    allele_matrix, phasing_matrix, p, missing_rep_val, na_rep_val = raw_block
                    binarization_id,codec_id,axis,sort_rows,sort_cols,transpose = ps_params
                    dist_f_name, solver_name, solver_profile = self.tsp_params

                    # Execute part 4.2 - binarization of allele matrix
                    bin_allele_matrices, additional_info = binarize_allele_matrix(
                        allele_matrix, 
                        binarization_id, 
                        axis=axis
                    )
                                
                    #? Create parameter based on binarization and encoder parameter
                    param_set = create_parameter_set(
                        missing_rep_val,
                        na_rep_val,
                        p,
                        phasing_matrix,
                        additional_info,
                        binarization_id,
                        codec_id,
                        axis,
                        sort_rows,
                        sort_cols,
                        transpose=transpose
                    )

                    # Execute part 4.3 - sorting
                    sorted_data = sort(
                        param_set, 
                        bin_allele_matrices, 
                        phasing_matrix, 
                        dist_f_name=dist_f_name, 
                        solver_name=solver_name,
                        solver_profile=solver_profile,
                    )

                    # Execute part 4.4 - entropy coding
                    data_bytes = encode(param_set, additional_info, *sorted_data)

                    # Initialize EncodedVariant, ParameterSet is not stored internally in EncodedVariants
                    # Order of arguments here is important (See EncodedVariants)
                    enc_variant = ds.GenotypePayload(param_set, *data_bytes, missing_rep_val, na_rep_val)

                    # Create new Block and store
                    block = ds.Block.from_encoded_variant(enc_variant)
                    
                    encoded_variants = block.block_payload
                    
                    for i_bin_mat in range(param_set.num_variants_flags):
                    
                        bin_mat, row_ids, col_ids = decode(
                            encoded_variants.variants_payloads[i_bin_mat],
                            encoded_variants.variants_row_ids_payloads[i_bin_mat],
                            encoded_variants.variants_col_ids_payloads[i_bin_mat],
                            param_set.variants_coder_ids[i_bin_mat],
                            unsort=False
                        )
                        
                        orig_row_ids = sorted_data[1][i_bin_mat]
                        if orig_row_ids is not None:
                            self.assertTrue(
                                np.array_equal(orig_row_ids, row_ids)
                            )

                        orig_col_ids = sorted_data[2][i_bin_mat]
                        if orig_col_ids is not None:
                            self.assertTrue(
                                np.array_equal(orig_col_ids, col_ids)
                            )
                        
                        orig_bin_mat = bin_allele_matrices[i_bin_mat]
                    
                        if row_ids is not None:
                            bin_mat = bin_mat[row_ids, :]
                        if col_ids is not None:
                            bin_mat = bin_mat[:, col_ids]
                            
                        self.assertTrue(
                            np.array_equal(bin_mat, orig_bin_mat)
                        )
                    
    def test_roundtrip_vcf01(self):
        
        fpath = join(getcwd(), 'tests', self.vcf01_fpath)
        for binarization_id in [BinarizationID.ROW_BIN_SPLIT, BinarizationID.BIT_PLANE]:
            for [sort_rows, sort_cols] in it.product([False, True], repeat=2):                
                ps_params = [
                    binarization_id,
                    self.codec_id,
                    0,
                    sort_rows,
                    sort_cols,
                    self.transpose
                ]
                
                f1_iter = reader.vcf_genotypes_reader(fpath, None, self.block_size)
                for block_ID, raw_block in enumerate(f1_iter):
                    
                    allele_mat, phase_mat, p, missing_repval, na_repval = raw_block                    
                    block, new_param_set = run_core(raw_block, ps_params, self.tsp_params)
                    
                    encoded_variants = block.block_payload
                    recon_allele_mat, recon_phase_mat = decode_encoded_variants(
                        new_param_set, encoded_variants, ret_gt=False)
                    
                    self.assertTrue(
                        np.array_equal(allele_mat, recon_allele_mat),
                        msg=f"{[sort_rows, sort_cols]}"
                    )
                
    def test_roundtrip_dummy01(self):
        
        dummy_out = """
        0/0\t0|1\t3/0\n
        3|1\t0/2\t1|1\n
        3|1\t0/2\t1|2\n
        0/0\t0|1\t3/1\n
        0/0\t.|.\t3/1\n
        """
        
        target_allele_mat = np.array([
            [0, 0, 0, 1, 3, 0],
            [3, 1, 0, 2, 1, 1],
            [3, 1, 0, 2, 1, 2],
            [0, 0, 0, 1, 3, 1],
            [0, 0, -1, -1, 3, 1]
        ], dtype=np.int8)
        
        allele_mat = target_allele_mat.copy()
        missing_repval = target_allele_mat.max() + 1
        allele_mat[allele_mat == -1] = missing_repval
        
        phasing_mat = np.array([
            [1, 0, 1],
            [0, 1, 0],
            [0, 1, 0],
            [1, 0, 1],
            [1, 0, 1],
        ], dtype=bool)
        
        raw_block = allele_mat, phasing_mat, 2, missing_repval, None
        
        for binarization_id in [BinarizationID.ROW_BIN_SPLIT, BinarizationID.BIT_PLANE]:
            for [sort_rows, sort_cols] in it.product([False, True], repeat=2):                
                ps_params = [
                    binarization_id,
                    self.codec_id,
                    0,
                    sort_rows,
                    sort_cols,
                    self.transpose
                ]
        
                allele_mat, phase_mat, p, missing_repval, na_repval = raw_block  
                block, new_param_set = run_core(raw_block, ps_params, self.tsp_params)
                
                encoded_variants = block.block_payload
                recon_allele_mat, recon_phase_mat = decode_encoded_variants(
                    new_param_set, encoded_variants, ret_gt=False)
                
                self.assertTrue(np.array_equal(target_allele_mat, recon_allele_mat))
                
                gt_mat = decode_encoded_variants(
                    new_param_set, encoded_variants, ret_gt=True)
                
                pass
                
    def test_roundtrip_dummy01_p_val(self):
        
        dummy_out = """
        0/0\t0/1\t3/0\n
        3/1\t0/2\t1/1\n
        3/1\t0/2\t1/2\n
        0/0\t0/1\t3/1\n
        0/0\t./.\t3/1\n
        """
        
        target_allele_mat = np.array([
            [0, 0, 0, 1, 3, 0],
            [3, 1, 0, 2, 1, 1],
            [3, 1, 0, 2, 1, 2],
            [0, 0, 0, 1, 3, 1],
            [0, 0, -1, -1, 3, 1]
        ], dtype=np.int8)
        
        allele_mat = target_allele_mat.copy()
        missing_repval = target_allele_mat.max() + 1
        allele_mat[allele_mat == -1] = missing_repval
        
        phasing_mat = np.ones((5,3), dtype=bool)
        
        raw_block = allele_mat, phasing_mat, 2, missing_repval, None
        
        for binarization_id in [BinarizationID.ROW_BIN_SPLIT, BinarizationID.BIT_PLANE]:
            for [sort_rows, sort_cols] in it.product([False, True], repeat=2):                
                ps_params = [
                    binarization_id,
                    self.codec_id,
                    0,
                    sort_rows,
                    sort_cols,
                    self.transpose
                ]
        
                allele_mat, phasing_mat, p, missing_repval, na_repval = raw_block  
                block, new_param_set = run_core(raw_block, ps_params, self.tsp_params)
                
                encoded_variants = block.block_payload
                recon_allele_mat, recon_phase_mat = decode_encoded_variants(
                    new_param_set, encoded_variants, ret_gt=False)
                
                self.assertTrue(np.array_equal(target_allele_mat, recon_allele_mat))
                
                gt_mat = decode_encoded_variants(
                    new_param_set, encoded_variants, ret_gt=True)
                
                pass
        