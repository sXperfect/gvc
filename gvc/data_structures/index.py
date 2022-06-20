from os.path import join
import numpy as np

class Index():
    """Index
    Consists of 2 Stages of Index:
        - Root Index (main.bin) contains the start POS and end POS of each block
        - Block Index ({block_id}.bin) contains the mapping from POS to row_ids
    """

    def __init__(self, index_dpath, decoder_context):
        self.index_dpath = index_dpath        

        # with open(join(self.index_dpath, 'main.bin'), 'rb') as f:
        #     payload = f.read()

        # self.root_idx = np.frombuffer(payload, dtype=int).reshape(-1, 2)
        # self.num_blocks = self.root_idx.shape[0]
        
        main_idx_fpath = join(self.index_dpath, 'main.npy')
        self.root_idx = np.load(main_idx_fpath)
        self.num_blocks = self.root_idx.shape[0]

        block_ptrs = []
        param_set_id_ptr = []
        for k in range(len(decoder_context.access_units)):
            blocks = decoder_context.access_units[k].blocks
            num_blocks = len(blocks)
            # param_set_id = decoder_context.access_units[k].get_param_set_id()
            param_set_id = decoder_context.access_units[k].header.parameter_set_id

            block_ptrs.extend(
                blocks
            )

            param_set_id_ptr.extend(
                [param_set_id]*num_blocks
            )

        block_ptrs = np.array(block_ptrs)
        param_set_id_ptr = np.array(param_set_id_ptr)

        # Root lookup contains 2 columns: block_pointer and parameter_set_id
        self.root_lookup = np.stack(
            (np.arange(self.num_blocks), block_ptrs, param_set_id_ptr), 
            axis=1
        )

        # Allocate memory for block_idx
        self.block_idx = [None] * self.num_blocks

    @classmethod
    def from_gvc_fpath(cls, input_fpath, decoder_context):
        input_dpath = input_fpath + '.metadata'
        try:
            return cls(input_dpath, decoder_context)
        except FileNotFoundError:
            return None

    def get_block_param_set_id_pairs(self, start_pos, end_pos):
        # NOTE: Slower version
        # block_mask = (self.root_idx[:, 0] <= start_pos) & (self.root_idx[:, 1] >= end_pos)
        # return self.root_lookup[block_mask]

        # ? NOTE: Faster version using binary search
        start_row = np.searchsorted(self.root_idx[:,0], start_pos, side='right')-1
        end_row = np.searchsorted(self.root_idx[:,1], end_pos, side='left')+1

        return self.root_lookup[start_row:end_row, :]

    def get_row_mask(self, block_id, start_pos, end_pos):

        curr_block_idx = self.block_idx[block_id]

        if curr_block_idx is None:
            # blk_idx_fpath = join(self.index_dpath, '.bin'.format(block_id))

            # with open(blk_idx_fpath, 'rb') as f:
            #     payload = f.read()

            # curr_block_idx = np.frombuffer(payload, dtype=int)
            blk_idx_fpath = join(self.index_dpath, f"{block_id}.npy")
            curr_block_idx = np.load(blk_idx_fpath)

            self.block_idx[block_id] = curr_block_idx

        # NOTE: Slower version
        # row_mask = (start_pos <= curr_block_idx) & (curr_block_idx <= end_pos)
        # return row_mask

        # NOTE: Faster version using binary search
        start_row = np.searchsorted(curr_block_idx, start_pos, side='left')
        end_row = np.searchsorted(curr_block_idx, end_pos, side='right')

        return slice(start_row, end_row, None)