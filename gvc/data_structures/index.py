from pathlib import Path
from os.path import join

import numpy as np


class Index:
    """Random-access sidecar index for an encoded GVC stream."""

    def __init__(self, index_dpath, decoder_context):
        self.index_dpath = index_dpath

        self.root_idx = np.load(join(self.index_dpath, "main.npy"), allow_pickle=False)
        self.samples = np.load(join(self.index_dpath, "samples.npy"), allow_pickle=False)

        if self.root_idx.ndim != 2 or self.root_idx.shape[1] != 2:
            raise ValueError("root index must have shape (n_blocks, 2)")
        if np.any(self.root_idx[:, 0] > self.root_idx[:, 1]):
            raise ValueError("root index contains a block with start > end")
        if self.root_idx.shape[0] > 1 and np.any(
            self.root_idx[1:, 0] < self.root_idx[:-1, 0]
        ):
            raise ValueError("root index must be sorted by block start position")
        if self.samples.ndim != 1:
            raise ValueError("sample index must be one-dimensional")
        if len(set(self.samples.tolist())) != len(self.samples):
            raise ValueError("sample index contains duplicate sample IDs")
        encoded_sample_count = getattr(decoder_context, "ncols", None)
        if (
            encoded_sample_count is not None
            and len(self.samples) != encoded_sample_count
        ):
            raise ValueError(
                "metadata sample count does not match encoded data: "
                "expected {}, got {}".format(
                    encoded_sample_count,
                    len(self.samples),
                )
            )

        block_ptrs = []
        param_set_ids = []
        for access_unit_id in sorted(decoder_context.access_units):
            access_unit = decoder_context.access_units[access_unit_id]
            block_ptrs.extend(access_unit.blocks)
            param_set_ids.extend(
                [access_unit.header.parameter_set_id] * len(access_unit.blocks)
            )

        if len(block_ptrs) != self.num_blocks:
            raise ValueError(
                "metadata block count does not match encoded access units"
            )

        self.root_lookup = np.empty((self.num_blocks, 3), dtype=object)
        if self.num_blocks:
            self.root_lookup[:, 0] = np.arange(self.num_blocks)
            self.root_lookup[:, 1] = block_ptrs
            self.root_lookup[:, 2] = param_set_ids

        self.block_idx = [None] * self.num_blocks

    @property
    def num_blocks(self):
        return self.root_idx.shape[0]

    @property
    def num_samples(self):
        return len(self.samples)

    @classmethod
    def from_gvc_fpath(cls, input_fpath, decoder_context):
        index_path = Path(input_fpath + ".metadata")
        if not index_path.exists():
            return None
        if not index_path.is_dir():
            raise NotADirectoryError(
                "metadata sidecar is not a directory: {}".format(index_path)
            )
        return cls(str(index_path), decoder_context)

    @staticmethod
    def _validate_interval(start_pos, end_pos):
        if start_pos > end_pos:
            raise ValueError("start position must be <= end position")

    def query_blk(self, start_pos, end_pos):
        self._validate_interval(start_pos, end_pos)
        if self.num_blocks == 0:
            return self.root_lookup

        # A block overlaps a query iff block_start <= query_end and
        # block_end >= query_start. This also handles queries before the first
        # block and gaps between blocks without negative searchsorted indices.
        overlaps = (
            (self.root_idx[:, 0] <= end_pos)
            & (self.root_idx[:, 1] >= start_pos)
        )
        return self.root_lookup[overlaps]

    def get_row_mask(self, block_id, start_pos, end_pos):
        self._validate_interval(start_pos, end_pos)
        if not isinstance(block_id, (int, np.integer)):
            raise TypeError("block_id must be an integer")
        block_id = int(block_id)
        if not 0 <= block_id < self.num_blocks:
            raise IndexError("block_id is outside the metadata index")

        curr_block_idx = self.block_idx[block_id]
        if curr_block_idx is None:
            curr_block_idx = np.load(
                join(self.index_dpath, "{}.npy".format(block_id)),
                allow_pickle=False,
            )
            if curr_block_idx.ndim != 1:
                raise ValueError("block position index must be one-dimensional")
            if curr_block_idx.size > 1 and np.any(
                curr_block_idx[1:] < curr_block_idx[:-1]
            ):
                raise ValueError("block position index must be sorted")
            self.block_idx[block_id] = curr_block_idx

        start_row = np.searchsorted(curr_block_idx, start_pos, side="left")
        end_row = np.searchsorted(curr_block_idx, end_pos, side="right")
        return slice(start_row, end_row, None)

    def query_columns(self, sample_ids):
        if sample_ids is None:
            return None
        if not isinstance(sample_ids, str):
            raise TypeError("sample_ids must be a semicolon-separated string or None")

        requested = [value for value in sample_ids.strip().split(";") if value]
        if not requested:
            raise ValueError("at least one sample ID must be provided")

        sample_to_index = {
            str(sample_id): index for index, sample_id in enumerate(self.samples)
        }
        unknown = [sample_id for sample_id in requested if sample_id not in sample_to_index]
        if unknown:
            raise ValueError(
                "unknown sample ID(s): {}".format(", ".join(unknown))
            )

        return np.asarray(
            [sample_to_index[sample_id] for sample_id in requested],
            dtype=np.uint32,
        )
