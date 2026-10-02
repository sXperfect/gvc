from io import BytesIO
from types import SimpleNamespace

import numpy as np
import pytest

from gvc import decoder
from gvc.data_structures import GenotypePayload, ParameterSet
from gvc.data_structures.consts import BinarizationID, CodecID


def _row_split_parameter_set():
    return ParameterSet(
        parameter_set_id=0,
        any_missing_flag=False,
        not_available_flag=False,
        p=2,
        binarization_id=BinarizationID.ROW_BIN_SPLIT,
        num_bin_mat=1,
        concat_axis=0,
        sort_variants_row_flags=[False],
        sort_variants_col_flags=[False],
        transpose_variants_mat_flags=[False],
        variants_coder_ids=[CodecID.JBIG1],
        encode_phase_data=False,
        phase_value=False,
    )


def test_amax_read_error_is_not_swallowed(monkeypatch):
    class BrokenAMax:
        def __len__(self):
            return 1

        def read(self):
            raise RuntimeError("synthetic AMax read failure")

    param_set = _row_split_parameter_set()
    payload = GenotypePayload(
        param_set,
        variants_payloads=[b"x"],
        variants_row_ids_payloads=[None],
        variants_col_ids_payloads=[None],
        variants_amax_payload=BrokenAMax(),
    )

    monkeypatch.setattr(
        decoder.codec,
        "decode",
        lambda *args, **kwargs: np.zeros((1, 2), dtype=bool),
    )

    with pytest.raises(RuntimeError, match="synthetic AMax read failure"):
        decoder.decode_encoded_variants(param_set, payload)


def test_decoder_context_missing_access_unit_raises_value_error():
    context = decoder.DecoderContext()

    with pytest.raises(ValueError, match="No access unit with id 7"):
        context.set_access_unit(7)


def test_cache_access_unit_rejects_inconsistent_sample_count(monkeypatch):
    instance = decoder.Decoder.__new__(decoder.Decoder)
    instance._f = BytesIO(b"")
    instance._bitstream_reader = object()
    instance.decoder_context = decoder.DecoderContext()
    instance.decoder_context.ncols = 2
    instance.decoder_context.parameter_sets[0] = SimpleNamespace(p=2)

    fake_block = SimpleNamespace(block_payload=object())
    fake_access_unit = SimpleNamespace(
        header=SimpleNamespace(access_unit_id=0, parameter_set_id=0),
        blocks=[fake_block],
    )

    monkeypatch.setattr(
        decoder.ds.AccessUnit,
        "from_bitstream",
        lambda *args, **kwargs: fake_access_unit,
    )
    monkeypatch.setattr(
        decoder,
        "_get_tensor_shape",
        lambda *args, **kwargs: (1, 3, 2),
    )

    with pytest.raises(ValueError, match="tensor sample count is inconsistent"):
        instance._cache_access_unit()



def test_amax_internal_attribute_error_is_not_swallowed(monkeypatch):
    class BrokenAMax:
        def __len__(self):
            return 1

        def read(self):
            raise AttributeError("synthetic internal attribute failure")

    param_set = _row_split_parameter_set()
    payload = GenotypePayload(
        param_set,
        variants_payloads=[b"x"],
        variants_row_ids_payloads=[None],
        variants_col_ids_payloads=[None],
        variants_amax_payload=BrokenAMax(),
    )
    monkeypatch.setattr(
        decoder.codec,
        "decode",
        lambda *args, **kwargs: np.zeros((1, 2), dtype=bool),
    )

    with pytest.raises(AttributeError, match="synthetic internal attribute failure"):
        decoder.decode_encoded_variants(param_set, payload)
