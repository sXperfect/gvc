import pytest

from gvc.encoder import _next_parameter_set_id


def test_parameter_set_id_limit_is_enforced_before_mutation():
    assert _next_parameter_set_id([None] * 255) == 255

    with pytest.raises(ValueError, match="serialized ID limit"):
        _next_parameter_set_id([None] * 256)
