import multiprocessing as mp
from queue import Queue
from types import SimpleNamespace

import pytest

from gvc.codec import MAT_CODECS
from gvc.data_structures.consts import BinarizationID, CodecID
from gvc.decoder import Decoder
from gvc.encoder import Encoder, run_multiprocessing, worker_writer
from gvc.multiprocessing import MultiprocessingEncodeError

from tests.mp_test_codec import decode as mp_decode
from tests.mp_test_codec import encode as mp_encode
from tests.mp_test_codec import install as install_mp_codec
from tests.test_file_roundtrip import (
    EXPECTED_GT,
    VCF_FIXTURE,
    _close_decoder,
)


class DummyEvent:
    def __init__(self):
        self._set = False

    def is_set(self):
        return self._set

    def set(self):
        self._set = True


def _status_queue():
    return Queue()


def test_multiprocessing_encoder_matches_single_process(monkeypatch, tmp_path):
    if mp.get_start_method() != "fork":
        pytest.skip("real multiprocessing codec test requires fork start method")

    codec = MAT_CODECS[CodecID.JBIG1]
    monkeypatch.setitem(codec, "encoder", mp_encode)
    monkeypatch.setitem(codec, "decoder", mp_decode)

    single = tmp_path / "single.gvc"
    multi = tmp_path / "multi.gvc"

    common = dict(
        binarization_name="bit_plane",
        axis=2,
        sort_rows=False,
        sort_cols=False,
        transpose=False,
        block_size=1,
        codec_name="jbig",
        preset_mode=0,
    )
    Encoder(str(VCF_FIXTURE), str(single), num_threads=0, **common).run()
    Encoder(str(VCF_FIXTURE), str(multi), num_threads=2, **common).run()

    assert multi.read_bytes() == single.read_bytes()

    decoder = Decoder(str(multi), str(tmp_path / "decoded.txt"))
    try:
        decoder.decode()
    finally:
        _close_decoder(decoder)
    assert (tmp_path / "decoded.txt").read_text() == EXPECTED_GT


def test_writer_orders_out_of_order_worker_results(monkeypatch, tmp_path):
    from gvc import encoder as encoder_module
    from gvc.multiprocessing import EncodedBlock, WorkerDone, WriterDone

    queue = Queue()
    status = _status_queue()
    stop = DummyEvent()
    stored = []

    class FakeParameterSet:
        parameter_set_id = 0

        def __eq__(self, other):
            return isinstance(other, FakeParameterSet)

        def to_bytes(self):
            return b"P"

    class FakeBlock:
        def __init__(self, value):
            self.value = value

        def __len__(self):
            return 1

    parameter_set = FakeParameterSet()

    def fake_store_access_unit(output_f, access_unit_id, param_set, blocks):
        values = [block.value for block in blocks]
        stored.extend(values)
        output_f.write(bytes(values))

    monkeypatch.setattr(
        encoder_module.gvc.common,
        "store_access_unit",
        fake_store_access_unit,
    )

    queue.put(EncodedBlock(2, parameter_set, FakeBlock(2)))
    queue.put(EncodedBlock(0, parameter_set, FakeBlock(0)))
    queue.put(EncodedBlock(1, parameter_set, FakeBlock(1)))
    queue.put(WorkerDone(0))
    queue.put(WorkerDone(1))

    output = tmp_path / "ordered.gvc"
    worker_writer(
        queue,
        status,
        stop,
        str(output),
        num_processes=2,
    )

    assert stored == [0, 1, 2]
    assert output.read_bytes() == b"P\x00\x01\x02"
    statuses = []
    while not status.empty():
        statuses.append(status.get_nowait())
    done = [item for item in statuses if isinstance(item, WriterDone)]
    assert len(done) == 1
    assert done[0].total_blocks == 3


def test_writer_rejects_missing_block(tmp_path):
    from gvc.multiprocessing import EncodedBlock, WorkerDone

    queue = Queue()
    status = _status_queue()
    stop = DummyEvent()
    parameter_set = type(
        "FakeParameterSet",
        (),
        {
            "parameter_set_id": 0,
            "to_bytes": lambda self: b"P",
            "__eq__": lambda self, other: True,
        },
    )()
    block = type("FakeBlock", (), {"__len__": lambda self: 1})()

    queue.put(EncodedBlock(1, parameter_set, block))
    queue.put(WorkerDone(0))

    with pytest.raises(RuntimeError, match="missing encoded block"):
        worker_writer(
            queue,
            status,
            stop,
            str(tmp_path / "missing.gvc"),
            num_processes=1,
        )


def test_reader_failure_preserves_existing_final_artifacts(tmp_path):
    final = tmp_path / "existing.gvc"
    final.write_bytes(b"ORIGINAL")
    metadata = tmp_path / "existing.gvc.metadata"
    metadata.mkdir()
    (metadata / "marker.txt").write_text("ORIGINAL-METADATA")

    before_children = {proc.pid for proc in mp.active_children()}
    with pytest.raises(MultiprocessingEncodeError, match="reader"):
        run_multiprocessing(
            str(tmp_path / "not-a-vcf.txt"),
            str(final),
            block_size=2,
            ps_params=[
                BinarizationID.BIT_PLANE,
                CodecID.JBIG1,
                2,
                False,
                False,
                False,
            ],
            tsp_params=["ham", "nn", 0],
            num_processes=2,
            stall_timeout=5,
        )

    assert final.read_bytes() == b"ORIGINAL"
    assert (metadata / "marker.txt").read_text() == "ORIGINAL-METADATA"
    assert not list(tmp_path.glob("existing.gvc.tmp.*"))
    assert not list(tmp_path.glob("existing.gvc.tmp.*.metadata"))

    after_children = {proc.pid for proc in mp.active_children()}
    assert after_children <= before_children


def test_parallel_encoder_rejects_invalid_worker_count(tmp_path):
    with pytest.raises(ValueError, match="positive"):
        run_multiprocessing(
            str(VCF_FIXTURE),
            str(tmp_path / "bad.gvc"),
            2,
            [],
            [],
            0,
        )


def test_supervisor_error_carries_stage_information(tmp_path):
    with pytest.raises(MultiprocessingEncodeError) as exc_info:
        run_multiprocessing(
            str(tmp_path / "invalid.input"),
            str(tmp_path / "never.gvc"),
            block_size=1,
            ps_params=[],
            tsp_params=[],
            num_processes=1,
            stall_timeout=5,
        )

    assert exc_info.value.errors
    assert any(error.stage == "reader" for error in exc_info.value.errors)
    assert "Invalid Format" in str(exc_info.value)



def test_worker_failure_is_atomic_and_leaves_no_children(tmp_path):
    final = tmp_path / "worker-failure.gvc"
    final.write_bytes(b"ORIGINAL")
    metadata = tmp_path / "worker-failure.gvc.metadata"
    metadata.mkdir()
    (metadata / "marker.txt").write_text("ORIGINAL-METADATA")

    before = {proc.pid for proc in mp.active_children()}
    with pytest.raises(MultiprocessingEncodeError) as exc_info:
        run_multiprocessing(
            str(VCF_FIXTURE),
            str(final),
            block_size=1,
            # Intentionally malformed: reader succeeds and encoder fails once
            # the first WorkItem is consumed.
            ps_params=[],
            tsp_params=[],
            num_processes=2,
            stall_timeout=10,
        )

    assert any(error.stage == "encoder" for error in exc_info.value.errors)
    assert final.read_bytes() == b"ORIGINAL"
    assert (metadata / "marker.txt").read_text() == "ORIGINAL-METADATA"
    assert not list(tmp_path.glob("worker-failure.gvc.tmp.*"))
    assert not list(tmp_path.glob("worker-failure.gvc.tmp.*.metadata"))
    assert {proc.pid for proc in mp.active_children()} <= before


def test_spawn_mode_failure_supervision_is_picklable(tmp_path):
    if "spawn" not in mp.get_all_start_methods():
        pytest.skip("spawn start method is unavailable")

    with pytest.raises(MultiprocessingEncodeError):
        run_multiprocessing(
            str(VCF_FIXTURE),
            str(tmp_path / "spawn-failure.gvc"),
            block_size=1,
            ps_params=[],
            tsp_params=[],
            num_processes=1,
            start_method="spawn",
            stall_timeout=15,
        )

    assert not (tmp_path / "spawn-failure.gvc").exists()
    assert not list(tmp_path.glob("spawn-failure.gvc.tmp.*"))


def test_transaction_commit_rolls_back_file_and_metadata(monkeypatch, tmp_path):
    import queue as local_queue

    from gvc.multiprocessing.supervisor import EncodeProcessSupervisor
    import gvc.multiprocessing.supervisor as supervisor_module

    final = tmp_path / "atomic.gvc"
    final.write_bytes(b"OLD-FILE")
    final_metadata = tmp_path / "atomic.gvc.metadata"
    final_metadata.mkdir()
    (final_metadata / "marker").write_text("OLD-METADATA")

    temp = tmp_path / "atomic.gvc.tmp.123"
    temp.write_bytes(b"NEW-FILE")
    temp_metadata = tmp_path / "atomic.gvc.tmp.123.metadata"
    temp_metadata.mkdir()
    (temp_metadata / "marker").write_text("NEW-METADATA")

    supervisor = EncodeProcessSupervisor(
        processes=[],
        error_q=local_queue.Queue(),
        status_q=local_queue.Queue(),
        stop_event=DummyEvent(),
        queues=[],
        temp_output=temp,
        final_output=final,
    )

    real_replace = supervisor_module.os.replace

    def fail_metadata_commit(src, dst):
        if str(src) == str(temp_metadata) and str(dst) == str(final_metadata):
            raise OSError("simulated metadata commit failure")
        return real_replace(src, dst)

    monkeypatch.setattr(supervisor_module.os, "replace", fail_metadata_commit)

    with pytest.raises(OSError, match="metadata commit failure"):
        supervisor.commit()

    assert final.read_bytes() == b"OLD-FILE"
    assert (final_metadata / "marker").read_text() == "OLD-METADATA"
    assert not temp.exists()
    assert not temp_metadata.exists()



def test_encoder_exposes_supervisor_configuration(tmp_path):
    encoder = Encoder(
        str(VCF_FIXTURE),
        str(tmp_path / "config.gvc"),
        num_threads=2,
        multiprocessing_start_method="spawn",
        multiprocessing_stall_timeout=30,
    )
    assert encoder.multiprocessing_start_method == "spawn"
    assert encoder.multiprocessing_stall_timeout == 30
    assert encoder.multiprocessing_initializer is None
    assert encoder.multiprocessing_initializer_args == ()



def _sleep_forever():
    import time

    time.sleep(60)


def test_supervisor_watchdog_terminates_stalled_child(tmp_path):
    import multiprocessing as local_mp

    from gvc.multiprocessing.supervisor import EncodeProcessSupervisor

    context = local_mp.get_context()
    error_q = context.Queue()
    status_q = context.Queue()
    stop_event = context.Event()
    stalled = context.Process(name="GVC-Stalled", target=_sleep_forever)
    stalled._gvc_stage = "test"
    stalled._gvc_worker_id = 0

    supervisor = EncodeProcessSupervisor(
        processes=[stalled],
        error_q=error_q,
        status_q=status_q,
        stop_event=stop_event,
        queues=[],
        temp_output=tmp_path / "stalled.tmp",
        final_output=tmp_path / "stalled.gvc",
        stall_timeout=0.2,
        poll_interval=0.05,
        graceful_timeout=0.1,
    )

    with pytest.raises(MultiprocessingEncodeError, match="no multiprocessing progress"):
        supervisor.run()

    assert not stalled.is_alive()
    assert not (tmp_path / "stalled.gvc").exists()



def test_parallel_encoder_rejects_unsafe_output_paths(tmp_path):
    output_dir = tmp_path / "directory.gvc"
    output_dir.mkdir()
    with pytest.raises(IsADirectoryError):
        run_multiprocessing(
            str(VCF_FIXTURE),
            str(output_dir),
            2,
            [],
            [],
            1,
        )

    missing_parent = tmp_path / "does-not-exist" / "out.gvc"
    with pytest.raises(FileNotFoundError, match="output directory"):
        run_multiprocessing(
            str(VCF_FIXTURE),
            str(missing_parent),
            2,
            [],
            [],
            1,
        )

    final = tmp_path / "metadata-file.gvc"
    sidecar = tmp_path / "metadata-file.gvc.metadata"
    sidecar.write_text("not-a-directory")
    with pytest.raises(NotADirectoryError, match="metadata"):
        run_multiprocessing(
            str(VCF_FIXTURE),
            str(final),
            2,
            [],
            [],
            1,
        )



def test_spawn_mode_successful_encode_roundtrip(monkeypatch, tmp_path):
    if "spawn" not in mp.get_all_start_methods():
        pytest.skip("spawn start method is unavailable")

    # Parent needs the decoder; spawned encoder children rebuild their own
    # registry using the explicit top-level initializer.
    codec = MAT_CODECS[CodecID.JBIG1]
    monkeypatch.setitem(codec, "encoder", mp_encode)
    monkeypatch.setitem(codec, "decoder", mp_decode)

    encoded = tmp_path / "spawn.gvc"
    decoded = tmp_path / "spawn.txt"

    Encoder(
        str(VCF_FIXTURE),
        str(encoded),
        binarization_name="bit_plane",
        axis=2,
        sort_rows=False,
        sort_cols=False,
        block_size=1,
        codec_name="jbig",
        preset_mode=0,
        num_threads=2,
        multiprocessing_start_method="spawn",
        multiprocessing_stall_timeout=20,
        multiprocessing_initializer=install_mp_codec,
    ).run()

    decoder = Decoder(str(encoded), str(decoded))
    try:
        decoder.decode()
    finally:
        _close_decoder(decoder)

    assert decoded.read_text() == EXPECTED_GT
    assert (tmp_path / "spawn.gvc.metadata" / "main.npy").is_file()


@pytest.mark.parametrize(
    "start_method",
    [
        method
        for method in ("spawn", "forkserver")
        if method in mp.get_all_start_methods()
    ],
)
def test_nonfork_initializer_must_be_picklable(tmp_path, start_method):
    def local_initializer():
        pass

    # A nested initializer is intentionally not picklable. The parent should
    # reject it before launching any child under contexts that serialize
    # process state.
    output = tmp_path / ("unpicklable-{}.gvc".format(start_method))
    with pytest.raises(TypeError, match="picklable"):
        run_multiprocessing(
            str(VCF_FIXTURE),
            str(output),
            block_size=1,
            ps_params=[
                BinarizationID.BIT_PLANE,
                CodecID.JBIG1,
                2,
                False,
                False,
                False,
            ],
            tsp_params=["ham", "nn", 0],
            num_processes=1,
            start_method=start_method,
            process_initializer=local_initializer,
            stall_timeout=10,
        )

    assert not output.exists()



def test_writer_io_failure_is_structured(monkeypatch, tmp_path):
    from gvc import encoder as encoder_module
    from gvc.multiprocessing import EncodedBlock, WorkerDone
    from gvc.multiprocessing.supervisor import child_entry

    queue = Queue()
    status = _status_queue()
    stop = DummyEvent()
    error_q = Queue()

    class FakeParameterSet:
        parameter_set_id = 0

        def __eq__(self, other):
            return isinstance(other, FakeParameterSet)

        def to_bytes(self):
            return b"P"

    class FakeBlock:
        def __len__(self):
            return 1

    def failing_store(*args, **kwargs):
        raise OSError("synthetic disk write failure")

    monkeypatch.setattr(
        encoder_module.gvc.common,
        "store_access_unit",
        failing_store,
    )

    queue.put(EncodedBlock(0, FakeParameterSet(), FakeBlock()))
    queue.put(WorkerDone(0))

    with pytest.raises(OSError, match="disk write failure"):
        child_entry(
            "writer",
            None,
            error_q,
            stop,
            worker_writer,
            (
                queue,
                status,
                stop,
                str(tmp_path / "writer-failure.gvc"),
                1,
            ),
        )

    error = error_q.get_nowait()
    assert error.stage == "writer"
    assert error.error_type == "OSError"
    assert "disk write failure" in error.message
    assert stop.is_set()


def test_supervisor_sigterm_handler_is_structured_and_restored(monkeypatch, tmp_path):
    import signal

    from gvc.multiprocessing.supervisor import EncodeProcessSupervisor

    error_q = Queue()
    status_q = Queue()
    stop_event = DummyEvent()
    supervisor = EncodeProcessSupervisor(
        processes=[],
        error_q=error_q,
        status_q=status_q,
        stop_event=stop_event,
        queues=[],
        temp_output=tmp_path / "signal.tmp",
        final_output=tmp_path / "signal.gvc",
    )

    previous_handler = object()
    installed = {}

    monkeypatch.setattr(signal, "getsignal", lambda signum: previous_handler)

    def fake_signal(signum, handler):
        installed[signum] = handler

    monkeypatch.setattr(signal, "signal", fake_signal)

    previous = supervisor._install_parent_signal_handlers()
    assert signal.SIGTERM in previous
    assert previous[signal.SIGTERM] is previous_handler

    with pytest.raises(MultiprocessingEncodeError, match="SignalTermination") as exc_info:
        installed[signal.SIGTERM](signal.SIGTERM, None)

    assert stop_event.is_set()
    assert exc_info.value.errors[0].stage == "supervisor"
    assert exc_info.value.errors[0].error_type == "SignalTermination"

    supervisor._restore_parent_signal_handlers(previous)
    assert installed[signal.SIGTERM] is previous_handler


def test_supervisor_run_restores_signal_handler_after_failure(monkeypatch, tmp_path):
    from gvc.multiprocessing.supervisor import EncodeProcessSupervisor

    supervisor = EncodeProcessSupervisor(
        processes=[],
        error_q=Queue(),
        status_q=Queue(),
        stop_event=DummyEvent(),
        queues=[],
        temp_output=tmp_path / "restore.tmp",
        final_output=tmp_path / "restore.gvc",
    )

    restored = []
    monkeypatch.setattr(supervisor, "start", lambda: None)
    monkeypatch.setattr(
        supervisor,
        "_install_parent_signal_handlers",
        lambda: {15: "previous"},
    )
    monkeypatch.setattr(
        supervisor,
        "_restore_parent_signal_handlers",
        lambda previous: restored.append(previous),
    )
    monkeypatch.setattr(
        supervisor,
        "wait",
        lambda: (_ for _ in ()).throw(RuntimeError("synthetic failure")),
    )
    monkeypatch.setattr(supervisor, "cancel", lambda: None)
    monkeypatch.setattr(supervisor, "cleanup_temp", lambda: None)
    monkeypatch.setattr(supervisor, "cleanup_ipc", lambda: None)

    with pytest.raises(RuntimeError, match="synthetic failure"):
        supervisor.run()

    assert restored == [{15: "previous"}]


def test_queue_put_stops_under_backpressure_when_cancelled():
    from gvc.encoder import _queue_put

    class FullQueue:
        def __init__(self, stop_event):
            self.calls = 0
            self.stop_event = stop_event

        def put(self, item, timeout=None):
            import queue as local_queue
            self.calls += 1
            self.stop_event.set()
            raise local_queue.Full

    stop = DummyEvent()
    target = FullQueue(stop)

    assert _queue_put(target, object(), stop, timeout=0.01) is False
    assert target.calls == 1


def test_queue_get_stops_when_cancelled_after_empty_poll():
    from gvc.encoder import _queue_get

    class EmptyQueue:
        def __init__(self, stop_event):
            self.calls = 0
            self.stop_event = stop_event

        def get(self, timeout=None):
            import queue as local_queue
            self.calls += 1
            self.stop_event.set()
            raise local_queue.Empty

    stop = DummyEvent()
    source = EmptyQueue(stop)

    assert _queue_get(source, stop, timeout=0.01) is None
    assert source.calls == 1


def test_supervisor_cancel_escalates_to_terminate_and_kill(tmp_path):
    from gvc.multiprocessing.supervisor import EncodeProcessSupervisor

    class FakeProcess:
        def __init__(self):
            self.pid = 1234
            self.alive = True
            self.join_calls = []
            self.terminate_calls = 0
            self.kill_calls = 0

        def join(self, timeout=None):
            self.join_calls.append(timeout)

        def is_alive(self):
            return self.alive

        def terminate(self):
            self.terminate_calls += 1

        def kill(self):
            self.kill_calls += 1
            self.alive = False

    proc = FakeProcess()
    stop = DummyEvent()
    supervisor = EncodeProcessSupervisor(
        processes=[proc],
        error_q=Queue(),
        status_q=Queue(),
        stop_event=stop,
        queues=[],
        temp_output=tmp_path / "cancel.tmp",
        final_output=tmp_path / "cancel.gvc",
        graceful_timeout=0,
    )

    supervisor.cancel()

    assert stop.is_set()
    assert proc.terminate_calls == 1
    assert proc.kill_calls == 1
    assert proc.alive is False
    assert len(proc.join_calls) >= 2


def _abrupt_exit():
    import os
    os._exit(7)


def test_supervisor_synthesizes_abrupt_child_exit(tmp_path):
    import multiprocessing as local_mp

    from gvc.multiprocessing.supervisor import EncodeProcessSupervisor

    context = local_mp.get_context()
    proc = context.Process(name="GVC-Abrupt", target=_abrupt_exit)
    proc._gvc_stage = "encoder"
    proc._gvc_worker_id = 3

    supervisor = EncodeProcessSupervisor(
        processes=[proc],
        error_q=context.Queue(),
        status_q=context.Queue(),
        stop_event=context.Event(),
        queues=[],
        temp_output=tmp_path / "abrupt.tmp",
        final_output=tmp_path / "abrupt.gvc",
        poll_interval=0.02,
        graceful_timeout=0.05,
    )

    with pytest.raises(MultiprocessingEncodeError) as exc_info:
        supervisor.run()

    assert any(
        error.stage == "encoder"
        and error.worker_id == 3
        and error.error_type == "ProcessExit"
        and "status 7" in error.message
        for error in exc_info.value.errors
    )
    assert not proc.is_alive()
    assert not (tmp_path / "abrupt.gvc").exists()


@pytest.mark.skipif(
    not hasattr(__import__("signal"), "SIGTERM"),
    reason="SIGTERM is unavailable",
)
def test_real_parent_sigterm_cleans_child_processes_and_temp_artifacts(tmp_path):
    import os
    import signal
    import subprocess
    import sys
    import textwrap
    from pathlib import Path

    marker = tmp_path / "child.pid"
    temp_output = tmp_path / "signal-integration.tmp"
    final_output = tmp_path / "signal-integration.gvc"

    script = textwrap.dedent(
        """
        import multiprocessing as mp
        import os
        from pathlib import Path

        from gvc.multiprocessing.supervisor import EncodeProcessSupervisor
        from tests.mp_test_codec import signal_child

        marker = Path({marker!r})
        temp_output = Path({temp_output!r})
        final_output = Path({final_output!r})

        if __name__ == "__main__":
            context = mp.get_context()
            proc = context.Process(
                name="GVC-Signal-Child",
                target=signal_child,
                args=(str(marker),),
            )
            proc._gvc_stage = "encoder"
            proc._gvc_worker_id = 0
            supervisor = EncodeProcessSupervisor(
                processes=[proc],
                error_q=context.Queue(),
                status_q=context.Queue(),
                stop_event=context.Event(),
                queues=[],
                temp_output=temp_output,
                final_output=final_output,
                poll_interval=0.02,
                graceful_timeout=0.05,
            )
            try:
                supervisor.run()
            except BaseException:
                raise SystemExit(143)
        """.format(
            marker=str(marker),
            temp_output=str(temp_output),
            final_output=str(final_output),
        )
    )

    process = subprocess.Popen(
        [sys.executable, "-c", script],
        cwd=str(Path(__file__).resolve().parents[1]),
    )
    try:
        deadline = __import__("time").time() + 10
        while not marker.exists() and __import__("time").time() < deadline:
            __import__("time").sleep(0.05)
        assert marker.exists(), "child process did not start"

        child_pid = int(marker.read_text())
        os.kill(process.pid, signal.SIGTERM)
        assert process.wait(timeout=10) == 143

        with pytest.raises(ProcessLookupError):
            os.kill(child_pid, 0)
        assert not temp_output.exists()
        assert not final_output.exists()
    finally:
        if process.poll() is None:
            process.kill()
            process.wait(timeout=5)


def test_supervisor_rejects_missing_completion_messages(tmp_path):
    from gvc.multiprocessing.supervisor import EncodeProcessSupervisor

    supervisor = EncodeProcessSupervisor(
        processes=[],
        error_q=Queue(),
        status_q=Queue(),
        stop_event=DummyEvent(),
        queues=[],
        temp_output=tmp_path / "missing-completion.tmp",
        final_output=tmp_path / "missing-completion.gvc",
        poll_interval=0.01,
    )

    with pytest.raises(MultiprocessingEncodeError, match="reader completion count"):
        supervisor.wait()


def test_supervisor_rejects_reader_writer_block_count_mismatch(tmp_path):
    from gvc.multiprocessing import ReaderDone, WriterDone
    from gvc.multiprocessing.supervisor import EncodeProcessSupervisor

    status = Queue()
    status.put(ReaderDone(5))
    status.put(WriterDone(4))

    supervisor = EncodeProcessSupervisor(
        processes=[],
        error_q=Queue(),
        status_q=status,
        stop_event=DummyEvent(),
        queues=[],
        temp_output=tmp_path / "mismatch.tmp",
        final_output=tmp_path / "mismatch.gvc",
        poll_interval=0.01,
    )

    with pytest.raises(MultiprocessingEncodeError, match="reader produced 5, writer committed 4"):
        supervisor.wait()


def test_transaction_commit_replaces_file_and_metadata_without_backup_leaks(tmp_path):
    import queue as local_queue

    from gvc.multiprocessing.supervisor import EncodeProcessSupervisor

    final = tmp_path / "replace.gvc"
    final.write_bytes(b"OLD-FILE")
    final_metadata = tmp_path / "replace.gvc.metadata"
    final_metadata.mkdir()
    (final_metadata / "marker").write_text("OLD-METADATA")

    temp = tmp_path / "replace.gvc.tmp.123"
    temp.write_bytes(b"NEW-FILE")
    temp_metadata = tmp_path / "replace.gvc.tmp.123.metadata"
    temp_metadata.mkdir()
    (temp_metadata / "marker").write_text("NEW-METADATA")

    supervisor = EncodeProcessSupervisor(
        processes=[],
        error_q=local_queue.Queue(),
        status_q=local_queue.Queue(),
        stop_event=DummyEvent(),
        queues=[],
        temp_output=temp,
        final_output=final,
    )

    supervisor.commit()

    assert final.read_bytes() == b"NEW-FILE"
    assert (final_metadata / "marker").read_text() == "NEW-METADATA"
    assert not temp.exists()
    assert not temp_metadata.exists()
    assert not list(tmp_path.glob("replace.gvc.gvc-backup-*"))
    assert not list(tmp_path.glob("replace.gvc.metadata.gvc-backup-*"))


def test_cleanup_ipc_closes_all_queues_even_when_one_cleanup_fails(tmp_path):
    from gvc.multiprocessing.supervisor import EncodeProcessSupervisor

    class FakeQueue:
        def __init__(self, fail_cancel=False, fail_close=False):
            self.cancel_calls = 0
            self.close_calls = 0
            self.fail_cancel = fail_cancel
            self.fail_close = fail_close

        def cancel_join_thread(self):
            self.cancel_calls += 1
            if self.fail_cancel:
                raise RuntimeError("cancel failed")

        def close(self):
            self.close_calls += 1
            if self.fail_close:
                raise RuntimeError("close failed")

    q1 = FakeQueue(fail_cancel=True)
    q2 = FakeQueue(fail_close=True)
    error_q = FakeQueue()
    status_q = FakeQueue()

    supervisor = EncodeProcessSupervisor(
        processes=[],
        error_q=error_q,
        status_q=status_q,
        stop_event=DummyEvent(),
        queues=[q1, q2],
        temp_output=tmp_path / "ipc.tmp",
        final_output=tmp_path / "ipc.gvc",
    )

    supervisor.cleanup_ipc()

    for q in (q1, q2, error_q, status_q):
        assert q.cancel_calls == 1
        assert q.close_calls == 1


def test_abrupt_child_exit_terminates_surviving_sibling(tmp_path):
    import multiprocessing as local_mp

    from gvc.multiprocessing.supervisor import EncodeProcessSupervisor

    context = local_mp.get_context()
    failed = context.Process(name="GVC-Failed", target=_abrupt_exit)
    failed._gvc_stage = "encoder"
    failed._gvc_worker_id = 0

    sibling = context.Process(name="GVC-Sibling", target=_sleep_forever)
    sibling._gvc_stage = "encoder"
    sibling._gvc_worker_id = 1

    supervisor = EncodeProcessSupervisor(
        processes=[failed, sibling],
        error_q=context.Queue(),
        status_q=context.Queue(),
        stop_event=context.Event(),
        queues=[],
        temp_output=tmp_path / "sibling.tmp",
        final_output=tmp_path / "sibling.gvc",
        poll_interval=0.02,
        graceful_timeout=0.05,
    )

    with pytest.raises(MultiprocessingEncodeError) as exc_info:
        supervisor.run()

    assert any(
        error.worker_id == 0 and error.error_type == "ProcessExit"
        for error in exc_info.value.errors
    )
    assert not failed.is_alive()
    assert not sibling.is_alive()
    assert not (tmp_path / "sibling.gvc").exists()
