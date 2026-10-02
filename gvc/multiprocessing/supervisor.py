"""Parent-owned lifecycle management for parallel GVC encoding."""

import os
import queue
import shutil
import signal
import threading
import time
import traceback
import uuid
from pathlib import Path

from .messages import Progress, ReaderDone, WorkerError, WriterDone


class MultiprocessingEncodeError(RuntimeError):
    """Raised when any child in the encode pipeline fails."""

    def __init__(self, errors):
        self.errors = list(errors)
        lines = ["multiprocessing encoder failed"]
        for error in self.errors:
            location = error.stage
            if error.worker_id is not None:
                location += "[{}]".format(error.worker_id)
            if error.block_id is not None:
                location += " block {}".format(error.block_id)
            lines.append(
                "{}: {}: {}".format(
                    location, error.error_type, error.message
                )
            )
        super().__init__("; ".join(lines))


def child_entry(
    stage,
    worker_id,
    error_q,
    stop_event,
    target,
    args,
    initializer=None,
    initializer_args=(),
):
    """Run one supervised child with optional process-local initialization."""
    try:
        if initializer is not None:
            initializer(*initializer_args)
        target(*args)
    except BaseException as exc:
        try:
            error_q.put(
                WorkerError(
                    stage=stage,
                    worker_id=worker_id,
                    block_id=getattr(exc, "block_id", None),
                    error_type=type(exc).__name__,
                    message=str(exc),
                    traceback=traceback.format_exc(),
                ),
                timeout=1.0,
            )
        except Exception:
            pass
        stop_event.set()
        raise


class EncodeProcessSupervisor:
    """Own all process lifecycle, failure cancellation, and output commit."""

    def __init__(
        self,
        processes,
        error_q,
        status_q,
        stop_event,
        queues,
        temp_output,
        final_output,
        stall_timeout=None,
        poll_interval=0.1,
        graceful_timeout=2.0,
    ):
        self.processes = list(processes)
        self.error_q = error_q
        self.status_q = status_q
        self.stop_event = stop_event
        self.queues = list(queues)
        self.temp_output = Path(temp_output)
        self.final_output = Path(final_output)
        self.stall_timeout = stall_timeout
        self.poll_interval = poll_interval
        self.graceful_timeout = graceful_timeout
        self._errors = []
        self._reader_total = None
        self._writer_total = None
        self._last_progress = time.monotonic()

    @property
    def temp_metadata(self):
        return Path(str(self.temp_output) + ".metadata")

    @property
    def final_metadata(self):
        return Path(str(self.final_output) + ".metadata")

    def start(self):
        for proc in self.processes:
            proc.start()

    def _drain_messages(self):
        progressed = False
        while True:
            try:
                error = self.error_q.get_nowait()
            except queue.Empty:
                break
            else:
                self._errors.append(error)
                progressed = True

        while True:
            try:
                status = self.status_q.get_nowait()
            except queue.Empty:
                break
            else:
                progressed = True
                if isinstance(status, ReaderDone):
                    self._reader_total = status.total_blocks
                elif isinstance(status, WriterDone):
                    self._writer_total = status.total_blocks
                elif isinstance(status, Progress):
                    pass

        if progressed:
            self._last_progress = time.monotonic()

    def _drain_completion_messages(self, timeout=1.0):
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            self._drain_messages()
            if self._reader_total is not None and self._writer_total is not None:
                return
            time.sleep(min(self.poll_interval, 0.05))

    def _unexpected_exits(self):
        known = {
            (error.stage, error.worker_id)
            for error in self._errors
        }
        synthetic = []
        for proc in self.processes:
            if proc.exitcode in (None, 0):
                continue
            stage = getattr(proc, "_gvc_stage", proc.name)
            worker_id = getattr(proc, "_gvc_worker_id", None)
            if (stage, worker_id) not in known:
                synthetic.append(
                    WorkerError(
                        stage=stage,
                        worker_id=worker_id,
                        block_id=None,
                        error_type="ProcessExit",
                        message="child exited with status {}".format(
                            proc.exitcode
                        ),
                        traceback="",
                    )
                )
        return synthetic

    def wait(self):
        while True:
            self._drain_messages()
            synthetic = self._unexpected_exits()
            if synthetic:
                self._errors.extend(synthetic)

            if self._errors:
                self.cancel()
                raise MultiprocessingEncodeError(self._errors)

            if all(proc.exitcode is not None for proc in self.processes):
                self._drain_completion_messages()
                synthetic = self._unexpected_exits()
                if synthetic:
                    self._errors.extend(synthetic)
                    raise MultiprocessingEncodeError(self._errors)
                break

            if (
                self.stall_timeout is not None
                and time.monotonic() - self._last_progress
                > self.stall_timeout
            ):
                self._errors.append(
                    WorkerError(
                        stage="supervisor",
                        worker_id=None,
                        block_id=None,
                        error_type="TimeoutError",
                        message="no multiprocessing progress for {:.1f}s".format(
                            self.stall_timeout
                        ),
                        traceback="",
                    )
                )
                self.cancel()
                raise MultiprocessingEncodeError(self._errors)

            time.sleep(self.poll_interval)

        if self._reader_total is None:
            raise MultiprocessingEncodeError(
                [
                    WorkerError(
                        "reader",
                        None,
                        None,
                        "ProtocolError",
                        "reader completion count was not reported",
                        "",
                    )
                ]
            )
        if self._writer_total is None:
            raise MultiprocessingEncodeError(
                [
                    WorkerError(
                        "writer",
                        None,
                        None,
                        "ProtocolError",
                        "writer completion count was not reported",
                        "",
                    )
                ]
            )
        if self._reader_total != self._writer_total:
            raise MultiprocessingEncodeError(
                [
                    WorkerError(
                        "supervisor",
                        None,
                        None,
                        "BlockCountMismatch",
                        "reader produced {}, writer committed {}".format(
                            self._reader_total, self._writer_total
                        ),
                        "",
                    )
                ]
            )
        return self._writer_total

    @staticmethod
    def _started(proc):
        return getattr(proc, "pid", None) is not None

    def cancel(self):
        self.stop_event.set()
        deadline = time.monotonic() + self.graceful_timeout
        for proc in self.processes:
            if not self._started(proc):
                continue
            remaining = max(0.0, deadline - time.monotonic())
            proc.join(timeout=remaining)

        for proc in self.processes:
            if self._started(proc) and proc.is_alive():
                proc.terminate()
        for proc in self.processes:
            if self._started(proc):
                proc.join(timeout=1.0)

        for proc in self.processes:
            if (
                self._started(proc)
                and proc.is_alive()
                and hasattr(proc, "kill")
            ):
                proc.kill()
                proc.join(timeout=1.0)

    def join(self):
        for proc in self.processes:
            if self._started(proc):
                proc.join()

    def cleanup_ipc(self):
        for q in self.queues + [self.error_q, self.status_q]:
            # A force-terminated child may leave buffered queue data whose
            # reader no longer exists. Do not let parent shutdown block while
            # waiting for feeder threads to flush unreachable data.
            try:
                q.cancel_join_thread()
            except Exception:
                pass
            try:
                q.close()
            except Exception:
                pass

    def cleanup_temp(self):
        try:
            self.temp_output.unlink()
        except FileNotFoundError:
            pass
        if self.temp_metadata.exists():
            shutil.rmtree(str(self.temp_metadata), ignore_errors=True)

    @staticmethod
    def _backup_path(path):
        return Path(
            str(path)
            + ".gvc-backup-{}-{}".format(
                os.getpid(),
                uuid.uuid4().hex,
            )
        )

    def commit(self):
        """Replace final file + sidecar with rollback on partial failure."""
        self.final_output.parent.mkdir(parents=True, exist_ok=True)
        file_backup = self._backup_path(self.final_output)
        metadata_backup = self._backup_path(self.final_metadata)

        for backup in (file_backup, metadata_backup):
            if backup.is_dir():
                shutil.rmtree(str(backup))
            elif backup.exists():
                backup.unlink()

        file_had_old = self.final_output.exists()
        metadata_had_old = self.final_metadata.exists()

        try:
            if file_had_old:
                os.replace(str(self.final_output), str(file_backup))
            if metadata_had_old:
                os.replace(str(self.final_metadata), str(metadata_backup))

            os.replace(str(self.temp_output), str(self.final_output))
            if self.temp_metadata.exists():
                os.replace(str(self.temp_metadata), str(self.final_metadata))

        except BaseException:
            try:
                if self.final_output.exists():
                    self.final_output.unlink()
                if self.final_metadata.exists():
                    shutil.rmtree(str(self.final_metadata))
                if file_had_old and file_backup.exists():
                    os.replace(str(file_backup), str(self.final_output))
                if metadata_had_old and metadata_backup.exists():
                    os.replace(str(metadata_backup), str(self.final_metadata))
            finally:
                self.cleanup_temp()
            raise
        else:
            if file_backup.exists():
                file_backup.unlink()
            if metadata_backup.exists():
                shutil.rmtree(str(metadata_backup))

    def _install_parent_signal_handlers(self):
        if threading.current_thread() is not threading.main_thread():
            return {}

        previous = {}
        for signum in (getattr(signal, "SIGTERM", None),):
            if signum is None:
                continue
            previous[signum] = signal.getsignal(signum)

            def _handler(received, frame, self=self):
                self.stop_event.set()
                raise MultiprocessingEncodeError(
                    [
                        WorkerError(
                            stage="supervisor",
                            worker_id=None,
                            block_id=None,
                            error_type="SignalTermination",
                            message="parent received signal {}".format(received),
                            traceback="",
                        )
                    ]
                )

            signal.signal(signum, _handler)
        return previous

    @staticmethod
    def _restore_parent_signal_handlers(previous):
        for signum, handler in previous.items():
            signal.signal(signum, handler)

    def run(self):
        self.cleanup_temp()
        previous_handlers = {}
        try:
            self.start()
            previous_handlers = self._install_parent_signal_handlers()
            count = self.wait()
            self.join()
            self.commit()
            return count
        except BaseException:
            self.cancel()
            self.cleanup_temp()
            raise
        finally:
            self._restore_parent_signal_handlers(previous_handlers)
            self.cleanup_ipc()
