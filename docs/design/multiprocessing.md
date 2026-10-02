# Multiprocessing design

GVC parallel encoding follows a parent-supervised reader/encoder/writer
pipeline. The compatibility requirement for the 1.0.x line is:

> Parallel encoding must either produce the same ordered GVC stream as
> sequential encoding, or fail atomically without publishing a partial final
> GVC file, leaving stale temporary artifacts, or leaking child processes.

## Ownership model

The parent process is the only lifecycle supervisor. Children process work but
do not decide global shutdown.

The pipeline is:

```text
Parent / EncodeProcessSupervisor
        |
        +-- Reader
        |     |
        |   work_q
        |     |
        +-- Encoder 0 --+
        +-- Encoder 1 --+--> result_q --> Writer
        +-- ...         |
        |
        +-- error_q
        +-- status_q / progress heartbeats
        +-- stop_event
```

The writer is the only child allowed to mutate the temporary GVC stream.
Encoder workers return complete encoded blocks tagged with their source block
ID. The writer buffers out-of-order results and commits blocks strictly in
ascending block-ID order.

## Message protocol

The queues use explicit picklable messages rather than sentinel tuples:

- `WorkItem`
- `EncodedBlock`
- `StopWork`
- `WorkerDone`
- `ReaderDone`
- `WriterDone`
- `Progress`
- `WorkerError`

No synchronization relies on `Queue.qsize()` or `Queue.empty()`.

## Backpressure and cancellation

Work and result queues are bounded. Queue operations use finite timeouts so a
producer can observe `stop_event` when a downstream consumer fails.

Any supervised child exception is converted into a structured `WorkerError`,
the shared stop event is set, and the parent performs:

1. graceful join for a bounded interval;
2. `terminate()` for remaining children;
3. `kill()` where available if termination does not complete;
4. temporary artifact cleanup.

Non-zero child exits without a reported exception are synthesized into
structured process-exit errors.

## Completion invariants

The reader reports its total number of source blocks and the writer reports its
total committed blocks. The parent commits output only when both values are
available and equal.

Progress heartbeats from reader, encoders, and writer feed an optional
stall-timeout watchdog. The watchdog is disabled by default because legitimate
large blocks may take a long time; CI enables it for deadlock/failure tests.

## Transactional output

Parallel encoding writes to unique temporary paths:

```text
<output>.tmp.<pid>.<uuid>
<output>.tmp.<pid>.<uuid>.metadata/
```

Only after all children exit successfully and the block counts match does the
parent replace the final GVC file and sidecar. Existing output is backed up
during commit and restored if either replacement fails.

Thus a failed encode must leave the previous final artifact intact, if one
existed.

## Fork and spawn

The process context can be selected with
`multiprocessing_start_method`.

With `fork`, inherited process state is naturally visible in workers. With
`spawn`, children import a clean interpreter and therefore do not inherit
runtime mutations such as dynamically registered codec functions.

For process-local state, callers may provide:

```python
Encoder(
    ...,
    multiprocessing_start_method="spawn",
    multiprocessing_initializer=my_initializer,
    multiprocessing_initializer_args=(...),
)
```

The initializer runs once in each encoder child before it receives work. For
`spawn`, the initializer and every value in its argument tuple must be
picklable; in practice the initializer should be a module-level function, not
a lambda, closure, or nested function.

The reader and writer do not run this initializer because codec execution is
isolated to encoder workers.

## Test requirements

The maintained suite covers:

- sequential vs two-worker byte-for-byte equivalence;
- deliberately out-of-order worker completion;
- reader and encoder failures;
- missing/duplicate block detection;
- structured error propagation;
- no child leakage after failure;
- preservation of pre-existing output on failure;
- transactional rollback if sidecar commit fails;
- watchdog termination of a stalled child;
- spawn-mode failure supervision;
- successful spawn-mode encode/decode with explicit child initialization;
- unsafe output path rejection.

Large production-scale queue pressure, real external JBIG subprocess behavior,
and non-Linux native multiprocessing remain release/offline verification items.


## Parent termination handling

On platforms that provide `SIGTERM`, the parent supervisor temporarily installs
a parent-only termination handler after child processes have started. A received
`SIGTERM` is converted into a structured `SignalTermination` supervisor
error. Normal supervisor unwinding then:

1. sets the shared stop event;
2. gives children a bounded graceful-exit window;
3. terminates/kills any remaining children;
4. removes temporary GVC and metadata artifacts;
5. restores the parent's previous signal handler;
6. closes IPC queues without waiting indefinitely for dead feeder threads.

This is intended to make scheduler cancellation and ordinary process
termination behave like any other supervised failure rather than publishing a
partial output.
