"""Thin ctypes wrapper for the fastboss C ABI."""

from __future__ import annotations

import ctypes as ct
import os
import pathlib
import sys
from dataclasses import dataclass
from typing import Optional

import numpy as np

FASTBOSS_ABI_VERSION = 1


class FastBossError(RuntimeError):
    pass


class _Options(ct.Structure):
    _fields_ = [
        ("struct_size", ct.c_size_t),
        ("abi_version", ct.c_int),
        ("alpha", ct.c_double),
        ("max_sweeps", ct.c_int),
        ("use_gst", ct.c_int),
        ("random_start", ct.c_int),
        ("seed", ct.c_uint),
        ("min_gain", ct.c_double),
        ("reserved", ct.c_int * 8),
    ]


class _Statistics(ct.Structure):
    _fields_ = [
        ("struct_size", ct.c_size_t),
        ("abi_version", ct.c_int),
        ("runtime_s", ct.c_double),
        ("sweeps", ct.c_int),
        ("accepted_moves", ct.c_int),
        ("score", ct.c_double),
        ("local_score_calls", ct.c_double),
        ("gst_trace_calls", ct.c_double),
        ("gst_cached_nodes_allocated", ct.c_double),
        ("gst_cached_distinct_node_visits", ct.c_double),
        ("gst_cached_child_materializations", ct.c_double),
        ("gst_cached_branch_entries", ct.c_double),
        ("gst_cached_branch_considerations", ct.c_double),
        ("gst_cached_branch_prefix_hits", ct.c_double),
        ("gst_cached_grow_candidate_checks", ct.c_double),
        ("gst_nocache_trace_calls", ct.c_double),
        ("reserved", ct.c_double * 8),
    ]


class _Result(ct.Structure):
    pass


@dataclass
class FastBossResult:
    score: float
    order: np.ndarray
    initial_order: np.ndarray
    adjacency: np.ndarray
    statistics: dict


def _library_candidates() -> list[pathlib.Path]:
    env = os.environ.get("FASTBOSS_LIBRARY")
    if env:
        return [pathlib.Path(env)]
    suffix = {
        "darwin": ".dylib",
        "win32": ".dll",
    }.get(sys.platform, ".so")
    here = pathlib.Path(__file__).resolve()
    roots = [here.parents[i] for i in range(min(6, len(here.parents)))]
    names = [f"libfastboss{suffix}", f"fastboss{suffix}"]
    return [root / "build" / name for root in roots for name in names]


def _load_library() -> ct.CDLL:
    for path in _library_candidates():
        if path.exists():
            lib = ct.CDLL(str(path))
            break
    else:
        searched = "\n".join(str(p) for p in _library_candidates())
        raise FastBossError(
            "Could not find fastboss shared library. Set FASTBOSS_LIBRARY.\n"
            f"Searched:\n{searched}"
        )

    lib.fastboss_default_options.restype = _Options
    lib.fastboss_fit_column_major.argtypes = [
        ct.POINTER(ct.c_double),
        ct.c_int,
        ct.c_int,
        ct.POINTER(_Options),
        ct.POINTER(ct.c_int),
        ct.POINTER(ct.POINTER(_Result)),
    ]
    lib.fastboss_fit_column_major.restype = ct.c_int
    lib.fastboss_score_order_column_major.argtypes = [
        ct.POINTER(ct.c_double),
        ct.c_int,
        ct.c_int,
        ct.POINTER(_Options),
        ct.POINTER(ct.c_int),
        ct.POINTER(ct.POINTER(_Result)),
    ]
    lib.fastboss_score_order_column_major.restype = ct.c_int
    lib.fastboss_free_result.argtypes = [ct.POINTER(_Result)]
    lib.fastboss_last_error.restype = ct.c_char_p
    lib.fastboss_result_score.argtypes = [ct.POINTER(_Result)]
    lib.fastboss_result_score.restype = ct.c_double
    lib.fastboss_result_num_variables.argtypes = [ct.POINTER(_Result)]
    lib.fastboss_result_num_variables.restype = ct.c_int
    lib.fastboss_result_order.argtypes = [ct.POINTER(_Result)]
    lib.fastboss_result_order.restype = ct.POINTER(ct.c_int)
    lib.fastboss_result_initial_order.argtypes = [ct.POINTER(_Result)]
    lib.fastboss_result_initial_order.restype = ct.POINTER(ct.c_int)
    lib.fastboss_result_adjacency.argtypes = [ct.POINTER(_Result)]
    lib.fastboss_result_adjacency.restype = ct.POINTER(ct.c_int)
    lib.fastboss_result_statistics.argtypes = [ct.POINTER(_Result)]
    lib.fastboss_result_statistics.restype = ct.POINTER(_Statistics)
    opt = lib.fastboss_default_options()
    if opt.struct_size != ct.sizeof(_Options) or opt.abi_version != FASTBOSS_ABI_VERSION:
        raise FastBossError(
            "fastboss option ABI mismatch: "
            f"C struct_size={opt.struct_size}, Python struct_size={ct.sizeof(_Options)}, "
            f"abi_version={opt.abi_version}"
        )
    return lib


_LIB: Optional[ct.CDLL] = None


def _lib() -> ct.CDLL:
    global _LIB
    if _LIB is None:
        _LIB = _load_library()
    return _LIB


def _options(
    alpha: float,
    sweeps: int,
    use_gst: bool,
    random_start: bool,
    seed: int,
    min_gain: float,
) -> _Options:
    if not isinstance(seed, (int, np.integer)) or not (0 <= seed <= 0xFFFFFFFF):
        raise ValueError("seed must be an integer in [0, 2**32 - 1]")
    opt = _lib().fastboss_default_options()
    opt.alpha = float(alpha)
    opt.max_sweeps = int(sweeps)
    opt.use_gst = int(bool(use_gst))
    opt.random_start = int(bool(random_start))
    opt.seed = int(seed)
    opt.min_gain = float(min_gain)
    return opt


def _matrix(data) -> np.ndarray:
    raw = np.asarray(data)
    if raw.ndim != 2:
        raise ValueError("data must be a two-dimensional array")
    if np.iscomplexobj(raw):
        raise ValueError("data must be real-valued")
    arr = np.asarray(raw, dtype=np.float64, order="F")
    if not np.all(np.isfinite(arr)):
        raise ValueError("data must contain only finite values")
    return arr


def _order(order, p: int) -> Optional[np.ndarray]:
    if order is None:
        return None
    raw = np.asarray(order)
    if raw.shape != (p,):
        raise ValueError("order must be a zero-based permutation of range(p)")

    if not np.issubdtype(raw.dtype, np.number):
        raise ValueError("order must be a zero-based permutation of range(p)")
    if np.iscomplexobj(raw):
        raise ValueError("order must contain real integer indices")

    numeric = raw.astype(np.float64, copy=False)
    if not np.all(np.isfinite(numeric)) or not np.all(numeric == np.rint(numeric)):
        raise ValueError("order must be a zero-based permutation of range(p)")

    values = numeric.astype(np.int64)
    if not np.array_equal(np.sort(values), np.arange(p, dtype=np.int64)):
        raise ValueError("order must be a zero-based permutation of range(p)")

    arr = values.astype(np.int32)
    return np.asfortranarray(arr)


def _collect(result_ptr: ct.POINTER(_Result)) -> FastBossResult:
    lib = _lib()
    p = lib.fastboss_result_num_variables(result_ptr)
    score = lib.fastboss_result_score(result_ptr)
    order = np.ctypeslib.as_array(lib.fastboss_result_order(result_ptr), shape=(p,)).copy()
    initial_order = np.ctypeslib.as_array(
        lib.fastboss_result_initial_order(result_ptr), shape=(p,)
    ).copy()
    adjacency = np.ctypeslib.as_array(
        lib.fastboss_result_adjacency(result_ptr), shape=(p * p,)
    ).reshape((p, p)).copy()
    stats = lib.fastboss_result_statistics(result_ptr).contents
    statistics = {
        name: getattr(stats, name)
        for name, _ctype in _Statistics._fields_
        if name not in {"reserved", "struct_size", "abi_version"}
    }
    return FastBossResult(score, order, initial_order, adjacency, statistics)


def _raise_last_error() -> None:
    msg = _lib().fastboss_last_error()
    raise FastBossError(msg.decode("utf-8") if msg else "fastboss failed")


def fit(
    data,
    *,
    alpha=2.0,
    sweeps=0,
    order=None,
    random_start=False,
    seed=20260525,
    use_gst=True,
    min_gain=1e-6,
):
    arr = _matrix(data)
    n, p = arr.shape
    order_arr = _order(order, p)
    opt = _options(alpha, sweeps, use_gst, random_start, seed, min_gain)
    result = ct.POINTER(_Result)()
    status = _lib().fastboss_fit_column_major(
        arr.ctypes.data_as(ct.POINTER(ct.c_double)),
        n,
        p,
        ct.byref(opt),
        None if order_arr is None else order_arr.ctypes.data_as(ct.POINTER(ct.c_int)),
        ct.byref(result),
    )
    if status != 0:
        _raise_last_error()
    try:
        return _collect(result)
    finally:
        _lib().fastboss_free_result(result)


def score_order(data, order, *, alpha=2.0, use_gst=True):
    arr = _matrix(data)
    n, p = arr.shape
    order_arr = _order(order, p)
    if order_arr is None:
        raise ValueError("order is required")
    opt = _options(alpha, 0, use_gst, False, 1, 1e-6)
    result = ct.POINTER(_Result)()
    status = _lib().fastboss_score_order_column_major(
        arr.ctypes.data_as(ct.POINTER(ct.c_double)),
        n,
        p,
        ct.byref(opt),
        order_arr.ctypes.data_as(ct.POINTER(ct.c_int)),
        ct.byref(result),
    )
    if status != 0:
        _raise_last_error()
    try:
        return _collect(result)
    finally:
        _lib().fastboss_free_result(result)
