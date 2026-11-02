"""Conservative worker allocation within declared and detected job resources."""
import math
import os
from pathlib import Path


def detected_cpu_limit():
    limits = []
    try:
        limits.append(len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        limits.append(os.cpu_count() or 1)
    try:
        quota, period = Path('/sys/fs/cgroup/cpu.max').read_text().split()
        if quota != 'max' and int(quota) > 0 and int(period) > 0:
            limits.append(max(1, int(quota) // int(period)))
    except (OSError, ValueError):
        pass
    for controller in ('cpu', 'cpu,cpuacct'):
        try:
            base = Path('/sys/fs/cgroup') / controller
            quota = int((base / 'cpu.cfs_quota_us').read_text())
            period = int((base / 'cpu.cfs_period_us').read_text())
            if quota > 0 and period > 0:
                limits.append(max(1, quota // period))
        except (OSError, ValueError):
            pass
    return max(1, min(limits))


def detected_memory_gib():
    for name in ('/sys/fs/cgroup/memory.max', '/sys/fs/cgroup/memory/memory.limit_in_bytes'):
        try:
            value = int(Path(name).read_text().strip())
            if value < 2**60:
                return value / 2**30
        except (OSError, ValueError):
            pass
    return None


def worker_allocation(requested, environment, *, cpu_limit=None, memory_gib=None):
    if isinstance(requested, bool) or not isinstance(requested, int) or not 1 <= requested <= 8:
        raise ValueError('Requested CPU workers must be an integer between 1 and 8.')
    declared_cpu = int(environment.get('ZINCSIGHT_CPU_BUDGET', '4'))
    declared_memory = float(environment.get('ZINCSIGHT_MEMORY_GIB', '12'))
    if declared_cpu < 1 or not math.isfinite(declared_memory) or declared_memory < 3:
        raise ValueError('Allocate at least one CPU and 3 GiB of memory for ZincSight.')
    cpu_limit = detected_cpu_limit() if cpu_limit is None else cpu_limit
    memory_gib = detected_memory_gib() if memory_gib is None else memory_gib
    if isinstance(cpu_limit, bool) or not isinstance(cpu_limit, int) or cpu_limit < 1:
        raise ValueError('The detected CPU limit must be a positive integer.')
    if memory_gib is not None and (isinstance(memory_gib, bool)
            or not isinstance(memory_gib, (int, float))
            or not math.isfinite(memory_gib) or memory_gib < 0):
        raise ValueError('The detected memory limit must be finite and nonnegative.')
    effective_memory = min(declared_memory, memory_gib) if memory_gib is not None else declared_memory
    if effective_memory < 3:
        raise ValueError('The detected job memory is below the 3 GiB worker budget.')
    workers = max(1, min(requested, declared_cpu, cpu_limit, int(effective_memory // 3), 8))
    return {'requested_cpu_workers': requested, 'cpu_workers': workers,
            'declared_cpu_budget': declared_cpu, 'declared_memory_gib': declared_memory,
            'detected_cpu_limit': cpu_limit, 'effective_memory_gib': effective_memory,
            'memory_budget_gib_per_worker': 3,
            'worker_limit_applied': workers < requested}
