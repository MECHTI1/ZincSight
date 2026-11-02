"""Worker budgets respect affinity, container quotas, and declared memory."""
import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('tamarind_resources', ROOT / 'integrations/tamarind/resources.py')
resources = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(resources)


class ResourceTests(unittest.TestCase):
    def fake_files(self, values):
        def read(path, *args, **kwargs):
            if str(path) not in values:
                raise OSError('missing controller')
            return values[str(path)]
        return patch.object(Path, 'read_text', read)

    def test_cpu_v2_quota_and_affinity_are_both_respected(self):
        with patch.object(resources.os, 'sched_getaffinity', return_value=set(range(8))), self.fake_files({
                '/sys/fs/cgroup/cpu.max': '250000 100000'}):
            self.assertEqual(resources.detected_cpu_limit(), 2)
        with patch.object(resources.os, 'sched_getaffinity', return_value={0, 1}), self.fake_files({
                '/sys/fs/cgroup/cpu.max': '800000 100000'}):
            self.assertEqual(resources.detected_cpu_limit(), 2)

    def test_cpu_v1_separate_and_combined_controllers_are_supported(self):
        for controller in ('cpu', 'cpu,cpuacct'):
            values = {f'/sys/fs/cgroup/{controller}/cpu.cfs_quota_us': '300000',
                      f'/sys/fs/cgroup/{controller}/cpu.cfs_period_us': '100000'}
            with self.subTest(controller=controller), patch.object(resources.os, 'sched_getaffinity', return_value=set(range(8))), self.fake_files(values):
                self.assertEqual(resources.detected_cpu_limit(), 3)

    def test_unlimited_or_malformed_cpu_quota_does_not_crash(self):
        for quota in ('max 100000', '0 100000', '100000 0', 'invalid', '-1 100000'):
            with self.subTest(quota=quota), patch.object(resources.os, 'sched_getaffinity', return_value=set(range(4))), self.fake_files({'/sys/fs/cgroup/cpu.max': quota}):
                self.assertEqual(resources.detected_cpu_limit(), 4)

    def test_cpu_count_fallback_has_at_least_one_worker(self):
        with patch.object(resources.os, 'sched_getaffinity', side_effect=OSError), patch.object(resources.os, 'cpu_count', return_value=None), self.fake_files({}):
            self.assertEqual(resources.detected_cpu_limit(), 1)

    def test_memory_v2_v1_and_unlimited_limits(self):
        with self.fake_files({'/sys/fs/cgroup/memory.max': str(12 * 2**30)}):
            self.assertEqual(resources.detected_memory_gib(), 12)
        with self.fake_files({'/sys/fs/cgroup/memory.max': 'max',
                              '/sys/fs/cgroup/memory/memory.limit_in_bytes': str(24 * 2**30)}):
            self.assertEqual(resources.detected_memory_gib(), 24)
        with self.fake_files({'/sys/fs/cgroup/memory.max': 'max',
                              '/sys/fs/cgroup/memory/memory.limit_in_bytes': str(2**63 - 1)}):
            self.assertIsNone(resources.detected_memory_gib())
        with self.fake_files({}):
            self.assertIsNone(resources.detected_memory_gib())

    def test_unlimited_container_uses_available_host_memory_before_total(self):
        with self.fake_files({'/sys/fs/cgroup/memory.max': 'max',
                              '/proc/meminfo': 'MemTotal:       50331648 kB\nMemAvailable:   12582912 kB\n'}):
            self.assertEqual(resources.detected_memory_gib(), 12)
            result = resources.worker_allocation(8, {}, cpu_limit=8)
            self.assertEqual(result['cpu_workers'], 4)
            self.assertEqual(result['effective_memory_gib'], 12)
        with self.fake_files({'/proc/meminfo': 'MemTotal: 25165824 kB\n'}):
            self.assertEqual(resources.detected_memory_gib(), 24)

    def test_worker_count_is_clamped_by_all_resource_constraints(self):
        environment = {'ZINCSIGHT_CPU_BUDGET': '8', 'ZINCSIGHT_MEMORY_GIB': '24'}
        cases = [(8, 8, 12, 4), (8, 2, 24, 2), (8, 8, 5.9, 1), (1, 8, 24, 1)]
        for requested, cpu, memory, expected in cases:
            with self.subTest(requested=requested, cpu=cpu, memory=memory):
                result = resources.worker_allocation(requested, environment, cpu_limit=cpu, memory_gib=memory)
                self.assertEqual(result['cpu_workers'], expected)
                self.assertEqual(result['memory_budget_gib_per_worker'], 3)
                self.assertEqual(result['effective_memory_gib'], memory)
                self.assertEqual(result['worker_limit_applied'], expected < requested)
        result = resources.worker_allocation(8, {'ZINCSIGHT_CPU_BUDGET': '2', 'ZINCSIGHT_MEMORY_GIB': '24'},
                                            cpu_limit=8, memory_gib=24)
        self.assertEqual(result['cpu_workers'], 2)

    def test_declared_defaults_apply_if_memory_detection_is_unavailable(self):
        with patch.object(resources, 'detected_cpu_limit', return_value=8), patch.object(resources, 'detected_memory_gib', return_value=None):
            result = resources.worker_allocation(8, {})
        self.assertEqual(result['cpu_workers'], 8)
        self.assertEqual(result['effective_memory_gib'], 24)
        self.assertEqual(result['declared_cpu_budget'], 8)

    def test_invalid_worker_requests_and_detected_resources_are_rejected(self):
        for request in (0, 9, True, 2.5, float('nan'), '4'):
            with self.subTest(request=request), self.assertRaises(ValueError):
                resources.worker_allocation(request, {}, cpu_limit=4, memory_gib=12)
        for cpu in (0, -1, 2.5, True):
            with self.subTest(cpu=cpu), self.assertRaises(ValueError):
                resources.worker_allocation(4, {}, cpu_limit=cpu, memory_gib=12)
        for memory in (0, 2.9, -1, float('nan'), float('inf'), '12', True):
            with self.subTest(memory=memory), self.assertRaises(ValueError):
                resources.worker_allocation(4, {}, cpu_limit=4, memory_gib=memory)

    def test_invalid_declared_cpu_or_memory_budget_is_rejected(self):
        for environment in ({'ZINCSIGHT_CPU_BUDGET': '0'}, {'ZINCSIGHT_CPU_BUDGET': '4.5'},
                            {'ZINCSIGHT_MEMORY_GIB': '2.9'}, {'ZINCSIGHT_MEMORY_GIB': 'nan'},
                            {'ZINCSIGHT_MEMORY_GIB': 'inf'}):
            with self.subTest(environment=environment), self.assertRaises(ValueError):
                resources.worker_allocation(4, environment, cpu_limit=4, memory_gib=12)


if __name__ == '__main__':
    unittest.main()
