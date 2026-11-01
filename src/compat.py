"""Compatibility with current versions of ZincSight's dependencies.

ZincSight was written for ``pgserver`` and for Biopython releases that still
ship ``Bio.PDB.QCPSuperimposer``. Current environments (Python 3.12+, Google
Colab) install ``pixeltable-pgserver`` instead, and Biopython 1.85+ no longer
has that module. Importing this module provides both names, exactly as the
Colab notebooks do, so a plain ``pip install -r requirements.txt`` works.
When the original packages are installed they are used unchanged.
"""
import sys
import types


def _provide_pgserver():
    try:
        import pgserver  # noqa: F401  (original package, or already provided)
        return
    except ImportError:
        pass
    try:
        import pixeltable_pgserver as pgserver
    except ImportError:
        return  # neither is installed; the caller's own import reports it
    original_get_server = pgserver.get_server

    def get_server(pgdata, cleanup_mode='stop', *args, **kwargs):
        return original_get_server(pgdata, *args, cleanup_mode=cleanup_mode, **kwargs)

    pgserver.get_server = get_server
    sys.modules['pgserver'] = pgserver


def _provide_qcp_superimposer():
    try:
        from Bio.PDB.QCPSuperimposer import QCPSuperimposer  # noqa: F401
        return
    except ImportError:
        pass
    import numpy as np
    from scipy.spatial.transform import Rotation

    class QCPSuperimposer:
        """Minimal-RMSD superposition with the API ZincSight uses (same as the Colab notebooks)."""

        def __init__(self):
            self.coords1 = None
            self.coords2 = None
            self.rms = 0.0

        def set(self, coords1, coords2):
            self.coords1 = coords1
            self.coords2 = coords2

        def run(self):
            centered1 = self.coords1 - self.coords1.mean(axis=0)
            centered2 = self.coords2 - self.coords2.mean(axis=0)
            rotation, _ = Rotation.align_vectors(centered1, centered2)
            diff = centered1 - rotation.apply(centered2)
            self.rms = float(np.sqrt(np.mean(np.sum(diff ** 2, axis=-1))))

        def get_rms(self):
            return self.rms

    module = types.ModuleType('Bio.PDB.QCPSuperimposer')
    module.QCPSuperimposer = QCPSuperimposer
    sys.modules['Bio.PDB.QCPSuperimposer'] = module


_provide_pgserver()
_provide_qcp_superimposer()
