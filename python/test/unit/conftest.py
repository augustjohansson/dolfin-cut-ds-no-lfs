import warnings
import pytest
import gc
from dolfin import MPI

try:
    from ffc.quadrature.deprecation import QuadratureRepresentationDeprecationWarning
    warnings.simplefilter("ignore", QuadratureRepresentationDeprecationWarning)
except ImportError:
    pass

def pytest_runtest_teardown(item):
    """Collect garbage after every test to force calling
    destructors which might be collective"""

    # Do the normal teardown
    item.teardown()

    # Collect the garbage (call destructors collectively)
    del item
    # NOTE: How are we sure that 'item' does not hold references
    #       to temporaries and someone else does not hold a reference
    #       to 'item'?! Well, it seems that it works...
    gc.collect()
    MPI.barrier(MPI.comm_world)
