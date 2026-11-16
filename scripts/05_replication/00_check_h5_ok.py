# ==============================================================================
# 00_check_h5_ok.py
# Check that a GEO HDF5 file opens (replication cohort).
# Produces: data check
# Run     : python3 scripts/05_replication/00_check_h5_ok.py   (from the repository root)
# ==============================================================================
import sys, h5py
try:
    h5py.File(sys.argv[1],'r'); sys.exit(0)
except Exception:
    sys.exit(1)
