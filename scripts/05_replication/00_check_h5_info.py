# ==============================================================================
# 00_check_h5_info.py
# List the contents of a GEO HDF5 file (replication cohort).
# Produces: data check
# Run     : python3 scripts/05_replication/00_check_h5_info.py   (from the repository root)
# ==============================================================================
import sys, os, re, h5py
f=sys.argv[1]
actual=os.path.getsize(f) if os.path.exists(f) else 0
stored=actual; ok=1
try:
    h5py.File(f,'r')
except Exception as e:
    ok=0
    m=re.search(r'stored_eof = (\d+)', str(e))
    if m: stored=int(m.group(1))
print(actual); print(stored); print(ok)
