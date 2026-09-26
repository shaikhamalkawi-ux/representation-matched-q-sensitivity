from pathlib import Path
import subprocess, sys

ROOT=Path(__file__).resolve().parent
jobs=[
    (ROOT/"case01_pinar_boran_thesis","reproduce_pinar_thesis.py"),
    (ROOT/"case02_khan_cosine_topsis","reproduce_khan2021.py"),
    (ROOT/"case03_liu_wang_qrofwa","reproduce_liu_wang.py"),
    (ROOT/"case04_seikh_mandal_archimedean","reproduce_seikh_mandal.py"),
    (ROOT/"case05_du2021_einstein","reproduce_du2021.py"),
    (ROOT/"case06_du2019","reproduce_du2019.py"),
    (ROOT/"case07_alkan_kahraman_2021","reproduce_alkan2021.py"),
    (ROOT/"theory","proof_checks.py"),
    (ROOT,"release_verify.py"),
]
for cwd,script in jobs:
    print("\n==>", cwd.name, script)
    subprocess.run([sys.executable,script],cwd=cwd,check=True)
print("\nALL EXECUTABLE CHECKS: PASS")
