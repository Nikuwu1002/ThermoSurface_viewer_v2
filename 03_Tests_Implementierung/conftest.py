import sys
from pathlib import Path

PROGRAM_DIR = Path(__file__).resolve().parents[1] / "01_Programm"
if str(PROGRAM_DIR) not in sys.path:
    sys.path.insert(0, str(PROGRAM_DIR))
