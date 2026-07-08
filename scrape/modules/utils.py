import logging
import sys
from pathlib import Path


class _LogStdout:
    """Route print() from modules to the log file instead of the terminal."""
    def __init__(self, logger):
        self._log = logger

    def write(self, msg):
        msg = msg.rstrip("\n")
        if msg:
            self._log.info(msg)

    def flush(self):
        pass

    def isatty(self):
        return False


def setup_logging(log_file, verbose=False):
    """Configure logging for a PepXPro script.

    Always writes INFO+ to log_file. Captures Python warnings (Biopython,
    PDBFixer, etc.) into the log. Redirects sys.stdout so that print() calls
    from modules also land in the log instead of the terminal.

    Returns (logger, original_stdout). Use original_stdout.write() for the
    two concise verbose lines — do not use print() for those.
    """
    Path(log_file).parent.mkdir(parents=True, exist_ok=True)
    root = logging.getLogger()
    root.handlers.clear()
    fh = logging.FileHandler(log_file, mode="w", encoding="ascii")
    fh.setFormatter(logging.Formatter("%(message)s"))
    root.addHandler(fh)
    root.setLevel(logging.INFO)
    logging.captureWarnings(True)

    logger = logging.getLogger("pepxpro")
    original_stdout = sys.stdout
    sys.stdout = _LogStdout(logger)
    return logger, original_stdout


def load_tsv(path):
    rows = []
    with open(path, encoding="ascii") as fh:
        header = fh.readline().strip().split("\t")
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            rows.append(dict(zip(header, parts)))
    return rows


def save_tsv(rows, path, columns):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="ascii") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c, "")) for c in columns) + "\n")
