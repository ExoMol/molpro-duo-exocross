"""Explicit ExoMol states-column mapping for one singlet Sigma+ state."""
import bz2
import math
from pathlib import Path

from .model import FitError, Observation, clean, number, observations


def read_states(path, *, j_column=4, v_column=5, weight_scale=100.0, state="X"):
    """Read E/cm-1 without shifting; column indices are one based.

    The first three columns are ID, energy, total statistical weight. The
    statistical weight is validated but does not enter the fitting weight.
    Explicitly assumes e parity and Lambda=Sigma=Omega=0 for every record.
    """
    if (any(isinstance(c, bool) or not isinstance(c, int) or c < 4 for c in (j_column, v_column))
            or j_column == v_column):
        raise FitError("States J/v columns must be distinct integers >= 4 (one based)")
    if not isinstance(weight_scale, (int, float)) or isinstance(weight_scale, bool) or not math.isfinite(weight_scale) or weight_scale <= 0:
        raise FitError("States weight scale must be finite and positive")
    if not state or len(state) > 3 or not state.isalnum():
        raise FitError("States electronic-state label must be 1-3 alphanumeric characters")
    path = Path(path)
    opener = bz2.open if path.suffix.lower() == ".bz2" else open
    result, ids = [], set()
    with opener(path, "rt", encoding="utf-8-sig") as stream:
        for lineno, line in enumerate(stream, 1):
            fields = clean(line).split()
            if not fields:
                continue
            try:
                if len(fields) < max(j_column, v_column):
                    raise FitError("Missing configured J/v column")
                identifier, energy, degeneracy = map(number, fields[:3])
                j, v = number(fields[j_column-1]), number(fields[v_column-1])
                if identifier < 1 or int(identifier) != identifier or identifier in ids:
                    raise FitError("States ID must be a unique positive integer")
                if degeneracy <= 0 or int(degeneracy) != degeneracy:
                    raise FitError("Statistical weight must be a positive integer")
                if any(x < 0 or int(x) != x for x in (j, v)) or energy < 0:
                    raise FitError("Expected nonnegative energy and integer J/v")
                ids.add(identifier)
                result.append(Observation(int(j), "e", int(v)+1, energy, state, int(v),
                                          weight_scale/math.sqrt(v+1)/math.sqrt(j+1)))
            except FitError as exc:
                raise FitError(f"{path.name}:{lineno}: {exc}") from exc
    # Reuse the fitting-table checks, including duplicate (J,v) detection.
    observations([o.line() for o in result])
    return sorted(result, key=lambda o: o.key)
