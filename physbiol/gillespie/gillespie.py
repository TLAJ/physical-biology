from dataclasses import dataclass
import numpy as np


@dataclass
class GillespieResult:
    next_time: float
    next_reaction: int


def get_next_time_reaction(v):
    v_total = np.sum(v)
    next_time = np.log(1 / np.random.rand()) / v_total

    r2 = np.random.rand()
    for ii in range(len(v) - 1):
        v_range = np.sum(v[: ii + 1]) / v_total
        if r2 < v_range:
            return GillespieResult(next_time, ii)
    return GillespieResult(next_time, len(v) - 1)
