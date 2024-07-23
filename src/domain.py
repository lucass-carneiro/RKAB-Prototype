from dataclasses import dataclass
import numpy as np

import itertools


@dataclass
class Point:
    i: int
    j: int
    k: int
    x: float
    y: float
    z: float

    def idx(self):
        return (self.i, self.j, self.k)

    def pos(self):
        return (self.x, self.y, self.z)


class Domain:
    def __init__(self, start: float, end: float, n: int, ghosts: int):
        self.int_start = start
        self.int_end = end
        self.ghosts = ghosts

        self.n_int = n
        self.n_all = self.n_int + 2 * self.ghosts
        self.delta = (self.int_end - self.int_start) / (self.n_int - 1)

        self.all_start = self.int_start - self.ghosts * self.delta
        self.all_end = self.int_end + self.ghosts * self.delta

    def all_points(self):
        for k in range(self.n_all):
            for j in range(self.n_all):
                for i in range(self.n_all):
                    x = self.all_start + i * self.delta
                    y = self.all_start + j * self.delta
                    z = self.all_start + k * self.delta
                    yield Point(i, j, k, x, y, z)

    def interior_points(self):
        for k in range(self.ghosts, self.n_all - self.ghosts):
            for j in range(self.ghosts, self.n_all - self.ghosts):
                for i in range(self.ghosts, self.n_all - self.ghosts):
                    x = self.all_start + i * self.delta
                    y = self.all_start + j * self.delta
                    z = self.all_start + k * self.delta
                    yield Point(i, j, k, x, y, z)

    def ghost_points(self):
        ghost_idx = *range(self.ghosts), *range(self.n_int + 1, self.n)
        for k in ghost_idx:
            for j in ghost_idx:
                for i in ghost_idx:
                    x = self.all_start + i * self.delta
                    y = self.all_start + j * self.delta
                    z = self.all_start + k * self.delta
                    yield Point(i, j, k, x, y, z)
