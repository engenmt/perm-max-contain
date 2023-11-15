#!/usr/bin/env python

import permpy as pp

from dataclasses import dataclass
from math import comb, factorial

import os
import sys
import time


def catalan(n):
    return comb(2 * n + 1, n) // (2 * n + 1)


@dataclass
class MaxContainSearch:
    n: int
    basis: list[pp.Permutation]
    num_blocks: int = (
        80  # This is the approximate width of the progress bar in the terminal.
    )

    def __post_init__(self):
        self.results: list[int] = [0 for _ in range(self.n + 1)]  # These are the maxes
        self.min_m = self.n // 3
        self.containment_totals: list[int] = [
            0 for _ in range(self.n + 1)
        ]  # These are the totals for computing the means
        self.containment_totals: list[int] = [
            0 for _ in range(self.n + 1)
        ]  # These are the totals for computing the means
        self.total_perms = factorial(self.n)
        self.num_perms_per_block = self.total_perms // self.num_blocks

    def permutations(self):
        for perm_idx, perm in enumerate(pp.Permutation.gen_all(self.n), start=1):
            # We 1-index so that perm_idx is human-readable for later.
            yield perm
            if perm_idx % self.num_perms_per_block == 0:
                num_blocks_filled = perm_idx // self.num_perms_per_block
                filled = "X" * num_blocks_filled
                # We initially print a carriage return '\r' to "print in place".
                print(
                    f"\r[{filled:<{self.num_blocks}}] {num_blocks_filled:>2}/{self.num_blocks}",
                    end="",
                )

    def compute_with_means(self):
        print(f"n = {self.n}; m = {self.min_m}, ..., {self.n};")
        # Lengths range from self.min_m to self.n

        for perm in self.permutations():
            downset = perm.downset(min_length=self.min_m)
            for m, patterns in enumerate(downset, start=self.min_m):
                # start=self.min_m sets the inital value for the index (m)
                # We still must filter with the basis
                num_proper_patterns = sum(p.avoids(B=self.basis) for p in patterns)
                self.containment_totals[m] += num_proper_patterns
                if num_proper_patterns > self.results[m]:
                    self.results[m] = num_proper_patterns

        self.means = [
            containment_total / self.total_perms
            for containment_total in self.containment_totals
        ]
        return self.results


def main(perm, N):
    basis = [pp.Permutation(perm)]
    time_start = time.time()

    S = MaxContainSearch(N, basis)
    maxes = S.compute_with_means()
    means = S.means

    time_end = time.time()
    seconds_elapsed = time_end - time_start
    print(f"{seconds_elapsed:8.1f}s elapsed.")

    log_to_file(perm, N, maxes, means, seconds_elapsed)


def log_to_file(perm, N, maxes, means, seconds_elapsed):
    dir = "data"
    if not os.path.exists(dir):
        os.makedirs(dir)

    log_file = os.path.join(dir, f"{perm}_{N}.txt")
    print(f"Logging to {log_file}")

    with open(log_file, "w") as f:
        # Herein, we print _to the file_!
        print(
            f"Seeking the maximum number of {perm}-avoiding patterns contained among all {N}-perms.",
            file=f,
        )

        print(f"{'m':>2} | {'MaxVal':>6} | {'MeanVal':>7}", file=f)
        for m, (max_val, mean_val) in enumerate(zip(maxes, means)):
            print(f"{m:>2} | {max_val:>6} | {mean_val:>7.2f}", file=f)

        print(f"{seconds_elapsed:8.1f}s elapsed.", file=f)


if __name__ == "__main__":
    try:
        perm = sys.argv[1]
        N = int(sys.argv[2])
    except IndexError:
        raise Exception(
            "Example usage for searching for 231-patterns in 8-perms: ./search.py 231 8"
        )

    main(perm, N)
