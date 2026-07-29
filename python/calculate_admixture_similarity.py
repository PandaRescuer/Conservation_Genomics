#!/usr/bin/env python3

import argparse
import itertools
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate average pairwise similarity among ADMIXTURE runs."
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing ADMIXTURE Q files"
    )
    parser.add_argument(
        "--prefix",
        default="admixture",
        help="Q-file prefix"
    )
    parser.add_argument(
        "--k-min",
        type=int,
        default=2,
        help="Minimum K value"
    )
    parser.add_argument(
        "--k-max",
        type=int,
        default=5,
        help="Maximum K value"
    )
    parser.add_argument(
        "--runs",
        type=int,
        default=3,
        help="Number of runs for each K"
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Output file prefix"
    )
    return parser.parse_args()


def clumpp_similarity(q1, q2):
    """Calculate the CLUMPP similarity coefficient."""
    k = q1.shape[1]
    uniform = np.full(q1.shape, 1.0 / k)

    numerator = np.linalg.norm(q1 - q2)
    denominator = np.sqrt(
        np.linalg.norm(q1 - uniform)
        * np.linalg.norm(q2 - uniform)
    )

    if denominator == 0:
        return 1.0 if np.allclose(q1, q2) else np.nan

    return 1.0 - numerator / denominator


def aligned_similarity(q1, q2):
    """Find the maximum similarity across all cluster permutations."""
    if q1.shape != q2.shape:
        raise ValueError(
            f"Q-matrix dimensions differ: {q1.shape} and {q2.shape}"
        )

    k = q1.shape[1]
    best_similarity = -np.inf
    best_permutation = None

    for permutation in itertools.permutations(range(k)):
        permuted_q2 = q2[:, permutation]
        similarity = clumpp_similarity(q1, permuted_q2)

        if similarity > best_similarity:
            best_similarity = similarity
            best_permutation = permutation

    return best_similarity, best_permutation


def main():
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_prefix = Path(args.output_prefix)

    summary_rows = []
    pairwise_rows = []

    for k in range(args.k_min, args.k_max + 1):
        q_matrices = []

        for run in range(1, args.runs + 1):
            q_file = input_dir / (
                f"{args.prefix}.K{k}.seed{run}.Q"
            )

            if not q_file.exists():
                raise FileNotFoundError(f"Missing Q file: {q_file}")

            q_matrix = np.loadtxt(q_file)

            if q_matrix.ndim == 1:
                q_matrix = q_matrix.reshape(-1, k)

            if q_matrix.shape[1] != k:
                raise ValueError(
                    f"{q_file} has {q_matrix.shape[1]} columns, expected {k}"
                )

            q_matrices.append(q_matrix)

        similarities = []

        for run1, run2 in itertools.combinations(range(args.runs), 2):
            similarity, permutation = aligned_similarity(
                q_matrices[run1],
                q_matrices[run2]
            )
            similarities.append(similarity)

            pairwise_rows.append(
                [
                    k,
                    run1 + 1,
                    run2 + 1,
                    similarity,
                    ",".join(str(value + 1) for value in permutation),
                ]
            )

        summary_rows.append(
            [
                k,
                len(similarities),
                np.mean(similarities),
                np.min(similarities),
                np.max(similarities),
                np.std(similarities, ddof=1)
                if len(similarities) > 1 else 0.0,
            ]
        )

    summary_file = Path(f"{output_prefix}.summary.tsv")
    pairwise_file = Path(f"{output_prefix}.pairwise.tsv")
    plot_file = Path(f"{output_prefix}.pdf")

    with open(summary_file, "w") as output:
        output.write(
            "K\tNumber_of_comparisons\t"
            "Average_pairwise_similarity\t"
            "Minimum_similarity\tMaximum_similarity\tSD\n"
        )

        for row in summary_rows:
            output.write(
                f"{row[0]}\t{row[1]}\t{row[2]:.6f}\t"
                f"{row[3]:.6f}\t{row[4]:.6f}\t{row[5]:.6f}\n"
            )

    with open(pairwise_file, "w") as output:
        output.write(
            "K\tRun_1\tRun_2\tSimilarity\t"
            "Best_column_permutation_for_run_2\n"
        )

        for row in pairwise_rows:
            output.write(
                f"{row[0]}\t{row[1]}\t{row[2]}\t"
                f"{row[3]:.6f}\t{row[4]}\n"
            )

    k_values = [row[0] for row in summary_rows]
    means = [row[2] for row in summary_rows]
    standard_deviations = [row[5] for row in summary_rows]

    plt.figure(figsize=(5, 4))
    plt.errorbar(
        k_values,
        means,
        yerr=standard_deviations,
        fmt="o-",
        color="black",
        capsize=4
    )
    plt.xticks(k_values)
    plt.ylim(0, 1.02)
    plt.xlabel("K")
    plt.ylabel("Average pairwise similarity")
    plt.tight_layout()
    plt.savefig(plot_file)
    plt.close()


if __name__ == "__main__":
    main()