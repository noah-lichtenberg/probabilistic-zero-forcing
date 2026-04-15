import math

from readrda import BEAinputoutputdata
from rzf_core import (
    graph_gen,
    tm_generation_directed_sparse,
    propagation_time_solver_sparse,
    nums_with_bitcount,
)
from rzf_plotting import plot_graph_with_zfs_heatmap


CATEGORY_KEY = {
    0: "Agriculture, forestry, fishing, and hunting",
    1: "Mining",
    2: "Utilities",
    3: "Construction",
    4: "Manufacturing",
    5: "Wholesale trade",
    6: "Retail trade",
    7: "Transportation and warehousing",
    8: "Information",
    9: "Finance, insurance, real estate, rental, and leasing",
    10: "Professional and business services",
    11: "Educational services, health care, and social assistance",
    12: "Arts, entertainment, recreation, accommodation, and food services",
    13: "Other services, except government",
    14: "Government",
}


def main():
    G = BEAinputoutputdata
    n = G.number_of_nodes()

    graphs = graph_gen(G)
    tm, _ = tm_generation_directed_sparse(graphs)

    print("solving")
    times, solve_time = propagation_time_solver_sparse(tm)

    singleton_indices = nums_with_bitcount(2 ** n, 1)
    zfs_size_one_times = [(int(math.log2(i)), times[i]) for i in singleton_indices]

    best_idx = min(range(len(zfs_size_one_times)), key=lambda idx: zfs_size_one_times[idx][1])
    print(zfs_size_one_times)
    print("Best singleton start:", zfs_size_one_times[best_idx])
    print("Solve time:", solve_time)

    plot_graph_with_zfs_heatmap(
        G,
        zfs_size_one_times,
        id_to_sector=CATEGORY_KEY,
        layout="circular",
        shift_labels_by_1=True,
        wrap_width=30,
        key_fontsize=11,
    )


if __name__ == "__main__":
    main()