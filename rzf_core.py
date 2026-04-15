import itertools
import random
import time
from typing import Iterable, Optional

import networkx as nx
import numpy as np
import scipy.linalg as la
import sympy as sp
from scipy.sparse import lil_matrix, csr_matrix, hstack, vstack


# -------------------------
# Binary / state utilities
# -------------------------

def d2b(i: int, n: int) -> str:
    return bin(i)[2:].zfill(n)


def b2d(binary_string: str) -> int:
    return int(binary_string, 2)


def extract(lst):
    return [item[0] for item in lst]


def same_ones(num: int, n: int) -> list[int]:
    ones = bin(num).count("1")
    return [i for i in range(2**n) if bin(i).count("1") == ones]


def same_ones_larger(num: int, n: int) -> list[int]:
    ones = bin(num).count("1")
    return [i for i in range(num, 2**n) if bin(i).count("1") == ones]


def nums_with_bitcount(max_num: int, ones: int) -> list[int]:
    return [i for i in range(max_num) if bin(i).count("1") == ones]


def transition_possible(i: int, j: int, n: int) -> bool:
    bi = d2b(i, n)
    bj = d2b(j, n)
    for k in range(n):
        if bi[k] == "1" and bj[k] == "0":
            return False
    return True


def differing_digits(i: int, j: int, n: int) -> list[int]:
    num1_rev = d2b(i, n)[::-1]
    num2_rev = d2b(j, n)[::-1]
    return [idx for idx in range(n) if num1_rev[idx] != num2_rev[idx]]


def same_zeros(i: int, j: int, n: int) -> list[int]:
    num1_rev = d2b(i, n)[::-1]
    num2_rev = d2b(j, n)[::-1]
    return [idx for idx in range(n) if num1_rev[idx] == num2_rev[idx] == "0"]


# -------------------------
# Color / neighbor helpers
# -------------------------

def _node_order(graph: nx.Graph) -> list:
    return list(graph.nodes())


def blue_neighbors(graph: nx.Graph, node) -> list:
    return [nbr for nbr in graph.neighbors(node) if graph.nodes[nbr].get("color") == "blue"]


def blue_suppliers(graph: nx.DiGraph, node) -> list:
    return [pred for pred in graph.predecessors(node) if graph.nodes[pred].get("color") == "blue"]


def num_blue_neighbors(graph: nx.Graph) -> dict:
    return {node: len(blue_neighbors(graph, node)) for node in graph.nodes()}


def num_blue_suppliers(graph: nx.DiGraph) -> dict:
    return {node: len(blue_suppliers(graph, node)) for node in graph.nodes()}


# -------------------------
# Forcing probabilities
# -------------------------

def forced_prob(graph: nx.Graph, node, blue_neighbor_counts: dict) -> float:
    temp = 1.0
    for u in blue_neighbors(graph, node):
        deg = graph.degree[u]
        if deg == 0:
            continue
        u_forces_node_prob = (blue_neighbor_counts[u] + 1) / deg
        u_doesnt_force_node_prob = 1 - u_forces_node_prob
        temp *= u_doesnt_force_node_prob
    return 1 - temp


def forced_prob_directed(graph: nx.DiGraph, node, weight_attr: str = "weight") -> float:
    suppliers = list(graph.predecessors(node))
    if not suppliers:
        return 0.0

    weighted = all(weight_attr in graph.edges[e] for e in graph.in_edges(node))
    if not weighted:
        return len(blue_suppliers(graph, node)) / len(suppliers)

    total_weight = 0.0
    blue_weight = 0.0
    for u in suppliers:
        w = graph[u][node].get(weight_attr, 1.0)
        total_weight += w
        if graph.nodes[u].get("color") == "blue":
            blue_weight += w

    if total_weight == 0:
        return 0.0
    return blue_weight / total_weight


# -------------------------
# State graph generation
# -------------------------

def graph_gen(target_graph: nx.Graph) -> list[nx.Graph]:
    nodes = list(target_graph.nodes())
    n = len(nodes)
    node_to_pos = {node: idx for idx, node in enumerate(nodes)}

    graphs = [None] * (2**n)
    for state in range(2**n):
        binary = d2b(state, n)
        g = target_graph.copy()
        for node in nodes:
            pos = node_to_pos[node]
            g.nodes[node]["color"] = "blue" if binary[n - pos - 1] == "1" else "white"
        graphs[state] = g
    return graphs


def graph_gen_sub(target_graph: nx.Graph, states: Iterable[int]) -> list[nx.Graph]:
    nodes = list(target_graph.nodes())
    n = len(nodes)
    node_to_pos = {node: idx for idx, node in enumerate(nodes)}

    out = []
    for state in states:
        binary = d2b(state, n)
        g = target_graph.copy()
        for node in nodes:
            pos = node_to_pos[node]
            g.nodes[node]["color"] = "blue" if binary[n - pos - 1] == "1" else "white"
        out.append(g)
    return out


# -------------------------
# Reachable states
# -------------------------

def larger_bin_old(num: int, n: int) -> list[int]:
    nums = []
    binary = d2b(num, n)
    zero_indices = [n - i - 1 for i in range(n) if binary[i] == "0"]
    for r in range(len(zero_indices) + 1):
        for subset in itertools.combinations(zero_indices, r):
            temp = num
            for j in subset:
                temp += 2**j
            nums.append(temp)
    return sorted(nums)


def larger_bin(num: int, graphs: list[nx.Graph]) -> list[int]:
    g = graphs[num]
    nodes = list(g.nodes())
    n = len(nodes)

    binary = d2b(num, n)
    blue_positions = [n - i - 1 for i in range(n) if binary[i] == "1"]
    white_positions = [n - i - 1 for i in range(n) if binary[i] == "0"]

    possible_positions = []
    for pos in white_positions:
        node = nodes[pos]
        blue_nbrs = blue_neighbors(g, node)
        if blue_nbrs:
            possible_positions.append(pos)

    nums = []
    for r in range(len(possible_positions) + 1):
        for subset in itertools.combinations(possible_positions, r):
            temp = num
            for j in subset:
                temp += 2**j
            nums.append(temp)
    return sorted(nums)


def larger_bin_directed_new_rule(num: int, graphs: list[nx.DiGraph]) -> list[int]:
    g = graphs[num]
    nodes = list(g.nodes())
    n = len(nodes)

    binary = d2b(num, n)
    blue_positions = [n - i - 1 for i in range(n) if binary[i] == "1"]
    blue_nodes = {nodes[pos] for pos in blue_positions}

    possible_positions = []
    seen = set()
    for pos in blue_positions:
        u = nodes[pos]
        for v in g.successors(u):
            if v not in blue_nodes and v not in seen:
                seen.add(v)
                possible_positions.append(nodes.index(v))

    nums = []
    for r in range(len(possible_positions) + 1):
        for subset in itertools.combinations(possible_positions, r):
            temp = num
            for j in subset:
                temp += 2**j
            nums.append(temp)
    return sorted(nums)


def path_binary_states(n: int) -> list[int]:
    results = []
    center = n // 2 if n % 2 == 1 else n // 2 - 1
    for i in range(n):
        for j in range(i, n):
            if i <= center <= j:
                s = ["0"] * n
                for k in range(i, j + 1):
                    s[k] = "1"
                results.append("".join(s))
    results.sort()
    return [int(x, 2) for x in results]


# -------------------------
# Transition matrices
# -------------------------

def tm_generation_old(graphs: list[nx.Graph]):
    g0 = graphs[0]
    n = g0.number_of_nodes()
    tm = np.zeros((2**n, 2**n))
    start = time.process_time()

    for i in range(2**n):
        state_graph = graphs[i]
        blue_counts = num_blue_neighbors(state_graph)
        for j in range(2**n):
            if not transition_possible(i, j, n):
                continue
            prob = 1.0
            for k in differing_digits(i, j, n):
                node = list(state_graph.nodes())[k]
                prob *= forced_prob(state_graph, node, blue_counts)
            for m in same_zeros(i, j, n):
                node = list(state_graph.nodes())[m]
                prob *= 1 - forced_prob(state_graph, node, blue_counts)
            tm[i, j] = prob

    tm[2**n - 1, 2**n - 1] = 1.0
    return tm, time.process_time() - start


def tm_generation(graphs: list[nx.Graph]):
    if graphs[0].is_directed():
        raise ValueError("Graphs must be undirected.")
    g0 = graphs[0]
    n = g0.number_of_nodes()
    tm = np.zeros((2**n, 2**n))
    start = time.process_time()

    for i in range(2**n):
        state_graph = graphs[i]
        blue_counts = num_blue_neighbors(state_graph)
        nodes = list(state_graph.nodes())
        for j in larger_bin(i, graphs):
            prob = 1.0
            for k in differing_digits(i, j, n):
                prob *= forced_prob(state_graph, nodes[k], blue_counts)
            for m in same_zeros(i, j, n):
                prob *= 1 - forced_prob(state_graph, nodes[m], blue_counts)
            tm[i, j] = prob

    tm[2**n - 1, 2**n - 1] = 1.0
    return tm, time.process_time() - start


def tm_generation_directed(graphs: list[nx.DiGraph], weight_attr: str = "weight"):
    if not graphs[0].is_directed():
        raise ValueError("Graphs must be directed.")
    g0 = graphs[0]
    n = g0.number_of_nodes()
    tm = np.zeros((2**n, 2**n))
    start = time.process_time()

    for i in range(2**n):
        state_graph = graphs[i]
        nodes = list(state_graph.nodes())
        for j in larger_bin_directed_new_rule(i, graphs):
            prob = 1.0
            for k in differing_digits(i, j, n):
                prob *= forced_prob_directed(state_graph, nodes[k], weight_attr)
            for m in same_zeros(i, j, n):
                prob *= 1 - forced_prob_directed(state_graph, nodes[m], weight_attr)
            tm[i, j] = prob

    tm[2**n - 1, 2**n - 1] = 1.0
    return tm, time.process_time() - start


def tm_generation_directed_sparse(graphs: list[nx.DiGraph], weight_attr: str = "weight"):
    if not graphs[0].is_directed():
        raise ValueError("Graphs must be directed.")
    g0 = graphs[0]
    n = g0.number_of_nodes()
    tm = lil_matrix((2**n, 2**n))
    start = time.process_time()

    for i in range(2**n):
        state_graph = graphs[i]
        nodes = list(state_graph.nodes())
        for j in larger_bin_directed_new_rule(i, graphs):
            prob = 1.0
            for k in differing_digits(i, j, n):
                prob *= forced_prob_directed(state_graph, nodes[k], weight_attr)
            for m in same_zeros(i, j, n):
                prob *= 1 - forced_prob_directed(state_graph, nodes[m], weight_attr)
            tm[i, j] = prob

    tm[2**n - 1, 2**n - 1] = 1.0
    return tm.tocsr(), time.process_time() - start


def tm_generation_sub(graphs: list[nx.Graph], starting_set: int):
    g0 = graphs[0]
    n = g0.number_of_nodes()
    tm = np.zeros((2**n, 2**n))
    start = time.process_time()
    states = larger_bin(starting_set, graphs)

    for i in states:
        state_graph = graphs[i]
        blue_counts = num_blue_neighbors(state_graph)
        nodes = list(state_graph.nodes())
        for j in larger_bin(i, graphs):
            prob = 1.0
            for k in differing_digits(i, j, n):
                prob *= forced_prob(state_graph, nodes[k], blue_counts)
            for m in same_zeros(i, j, n):
                prob *= 1 - forced_prob(state_graph, nodes[m], blue_counts)
            tm[i, j] = prob

    tm[2**n - 1, 2**n - 1] = 1.0
    subtm = tm[np.ix_(states, states)]
    return subtm, time.process_time() - start


def tm_generation_directed_sub(graphs: list[nx.DiGraph], states: list[int], weight_attr: str = "weight"):
    if not graphs[0].is_directed():
        raise ValueError("Graphs must be directed.")
    n = graphs[0].number_of_nodes()
    num_states = len(states)
    tm = np.zeros((num_states, num_states))
    start = time.process_time()

    for i in range(num_states):
        state_graph = graphs[i]
        nodes = list(state_graph.nodes())
        for j in range(i, num_states):
            if transition_possible(states[i], states[j], n):
                prob = 1.0
                for k in differing_digits(states[i], states[j], n):
                    prob *= forced_prob_directed(state_graph, nodes[k], weight_attr)
                for m in same_zeros(states[i], states[j], n):
                    prob *= 1 - forced_prob_directed(state_graph, nodes[m], weight_attr)
                tm[i, j] = prob

    tm[num_states - 1, num_states - 1] = 1.0
    return tm, time.process_time() - start


def tm_generation_directed_sub_sparse(graphs: list[nx.DiGraph], states: list[int], weight_attr: str = "weight"):
    if not graphs[0].is_directed():
        raise ValueError("Graphs must be directed.")
    n = graphs[0].number_of_nodes()
    num_states = len(states)
    tm = lil_matrix((num_states, num_states))
    start = time.process_time()

    for i in range(num_states):
        state_graph = graphs[i]
        nodes = list(state_graph.nodes())
        for j in range(i, num_states):
            if transition_possible(states[i], states[j], n):
                prob = 1.0
                for k in differing_digits(states[i], states[j], n):
                    prob *= forced_prob_directed(state_graph, nodes[k], weight_attr)
                for m in same_zeros(states[i], states[j], n):
                    prob *= 1 - forced_prob_directed(state_graph, nodes[m], weight_attr)
                tm[i, j] = prob

    tm[num_states - 1, num_states - 1] = 1.0
    return tm.tocsr(), time.process_time() - start


# -------------------------
# Expected propagation time
# -------------------------

def calc_expected(tm: np.ndarray, row: int, solutions: np.ndarray) -> float:
    size = len(tm)
    coeff = 1 - tm[row][row]
    if coeff == 0:
        return float("inf")
    temp = 0.0
    for i in range(row + 1, size):
        if tm[row][i] != 0:
            temp += tm[row][i] * (solutions[i] + 1)
    temp += 1 - coeff
    temp /= coeff
    return temp


def calc_expected_sparse(tm: csr_matrix, row: int, solutions: np.ndarray) -> float:
    coeff = 1 - tm[row, row]
    if coeff == 0:
        return float("inf")
    temp = 0.0
    row_data = tm.getrow(row)
    for col_idx, val in zip(row_data.indices, row_data.data):
        if col_idx != row and val != 0:
            temp += val * (solutions[col_idx] + 1)
    temp += 1 - coeff
    temp /= coeff
    return temp


def can_be_calculated(tm: np.ndarray, row: int, solutions: np.ndarray) -> bool:
    for i in range(len(tm) - 1):
        if tm[row][i] != 0 and solutions[i] == 0 and i != row:
            return False
    return True


def propagation_time_solver(tm: np.ndarray):
    start = time.process_time()
    size = len(tm)
    solutions = np.zeros(size)
    solutions[size - 1] = 0
    for i in range(size - 2, 0, -1):
        solutions[i] = calc_expected(tm, i, solutions)
    solutions[0] = float("inf")
    return solutions, time.process_time() - start


def propagation_time_solver_sparse(tm: csr_matrix):
    start = time.process_time()
    size = tm.shape[0]
    solutions = np.zeros(size)
    solutions[size - 1] = 0
    for i in range(size - 2, 0, -1):
        solutions[i] = calc_expected_sparse(tm, i, solutions)
    solutions[0] = float("inf")
    return solutions, time.process_time() - start


def propagation_time_solver_inverse(tm: np.ndarray):
    sub = tm[1:, 1:].copy()
    size = len(sub)
    for row in range(size):
        sub[row, -1] -= 1
    start = time.process_time()
    temp = la.solve_triangular(sub - np.identity(size), np.identity(size))
    vector = temp[:, size - 1] + 1
    vector = np.insert(vector, 0, float("inf"))
    return vector, time.process_time() - start


def propagation_time_solver_algebraic(tm: np.ndarray):
    size = len(tm)
    mus = sp.symbols(f"mu0:{size}")
    a = sp.Matrix(tm)
    eqns = []
    for i in range(size):
        if i != size - 1:
            eqns.append(a.row(i).dot(mus) + 1 - mus[i])
        else:
            eqns.append(a.row(i).dot(mus))
    return list(sp.linsolve(eqns[-(size - 1):], mus[-(size - 1):]))


# -------------------------
# Simulation + misc
# -------------------------

def all_nodes_blue(graph: nx.Graph) -> bool:
    return all(data.get("color") == "blue" for _, data in graph.nodes(data=True))


def rzf_simulation(target_graph: nx.DiGraph, starting_set: list, seed: Optional[int] = None) -> int:
    rng = random.Random(seed)
    graph = target_graph.copy()
    start_set = set(starting_set)

    for node in graph.nodes():
        graph.nodes[node]["color"] = "blue" if node in start_set else "white"

    iterations = 0
    while not all_nodes_blue(graph):
        iterations += 1
        turned_blue = []
        for node in graph.nodes():
            if graph.nodes[node]["color"] == "blue":
                continue
            indeg = graph.in_degree(node)
            if indeg > 0:
                prob = len(blue_suppliers(graph, node)) / indeg
                if rng.random() <= prob:
                    turned_blue.append(node)
        for node in turned_blue:
            graph.nodes[node]["color"] = "blue"
    return iterations


def rzf_simulation_weighted(graph: nx.DiGraph, starting_set: list, weight_attr: str = "weight",
                            seed: Optional[int] = None, max_steps: int = 1_000_000) -> int:
    rng = random.Random(seed)
    g = graph.copy()
    start_set = set(starting_set)

    for node in g.nodes():
        g.nodes[node]["color"] = "blue" if node in start_set else "white"

    steps = 0
    while not all_nodes_blue(g):
        steps += 1
        if steps > max_steps:
            raise RuntimeError("Reached max_steps without full propagation.")
        turned = []
        for v in g.nodes():
            if g.nodes[v]["color"] == "blue":
                continue
            total_w = 0.0
            blue_w = 0.0
            for u, _, data in g.in_edges(v, data=True):
                w = data.get(weight_attr, 1.0)
                total_w += w
                if g.nodes[u]["color"] == "blue":
                    blue_w += w
            if total_w > 0 and rng.random() <= blue_w / total_w:
                turned.append(v)
        for v in turned:
            g.nodes[v]["color"] = "blue"
    return steps


def rzf_simulation_all_nodes(graph: nx.DiGraph, trials_per_node: int) -> list[tuple]:
    results = []
    for node in graph.nodes():
        avg = 0.0
        for _ in range(trials_per_node):
            avg += rzf_simulation_weighted(graph, [node])
        results.append((node, avg / trials_per_node))
    return results


def optimal_zfs_by_size(propagation_times: np.ndarray, size: int, n: int):
    indices = nums_with_bitcount(2**n, size)
    best = min(propagation_times[i] for i in indices)
    return [(i, best) for i in indices if propagation_times[i] == best]


def return_transition_times(target_graph: nx.Graph):
    graphs = graph_gen(target_graph)
    if nx.is_directed(target_graph):
        tm, _ = tm_generation_directed(graphs)
    else:
        tm, _ = tm_generation(graphs)
    times, _ = propagation_time_solver(tm)
    return times


def pad_top_left_with_one(matrix: np.ndarray) -> np.ndarray:
    rows, cols = matrix.shape
    padded = np.zeros((rows + 1, cols + 1), dtype=matrix.dtype)
    padded[1:, 1:] = matrix
    padded[0, 0] = 1
    return padded


def pad_top_left_with_one_sparse(matrix: csr_matrix) -> csr_matrix:
    rows, cols = matrix.shape
    top_row = csr_matrix(([1.0], ([0], [0])), shape=(1, cols + 1))
    left_col = csr_matrix((rows, 1))
    bottom = hstack([left_col, matrix])
    return vstack([top_row, bottom]).tocsr()