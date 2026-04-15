import time
from collections import defaultdict

import networkx as nx
import numpy as np

from rzf_core import (
    larger_bin,
    differing_digits,
    same_zeros,
    forced_prob,
    num_blue_neighbors,
)


def color_match(n1, n2):
    return n1["color"] == n2["color"]


def blue_degree_sequence(graphs, index):
    graph = graphs[index]
    seq = [graph.degree[v] for v in graph.nodes() if graph.nodes[v].get("color") == "blue"]
    seq.sort(reverse=True)
    return seq


def white_degree_sequence(graphs, index):
    graph = graphs[index]
    seq = [graph.degree[v] for v in graph.nodes() if graph.nodes[v].get("color") == "white"]
    seq.sort(reverse=True)
    return seq


def blue_then_white_degree_sequence(graphs, index):
    return tuple(blue_degree_sequence(graphs, index)), tuple(white_degree_sequence(graphs, index))


def blue_degree_sequence_groups(graphs):
    groups = defaultdict(list)
    for i in range(len(graphs)):
        groups[tuple(blue_degree_sequence(graphs, i))].append(i)
    return groups


def white_degree_sequence_groups(graphs):
    groups = defaultdict(list)
    for i in range(len(graphs)):
        groups[tuple(white_degree_sequence(graphs, i))].append(i)
    return groups


def degree_sequence_groups(graphs):
    groups = defaultdict(list)
    for i in range(len(graphs)):
        groups[tuple(blue_then_white_degree_sequence(graphs, i))].append(i)
    return groups


def pre_automorphism_checking_dictionary(graphs):
    groups = defaultdict(list)
    for i in range(len(graphs)):
        key = tuple(blue_degree_sequence(graphs, i)) + tuple(sorted(num_blue_neighbors(graphs[i]).values()))
        groups[key].append(i)
    return groups


def return_automorphic_states(graphs, index, automorphic_keys):
    key = tuple(blue_degree_sequence(graphs, index)) + tuple(sorted(num_blue_neighbors(graphs[index]).values()))
    ans = []
    for state in automorphic_keys[key]:
        if nx.is_isomorphic(graphs[index], graphs[state], node_match=color_match):
            ans.append(state)
    return ans


def automorphism_groups(graphs):
    start = time.process_time()
    automorphic_keys = pre_automorphism_checking_dictionary(graphs)
    states = list(range(len(graphs)))
    groups = []

    while states:
        group = return_automorphic_states(graphs, states[0], automorphic_keys)
        groups.append(group)
        state_set = set(group)
        states = [x for x in states if x not in state_set]

    print(f"Time needed to generate automorphism groups {time.process_time() - start} seconds")
    return groups


def generate_automorphism_matrix(automorphism_groups_list, graphs):
    start = time.process_time()
    num_states = len(automorphism_groups_list)
    matrix = np.zeros((num_states, num_states))
    matrix[0][0] = 1

    g0 = graphs[0]
    n = g0.number_of_nodes()

    for i in range(num_states):
        starting_state = automorphism_groups_list[i][0]
        starting_graph = graphs[starting_state]
        possible_states = larger_bin(starting_state, graphs)
        blue_counts = num_blue_neighbors(starting_graph)
        nodes = list(starting_graph.nodes())

        for j in range(i, num_states):
            prob = 0.0
            target_group = automorphism_groups_list[j]
            transition_possible_set = sorted(set(possible_states).intersection(target_group))

            for graph_state in transition_possible_set:
                temp = 1.0
                for k in differing_digits(starting_state, graph_state, n):
                    temp *= forced_prob(starting_graph, nodes[k], blue_counts)
                for m in same_zeros(starting_state, graph_state, n):
                    temp *= 1 - forced_prob(starting_graph, nodes[m], blue_counts)
                prob += temp

            matrix[i][j] = prob

    print(f"Time needed to generate automorphism transition matrix {time.process_time() - start} seconds")
    return matrix