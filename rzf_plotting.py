import textwrap

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.cm import ScalarMappable, get_cmap
from matplotlib.colors import Normalize
from matplotlib.patches import FancyArrowPatch


def print_graph(G):
    pos = nx.spring_layout(G, seed=42)
    plt.figure(figsize=(6, 4))
    nx.draw_networkx(
        G,
        pos,
        with_labels=True,
        node_color="skyblue",
        edge_color="gray",
        node_size=1500,
        font_size=12,
    )
    edge_labels = {(u, v): d.get("weight", "") for u, v, d in G.edges(data=True)}
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color="red")
    plt.axis("off")
    plt.show()


def print_graph_bidirectional(G):
    if not G.is_directed():
        G = nx.DiGraph(G)

    pos = nx.spring_layout(G, seed=42, k=0.8)
    fig, ax = plt.subplots(figsize=(6, 4))

    nx.draw_networkx_nodes(
        G,
        pos,
        ax=ax,
        node_color="skyblue",
        edgecolors="black",
        node_size=1200,
        linewidths=0.8,
    )
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=12)

    def curved_arrow(u, v, rad):
        arrow = FancyArrowPatch(
            pos[u],
            pos[v],
            connectionstyle=f"arc3,rad={rad}",
            arrowstyle="-|>",
            mutation_scale=16,
            lw=1.8,
            color="gray",
            shrinkA=18,
            shrinkB=18,
            zorder=3,
            clip_on=False,
        )
        ax.add_patch(arrow)

    seen = set()
    singles = []
    for u, v in G.edges():
        if (u, v) in seen:
            continue
        if G.has_edge(v, u):
            curved_arrow(u, v, 0.25)
            curved_arrow(v, u, -0.25)
            seen.add((u, v))
            seen.add((v, u))
        else:
            singles.append((u, v))
            seen.add((u, v))

    for u, v in singles:
        curved_arrow(u, v, 0.0)

    ax.set_axis_off()
    fig.tight_layout()
    plt.show()


def plot_rzf_simulation_results(
    G,
    results,
    pos=None,
    node_size=350,
    cmap="viridis",
    with_labels=True,
    title=None,
    vmin=None,
    vmax=None,
    ax=None,
    edge_arrows=None,
    save_path=None,
):
    res_map = dict(results)
    node_list = list(G.nodes())
    values = np.array([res_map.get(n, np.nan) for n in node_list], dtype=float)

    finite_vals = values[np.isfinite(values)]
    if finite_vals.size == 0:
        raise ValueError("No finite values found.")

    if vmin is None:
        vmin = float(np.nanmin(finite_vals))
    if vmax is None:
        vmax = float(np.nanmax(finite_vals))
    if vmin == vmax:
        vmax = vmin + 1e-9

    if pos is None:
        pos = nx.spring_layout(G, seed=42)

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    else:
        fig = ax.figure

    if edge_arrows is None:
        edge_arrows = G.is_directed()

    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=node_list,
        node_size=node_size,
        node_color=values,
        cmap=get_cmap(cmap),
        vmin=vmin,
        vmax=vmax,
        ax=ax,
    )
    nx.draw_networkx_edges(G, pos, ax=ax, arrows=edge_arrows)

    if with_labels:
        nx.draw_networkx_labels(G, pos, font_size=9, ax=ax)

    norm = Normalize(vmin=vmin, vmax=vmax)
    sm = ScalarMappable(norm=norm, cmap=plt.cm.get_cmap(cmap))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.85, pad=0.02)
    cbar.set_label("Avg RZF iterations", rotation=90)

    if title:
        ax.set_title(title)
    ax.set_axis_off()
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig, ax


def plot_graph_with_zfs_heatmap(
    G,
    zfs_size_one_times,
    id_to_sector=None,
    layout="circular",
    shift_labels_by_1=True,
    wrap_width=32,
    key_fontsize=9,
):
    """
    Plot a graph with singleton-start ZFS/RZF propagation times shown as node colors.

    Parameters
    ----------
    G : networkx graph
        Graph to plot.
    zfs_size_one_times : list of (node, time)
        Example: [(0, 4.2), (1, 5.1), ...]
    id_to_sector : dict, optional
        Mapping from node id to label text for the legend/key.
    layout : str, optional
        One of {"circular", "spring", "kamada"}.
        For BEA, use "circular".
    shift_labels_by_1 : bool, optional
        If True, node labels are displayed as 1..n instead of 0..n-1.
    wrap_width : int, optional
        Width for wrapped legend text.
    key_fontsize : int, optional
        Font size for the legend/key.
    """
    times = {node: t for node, t in zfs_size_one_times}
    nodes = list(G.nodes())
    node_values = np.array([times[node] for node in nodes], dtype=float)

    finite_vals = node_values[np.isfinite(node_values)]
    if finite_vals.size == 0:
        raise ValueError("All node values are non-finite (inf/NaN).")
    vmin, vmax = float(finite_vals.min()), float(finite_vals.max())

    if layout == "circular":
        pos = nx.circular_layout(G)
    elif layout == "spring":
        pos = nx.spring_layout(G, seed=42)
    elif layout == "kamada":
        pos = nx.kamada_kawai_layout(G)
    else:
        raise ValueError("Unknown layout. Use 'circular', 'spring', or 'kamada'.")

    base = plt.cm.autumn
    compressed_colors = base(np.linspace(0.35, 0.65, 256))
    compressed_cmap = mcolors.LinearSegmentedColormap.from_list(
        "compressed_autumn", compressed_colors
    )

    fig, ax = plt.subplots(figsize=(14, 10))
    ax.set_axis_off()

    weights = [d.get("weight", 1.0) for _, _, d in G.edges(data=True)]
    w_min, w_max = (min(weights), max(weights)) if weights else (1.0, 1.0)

    def norm_w(w):
        if abs(w_max - w_min) < 1e-12:
            return 0.6
        return 0.15 + 0.85 * ((w - w_min) / (w_max - w_min))

    for u, v, d in G.edges(data=True):
        nx.draw_networkx_edges(
            G,
            pos,
            ax=ax,
            edgelist=[(u, v)],
            arrows=True,
            arrowstyle="->",
            width=1.5,
            alpha=norm_w(float(d.get("weight", 1.0))),
            connectionstyle="arc3,rad=0.05",
            edge_color="black",
        )

    node_collection = nx.draw_networkx_nodes(
        G,
        pos,
        ax=ax,
        nodelist=nodes,
        node_color=node_values,
        cmap=compressed_cmap,
        vmin=vmin,
        vmax=vmax,
        node_size=800,
    )

    labels = (
        {node: str(int(node) + 1) for node in nodes}
        if shift_labels_by_1
        else {node: str(int(node)) for node in nodes}
    )
    nx.draw_networkx_labels(G, pos, ax=ax, labels=labels, font_size=8)

    cbar = fig.colorbar(node_collection, ax=ax, fraction=0.046, pad=0.02)
    cbar.set_label("ZFS propagation time from singleton start")

    if id_to_sector is not None:
        def wrap_entry(prefix, text, width):
            wrapped = textwrap.fill(
                text,
                width=width,
                subsequent_indent=" " * (len(prefix) + 1),
                break_long_words=False,
                break_on_hyphens=False,
            )
            return f"{prefix} {wrapped}"

        key_lines = ["Node key", ""]
        for k in sorted(id_to_sector):
            shown_k = k + 1 if shift_labels_by_1 else k
            prefix = f"{shown_k:>2} →"
            key_lines.append(wrap_entry(prefix, id_to_sector[k], wrap_width))

        key_text = "\n".join(key_lines)

        fig.text(
            0.02,
            0.5,
            key_text,
            va="center",
            ha="left",
            fontsize=key_fontsize,
            family="monospace",
            bbox=dict(
                boxstyle="round,pad=0.4",
                facecolor="none",
                edgecolor="0.4",
                linewidth=0.8,
            ),
        )

    plt.subplots_adjust(left=0.22, right=0.98, top=0.98, bottom=0.02)
    plt.show()