"""Generate offline interactive websites from :class:`DiNetxify.visualization.Plot`."""

from __future__ import annotations

import json
import math
import re
import shutil
from datetime import datetime
from decimal import Decimal, InvalidOperation
from importlib import resources
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, TYPE_CHECKING

import networkx as nx
import numpy as np
import pandas as pd
from plotly.offline import get_plotlyjs

if TYPE_CHECKING:
    from .visualization import Plot


LAYOUT_SEED = 42
LOUVAIN_ATTEMPTS = 100
REFERENCE_2D_MEDIAN_DIAMETER_FRACTION = 0.024
MAX_2D_DISPLAY_SCALE = 4.0
REFERENCE_3D_MEDIAN_RADIUS = 3.75
MODULE_COLORS = [
    "#2563eb", "#dc2626", "#059669", "#7c3aed", "#d97706",
    "#0891b2", "#be185d", "#4d7c0f", "#9333ea", "#0f766e",
    "#c2410c", "#475569", "#0369a1", "#a21caf", "#15803d",
]


def _to_float(value: Any) -> Optional[float]:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _safe_exp(value: Any) -> Optional[float]:
    numeric = _to_float(value)
    if numeric is None:
        return None
    try:
        result = math.exp(numeric)
    except OverflowError:
        return None
    return result if math.isfinite(result) and result > 0 else None


def _normalize_phecode(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    if not text:
        return ""
    try:
        normalized = format(Decimal(text).normalize(), "f")
        if "." in normalized:
            normalized = normalized.rstrip("0").rstrip(".")
        return normalized
    except (InvalidOperation, ValueError):
        return text


def _as_bool(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if value is None or pd.isna(value):
        return False
    return str(value).strip().lower() in {"1", "true", "yes"}


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (pd.Timestamp, datetime)):
        return value.isoformat()
    if pd.isna(value):
        return None
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def _effect_measure(model_type: Any) -> str:
    model = str(model_type or "").lower()
    if "logistic" in model:
        return "Odds ratio"
    return "Hazard ratio"


def _system_record(plot: "Plot", system: str) -> Dict[str, str]:
    key = str(system or "").strip()
    return {
        "id": key,
        "name": plot.sys_dict.get(key, key.title() or "Unclassified"),
        "color": plot._system_color.get(key, "#94a3b8"),
    }


def _count_label(plot: "Plot") -> str:
    if plot._outcome_name:
        return "Outcome-group individuals with condition"
    if plot._exposure_name:
        return "Exposed individuals with condition"
    return "Individuals with condition"


def _count_description(plot: "Plot") -> str:
    if plot._outcome_name:
        return "the number of individuals in the outcome group with the condition"
    if plot._exposure_name:
        return "the number of exposed individuals with the condition"
    return "the number of individuals with the condition"


def _phewas_lookup(plot: "Plot") -> Dict[str, Dict[str, Any]]:
    lookup: Dict[str, Dict[str, Any]] = {}
    for _, row in plot._phewas.iterrows():
        phecode = _normalize_phecode(row.get(plot.phecode_col))
        if phecode:
            lookup[phecode] = row.to_dict()
    return lookup


def _parse_group(value: Any) -> Dict[str, Optional[float]]:
    match = re.search(r"([\d,]+)\s*/\s*([\d,.]+)\s*\(([\d.]+)\)", str(value or ""))
    if match is None:
        return {"events": None, "person_years": None, "rate": None}
    return {
        "events": _to_float(match.group(1).replace(",", "")),
        "person_years": _to_float(match.group(2).replace(",", "")),
        "rate": _to_float(match.group(3)),
    }


def _build_phewas_payload(plot: "Plot") -> Dict[str, Any]:
    rows = []
    systems = set()
    mode = "before_outcome" if plot._outcome_name else "after_exposure"
    subject_name = plot._outcome_name or plot._exposure_name or "DiNetxify analysis"
    analysis_name = subject_name
    group_one_col, group_two_col = (
        ("outcome_group", "no_outcome_group")
        if mode == "before_outcome"
        else ("exposed_group", "unexposed_group")
    )
    group_one_label, group_two_label = (
        ("Outcome group", "No-outcome group")
        if mode == "before_outcome"
        else ("Exposed group", "Unexposed group")
    )
    for _, source in plot._phewas.iterrows():
        phecode = _normalize_phecode(source.get(plot.phecode_col))
        if not phecode:
            continue
        system = str(source.get(plot.system_col) or "").strip()
        systems.add(system)
        beta = _to_float(source.get(plot.phewas_coef_col))
        se = _to_float(source.get(plot.phewas_se_col))
        ratio = _safe_exp(beta)
        lower = _safe_exp(beta - 1.96 * se) if beta is not None and se is not None else None
        upper = _safe_exp(beta + 1.96 * se) if beta is not None and se is not None else None
        p_value = _to_float(source.get("phewas_p"))
        q_value = _to_float(source.get("phewas_p_adjusted"))
        count = _to_float(source.get(plot.phewas_number_col))
        effect_measure = _effect_measure(source.get("model_type"))
        group_one = _parse_group(source.get(group_one_col))
        group_two = _parse_group(source.get(group_two_col))
        significant = _as_bool(source.get(plot.phewas_significance_col))
        rows.append(
            {
                "row_id": len(rows) + 1,
                "analysis_id": "analysis",
                "analysis_name": analysis_name,
                "cancer_id": "analysis",
                "cancer_name": analysis_name,
                "phecode": phecode,
                "phecode_sort": _to_float(phecode) if _to_float(phecode) is not None else 1e12,
                "disease_name": str(source.get(plot.disease_col) or phecode).strip(),
                "system": system,
                "system_display": plot.sys_dict.get(system, system.title() or "Unclassified"),
                "disease_system": system,
                "disease_system_display": plot.sys_dict.get(system, system.title() or "Unclassified"),
                "system_color": plot._system_color.get(system, "#94a3b8"),
                "sex": str(source.get("sex") or "Both"),
                "count": count,
                "beta": beta,
                "se": se,
                "effect_ratio": ratio,
                "hr_original": ratio,
                "hr_plot": min(ratio, 30.0) if ratio is not None else None,
                "ci_lower": lower,
                "ci_upper": upper,
                "p_value": p_value,
                "q_value": q_value,
                "significant": significant,
                "positive_significant": bool(significant and (ratio is None or ratio > 1.0)),
                "direction": "protective" if ratio is not None and ratio < 1.0 else "elevated",
                "model_type": str(source.get("model_type") or ""),
                "effect_measure": effect_measure,
                "events_group_one": group_one["events"],
                "person_years_group_one": group_one["person_years"],
                "incidence_rate_group_one": group_one["rate"],
                "events_group_two": group_two["events"],
                "person_years_group_two": group_two["person_years"],
                "incidence_rate_group_two": group_two["rate"],
                "group_one_label": group_one_label,
                "group_two_label": group_two_label,
                "incidence_rate_cancer": group_one["rate"],
                "incidence_rate_matched_unexposed": group_two["rate"],
                "comparison_label": group_two_label,
            }
        )

    rows.sort(
        key=lambda row: (
            row["system_display"],
            row["p_value"] if row["p_value"] is not None else math.inf,
            row["disease_name"],
        )
    )
    has_effects = any(row["effect_ratio"] is not None for row in rows)
    effect_measure = next(
        (row["effect_measure"] for row in rows if row["effect_ratio"] is not None),
        "Case count",
    )
    count_label = _count_label(plot)
    analysis = {
        "id": "analysis",
        "slug": "analysis",
        "name": analysis_name,
        "short_name": analysis_name,
        "design": mode,
        "sex": "All individuals",
        "comparison_label": group_two_label,
    }
    return {
        "metadata": {
            "generated_at": datetime.now().isoformat(timespec="seconds"),
            "mode": mode,
            "subject_name": subject_name,
            "count_column": plot.phewas_number_col,
            "count_label": count_label,
            "count_description": _count_description(plot),
            "row_count": len(rows),
            "has_effect_estimates": has_effects,
            "effect_measure": effect_measure,
            "effect_cap_for_plotting": 30,
        },
        "cancers": [analysis],
        "analyses": [analysis],
        "systems": [_system_record(plot, system) for system in sorted(systems)],
        "rows": rows,
    }


def _node_radius(count: Any, low: float = 0.0, high: float = 1.0) -> float:
    value = _to_float(count)
    if value is None:
        return 5.0
    return (1.0 / 2.4) * max(value, 1.0) ** 0.46600437333691014


def _network_display_scale(nodes: list[Dict[str, Any]], stats: Dict[str, Any]) -> float:
    radii = [float(node["radius"]) for node in nodes if _to_float(node.get("radius"))]
    width = _to_float(stats.get("max_x")) - _to_float(stats.get("min_x"))
    height = _to_float(stats.get("max_y")) - _to_float(stats.get("min_y"))
    if not radii or width <= 0 or height <= 0:
        return 1.0
    scale = REFERENCE_2D_MEDIAN_DIAMETER_FRACTION * min(width, height) / (2 * float(np.median(radii)))
    return float(np.clip(scale, 0.5, MAX_2D_DISPLAY_SCALE))


def _three_d_size_reduction(nodes: list[Dict[str, Any]]) -> float:
    radii = [
        (3 * float(node["cases"]) / (4 * math.pi)) ** (1 / 3)
        for node in nodes
        if _to_float(node.get("cases")) is not None and float(node["cases"]) > 0
    ]
    if not radii:
        return 1.0
    return float(np.clip(REFERENCE_3D_MEDIAN_RADIUS / float(np.median(radii)), 0.5, 2.0))


def _resolve_collisions(
    positions: Dict[Any, np.ndarray],
    radii: Dict[Any, float],
    padding: float,
    max_iterations: int = 250,
) -> None:
    nodes = sorted(positions)
    if len(nodes) < 2:
        return
    cell_size = max(radii.values(), default=8.0) * 2 + padding
    for iteration in range(max_iterations):
        buckets: Dict[tuple[int, int], list[str]] = {}
        for node in nodes:
            x, y = positions[node]
            key = (int(math.floor(x / cell_size)), int(math.floor(y / cell_size)))
            buckets.setdefault(key, []).append(node)

        moved = False
        checked = set()
        for key, bucket_nodes in buckets.items():
            candidates = []
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    candidates.extend(buckets.get((key[0] + dx, key[1] + dy), []))
            for left in bucket_nodes:
                for right in candidates:
                    if left == right:
                        continue
                    pair = tuple(sorted((left, right)))
                    if pair in checked:
                        continue
                    checked.add(pair)
                    delta = positions[right] - positions[left]
                    distance = float(np.linalg.norm(delta))
                    required = radii[left] + radii[right] + padding
                    if distance >= required:
                        continue
                    if distance < 1e-9:
                        angle = (sum(ord(char) for char in str(right)) + iteration * 37) % 360
                        direction = np.array(
                            [math.cos(math.radians(angle)), math.sin(math.radians(angle))]
                        )
                    else:
                        direction = delta / distance
                    shift = direction * ((required - distance) / 2 + 0.05)
                    positions[left] -= shift
                    positions[right] += shift
                    moved = True
        if not moved:
            break


def _module_aware_layout(
    graph: nx.Graph,
    module_by_node: Dict[str, int],
    radii: Dict[str, float],
) -> tuple[Dict[str, tuple[float, float]], list[Dict[str, Any]]]:
    module_nodes: Dict[int, list[str]] = {}
    for node in graph.nodes:
        module_nodes.setdefault(module_by_node[node], []).append(node)

    local_positions: Dict[int, Dict[str, np.ndarray]] = {}
    module_radii: Dict[int, float] = {}
    for module, nodes in sorted(module_nodes.items()):
        subgraph = graph.subgraph(nodes)
        if len(nodes) == 1:
            positions = {nodes[0]: np.zeros(2, dtype=float)}
        else:
            scale = max(45.0, 20.0 * math.sqrt(len(nodes)))
            raw = nx.spring_layout(
                subgraph,
                seed=LAYOUT_SEED + int(module),
                weight="layout_weight",
                scale=scale,
                iterations=150,
            )
            positions = {str(node): np.asarray(value, dtype=float) for node, value in raw.items()}
        _resolve_collisions(positions, {node: radii[node] for node in nodes}, padding=4.0)
        center = np.mean(list(positions.values()), axis=0)
        for node in positions:
            positions[node] -= center
        local_positions[module] = positions
        module_radii[module] = max(
            float(np.linalg.norm(positions[node])) + radii[node] for node in nodes
        ) + 24.0

    module_graph = nx.Graph()
    module_graph.add_nodes_from(module_nodes)
    for source, target, attrs in graph.edges(data=True):
        left = module_by_node[source]
        right = module_by_node[target]
        if left == right:
            continue
        weight = float(attrs.get("layout_weight", 1.0))
        previous = module_graph.get_edge_data(left, right, {}).get("weight", 0.0)
        module_graph.add_edge(left, right, weight=previous + weight)

    if len(module_graph) == 1:
        centers = {next(iter(module_graph.nodes)): np.zeros(2, dtype=float)}
    else:
        center_scale = max(280.0, sum(module_radii.values()) / math.pi)
        raw_centers = nx.spring_layout(
            module_graph,
            seed=LAYOUT_SEED,
            weight="weight",
            scale=center_scale,
            iterations=200,
        )
        centers = {
            int(module): np.asarray(position, dtype=float)
            for module, position in raw_centers.items()
        }
        _resolve_collisions(centers, module_radii, padding=65.0, max_iterations=400)

    final_positions: Dict[str, tuple[float, float]] = {}
    module_records = []
    for module, nodes in sorted(module_nodes.items()):
        points = []
        for node in nodes:
            position = local_positions[module][node] + centers[module]
            final_positions[node] = (float(position[0]), float(position[1]))
            points.append(position)
        xs = [float(point[0]) for point in points]
        ys = [float(point[1]) for point in points]
        pad = max(radii[node] for node in nodes) + 22.0
        module_records.append(
            {
                "display_number": module,
                "label": f"Module {module}",
                "color": MODULE_COLORS[(module - 1) % len(MODULE_COLORS)],
                "node_count": len(nodes),
                "min_x": min(xs) - pad,
                "max_x": max(xs) + pad,
                "min_y": min(ys) - pad,
                "max_y": max(ys) + pad,
            }
        )
    return final_positions, module_records


def _build_comorbidity_payload(
    plot: "Plot",
) -> tuple[Dict[str, Any], nx.Graph, Dict[str, int]]:
    frame = plot._base_comorbidity.copy()
    graph = nx.Graph()
    phewas = _phewas_lookup(plot)
    analysis_name = plot._outcome_name or plot._exposure_name or "DiNetxify analysis"
    if frame.empty:
        empty_network = {
            "id": "analysis", "name": analysis_name, "short_name": analysis_name,
            "nodes": [], "edges": [], "modules": [], "systems": [],
            "stats": {"node_count": 0, "edge_count": 0, "module_count": 0},
        }
        return {
            "metadata": {"generated_at": datetime.now().isoformat(timespec="seconds")},
            "networks": [empty_network],
        }, graph, {}

    raw_clusters = plot.get_louvain_clusters(max_attempts=LOUVAIN_ATTEMPTS)
    raw_cluster_by_id = {
        _normalize_phecode(node): cluster for node, cluster in raw_clusters.items()
    }
    ordered_raw = sorted(set(raw_cluster_by_id.values()), key=lambda value: str(value))
    display_by_raw = {raw: index + 1 for index, raw in enumerate(ordered_raw)}
    module_by_id = {
        node: display_by_raw[cluster] for node, cluster in raw_cluster_by_id.items()
    }

    edges = []
    for _, row in frame.iterrows():
        source = _normalize_phecode(row.get(plot.source_col))
        target = _normalize_phecode(row.get(plot.target_col))
        beta = _to_float(row.get(plot.comorbidity_beta_col))
        if not source or not target or beta is None:
            continue
        effect_ratio = _safe_exp(beta)
        if effect_ratio is None:
            continue
        layout_weight = min(max(effect_ratio, 1.0), 10.0)
        if graph.has_edge(source, target):
            if beta <= graph[source][target]["beta"]:
                continue
        graph.add_edge(
            source,
            target,
            beta=beta,
            effect_ratio=effect_ratio,
            layout_weight=layout_weight,
        )
        se = _to_float(row.get("comorbidity_se"))
        edges.append(
            {
                "id": f"edge-{len(edges)}",
                "source": source,
                "target": target,
                "beta": beta,
                "effect_ratio": effect_ratio,
                "filter_weight": layout_weight,
                "weight": layout_weight,
                "width": math.exp(-1.60276782) * layout_weight ** 1.1944912184039562,
                "or": effect_ratio,
                "hr": effect_ratio,
                "ci_lower": _safe_exp(beta - 1.96 * se) if se is not None else None,
                "ci_upper": _safe_exp(beta + 1.96 * se) if se is not None else None,
                "p_value": _to_float(row.get("comorbidity_p")),
                "q_value": _to_float(row.get("comorbidity_p_adjusted")),
                "n_total": _to_float(row.get("n_total")),
            }
        )

    radii = {
        node: _node_radius(phewas.get(node, {}).get(plot.phewas_number_col))
        for node in graph.nodes
    }
    graph_modules = {
        node: module_by_id.get(node, max(module_by_id.values(), default=0) + 1)
        for node in graph.nodes
    }
    positions, modules = _module_aware_layout(graph, graph_modules, radii)
    preliminary_xs = [position[0] for position in positions.values()] or [0.0]
    preliminary_ys = [position[1] for position in positions.values()] or [0.0]
    preliminary_max_radius = max(radii.values(), default=10.0)
    preliminary_stats = {
        "min_x": min(preliminary_xs) - preliminary_max_radius - 40,
        "max_x": max(preliminary_xs) + preliminary_max_radius + 40,
        "min_y": min(preliminary_ys) - preliminary_max_radius - 40,
        "max_y": max(preliminary_ys) + preliminary_max_radius + 40,
    }
    layout_scale = _network_display_scale(
        [{"radius": radius} for radius in radii.values()],
        preliminary_stats,
    )
    if layout_scale > 1.0:
        display_radii = {node: radius * layout_scale for node, radius in radii.items()}
        positions, modules = _module_aware_layout(graph, graph_modules, display_radii)
        adjusted_positions = {
            node: np.asarray(position, dtype=float) for node, position in positions.items()
        }
        maximum_display_radii = {
            node: radius * MAX_2D_DISPLAY_SCALE for node, radius in radii.items()
        }
        _resolve_collisions(
            adjusted_positions,
            maximum_display_radii,
            padding=4.0,
            max_iterations=400,
        )
        positions = {
            node: (float(position[0]), float(position[1]))
            for node, position in adjusted_positions.items()
        }
        for module in modules:
            members = [
                node for node, display_number in graph_modules.items()
                if display_number == module["display_number"]
            ]
            module_xs = [positions[node][0] for node in members]
            module_ys = [positions[node][1] for node in members]
            pad = max(maximum_display_radii[node] for node in members) + 22.0
            module.update(
                {
                    "min_x": min(module_xs) - pad,
                    "max_x": max(module_xs) + pad,
                    "min_y": min(module_ys) - pad,
                    "max_y": max(module_ys) + pad,
                }
            )

    nodes = []
    systems = set()
    for node in sorted(graph.nodes):
        row = phewas.get(node, {})
        system = str(row.get(plot.system_col) or "").strip()
        systems.add(system)
        count = _to_float(row.get(plot.phewas_number_col))
        beta = _to_float(row.get(plot.phewas_coef_col))
        se = _to_float(row.get(plot.phewas_se_col))
        module = graph_modules[node]
        x, y = positions[node]
        node_record = {
            "id": node,
            "phecode": node,
            "disease_name": str(row.get(plot.disease_col) or node).strip(),
            "system": system,
            "system_display": plot.sys_dict.get(system, system.title() or "Unclassified"),
            "disease_system": system,
            "disease_system_display": plot.sys_dict.get(system, system.title() or "Unclassified"),
            "system_color": plot._system_color.get(system, "#94a3b8"),
            "module": module,
            "module_raw_id": module,
            "module_display_number": module,
            "module_label": f"Module {module}",
            "module_clinical_name": None,
            "count": count,
            "cases": count,
            "phewas_effect_ratio": _safe_exp(beta),
            "phewas_hr": _safe_exp(beta),
            "phewas_ci_lower": _safe_exp(beta - 1.96 * se) if beta is not None and se is not None else None,
            "phewas_ci_upper": _safe_exp(beta + 1.96 * se) if beta is not None and se is not None else None,
            "phewas_p_value": _to_float(row.get("phewas_p")),
            "phewas_q_value": _to_float(row.get("phewas_p_adjusted")),
            "phewas_effect_measure": _effect_measure(row.get("model_type")),
            "x": x,
            "y": y,
            "radius": radii[node],
            "degree": int(graph.degree(node)),
            "weighted_degree": float(graph.degree(node, weight="layout_weight")),
        }
        nodes.append(node_record)
        graph.nodes[node].update(
            {
                "x": x,
                "y": y,
                "radius": radii[node],
                "module": module,
                "disease": node_record["disease_name"],
                "system": system,
                "color": node_record["system_color"],
                "count": count if count is not None else 0.0,
            }
        )

    xs = [node["x"] for node in nodes] or [0.0]
    ys = [node["y"] for node in nodes] or [0.0]
    max_radius = max((node["radius"] for node in nodes), default=10.0)
    weights = [edge["filter_weight"] for edge in edges] or [1.0]
    degrees = sorted(node["degree"] for node in nodes)
    for module in modules:
        module["raw_id"] = module["display_number"]
        module["clinical_name"] = None
    stats = {
        "node_count": len(nodes),
        "edge_count": len(edges),
        "module_count": len(modules),
        "min_x": min(xs) - max_radius - 40,
        "max_x": max(xs) + max_radius + 40,
        "min_y": min(ys) - max_radius - 40,
        "max_y": max(ys) + max_radius + 40,
        "min_edge_weight": min(weights),
        "max_edge_weight": max(weights),
        "major_degree_threshold": degrees[int(len(degrees) * 0.75)] if degrees else 0,
    }
    network = {
        "id": "analysis",
        "name": analysis_name,
        "short_name": analysis_name,
        "type": "analysis_specific",
        "is_consensus": False,
        "edge_strength_label": "Minimum displayed comorbidity effect ratio (capped at 10)",
        "node_count_label": _count_label(plot),
        "node_count_description": _count_description(plot),
        "nodes": nodes,
        "edges": edges,
        "modules": modules,
        "systems": [_system_record(plot, system) for system in sorted(systems)],
        "default_node_scale": _network_display_scale(nodes, stats),
        "three_d_size_reduction": _three_d_size_reduction(nodes),
        "stats": stats,
    }
    return {
        "metadata": {
            "generated_at": datetime.now().isoformat(timespec="seconds"),
            "network_count": 1,
            "default_network_id": "analysis",
            "edge_definition": "Statistically significant positive comorbidity associations retained by Plot filtering.",
            "effect_measure": "Exponentiated model coefficient.",
        },
        "networks": [network],
    }, graph, module_by_id


def _trajectory_positions(graph: nx.DiGraph) -> Dict[str, tuple[float, float]]:
    if len(graph) == 0:
        return {}
    condensation = nx.condensation(graph)
    depth: Dict[int, int] = {}
    for component in nx.topological_sort(condensation):
        predecessors = list(condensation.predecessors(component))
        depth[component] = 0 if not predecessors else max(depth[pred] + 1 for pred in predecessors)
    component_by_node = condensation.graph["mapping"]
    layers: Dict[int, list[str]] = {}
    for node in graph.nodes:
        layers.setdefault(depth[component_by_node[node]], []).append(node)
    max_depth = max(layers, default=0)
    positions = {}
    for layer, nodes in sorted(layers.items()):
        nodes.sort()
        x = 175.0 if max_depth == 0 else 175.0 + 650.0 * layer / max_depth
        for index, node in enumerate(nodes):
            y = 70.0 + 560.0 * (index + 1) / (len(nodes) + 1)
            positions[node] = (x, y)
    return positions


def _build_trajectory_payload(
    plot: "Plot",
    module_by_id: Dict[str, int],
) -> Dict[str, Any]:
    frame = plot._base_trajectory.copy()
    phewas = _phewas_lookup(plot)
    grouped: Dict[int, list[Dict[str, Any]]] = {}
    for _, row in frame.iterrows():
        source = _normalize_phecode(row.get(plot.source_col))
        target = _normalize_phecode(row.get(plot.target_col))
        module = module_by_id.get(source)
        if not source or not target or module is None or module != module_by_id.get(target):
            continue
        beta = _to_float(row.get(plot.trajectory_beta_col))
        effect_ratio = _safe_exp(beta)
        if beta is None or effect_ratio is None:
            continue
        se = _to_float(row.get("trajectory_se"))
        grouped.setdefault(module, []).append(
            {
                "source": source,
                "target": target,
                "beta": beta,
                "effect_ratio": effect_ratio,
                "or": effect_ratio,
                "filter_weight": min(max(effect_ratio, 1.0), 10.0),
                "ci_lower": _safe_exp(beta - 1.96 * se) if se is not None else None,
                "ci_upper": _safe_exp(beta + 1.96 * se) if se is not None else None,
                "p_value": _to_float(row.get("trajectory_p")),
                "q_value": _to_float(row.get("trajectory_p_adjusted")),
            }
        )

    modules = []
    for module, edges in sorted(grouped.items()):
        graph = nx.DiGraph()
        graph.add_edges_from((edge["source"], edge["target"]) for edge in edges)
        positions = _trajectory_positions(graph)
        endpoint = None
        if plot._exposure_name:
            endpoint = {
                "id": "__exposure__",
                "disease_name": plot._exposure_name,
                "is_endpoint": True,
                "endpoint_type": "exposure",
                "x": 25.0,
                "y": 350.0,
                "radius": float(plot._exposure_size or 15),
            }
            for node in sorted(node for node in graph if graph.in_degree(node) == 0):
                edges.append(
                    {
                        "id": f"endpoint-{len(edges)}",
                        "source": endpoint["id"],
                        "target": node,
                        "effect_ratio": None,
                        "filter_weight": 1.0,
                        "is_endpoint_edge": True,
                    }
                )
        elif plot._outcome_name:
            endpoint = {
                "id": "__outcome__",
                "disease_name": plot._outcome_name,
                "is_endpoint": True,
                "endpoint_type": "outcome",
                "x": 975.0,
                "y": 350.0,
                "radius": float(plot._outcome_size or 15),
            }
            for node in sorted(node for node in graph if graph.out_degree(node) == 0):
                edges.append(
                    {
                        "id": f"endpoint-{len(edges)}",
                        "source": node,
                        "target": endpoint["id"],
                        "effect_ratio": None,
                        "filter_weight": 1.0,
                        "is_endpoint_edge": True,
                    }
                )

        nodes = []
        systems = set()
        for node in sorted(graph.nodes):
            row = phewas.get(node, {})
            system = str(row.get(plot.system_col) or "").strip()
            systems.add(system)
            x, y = positions[node]
            nodes.append(
                {
                    "id": node,
                    "phecode": node,
                    "disease_name": str(row.get(plot.disease_col) or node).strip(),
                    "system": system,
                    "system_display": plot.sys_dict.get(system, system.title() or "Unclassified"),
                    "disease_system": system,
                    "disease_system_display": plot.sys_dict.get(system, system.title() or "Unclassified"),
                    "system_color": plot._system_color.get(system, "#94a3b8"),
                    "module": module,
                    "module_raw_id": module,
                    "module_display_number": module,
                    "module_label": f"Module {module}",
                    "count": _to_float(row.get(plot.phewas_number_col)),
                    "cases": _to_float(row.get(plot.phewas_number_col)),
                    "x": x,
                    "y": y,
                    "radius": _node_radius(row.get(plot.phewas_number_col)),
                    "in_degree": int(graph.in_degree(node)),
                    "out_degree": int(graph.out_degree(node)),
                    "is_endpoint": False,
                }
            )
        if endpoint is not None:
            endpoint.update(
                {
                    "phecode": "",
                    "system": "",
                    "system_display": endpoint["endpoint_type"].title(),
                    "disease_system": "",
                    "disease_system_display": endpoint["endpoint_type"].title(),
                    "system_color": "#111827",
                    "module": module,
                    "module_raw_id": module,
                    "module_display_number": module,
                    "module_label": f"Module {module}",
                    "count": None,
                    "cases": None,
                    "in_degree": 0,
                    "out_degree": 0,
                }
            )
            nodes.append(endpoint)
        for index, edge in enumerate(edges):
            edge.setdefault("id", f"trajectory-{module}-{index}")
        modules.append(
            {
                "raw_id": module,
                "display_number": module,
                "label": f"Module {module}",
                "nodes": nodes,
                "edges": edges,
                "systems": [_system_record(plot, system) for system in sorted(systems)],
                "stats": {
                    "node_count": len(nodes),
                    "edge_count": len(edges),
                    "min_edge_weight": min(
                        (edge["filter_weight"] for edge in edges if not edge.get("is_endpoint_edge")),
                        default=1.0,
                    ),
                    "max_edge_weight": max(
                        (edge["filter_weight"] for edge in edges if not edge.get("is_endpoint_edge")),
                        default=1.0,
                    ),
                    "major_degree_threshold": 2,
                    "min_x": 0,
                    "max_x": 960,
                    "min_y": 0,
                    "max_y": 700,
                },
            }
        )

    mode = "before_outcome" if plot._outcome_name else "after_exposure"
    subject_name = plot._outcome_name or plot._exposure_name or "DiNetxify analysis"
    analysis_name = subject_name
    stats = {
        "module_count": len(modules),
        "node_count": len({node["id"] for module in modules for node in module["nodes"]}),
        "edge_count": sum(len(module["edges"]) for module in modules),
    }
    return {
        "metadata": {
            "generated_at": datetime.now().isoformat(timespec="seconds"),
            "mode": mode,
            "subject_name": subject_name,
            "analysis_count": 1,
            "default_analysis_id": "analysis",
            "interpretation": "Arrows indicate recorded temporal ordering and do not establish causality.",
            "count_label": _count_label(plot),
            "count_description": _count_description(plot),
        },
        "analyses": [{
            "id": "analysis",
            "slug": "analysis",
            "name": analysis_name,
            "short_name": analysis_name,
            "design": mode,
            "modules": modules,
            "stats": stats,
        }],
    }


def _nav(root: str, active: str, sections: Iterable[str]) -> str:
    definitions = [
        ("home", "Home", "index.html"),
        ("phewas", "PheWAS", "pages/phewas/index.html"),
        ("comorbidity", "Comorbidity Network", "pages/comorbidity-network/index.html"),
        ("trajectory", "Disease Trajectory", "pages/disease-trajectory/index.html"),
        ("three-d", "3D Disease Network", "pages/3d-disease-network/index.html"),
    ]
    links = []
    available = {"home", *sections}
    for section, label, target in definitions:
        if section not in available:
            continue
        class_name = ' class="active"' if section == active else ""
        links.append(f'<a href="{root}/{target}"{class_name}>{label}</a>')
    return "\n".join(links)


def _shell(
    title: str,
    description: str,
    body: str,
    root: str,
    active: str,
    sections: Iterable[str],
    scripts: str = "",
) -> str:
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <meta name="description" content="{description}">
  <title>{title}</title>
  <link rel="stylesheet" href="{root}/assets/style.css">
</head>
<body>
  <header class="site-header">
    <span class="eyebrow">DiNetxify results</span>
    <h1>{title}</h1>
    <p>{description}</p>
  </header>
  <nav class="site-nav" aria-label="Visualization navigation">
    {_nav(root, active, sections)}
  </nav>
  <main>{body}</main>
  <footer>Interactive visualization generated by DiNetxify.</footer>
  {scripts}
</body>
</html>
"""


def _home_page(sections: list[str], phewas_payload: Dict[str, Any]) -> str:
    phewas_text = (
        "Explore disease-level effect estimates and PheWAS counts."
        if phewas_payload["metadata"]["has_effect_estimates"]
        else "Explore disease-level PheWAS counts."
    )
    cards = [
        ("phewas", "Association results", "PheWAS", phewas_text, "pages/phewas/index.html"),
        ("comorbidity", "Non-temporal relationships", "Comorbidity Network", "Explore positive comorbidity associations retained for visualization and detected modules.", "pages/comorbidity-network/index.html"),
        ("trajectory", "Temporal ordering", "Disease Trajectory", "Explore directed diagnosis ordering within disease modules.", "pages/disease-trajectory/index.html"),
        ("three-d", "Integrated network", "3D Disease Network", "Explore comorbidity structure and temporal ordering together.", "pages/3d-disease-network/index.html"),
    ]
    card_html = "".join(
        f"""<article class="card"><span class="card-kicker">{kicker}</span><h2>{title}</h2>
        <p>{text}</p><a class="button" href="./{link}">Open {title}</a></article>"""
        for section, kicker, title, text, link in cards if section in sections
    )
    section_labels = {
        "phewas": "PheWAS",
        "comorbidity": "comorbidity-network",
        "trajectory": "disease-trajectory",
        "three-d": "integrated 3D",
    }
    available = [section_labels[section] for section in sections]
    if len(available) == 1:
        available_text = available[0]
    elif len(available) == 2:
        available_text = " and ".join(available)
    else:
        available_text = ", ".join(available[:-1]) + ", and " + available[-1]
    body = f"""
    <section class="section-panel"><h2>Interactive disease-network results</h2>
      <p>This offline website presents the available {available_text} results produced by DiNetxify.</p>
    </section>
    <section class="card-grid">{card_html}</section>
    """
    return _shell(
        "DiNetxify Interactive Visualization",
        "Interactive visualization of available DiNetxify analysis results.",
        body, ".", "home", sections,
    )


def _phewas_landing(payload: Dict[str, Any], sections: list[str]) -> str:
    has_effects = payload["metadata"]["has_effect_estimates"]
    scatter_text = (
        "Explore effect estimates and PheWAS counts across all modeled conditions."
        if has_effects else
        "Explore PheWAS counts across all modeled conditions."
    )
    forest = (
        '<article class="card"><span class="card-kicker">Effect estimates</span><h2>Forest plot</h2>'
        '<p>Compare selected condition-specific effect estimates and confidence intervals.</p>'
        '<a class="button" href="./forest.html">Open forest plot</a></article>'
        if has_effects else ""
    )
    body = f"""
    <section class="section-panel"><h2>PheWAS explorer</h2>
      <p>Review modeled conditions by physiological system, significance, disease name, and Phecode.</p>
      <div class="metric-grid"><div class="metric-card"><strong>{payload['metadata']['row_count']}</strong><span>Modeled conditions</span></div></div>
    </section>
    <section class="card-grid">
      <article class="card"><span class="card-kicker">Overview</span><h2>Interactive scatter</h2>
        <p>{scatter_text}</p>
        <a class="button" href="./scatter.html">Open scatter plot</a></article>
      {forest}
    </section>
    """
    return _shell("Interactive PheWAS Results", "Interactive exploration of DiNetxify PheWAS results.", body, "../..", "phewas", sections)


def _phewas_scatter_page(payload: Dict[str, Any], sections: list[str]) -> str:
    if payload["metadata"]["has_effect_estimates"]:
        plot_note = (
            "Effect estimates above 30 are capped for plotting while uncapped values remain "
            "available in tooltips and data files. Darker markers denote statistically "
            "significant results. Search highlights matching conditions without removing "
            "other results."
        )
    else:
        plot_note = (
            f"Vertical position shows {payload['metadata']['count_description']}. Search "
            "highlights matching conditions without removing other results."
        )
    body = f"""
    <section class="filter-panel" aria-label="PheWAS scatter plot filters">
      <div><label for="cancerFilter">Analysis</label><select id="cancerFilter"></select></div>
      <div><label for="systemFilter">Physiological system</label><select id="systemFilter"></select></div>
      <div><label for="scatterSearch">Highlight condition or Phecode</label><input id="scatterSearch" type="search" placeholder="Example: diabetes or 250.2"></div>
    </section>
    <section class="summary-row" aria-label="Current scatter summary">
      <span><strong id="cancerName"></strong></span>
      <span><strong id="retainedCount">0</strong> available results</span>
      <span><strong id="significantCount">0</strong> significant results</span>
      <span><strong id="visibleCount">0</strong> displayed after system filter</span>
      <span><strong id="highlightedCount">0</strong> highlighted by search</span>
    </section>
    <section class="plot-card">
      <div class="plot-toolbar" aria-label="Scatter plot toolbar">
        <button id="scatterDownloadPng" type="button">PNG</button><button id="scatterDownloadSvg" type="button">SVG</button>
        <button id="scatterZoomIn" type="button">+</button><button id="scatterZoomOut" type="button">-</button><button id="scatterResetZoom" type="button">Reset</button>
      </div>
      <div id="phewasPlot"></div>
    </section>
    <p class="plot-note">{plot_note}</p>
    """
    scripts = '<script src="../../assets/plotly.min.js"></script><script src="../../data/phewas_data.js"></script><script src="../../assets/phewas_app.js"></script><script>window.PhewasApp.initScatter();</script>'
    return _shell("Interactive PheWAS Scatter Plot", "Associations grouped by physiological system and ordered by Phecode.", body, "../..", "phewas", sections, scripts)


def _phewas_forest_page(payload: Dict[str, Any], sections: list[str]) -> str:
    effect_measure = payload["metadata"]["effect_measure"]
    body = f"""
    <section class="settings-grid" aria-label="Forest plot settings">
      <div class="settings-card"><h2>Condition selection</h2>
        <div class="field-row"><div><label for="systemFilter">Physiological system</label><select id="systemFilter"></select></div>
        <div><label for="diseaseSearch">Search condition or Phecode</label><input id="diseaseSearch" type="search" placeholder="Example: diabetes or 250.2"></div></div>
        <div class="inline-actions"><div class="grow"><label for="diseaseSelect">Condition to add</label><select id="diseaseSelect"></select></div>
          <button id="addDiseaseButton" type="button">Add condition</button><button id="clearDiseases" type="button" class="secondary">Clear conditions</button></div>
        <div id="selectedDiseaseChips" class="chips"></div>
      </div>
      <div class="settings-card"><h2>Analyses</h2><div class="checkbox-list" id="cancerCheckboxes"></div><div class="chips" id="selectedCancerChips"></div>
        <div class="inline-actions"><button id="selectAllCancers" type="button" class="secondary">All analyses</button><button id="clearCancers" type="button" class="secondary">Clear analyses</button></div></div>
    </section>
    <section class="status-row"><span><strong id="diseaseCount">0</strong> conditions selected</span><span><strong id="cancerCount">0</strong> analyses selected</span><span>Maximum selected conditions: <strong>20</strong></span><span>Scale: <strong>{effect_measure} (log scale), 0.7 to 30</strong></span></section>
    <section class="plot-card"><div class="plot-toolbar"><button id="forestDownloadPng" type="button">PNG</button><button id="forestDownloadSvg" type="button">SVG</button><button id="forestZoomIn" type="button">+</button><button id="forestZoomOut" type="button">-</button><button id="forestResetZoom" type="button">Reset</button></div>
      <div id="plotContainer" class="plot-scroll"></div><div id="tooltip" class="tooltip"></div>
    </section>
    <p class="plot-note">Arrows indicate estimates or confidence intervals extending beyond the displayed range. Marker opacity indicates statistical significance.</p>
    """
    scripts = '<script src="../../data/phewas_data.js"></script><script src="../../assets/phewas_app.js"></script><script>window.PhewasApp.initForest();</script>'
    return _shell("Interactive PheWAS Forest Plot", "Comparison of selected condition-specific effect estimates.", body, "../..", "phewas", sections, scripts)


def _network_page(payload: Dict[str, Any], sections: list[str]) -> str:
    network = payload["networks"][0]
    if not network["nodes"]:
        content = '<section class="empty-state"><h2>No displayable comorbidity associations</h2><p>The supplied result contains no significant positive associations after Plot filtering.</p></section>'
        return _shell("Interactive Comorbidity Network", "Interactive visualization of comorbidity associations.", content, "../..", "comorbidity", sections)
    count_description = network["node_count_description"]
    body = f"""
    <section class="section-panel"><h2>Comorbidity Network Explorer</h2><p>Nodes represent modeled medical conditions and edges represent statistically significant positive comorbidity associations retained by Plot filtering.</p>
      <div class="metric-grid"><div class="metric-card"><strong>{network['stats']['node_count']}</strong><span>Network nodes</span></div><div class="metric-card"><strong>{network['stats']['edge_count']}</strong><span>Network edges</span></div><div class="metric-card"><strong>{network['stats']['module_count']}</strong><span>Detected modules</span></div></div></section>
    <section class="network-control-grid" aria-label="Comorbidity-network controls">
      <div class="settings-card"><h2>Network selection</h2><label for="networkSelect">Analysis</label><select id="networkSelect"></select>
        <div class="field-row"><div><label for="networkSystemFilter">Physiological system</label><select id="networkSystemFilter"></select></div><div><label for="networkModuleFilter">Comorbidity module</label><select id="networkModuleFilter"></select></div></div></div>
      <div class="settings-card"><h2>Search and labeling</h2><label for="networkSearch">Search condition or Phecode</label><input id="networkSearch" type="search" list="networkNodeOptions" placeholder="Example: chronic pain or 338.2"><datalist id="networkNodeOptions"></datalist>
        <div class="field-row"><div><label for="networkLabelMode">Labels</label><select id="networkLabelMode"><option value="selected">Selected and neighbors</option><option value="major">Major nodes</option><option value="all">All visible nodes</option><option value="none">No labels</option></select></div>
        <div><label for="networkNodeScale">Node size <span id="networkNodeScaleValue" class="network-range-value">100%</span></label><input id="networkNodeScale" type="range" min="0.5" max="1.75" step="0.05" value="1"></div></div></div>
      <div class="settings-card"><h2>Edge filtering</h2><label for="networkEdgeFilter"><span id="networkEdgeFilterLabel">Minimum displayed edge strength</span> <span id="networkEdgeFilterValue" class="network-range-value">0</span></label><input id="networkEdgeFilter" type="range">
        <div class="network-inline"><label for="networkHideIsolated"><input id="networkHideIsolated" type="checkbox" checked> Hide nodes without visible edges</label><button id="networkResetFilters" type="button" class="secondary">Reset filters</button></div></div>
    </section>
    <section class="status-row"><span><strong id="networkName"></strong></span><span><strong id="networkVisibleNodes">0</strong> visible nodes</span><span><strong id="networkVisibleEdges">0</strong> visible edges</span><span><strong id="networkModuleCount">0</strong> modules</span><span><strong id="networkHiddenIsolated">0</strong> hidden nodes</span><span>Selected condition: <strong id="networkSelectedNode">None</strong></span></section>
    <section class="plot-card"><div class="plot-toolbar"><button id="networkDownloadPng" type="button">PNG</button><button id="networkDownloadSvg" type="button">SVG</button><button id="networkZoomIn" type="button">+</button><button id="networkZoomOut" type="button">-</button><button id="networkResetView" type="button">Reset</button></div>
      <div class="network-viewer"><div class="network-canvas-wrap"><svg id="comorbidityNetworkSvg"><g id="comorbidityNetworkViewport"><g id="comorbidityNetworkEdges"></g><g id="comorbidityNetworkNodes"></g><g id="comorbidityNetworkLabels"></g></g></svg><div id="comorbidityNetworkTooltip" class="tooltip"></div></div>
      <aside class="network-side"><section class="network-panel"><h3>Selected Condition</h3><div id="networkDetails" class="empty-state">No condition selected.</div></section><section class="network-panel"><h3>Strongest Visible Neighbors</h3><div id="networkNeighbors" class="empty-state">No condition selected.</div></section><section class="network-panel"><h3>Physiological Systems</h3><div id="networkLegend" class="network-legend"></div></section></aside></div>
    </section>
    <p class="plot-note">Node size reflects {count_description}, and node color indicates physiological system. Edge thickness reflects the capped comorbidity effect ratio.</p>
    """
    scripts = '<script src="../../data/comorbidity_data.js"></script><script src="../../assets/network_app.js"></script><script>window.ComorbidityNetworkApp.init();</script>'
    return _shell("Interactive Comorbidity Network", "Interactive visualization of comorbidity associations.", body, "../..", "comorbidity", sections, scripts)


def _trajectory_page(payload: Dict[str, Any], sections: list[str]) -> str:
    analysis = payload["analyses"][0]
    if not analysis["modules"]:
        content = '<section class="empty-state"><h2>No displayable disease trajectories</h2><p>No significant positive within-module trajectory associations are available.</p></section>'
        return _shell("Interactive Disease Trajectory", "Interactive visualization of recorded disease ordering.", content, "../..", "trajectory", sections)
    count_description = payload["metadata"]["count_description"]
    body = f"""
    <section class="section-panel"><h2>Within-module Temporal Ordering of Recorded Diagnoses</h2><p>{payload['metadata']['interpretation']}</p>
      <div class="metric-grid"><div class="metric-card"><strong>{analysis['stats']['module_count']}</strong><span>Modules with temporal ordering</span></div><div class="metric-card"><strong>{analysis['stats']['node_count']}</strong><span>Trajectory nodes</span></div><div class="metric-card"><strong>{analysis['stats']['edge_count']}</strong><span>Directed temporal pairs</span></div></div></section>
    <section class="network-control-grid" aria-label="Disease-trajectory controls">
      <div class="settings-card"><h2>Analysis and module</h2><label for="trajectoryAnalysisSelect">Analysis</label><select id="trajectoryAnalysisSelect"></select><label for="trajectoryModuleSelect">Comorbidity module</label><select id="trajectoryModuleSelect"></select></div>
      <div class="settings-card"><h2>Search and labeling</h2><label for="trajectorySearch">Search condition or Phecode</label><input id="trajectorySearch" type="search" list="trajectoryNodeOptions" placeholder="Example: anxiety or 300.1"><datalist id="trajectoryNodeOptions"></datalist><label for="trajectoryLabelMode">Labels</label><select id="trajectoryLabelMode"><option value="all">All visible nodes</option><option value="major">Major nodes</option><option value="selected">Selected and connected</option><option value="none">No labels</option></select></div>
      <div class="settings-card"><h2>Display settings</h2><label for="trajectoryEdgeFilter">Minimum temporal-ordering OR <span id="trajectoryEdgeFilterValue" class="network-range-value">1</span></label><input id="trajectoryEdgeFilter" type="range"><label for="trajectoryNodeScale">Node size <span id="trajectoryNodeScaleValue" class="network-range-value">100%</span></label><input id="trajectoryNodeScale" type="range" min="0.45" max="1.6" step="0.05" value="1"><button id="trajectoryResetFilters" type="button" class="secondary">Reset settings</button></div>
    </section>
    <section class="status-row"><span><strong id="trajectoryAnalysisName"></strong></span><span><strong id="trajectoryModuleName"></strong></span><span><strong id="trajectoryVisibleNodes">0</strong> visible nodes</span><span><strong id="trajectoryVisibleEdges">0</strong> visible temporal pairs</span><span>Selected condition: <strong id="trajectorySelectedNode">None</strong></span></section>
    <section class="plot-card"><div class="plot-toolbar"><button id="trajectoryDownloadPng" type="button">PNG</button><button id="trajectoryDownloadSvg" type="button">SVG</button><button id="trajectoryZoomIn" type="button">+</button><button id="trajectoryZoomOut" type="button">-</button><button id="trajectoryResetView" type="button">Reset</button></div>
      <div class="network-viewer"><div class="network-canvas-wrap"><svg id="trajectoryNetworkSvg"><defs><marker id="trajectoryArrow" markerWidth="8" markerHeight="8" refX="7" refY="3" orient="auto" markerUnits="strokeWidth"><path d="M0,0 L0,6 L7,3 z" fill="#475467"></path></marker></defs><g id="trajectoryNetworkViewport"><g id="trajectoryEdges"></g><g id="trajectoryNodes"></g><g id="trajectoryLabels"></g></g></svg><div id="trajectoryTooltip" class="tooltip"></div></div>
      <aside class="network-side"><section class="network-panel"><h3>Selected Condition</h3><div id="trajectoryDetails" class="empty-state">No condition selected.</div></section><section class="network-panel"><h3>Earlier-recorded Conditions</h3><div id="trajectoryEarlier" class="empty-state">No condition selected.</div></section><section class="network-panel"><h3>Later-recorded Conditions</h3><div id="trajectoryLater" class="empty-state">No condition selected.</div></section><section class="network-panel"><h3>Physiological Systems</h3><div id="trajectoryLegend" class="network-legend"></div></section></aside></div>
    </section>
    <p class="plot-note">Arrows point from the earlier-recorded condition to the later-recorded condition. Node size reflects {count_description}; node color indicates physiological system.</p>
    """
    scripts = '<script src="../../data/trajectory_data.js"></script><script src="../../assets/trajectory_app.js"></script><script>window.TrajectoryApp.init();</script>'
    return _shell("Interactive Disease Trajectory", "Interactive visualization of recorded disease ordering.", body, "../..", "trajectory", sections, scripts)


def _three_d_page(sections: list[str]) -> str:
    body = """
    <section class="section-panel"><h2>Integrated 3D disease network</h2><p>This view combines the available comorbidity structure and temporal ordering using the existing DiNetxify 3D plot.</p></section>
    <section class="plot-card three-d-plot-card"><div class="three-d-status"><span>Interactive 3D view</span><a class="button" href="./figures/network_3d.html" target="_blank" rel="noopener">Open full size</a></div>
      <div class="three-d-frame-wrap"><iframe class="three-d-frame" src="./figures/network_3d.html" title="Interactive 3D disease network" loading="eager" allowfullscreen></iframe></div></section>
    """
    return _shell("Interactive 3D Disease Network", "Integrated visualization of comorbidity and temporal disease relationships.", body, "../..", "three-d", sections)


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, ensure_ascii=True, separators=(",", ":"), default=_json_default),
        encoding="utf-8",
    )


def _write_data_js(path: Path, variable: str, payload: Dict[str, Any]) -> None:
    encoded = json.dumps(payload, ensure_ascii=True, separators=(",", ":"), default=_json_default)
    path.write_text(f"window.{variable} = {encoded};\n", encoding="utf-8")


def _make_three_d_responsive(path: Path) -> None:
    html = path.read_text(encoding="utf-8")
    if "</head>" not in html:
        raise ValueError(f"Could not make generated 3D HTML responsive: {path}")
    responsive_style = (
        '<meta name="viewport" content="width=device-width, initial-scale=1.0" />'
        "<style>html,body{width:100%;height:100%;margin:0;overflow:hidden;background:#fff;}"
        ".plotly-graph-div{width:100%!important;height:100vh!important;min-height:760px;}</style>"
    )
    path.write_text(html.replace("</head>", f"{responsive_style}</head>", 1), encoding="utf-8")


def _copy_assets(output_dir: Path) -> None:
    asset_root = resources.files("DiNetxify").joinpath("interactive_assets")
    for name in ("style.css", "common.js", "phewas_app.js", "network_app.js", "trajectory_app.js"):
        output = output_dir / "assets" / name
        output.write_text(asset_root.joinpath(name).read_text(encoding="utf-8"), encoding="utf-8")
    (output_dir / "assets" / "plotly.min.js").write_text(get_plotlyjs(), encoding="utf-8")


def _write_readme(output_dir: Path, sections: list[str]) -> None:
    labels = {
        "phewas": "PheWAS",
        "comorbidity": "Comorbidity Network",
        "trajectory": "Disease Trajectory",
        "three-d": "3D Disease Network",
    }
    available = "\n".join(f"- {labels[section]}" for section in sections)
    temporal_note = (
        "\nTemporal arrows represent recorded diagnosis ordering and should not be interpreted as causal effects.\n"
        if "trajectory" in sections else ""
    )
    (output_dir / "README.md").write_text(
        "# DiNetxify Interactive Visualization\n\n"
        "Open `index.html` in a modern browser. The site is self-contained and does not require an internet connection.\n\n"
        "## Included sections\n\n"
        f"{available}\n"
        f"{temporal_note}",
        encoding="utf-8",
    )


def build_interactive_website(plot: "Plot", output_dir: str | Path) -> None:
    """Build a self-contained, multi-page interactive visualization website."""
    destination = Path(output_dir).expanduser()
    if destination.exists() and not destination.is_dir():
        raise ValueError(f"output_dir is not a directory: {destination}")
    destination.mkdir(parents=True, exist_ok=True)
    for generated_dir in ("assets", "data", "pages"):
        path = destination / generated_dir
        if path.exists():
            shutil.rmtree(path)
        path.mkdir(parents=True)

    phewas_payload = _build_phewas_payload(plot)
    sections = ["phewas"]
    has_comorbidity = plot._base_comorbidity is not None
    has_trajectory = has_comorbidity and plot._base_trajectory is not None

    network_payload = None
    network_graph = None
    module_by_id: Dict[str, int] = {}
    if has_comorbidity:
        network_payload, network_graph, module_by_id = _build_comorbidity_payload(plot)
        sections.append("comorbidity")

    trajectory_payload = None
    if has_trajectory:
        trajectory_payload = _build_trajectory_payload(plot, module_by_id)
        sections.append("trajectory")
        if (
            network_payload
            and network_payload["networks"][0]["nodes"]
            and trajectory_payload["analyses"][0]["modules"]
        ):
            sections.append("three-d")

    _copy_assets(destination)
    _write_json(destination / "data" / "phewas.json", phewas_payload)
    _write_data_js(destination / "data" / "phewas_data.js", "PHEWAS_DATA", phewas_payload)
    pd.DataFrame(phewas_payload["rows"]).to_csv(destination / "data" / "phewas.csv", index=False)

    if network_payload is not None and network_graph is not None:
        _write_json(destination / "data" / "comorbidity.json", network_payload)
        _write_data_js(destination / "data" / "comorbidity_data.js", "COMORBIDITY_NETWORK_DATA", network_payload)
        nx.write_gexf(network_graph, destination / "data" / "comorbidity.gexf")
    if trajectory_payload is not None:
        _write_json(destination / "data" / "trajectory.json", trajectory_payload)
        _write_data_js(destination / "data" / "trajectory_data.js", "TRAJECTORY_DATA", trajectory_payload)

    (destination / "index.html").write_text(_home_page(sections, phewas_payload), encoding="utf-8")
    phewas_dir = destination / "pages" / "phewas"
    phewas_dir.mkdir(parents=True)
    (phewas_dir / "index.html").write_text(_phewas_landing(phewas_payload, sections), encoding="utf-8")
    (phewas_dir / "scatter.html").write_text(_phewas_scatter_page(phewas_payload, sections), encoding="utf-8")
    if phewas_payload["metadata"]["has_effect_estimates"]:
        (phewas_dir / "forest.html").write_text(_phewas_forest_page(phewas_payload, sections), encoding="utf-8")

    if network_payload is not None:
        network_dir = destination / "pages" / "comorbidity-network"
        network_dir.mkdir(parents=True)
        (network_dir / "index.html").write_text(_network_page(network_payload, sections), encoding="utf-8")
    if trajectory_payload is not None:
        trajectory_dir = destination / "pages" / "disease-trajectory"
        trajectory_dir.mkdir(parents=True)
        (trajectory_dir / "index.html").write_text(_trajectory_page(trajectory_payload, sections), encoding="utf-8")
    if "three-d" in sections:
        three_d_dir = destination / "pages" / "3d-disease-network"
        figure_dir = three_d_dir / "figures"
        figure_dir.mkdir(parents=True)
        three_d_path = figure_dir / "network_3d.html"
        network = network_payload["networks"][0]
        plot.three_dimension_plot(
            str(three_d_path),
            size_reduction=network["three_d_size_reduction"],
            layout_width=None,
            layout_height=None,
        )
        _make_three_d_responsive(three_d_path)
        (three_d_dir / "index.html").write_text(_three_d_page(sections), encoding="utf-8")

    _write_readme(destination, sections)
