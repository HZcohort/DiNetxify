(function () {
  "use strict";

  const DATA = window.COMORBIDITY_NETWORK_DATA || { networks: [] };
  const NETWORKS = DATA.networks || [];
  const nodeById = new Map();
  const adjacency = new Map();
  let current = null;

  const state = {
    system: "all",
    module: "all",
    minEdge: 0,
    nodeScale: 1,
    labelMode: "major",
    hideIsolated: true,
    selected: null,
    viewBox: { x: 0, y: 0, w: 100, h: 100 }
  };

  let svg;
  let viewport;
  let edgeLayer;
  let nodeLayer;
  let labelLayer;
  let tooltip;

  function byId(id) {
    return document.getElementById(id);
  }

  function escapeHtml(value) {
    return String(value == null ? "" : value)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;")
      .replace(/'/g, "&#39;");
  }

  function fmt(value, digits) {
    if (value == null || !Number.isFinite(Number(value))) {
      return "NA";
    }
    return Number(value).toLocaleString(undefined, {
      minimumFractionDigits: digits,
      maximumFractionDigits: digits
    });
  }

  function fmtInt(value) {
    if (value == null || !Number.isFinite(Number(value))) {
      return "NA";
    }
    return Number(value).toLocaleString(undefined, { maximumFractionDigits: 0 });
  }

  function fmtP(value) {
    if (value == null || !Number.isFinite(Number(value))) {
      return "NA";
    }
    const number = Number(value);
    if (number === 0) {
      return "<1e-300";
    }
    if (Math.abs(number) < 0.001) {
      return number.toExponential(2);
    }
    return number.toPrecision(3);
  }

  function fmtHR(value) {
    if (value == null || !Number.isFinite(Number(value))) {
      return "NA";
    }
    const number = Number(value);
    if (number >= 10) {
      return number.toFixed(1);
    }
    return number.toFixed(2);
  }

  function fmtCI(node) {
    if (node.phewas_ci_lower == null || node.phewas_ci_upper == null) {
      return "95% CI unavailable";
    }
    return "95% CI: " + fmtHR(node.phewas_ci_lower) + "-" + fmtHR(node.phewas_ci_upper);
  }

  function phewasLabel(node) {
    return node.phewas_effect_measure || "Effect ratio";
  }

  function nodeCountLabel() {
    return current && current.node_count_label ? current.node_count_label : "Individuals with condition";
  }

  function displayName(node) {
    return (node.disease_name || node.phecode || node.id) + " (" + (node.phecode || node.id) + ")";
  }

  function initialViewBox(network) {
    const stats = network.stats || {};
    const width = Math.max(100, Number(stats.max_x) - Number(stats.min_x));
    const height = Math.max(100, Number(stats.max_y) - Number(stats.min_y));
    return {
      x: Number(stats.min_x) || 0,
      y: Number(stats.min_y) || 0,
      w: width,
      h: height
    };
  }

  function applyViewBox() {
    svg.setAttribute("viewBox", state.viewBox.x + " " + state.viewBox.y + " " + state.viewBox.w + " " + state.viewBox.h);
  }

  function scaledRadius(node) {
    const baseline = Number(current && current.default_node_scale) || 1;
    return Math.max(2.3, Number(node.radius || 5) * baseline * state.nodeScale);
  }

  function edgeStrength(edge) {
    return Number(edge.filter_weight == null ? edge.weight : edge.filter_weight) || 0;
  }

  function sliderDigits() {
    return current && current.is_consensus ? 2 : 1;
  }

  function resetStateForNetwork(network) {
    state.system = "all";
    state.module = "all";
    state.minEdge = Number(network.stats.min_edge_weight) || 0;
    state.nodeScale = 1;
    state.labelMode = "major";
    state.hideIsolated = true;
    state.selected = null;
    state.viewBox = initialViewBox(network);
  }

  function setupNetworkSelect() {
    const select = byId("networkSelect");
    select.innerHTML = "";
    NETWORKS.forEach((network) => {
      select.append(new Option(network.name, network.id));
    });
    select.addEventListener("change", () => {
      setCurrentNetwork(select.value);
      const url = new URL(window.location.href);
      url.searchParams.set("analysis", select.value);
      history.replaceState(null, "", url);
    });
  }

  function populateFilters() {
    const systemSelect = byId("networkSystemFilter");
    systemSelect.innerHTML = "";
    systemSelect.append(new Option("All physiological systems", "all"));
    (current.systems || []).forEach((system) => {
      systemSelect.append(new Option(system.name, system.id));
    });
    systemSelect.value = state.system;

    const moduleSelect = byId("networkModuleFilter");
    moduleSelect.innerHTML = "";
    moduleSelect.append(new Option("All modules", "all"));
    (current.modules || []).forEach((module) => {
      moduleSelect.append(new Option(module.label, String(module.display_number)));
    });
    moduleSelect.value = state.module;

    const nodeOptions = byId("networkNodeOptions");
    nodeOptions.innerHTML = "";
    (current.nodes || [])
      .slice()
      .sort((a, b) => displayName(a).localeCompare(displayName(b)))
      .forEach((node) => {
        const option = document.createElement("option");
        option.value = displayName(node);
        nodeOptions.appendChild(option);
      });

    const edgeFilter = byId("networkEdgeFilter");
    const minEdge = Number(current.stats.min_edge_weight) || 0;
    const maxEdge = Number(current.stats.max_edge_weight) || minEdge;
    const step = current.is_consensus ? 0.05 : 0.1;
    edgeFilter.min = String(minEdge);
    edgeFilter.max = String(maxEdge > minEdge ? maxEdge : minEdge + step);
    edgeFilter.step = String(step);
    edgeFilter.value = String(state.minEdge);
    byId("networkEdgeFilterLabel").textContent = current.edge_strength_label || "Minimum edge strength";
    byId("networkEdgeFilterValue").textContent = fmt(state.minEdge, sliderDigits());

    byId("networkNodeScale").value = String(state.nodeScale);
    byId("networkNodeScaleValue").textContent = Math.round(state.nodeScale * 100) + "%";
    byId("networkLabelMode").value = state.labelMode;
    byId("networkHideIsolated").checked = state.hideIsolated;
    byId("networkSearch").value = "";
  }

  function setCurrentNetwork(networkId) {
    current = NETWORKS.find((network) => network.id === networkId) || NETWORKS[0];
    if (!current) {
      return;
    }
    byId("networkSelect").value = current.id;
    resetStateForNetwork(current);
    populateFilters();
    buildGraph();
    buildLegend();
    applyViewBox();
    updateVisibility();
  }

  function clearGraph() {
    edgeLayer.innerHTML = "";
    nodeLayer.innerHTML = "";
    labelLayer.innerHTML = "";
    nodeById.clear();
    adjacency.clear();
  }

  function buildGraph() {
    clearGraph();
    (current.nodes || []).forEach((node) => {
      nodeById.set(node.id, node);
      adjacency.set(node.id, []);
    });

    (current.edges || []).forEach((edge, index) => {
      edge.id = edge.id || "edge-" + index;
      const source = nodeById.get(edge.source);
      const target = nodeById.get(edge.target);
      if (!source || !target) {
        return;
      }
      adjacency.get(edge.source).push({ id: edge.target, edge });
      adjacency.get(edge.target).push({ id: edge.source, edge });

      const line = document.createElementNS("http://www.w3.org/2000/svg", "line");
      line.setAttribute("x1", source.x);
      line.setAttribute("y1", source.y);
      line.setAttribute("x2", target.x);
      line.setAttribute("y2", target.y);
      line.setAttribute("stroke-width", Math.max(0.45, Number(edge.width || 0.8) * 1.18));
      line.classList.add("network-edge");
      line.dataset.edgeId = edge.id;
      line.addEventListener("mouseenter", (event) => showEdgeTooltip(event, edge));
      line.addEventListener("mousemove", moveTooltip);
      line.addEventListener("mouseleave", hideTooltip);
      edgeLayer.appendChild(line);
      edge.el = line;
    });

    (current.nodes || []).forEach((node) => {
      const circle = document.createElementNS("http://www.w3.org/2000/svg", "circle");
      circle.setAttribute("cx", node.x);
      circle.setAttribute("cy", node.y);
      circle.setAttribute("r", scaledRadius(node));
      circle.setAttribute("fill", node.system_color || "#999999");
      circle.classList.add("network-node");
      circle.dataset.nodeId = node.id;
      circle.addEventListener("mouseenter", (event) => showNodeTooltip(event, node));
      circle.addEventListener("mousemove", moveTooltip);
      circle.addEventListener("mouseleave", hideTooltip);
      circle.addEventListener("click", (event) => {
        event.stopPropagation();
        selectNode(node.id, false);
      });
      nodeLayer.appendChild(circle);
      node.el = circle;

      const label = document.createElementNS("http://www.w3.org/2000/svg", "text");
      label.setAttribute("x", node.x);
      label.setAttribute("y", node.y - scaledRadius(node) - 5);
      label.textContent = node.phecode || node.id;
      label.classList.add("network-label");
      labelLayer.appendChild(label);
      node.labelEl = label;
    });
  }

  function buildLegend() {
    const legend = byId("networkLegend");
    legend.innerHTML = "";
    (current.systems || []).forEach((system) => {
      const row = document.createElement("div");
      row.className = "network-legend-row";
      row.innerHTML = '<span class="network-swatch" style="background:' + escapeHtml(system.color) + '"></span><span>' + escapeHtml(system.name) + '</span>';
      row.addEventListener("click", () => {
        state.system = system.id;
        byId("networkSystemFilter").value = system.id;
        updateVisibility();
      });
      legend.appendChild(row);
    });
  }

  function nodePassesFilters(node) {
    const moduleValue = node.module_display_number == null ? "" : String(node.module_display_number);
    return (state.system === "all" || node.disease_system === state.system) &&
      (state.module === "all" || moduleValue === state.module);
  }

  function computeVisible() {
    const nodePass = new Set((current.nodes || []).filter(nodePassesFilters).map((node) => node.id));
    const visibleEdges = new Set();
    const connectedNodes = new Set();
    (current.edges || []).forEach((edge) => {
      const pass = nodePass.has(edge.source) &&
        nodePass.has(edge.target) &&
        edgeStrength(edge) >= state.minEdge;
      if (pass) {
        visibleEdges.add(edge.id);
        connectedNodes.add(edge.source);
        connectedNodes.add(edge.target);
      }
    });

    const visibleNodes = new Set();
    (current.nodes || []).forEach((node) => {
      if (!nodePass.has(node.id)) {
        return;
      }
      if (!state.hideIsolated || connectedNodes.has(node.id)) {
        visibleNodes.add(node.id);
      }
    });
    const isolatedAfterFilters = Math.max(0, nodePass.size - connectedNodes.size);
    const hiddenIsolated = state.hideIsolated ? isolatedAfterFilters : 0;
    return { visibleNodes, visibleEdges, isolatedAfterFilters, hiddenIsolated };
  }

  function updateVisibility() {
    if (!current) {
      return;
    }
    const visible = computeVisible();
    const selectedNeighbors = new Set();
    if (state.selected) {
      (adjacency.get(state.selected) || []).forEach((item) => selectedNeighbors.add(item.id));
    }

    (current.edges || []).forEach((edge) => {
      const isVisible = visible.visibleEdges.has(edge.id);
      const selectedEdge = state.selected && (edge.source === state.selected || edge.target === state.selected);
      edge.el.style.display = isVisible ? "" : "none";
      edge.el.classList.toggle("highlight", Boolean(selectedEdge));
      edge.el.classList.toggle("dimmed", Boolean(state.selected && !selectedEdge));
    });

    (current.nodes || []).forEach((node) => {
      const isVisible = visible.visibleNodes.has(node.id);
      const selectedNode = node.id === state.selected;
      const neighborNode = selectedNeighbors.has(node.id);
      const radius = scaledRadius(node);
      node.el.setAttribute("r", radius);
      node.labelEl.setAttribute("y", node.y - radius - 5);
      node.el.style.display = isVisible ? "" : "none";
      node.labelEl.style.display = isVisible ? "" : "none";
      node.el.classList.toggle("highlight", selectedNode);
      node.el.classList.toggle("dimmed", Boolean(state.selected && !selectedNode && !neighborNode));
      node.labelEl.classList.toggle("visible", shouldShowLabel(node, selectedNode, neighborNode));
    });

    renderStatus(visible);
    renderDetails(visible);
  }

  function shouldShowLabel(node, selectedNode, neighborNode) {
    if (state.labelMode === "none") {
      return false;
    }
    if (state.labelMode === "all") {
      return true;
    }
    if (state.labelMode === "selected") {
      return selectedNode || neighborNode;
    }
    return selectedNode || neighborNode || Number(node.degree || 0) >= Number(current.stats.major_degree_threshold || 0);
  }

  function renderStatus(visible) {
    byId("networkName").textContent = current.name;
    byId("networkVisibleNodes").textContent = visible.visibleNodes.size.toLocaleString();
    byId("networkVisibleEdges").textContent = visible.visibleEdges.size.toLocaleString();
    byId("networkModuleCount").textContent = Number(current.stats.module_count || 0).toLocaleString();
    byId("networkHiddenIsolated").textContent = Number(visible.hiddenIsolated || 0).toLocaleString();
    byId("networkSelectedNode").textContent = state.selected ? displayName(nodeById.get(state.selected)) : "None";
  }

  function renderDetails(visible) {
    const details = byId("networkDetails");
    const neighbors = byId("networkNeighbors");
    if (!state.selected || !nodeById.has(state.selected)) {
      details.className = "empty-state";
      details.textContent = "No condition selected.";
      neighbors.className = "empty-state";
      neighbors.textContent = "No condition selected.";
      return;
    }

    const node = nodeById.get(state.selected);
    details.className = "";
    const phewasBlock = current.is_consensus ? "" : `
      <dt>PheWAS ${escapeHtml(phewasLabel(node))}</dt><dd>${escapeHtml(fmtHR(node.phewas_hr))} (${escapeHtml(fmtCI(node))})</dd>
      <dt>p-value</dt><dd>${escapeHtml(fmtP(node.phewas_p_value))}</dd>
      <dt>q-value</dt><dd>${escapeHtml(fmtP(node.phewas_q_value))}</dd>`;
    const consensusBlock = current.is_consensus ? `
      <dt>Within-module z</dt><dd>${escapeHtml(fmt(node.within_module_zscore, 2))}</dd>` : "";
    details.innerHTML = `
      <div class="network-selected-title">${escapeHtml(node.disease_name || node.id)}</div>
      <div class="network-selected-code">Phecode ${escapeHtml(node.phecode || node.id)}</div>
      <dl class="network-kv">
        <dt>Physiological system</dt><dd>${escapeHtml(node.disease_system_display || "Unclassified")}</dd>
        <dt>Module</dt><dd>${escapeHtml(node.module_label || "Unassigned")}</dd>
        <dt>${escapeHtml(nodeCountLabel())}</dt><dd>${escapeHtml(fmtInt(node.cases))}</dd>
        ${phewasBlock}
        ${consensusBlock}
        <dt>Degree</dt><dd>${escapeHtml(fmtInt(node.degree))}</dd>
        <dt>Weighted degree (capped)</dt><dd>${escapeHtml(fmt(node.weighted_degree, 2))}</dd>
      </dl>`;

    const rows = (adjacency.get(node.id) || [])
      .filter((item) => visible.visibleNodes.has(item.id) && visible.visibleEdges.has(item.edge.id))
      .sort((a, b) => edgeStrength(b.edge) - edgeStrength(a.edge))
      .slice(0, 12);

    if (!rows.length) {
      neighbors.className = "empty-state";
      neighbors.textContent = "No visible neighbors.";
      return;
    }
    neighbors.className = "network-neighbor-list";
    neighbors.innerHTML = "";
    rows.forEach((item) => {
      const neighbor = nodeById.get(item.id);
      const row = document.createElement("div");
      row.className = "network-neighbor";
      row.innerHTML = `
        <div class="network-neighbor-name">${escapeHtml(displayName(neighbor))}</div>
        <div class="network-neighbor-weight">${escapeHtml(fmt(edgeStrength(item.edge), sliderDigits()))}</div>`;
      row.addEventListener("click", () => selectNode(neighbor.id, true));
      neighbors.appendChild(row);
    });
  }

  function selectNode(nodeId, center) {
    state.selected = nodeId;
    if (center) {
      centerNode(nodeById.get(nodeId));
    }
    updateVisibility();
  }

  function centerNode(node) {
    if (!node) {
      return;
    }
    state.viewBox.x = Number(node.x) - state.viewBox.w / 2;
    state.viewBox.y = Number(node.y) - state.viewBox.h / 2;
    applyViewBox();
  }

  function showNodeTooltip(event, node) {
    const phewasLine = current.is_consensus ? "" :
      "<div>PheWAS " + escapeHtml(phewasLabel(node)) + ": " + escapeHtml(fmtHR(node.phewas_hr)) + " (" + escapeHtml(fmtCI(node)) + ")</div>" +
      "<div>p-value: " + escapeHtml(fmtP(node.phewas_p_value)) + "; q-value: " + escapeHtml(fmtP(node.phewas_q_value)) + "</div>";
    const zLine = current.is_consensus ? "<div>Within-module z-score: " + escapeHtml(fmt(node.within_module_zscore, 2)) + "</div>" : "";
    tooltip.innerHTML = `
      <strong>${escapeHtml(displayName(node))}</strong>
      <div>Physiological system: ${escapeHtml(node.disease_system_display || "Unclassified")}</div>
      <div>Module: ${escapeHtml(node.module_label || "Unassigned")}</div>
      <div>${escapeHtml(nodeCountLabel())}: ${escapeHtml(fmtInt(node.cases))}</div>
      ${phewasLine}
      ${zLine}
      <div>Degree: ${escapeHtml(fmtInt(node.degree))}; weighted degree: ${escapeHtml(fmt(node.weighted_degree, 2))}</div>`;
    tooltip.style.display = "block";
    moveTooltip(event);
  }

  function showEdgeTooltip(event, edge) {
    const source = nodeById.get(edge.source);
    const target = nodeById.get(edge.target);
    if (current.is_consensus) {
      tooltip.innerHTML = `
        <strong>${escapeHtml(source.phecode)} - ${escapeHtml(target.phecode)}</strong>
        <div>${escapeHtml(source.disease_name || source.id)}</div>
        <div>${escapeHtml(target.disease_name || target.id)}</div>
        <div>Recurrence proportion: ${escapeHtml(fmt(edge.recurrence_proportion, 2))}</div>
        <div>Approximate supporting networks: ${escapeHtml(fmtInt(edge.support_count))} of ${escapeHtml(fmtInt(edge.support_denominator))}</div>`;
    } else {
      tooltip.innerHTML = `
        <strong>${escapeHtml(source.phecode)} - ${escapeHtml(target.phecode)}</strong>
        <div>${escapeHtml(source.disease_name || source.id)}</div>
        <div>${escapeHtml(target.disease_name || target.id)}</div>
        <div>Comorbidity effect ratio: ${escapeHtml(fmtHR(edge.or || edge.hr))}</div>
        <div>95% CI: ${escapeHtml(fmtHR(edge.ci_lower))}-${escapeHtml(fmtHR(edge.ci_upper))}</div>
        <div>p-value: ${escapeHtml(fmtP(edge.p_value))}; q-value: ${escapeHtml(fmtP(edge.q_value))}</div>
        <div>Individuals in edge model: ${escapeHtml(fmtInt(edge.n_total))}</div>`;
    }
    tooltip.style.display = "block";
    moveTooltip(event);
  }

  function moveTooltip(event) {
    const pad = 14;
    const rect = tooltip.getBoundingClientRect();
    let left = event.clientX + pad;
    let top = event.clientY + pad;
    if (left + rect.width > window.innerWidth - 8) {
      left = event.clientX - rect.width - pad;
    }
    if (top + rect.height > window.innerHeight - 8) {
      top = event.clientY - rect.height - pad;
    }
    tooltip.style.left = left + "px";
    tooltip.style.top = top + "px";
  }

  function hideTooltip() {
    tooltip.style.display = "none";
  }

  function clientToSvg(clientX, clientY) {
    const rect = svg.getBoundingClientRect();
    return {
      x: state.viewBox.x + (clientX - rect.left) * state.viewBox.w / rect.width,
      y: state.viewBox.y + (clientY - rect.top) * state.viewBox.h / rect.height
    };
  }

  function zoomBy(factor) {
    const center = {
      x: state.viewBox.x + state.viewBox.w / 2,
      y: state.viewBox.y + state.viewBox.h / 2
    };
    state.viewBox.x = center.x - (state.viewBox.w * factor) / 2;
    state.viewBox.y = center.y - (state.viewBox.h * factor) / 2;
    state.viewBox.w *= factor;
    state.viewBox.h *= factor;
    applyViewBox();
  }

  function setupPanZoom() {
    svg.addEventListener("click", (event) => {
      if (event.target === svg) {
        state.selected = null;
        updateVisibility();
      }
    });

    svg.addEventListener("wheel", (event) => {
      event.preventDefault();
      const point = clientToSvg(event.clientX, event.clientY);
      const factor = event.deltaY > 0 ? 1.12 : 0.88;
      state.viewBox.x = point.x - (point.x - state.viewBox.x) * factor;
      state.viewBox.y = point.y - (point.y - state.viewBox.y) * factor;
      state.viewBox.w *= factor;
      state.viewBox.h *= factor;
      applyViewBox();
    }, { passive: false });

    let dragStart = null;
    svg.addEventListener("pointerdown", (event) => {
      if (event.target !== svg) {
        return;
      }
      svg.setPointerCapture(event.pointerId);
      svg.classList.add("dragging");
      dragStart = {
        x: event.clientX,
        y: event.clientY,
        viewBox: { ...state.viewBox }
      };
    });
    svg.addEventListener("pointermove", (event) => {
      if (!dragStart) {
        return;
      }
      const rect = svg.getBoundingClientRect();
      const dx = (event.clientX - dragStart.x) * dragStart.viewBox.w / rect.width;
      const dy = (event.clientY - dragStart.y) * dragStart.viewBox.h / rect.height;
      state.viewBox.x = dragStart.viewBox.x - dx;
      state.viewBox.y = dragStart.viewBox.y - dy;
      applyViewBox();
    });
    svg.addEventListener("pointerup", (event) => {
      dragStart = null;
      svg.classList.remove("dragging");
      if (svg.hasPointerCapture(event.pointerId)) {
        svg.releasePointerCapture(event.pointerId);
      }
    });
  }

  function findSearchMatch(query) {
    const normalized = query.trim().toLowerCase();
    if (!normalized) {
      return null;
    }
    return (current.nodes || []).find((node) =>
      String(node.phecode || node.id).toLowerCase() === normalized ||
      displayName(node).toLowerCase() === normalized ||
      displayName(node).toLowerCase().includes(normalized)
    ) || null;
  }

  function resetFilters() {
    resetStateForNetwork(current);
    populateFilters();
    applyViewBox();
    updateVisibility();
  }

  function downloadBlob(blob, filename) {
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    link.remove();
    URL.revokeObjectURL(url);
  }

  function exportSvgText() {
    const clone = svg.cloneNode(true);
    clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
    clone.setAttribute("width", "1600");
    clone.setAttribute("height", "1050");
    const style = document.createElementNS("http://www.w3.org/2000/svg", "style");
    style.textContent = `
      .network-edge{stroke:#1f2937;stroke-opacity:.22;vector-effect:non-scaling-stroke}
      .network-edge.dimmed{stroke-opacity:.04}.network-edge.highlight{stroke:#111827;stroke-opacity:.76}
      .network-node{stroke:#fff;stroke-width:1.35;vector-effect:non-scaling-stroke}
      .network-node.dimmed{opacity:.16}.network-node.highlight{stroke:#111827;stroke-width:2.4}
      .network-label{display:none;fill:#263341;stroke:#fff;stroke-width:2.2px;paint-order:stroke;font:700 12px Arial;text-anchor:middle}
      .network-label.visible{display:block}`;
    clone.insertBefore(style, clone.firstChild);
    return new XMLSerializer().serializeToString(clone);
  }

  function downloadSvg() {
    downloadBlob(new Blob([exportSvgText()], { type: "image/svg+xml;charset=utf-8" }), "comorbidity_network.svg");
  }

  function downloadPng() {
    const text = exportSvgText();
    const image = new Image();
    const blob = new Blob([text], { type: "image/svg+xml;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    image.onload = () => {
      const canvas = document.createElement("canvas");
      canvas.width = 1600;
      canvas.height = 1050;
      const context = canvas.getContext("2d");
      context.fillStyle = "#ffffff";
      context.fillRect(0, 0, canvas.width, canvas.height);
      context.drawImage(image, 0, 0, canvas.width, canvas.height);
      URL.revokeObjectURL(url);
      canvas.toBlob((pngBlob) => {
        if (pngBlob) {
          downloadBlob(pngBlob, "comorbidity_network.png");
        }
      });
    };
    image.src = url;
  }

  function setupControls() {
    setupNetworkSelect();
    byId("networkSystemFilter").addEventListener("change", (event) => {
      state.system = event.target.value;
      state.selected = null;
      updateVisibility();
    });
    byId("networkModuleFilter").addEventListener("change", (event) => {
      state.module = event.target.value;
      state.selected = null;
      updateVisibility();
    });
    byId("networkEdgeFilter").addEventListener("input", (event) => {
      state.minEdge = Number(event.target.value);
      byId("networkEdgeFilterValue").textContent = fmt(state.minEdge, sliderDigits());
      updateVisibility();
    });
    byId("networkNodeScale").addEventListener("input", (event) => {
      state.nodeScale = Number(event.target.value);
      byId("networkNodeScaleValue").textContent = Math.round(state.nodeScale * 100) + "%";
      updateVisibility();
    });
    byId("networkLabelMode").addEventListener("change", (event) => {
      state.labelMode = event.target.value;
      updateVisibility();
    });
    byId("networkHideIsolated").addEventListener("change", (event) => {
      state.hideIsolated = event.target.checked;
      updateVisibility();
    });
    byId("networkSearch").addEventListener("change", (event) => {
      const searchValue = event.target.value;
      const match = findSearchMatch(searchValue);
      if (match) {
        resetFilters();
        byId("networkSearch").value = searchValue;
        selectNode(match.id, true);
      }
    });
    byId("networkResetFilters").addEventListener("click", resetFilters);
    byId("networkDownloadSvg").addEventListener("click", downloadSvg);
    byId("networkDownloadPng").addEventListener("click", downloadPng);
    byId("networkZoomIn").addEventListener("click", () => zoomBy(0.82));
    byId("networkZoomOut").addEventListener("click", () => zoomBy(1.18));
    byId("networkResetView").addEventListener("click", () => {
      state.viewBox = initialViewBox(current);
      applyViewBox();
    });
  }

  function init() {
    svg = byId("comorbidityNetworkSvg");
    viewport = byId("comorbidityNetworkViewport");
    edgeLayer = byId("comorbidityNetworkEdges");
    nodeLayer = byId("comorbidityNetworkNodes");
    labelLayer = byId("comorbidityNetworkLabels");
    tooltip = byId("comorbidityNetworkTooltip");
    if (!svg || !viewport || !edgeLayer || !nodeLayer || !labelLayer || !NETWORKS.length) {
      return;
    }
    setupControls();
    setupPanZoom();
    const requested = new URLSearchParams(window.location.search).get("analysis");
    const initial = NETWORKS.some((network) => network.id === requested) ? requested : NETWORKS[0].id;
    setCurrentNetwork(initial);
  }

  window.ComorbidityNetworkApp = { init };
})();
