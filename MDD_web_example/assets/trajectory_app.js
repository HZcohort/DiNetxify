(function () {
  "use strict";

  const DATA = window.TRAJECTORY_DATA || { analyses: [] };
  const ANALYSES = DATA.analyses || [];
  const COUNT_LABEL = DATA.metadata && DATA.metadata.count_label ? DATA.metadata.count_label : "Individuals with condition";
  const state = {
    minEdge: 1,
    nodeScale: 1,
    labelMode: "all",
    selected: null,
    viewBox: { x: 0, y: 0, w: 1000, h: 700 }
  };

  let currentAnalysis = null;
  let currentModule = null;
  let svg;
  let viewport;
  let edgeLayer;
  let nodeLayer;
  let labelLayer;
  let tooltip;
  const nodeById = new Map();
  const incoming = new Map();
  const outgoing = new Map();

  function byId(id) { return document.getElementById(id); }

  function escapeHtml(value) {
    return String(value == null ? "" : value)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;")
      .replace(/'/g, "&#39;");
  }

  function fmt(value, digits) {
    if (value == null || !Number.isFinite(Number(value))) return "NA";
    return Number(value).toLocaleString(undefined, { minimumFractionDigits: digits, maximumFractionDigits: digits });
  }

  function fmtInt(value) {
    if (value == null || !Number.isFinite(Number(value))) return "NA";
    return Number(value).toLocaleString(undefined, { maximumFractionDigits: 0 });
  }

  function fmtP(value) {
    if (value == null || !Number.isFinite(Number(value))) return "NA";
    const number = Number(value);
    if (number === 0) return "<1e-300";
    if (Math.abs(number) < .001) return number.toExponential(2);
    return number.toPrecision(3);
  }

  function fmtEffect(value) {
    if (value == null || !Number.isFinite(Number(value))) return "NA";
    const number = Number(value);
    if (number >= 100) return number.toFixed(0);
    if (number >= 10) return number.toFixed(1);
    return number.toFixed(2);
  }

  function displayName(node) {
    if (node.is_endpoint) return node.disease_name || node.endpoint_type || "Endpoint";
    return (node.disease_name || node.phecode) + " (" + node.phecode + ")";
  }

  function shortLabel(node) {
    if (node.is_endpoint) return node.disease_name || node.endpoint_type || "Endpoint";
    const name = String(node.disease_name || node.phecode);
    return (name.length > 28 ? name.slice(0, 26) + "…" : name) + " (" + node.phecode + ")";
  }

  function updateUrl() {
    if (!currentAnalysis || !currentModule) return;
    const url = new URL(window.location.href);
    url.searchParams.set("analysis", currentAnalysis.id);
    url.searchParams.set("module", String(currentModule.display_number));
    history.replaceState(null, "", url);
  }

  function applyViewBox() {
    svg.setAttribute("viewBox", state.viewBox.x + " " + state.viewBox.y + " " + state.viewBox.w + " " + state.viewBox.h);
  }

  function resetView() {
    state.viewBox = { x: 0, y: 0, w: 1000, h: 700 };
    applyViewBox();
  }

  function setupAnalysisSelect() {
    const select = byId("trajectoryAnalysisSelect");
    select.innerHTML = "";
    ANALYSES.forEach((analysis) => select.append(new Option(analysis.name, analysis.id)));
    select.addEventListener("change", () => setAnalysis(select.value, null));
  }

  function populateModuleSelect(requestedModule) {
    const select = byId("trajectoryModuleSelect");
    select.innerHTML = "";
    (currentAnalysis.modules || []).forEach((module) => {
      select.append(new Option(module.label + " (" + module.stats.edge_count + " temporal pairs)", String(module.display_number)));
    });
    const requested = Number(requestedModule);
    const selected = (currentAnalysis.modules || []).find((module) => module.display_number === requested) || currentAnalysis.modules[0];
    if (selected) select.value = String(selected.display_number);
  }

  function setAnalysis(analysisId, requestedModule) {
    currentAnalysis = ANALYSES.find((analysis) => analysis.id === analysisId) || ANALYSES[0];
    if (!currentAnalysis) return;
    byId("trajectoryAnalysisSelect").value = currentAnalysis.id;
    populateModuleSelect(requestedModule);
    const moduleNumber = Number(byId("trajectoryModuleSelect").value);
    setModule(moduleNumber);
  }

  function setModule(displayNumber) {
    currentModule = (currentAnalysis.modules || []).find((module) => module.display_number === Number(displayNumber)) || currentAnalysis.modules[0];
    if (!currentModule) return;
    byId("trajectoryModuleSelect").value = String(currentModule.display_number);
    state.selected = null;
    state.minEdge = Number(currentModule.stats.min_edge_weight || 1);
    state.nodeScale = 1;
    state.labelMode = "all";
    resetView();
    populateModuleControls();
    buildGraph();
    renderLegend();
    updateVisibility();
    updateUrl();
  }

  function populateModuleControls() {
    const edge = byId("trajectoryEdgeFilter");
    edge.min = String(currentModule.stats.min_edge_weight || 1);
    edge.max = String(currentModule.stats.max_edge_weight || 10);
    edge.step = ".1";
    edge.value = String(state.minEdge);
    byId("trajectoryEdgeFilterValue").textContent = fmt(state.minEdge, 1);
    byId("trajectoryNodeScale").value = "1";
    byId("trajectoryNodeScaleValue").textContent = "100%";
    byId("trajectoryLabelMode").value = state.labelMode;
    byId("trajectorySearch").value = "";
    const options = byId("trajectoryNodeOptions");
    options.innerHTML = "";
    currentModule.nodes.slice().sort((a, b) => displayName(a).localeCompare(displayName(b))).forEach((node) => {
      const option = document.createElement("option");
      option.value = displayName(node);
      options.appendChild(option);
    });
  }

  function clearGraph() {
    edgeLayer.innerHTML = "";
    nodeLayer.innerHTML = "";
    labelLayer.innerHTML = "";
    nodeById.clear();
    incoming.clear();
    outgoing.clear();
  }

  function pathFor(source, target) {
    const direction = Number(target.x) >= Number(source.x) ? 1 : -1;
    const sourceX = Number(source.x) + direction * (scaledRadius(source) + 2);
    const targetX = Number(target.x) - direction * (scaledRadius(target) + 4);
    const middle = (sourceX + targetX) / 2;
    return "M " + sourceX + " " + source.y + " C " + middle + " " + source.y + ", " + middle + " " + target.y + ", " + targetX + " " + target.y;
  }

  function buildGraph() {
    clearGraph();
    currentModule.nodes.forEach((node) => {
      nodeById.set(node.id, node);
      incoming.set(node.id, []);
      outgoing.set(node.id, []);
    });
    currentModule.edges.forEach((edge) => {
      const source = nodeById.get(edge.source);
      const target = nodeById.get(edge.target);
      if (!source || !target) return;
      outgoing.get(source.id).push({ id: target.id, edge: edge });
      incoming.get(target.id).push({ id: source.id, edge: edge });
      const path = document.createElementNS("http://www.w3.org/2000/svg", "path");
      path.setAttribute("d", pathFor(source, target));
      path.setAttribute("class", "trajectory-edge");
      path.setAttribute("data-edge-id", edge.id);
      path.setAttribute("marker-end", "url(#trajectoryArrow)");
      path.setAttribute("stroke-width", String(.8 + Math.log(Number(edge.filter_weight || 1)) * .85));
      path.addEventListener("mousemove", (event) => showEdgeTooltip(event, edge));
      path.addEventListener("mouseleave", hideTooltip);
      edgeLayer.appendChild(path);
    });
    currentModule.nodes.forEach((node) => {
      const circle = document.createElementNS("http://www.w3.org/2000/svg", "circle");
      circle.setAttribute("cx", node.x);
      circle.setAttribute("cy", node.y);
      circle.setAttribute("r", scaledRadius(node));
      circle.setAttribute("fill", node.system_color || "#999999");
      circle.setAttribute("class", "trajectory-node");
      circle.setAttribute("data-node-id", node.id);
      circle.addEventListener("click", () => selectNode(node.id, false));
      circle.addEventListener("mousemove", (event) => showNodeTooltip(event, node));
      circle.addEventListener("mouseleave", hideTooltip);
      nodeLayer.appendChild(circle);

      const label = document.createElementNS("http://www.w3.org/2000/svg", "text");
      label.setAttribute("x", node.x);
      label.setAttribute("y", Number(node.y) + scaledRadius(node) + 14);
      label.setAttribute("class", "trajectory-label");
      label.setAttribute("data-label-id", node.id);
      label.textContent = shortLabel(node);
      labelLayer.appendChild(label);
    });
  }

  function scaledRadius(node) {
    return Math.max(3.5, Number(node.radius || 7) * state.nodeScale);
  }

  function visibleEdges() {
    return new Set(currentModule.edges.filter((edge) => edge.is_endpoint_edge || Number(edge.filter_weight || 0) >= state.minEdge).map((edge) => edge.id));
  }

  function connectedNodes(edgeIds) {
    const nodes = new Set();
    currentModule.edges.forEach((edge) => {
      if (edgeIds.has(edge.id)) {
        nodes.add(edge.source);
        nodes.add(edge.target);
      }
    });
    return nodes;
  }

  function isAdjacent(nodeId) {
    if (!state.selected) return false;
    return (incoming.get(state.selected) || []).some((item) => item.id === nodeId) ||
      (outgoing.get(state.selected) || []).some((item) => item.id === nodeId);
  }

  function labelVisible(node, connected) {
    if (!connected.has(node.id)) return false;
    if (state.labelMode === "none") return false;
    if (state.labelMode === "all") return true;
    if (state.labelMode === "selected") return node.id === state.selected || isAdjacent(node.id);
    return node.id === state.selected || isAdjacent(node.id) || Number(node.in_degree + node.out_degree) >= 2;
  }

  function updateVisibility() {
    if (!currentModule) return;
    const edgeIds = visibleEdges();
    const connected = connectedNodes(edgeIds);
    edgeLayer.querySelectorAll(".trajectory-edge").forEach((element) => {
      const edge = currentModule.edges.find((item) => item.id === element.dataset.edgeId);
      const visible = edge && edgeIds.has(edge.id);
      element.style.display = visible ? "" : "none";
      if (edge) {
        const source = nodeById.get(edge.source);
        const target = nodeById.get(edge.target);
        if (source && target) element.setAttribute("d", pathFor(source, target));
      }
      const highlight = visible && state.selected && (edge.source === state.selected || edge.target === state.selected);
      element.classList.toggle("highlight", Boolean(highlight));
      element.classList.toggle("dimmed", Boolean(state.selected && !highlight));
    });
    nodeLayer.querySelectorAll(".trajectory-node").forEach((element) => {
      const nodeId = element.dataset.nodeId;
      const visible = connected.has(nodeId);
      element.style.display = visible ? "" : "none";
      const highlight = nodeId === state.selected || isAdjacent(nodeId);
      element.classList.toggle("highlight", Boolean(state.selected && highlight));
      element.classList.toggle("dimmed", Boolean(state.selected && !highlight));
      const node = nodeById.get(nodeId);
      if (node) element.setAttribute("r", scaledRadius(node));
    });
    labelLayer.querySelectorAll(".trajectory-label").forEach((element) => {
      const node = nodeById.get(element.dataset.labelId);
      const visible = node && labelVisible(node, connected);
      element.classList.toggle("visible", Boolean(visible));
      if (node) element.setAttribute("y", Number(node.y) + scaledRadius(node) + 14);
    });
    renderStatus(edgeIds, connected);
    renderDetails(edgeIds);
  }

  function renderStatus(edgeIds, connected) {
    byId("trajectoryAnalysisName").textContent = currentAnalysis.name;
    byId("trajectoryModuleName").textContent = currentModule.label;
    byId("trajectoryVisibleNodes").textContent = connected.size.toLocaleString();
    byId("trajectoryVisibleEdges").textContent = edgeIds.size.toLocaleString();
    byId("trajectorySelectedNode").textContent = state.selected ? displayName(nodeById.get(state.selected)) : "None";
  }

  function listDirection(container, rows, edgeIds) {
    container.innerHTML = "";
    const visible = rows.filter((item) => edgeIds.has(item.edge.id)).sort((a, b) => Number(b.edge.or) - Number(a.edge.or));
    if (!visible.length) {
      container.className = "empty-state";
      container.textContent = "No visible conditions.";
      return;
    }
    container.className = "direction-list";
    visible.forEach((item) => {
      const node = nodeById.get(item.id);
      const row = document.createElement("div");
      row.className = "direction-row";
      row.innerHTML = '<div class="direction-name">' + escapeHtml(displayName(node)) + '</div><div class="direction-weight">OR ' + escapeHtml(fmtEffect(item.edge.or)) + '</div>';
      row.addEventListener("click", () => selectNode(node.id, true));
      container.appendChild(row);
    });
  }

  function renderDetails(edgeIds) {
    const details = byId("trajectoryDetails");
    const earlier = byId("trajectoryEarlier");
    const later = byId("trajectoryLater");
    if (!state.selected || !nodeById.has(state.selected)) {
      details.className = "empty-state";
      details.textContent = "No condition selected.";
      earlier.className = "empty-state";
      earlier.textContent = "No condition selected.";
      later.className = "empty-state";
      later.textContent = "No condition selected.";
      return;
    }
    const node = nodeById.get(state.selected);
    details.className = "";
    details.innerHTML = '<div class="network-selected-title">' + escapeHtml(node.disease_name) + '</div>' +
      '<div class="network-selected-code">Phecode ' + escapeHtml(node.phecode) + '</div>' +
      '<dl class="network-kv"><dt>Physiological system</dt><dd>' + escapeHtml(node.disease_system_display) + '</dd>' +
      '<dt>Module</dt><dd>' + escapeHtml(node.module_label) + '</dd>' +
      '<dt>' + escapeHtml(COUNT_LABEL) + '</dt><dd>' + escapeHtml(fmtInt(node.cases)) + '</dd>' +
      '<dt>Earlier-recorded links</dt><dd>' + escapeHtml(fmtInt(node.in_degree)) + '</dd>' +
      '<dt>Later-recorded links</dt><dd>' + escapeHtml(fmtInt(node.out_degree)) + '</dd></dl>';
    listDirection(earlier, incoming.get(node.id) || [], edgeIds);
    listDirection(later, outgoing.get(node.id) || [], edgeIds);
  }

  function selectNode(nodeId, center) {
    state.selected = nodeId;
    if (center) centerNode(nodeById.get(nodeId));
    updateVisibility();
  }

  function centerNode(node) {
    if (!node) return;
    state.viewBox.x = Number(node.x) - state.viewBox.w / 2;
    state.viewBox.y = Number(node.y) - state.viewBox.h / 2;
    applyViewBox();
  }

  function showNodeTooltip(event, node) {
    tooltip.innerHTML = '<strong>' + escapeHtml(displayName(node)) + '</strong>' +
      '<div>Physiological system: ' + escapeHtml(node.disease_system_display) + '</div>' +
      '<div>' + escapeHtml(node.module_label) + '</div>' +
      '<div>' + escapeHtml(COUNT_LABEL) + ': ' + escapeHtml(fmtInt(node.cases)) + '</div>';
    tooltip.style.display = "block";
    moveTooltip(event);
  }

  function showEdgeTooltip(event, edge) {
    const source = nodeById.get(edge.source);
    const target = nodeById.get(edge.target);
    if (edge.is_endpoint_edge) {
      tooltip.innerHTML = '<strong>Analysis endpoint link</strong>' +
        '<div>' + escapeHtml(displayName(source)) + '</div><div>→ ' + escapeHtml(displayName(target)) + '</div>';
      tooltip.style.display = "block";
      moveTooltip(event);
      return;
    }
    const countLine = edge.earlier_later_count == null ? "" :
      '<div>Recorded order counts: ' + escapeHtml(fmtInt(edge.earlier_later_count)) + ' forward, ' + escapeHtml(fmtInt(edge.reverse_count)) + ' reverse</div>';
    tooltip.innerHTML = '<strong>Earlier recorded → later recorded</strong>' +
      '<div>' + escapeHtml(displayName(source)) + '</div><div>→ ' + escapeHtml(displayName(target)) + '</div>' +
      '<div>Adjusted temporal-ordering OR: ' + escapeHtml(fmtEffect(edge.or)) + '</div>' +
      '<div>95% CI: ' + escapeHtml(fmtEffect(edge.ci_lower)) + '–' + escapeHtml(fmtEffect(edge.ci_upper)) + '</div>' +
      '<div>p-value: ' + escapeHtml(fmtP(edge.p_value)) + '; q-value: ' + escapeHtml(fmtP(edge.q_value)) + '</div>' + countLine;
    tooltip.style.display = "block";
    moveTooltip(event);
  }

  function moveTooltip(event) {
    const pad = 14;
    const rect = tooltip.getBoundingClientRect();
    let left = event.clientX + pad;
    let top = event.clientY + pad;
    if (left + rect.width > window.innerWidth - 8) left = event.clientX - rect.width - pad;
    if (top + rect.height > window.innerHeight - 8) top = event.clientY - rect.height - pad;
    tooltip.style.left = left + "px";
    tooltip.style.top = top + "px";
  }

  function hideTooltip() { tooltip.style.display = "none"; }

  function renderLegend() {
    const legend = byId("trajectoryLegend");
    legend.innerHTML = "";
    currentModule.systems.forEach((system) => {
      const row = document.createElement("div");
      row.className = "network-legend-row";
      row.innerHTML = '<span class="network-swatch" style="background:' + escapeHtml(system.color) + '"></span><span>' + escapeHtml(system.name) + '</span>';
      legend.appendChild(row);
    });
  }

  function findSearchMatch(value) {
    const query = String(value || "").trim().toLowerCase();
    if (!query) return null;
    return currentModule.nodes.find((node) => displayName(node).toLowerCase() === query) ||
      currentModule.nodes.find((node) => displayName(node).toLowerCase().includes(query));
  }

  function zoomBy(factor) {
    const centerX = state.viewBox.x + state.viewBox.w / 2;
    const centerY = state.viewBox.y + state.viewBox.h / 2;
    state.viewBox.w *= factor;
    state.viewBox.h *= factor;
    state.viewBox.x = centerX - state.viewBox.w / 2;
    state.viewBox.y = centerY - state.viewBox.h / 2;
    applyViewBox();
  }

  function setupPanZoom() {
    let dragging = false;
    let start = null;
    svg.addEventListener("mousedown", (event) => {
      dragging = true;
      start = { x: event.clientX, y: event.clientY, viewX: state.viewBox.x, viewY: state.viewBox.y };
      svg.classList.add("dragging");
    });
    window.addEventListener("mousemove", (event) => {
      if (!dragging || !start) return;
      const rect = svg.getBoundingClientRect();
      state.viewBox.x = start.viewX - (event.clientX - start.x) * state.viewBox.w / rect.width;
      state.viewBox.y = start.viewY - (event.clientY - start.y) * state.viewBox.h / rect.height;
      applyViewBox();
    });
    window.addEventListener("mouseup", () => { dragging = false; start = null; svg.classList.remove("dragging"); });
    svg.addEventListener("wheel", (event) => { event.preventDefault(); zoomBy(event.deltaY > 0 ? 1.12 : .89); }, { passive: false });
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
    style.textContent = '.trajectory-edge{fill:none;stroke:#475467;stroke-opacity:.4}.trajectory-node{stroke:#fff;stroke-width:1.5}.trajectory-label{display:none;fill:#263341;stroke:#fff;stroke-width:3;paint-order:stroke;font:700 11px Arial;text-anchor:middle}.trajectory-label.visible{display:block}';
    clone.insertBefore(style, clone.firstChild);
    return new XMLSerializer().serializeToString(clone);
  }

  function downloadSvg() {
    downloadBlob(new Blob([exportSvgText()], { type: "image/svg+xml;charset=utf-8" }), "within_module_temporal_ordering.svg");
  }

  function downloadPng() {
    const text = exportSvgText();
    const image = new Image();
    const url = URL.createObjectURL(new Blob([text], { type: "image/svg+xml;charset=utf-8" }));
    image.onload = () => {
      const canvas = document.createElement("canvas");
      canvas.width = 1600;
      canvas.height = 1050;
      const context = canvas.getContext("2d");
      context.fillStyle = "#ffffff";
      context.fillRect(0, 0, canvas.width, canvas.height);
      context.drawImage(image, 0, 0, canvas.width, canvas.height);
      URL.revokeObjectURL(url);
      canvas.toBlob((blob) => { if (blob) downloadBlob(blob, "within_module_temporal_ordering.png"); }, "image/png");
    };
    image.src = url;
  }

  function setupControls() {
    setupAnalysisSelect();
    byId("trajectoryModuleSelect").addEventListener("change", (event) => setModule(Number(event.target.value)));
    byId("trajectoryEdgeFilter").addEventListener("input", (event) => {
      state.minEdge = Number(event.target.value);
      byId("trajectoryEdgeFilterValue").textContent = fmt(state.minEdge, 1);
      updateVisibility();
    });
    byId("trajectoryNodeScale").addEventListener("input", (event) => {
      state.nodeScale = Number(event.target.value);
      byId("trajectoryNodeScaleValue").textContent = Math.round(state.nodeScale * 100) + "%";
      updateVisibility();
    });
    byId("trajectoryLabelMode").addEventListener("change", (event) => { state.labelMode = event.target.value; updateVisibility(); });
    byId("trajectorySearch").addEventListener("change", (event) => {
      const match = findSearchMatch(event.target.value);
      if (match) {
        state.minEdge = Number(currentModule.stats.min_edge_weight || 1);
        byId("trajectoryEdgeFilter").value = String(state.minEdge);
        byId("trajectoryEdgeFilterValue").textContent = fmt(state.minEdge, 1);
        selectNode(match.id, true);
      }
    });
    byId("trajectoryResetFilters").addEventListener("click", () => setModule(currentModule.display_number));
    byId("trajectoryDownloadSvg").addEventListener("click", downloadSvg);
    byId("trajectoryDownloadPng").addEventListener("click", downloadPng);
    byId("trajectoryZoomIn").addEventListener("click", () => zoomBy(.82));
    byId("trajectoryZoomOut").addEventListener("click", () => zoomBy(1.18));
    byId("trajectoryResetView").addEventListener("click", resetView);
  }

  function init() {
    svg = byId("trajectoryNetworkSvg");
    viewport = byId("trajectoryNetworkViewport");
    edgeLayer = byId("trajectoryEdges");
    nodeLayer = byId("trajectoryNodes");
    labelLayer = byId("trajectoryLabels");
    tooltip = byId("trajectoryTooltip");
    if (!svg || !viewport || !edgeLayer || !nodeLayer || !labelLayer || !ANALYSES.length) return;
    setupControls();
    setupPanZoom();
    const params = new URLSearchParams(window.location.search);
    const requestedAnalysis = params.get("analysis");
    const requestedModule = params.get("module");
    const initial = ANALYSES.some((analysis) => analysis.id === requestedAnalysis) ? requestedAnalysis : ANALYSES[0].id;
    setAnalysis(initial, requestedModule);
  }

  window.TrajectoryApp = { init: init };
})();
