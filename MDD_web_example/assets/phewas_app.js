(function () {
  "use strict";

  const DATA = window.PHEWAS_DATA || { rows: [], cancers: [], systems: [] };
  const ROWS = DATA.rows || [];
  const CANCERS = DATA.cancers || [];
  const SYSTEMS = DATA.systems || [];
  const SYSTEM_ORDER = new Map(SYSTEMS.map((system, index) => [system.id, index]));
  const SCATTER_HR_CAP = 30;
  const MAX_DISEASES = 20;
  const FOREST_MIN = 0.7;
  const FOREST_MAX = 30;
  const FOREST_TICKS = [1, 2, 5, 10, 20, 30];
  const EFFECT_LABEL = DATA.metadata && DATA.metadata.effect_measure ? DATA.metadata.effect_measure : "Effect ratio";

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

  function formatNumber(value, digits) {
    if (value == null || !Number.isFinite(Number(value))) {
      return "NA";
    }
    return Number(value).toFixed(digits);
  }

  function formatP(value) {
    if (value == null || !Number.isFinite(Number(value))) {
      return "NA";
    }
    const numeric = Number(value);
    if (numeric === 0) {
      return "<1e-300";
    }
    if (Math.abs(numeric) < 0.001) {
      return numeric.toExponential(2);
    }
    return numeric.toPrecision(3);
  }

  function formatHR(value) {
    if (value == null || !Number.isFinite(Number(value))) {
      return "NA";
    }
    const numeric = Number(value);
    if (numeric > 100) {
      return ">100";
    }
    if (numeric >= 10) {
      return numeric.toFixed(1);
    }
    return numeric.toFixed(2);
  }

  function formatCI(row) {
    if (row.ci_lower == null || row.ci_upper == null) {
      return "95% CI unavailable";
    }
    return "95% CI: " + formatHR(row.ci_lower) + "-" + formatHR(row.ci_upper);
  }

  function formatRate(value) {
    if (value == null || !Number.isFinite(Number(value))) {
      return "NA";
    }
    return Number(value).toFixed(2) + " per 1,000 person-years";
  }

  function hexToRgba(hex, alpha) {
    const clean = String(hex || "#999999").replace("#", "");
    const r = parseInt(clean.slice(0, 2), 16);
    const g = parseInt(clean.slice(2, 4), 16);
    const b = parseInt(clean.slice(4, 6), 16);
    return "rgba(" + r + "," + g + "," + b + "," + alpha + ")";
  }

  function compareRows(a, b) {
    const systemDelta = (SYSTEM_ORDER.get(a.disease_system) ?? 999) - (SYSTEM_ORDER.get(b.disease_system) ?? 999);
    if (systemDelta !== 0) {
      return systemDelta;
    }
    const pheDelta = Number(a.phecode_sort) - Number(b.phecode_sort);
    if (pheDelta !== 0) {
      return pheDelta;
    }
    return String(a.disease_name).localeCompare(String(b.disease_name));
  }

  function rowTooltip(row) {
    return [
      "Medical condition: " + escapeHtml(row.disease_name),
      "Phecode: " + escapeHtml(row.phecode),
      "Physiological system: " + escapeHtml(row.disease_system_display),
      "Analysis: " + escapeHtml(row.analysis_name || row.cancer_name),
      escapeHtml(row.effect_measure || EFFECT_LABEL) + ": " + formatHR(row.hr_original) + " (" + formatCI(row) + ")",
      "p-value: " + formatP(row.p_value),
      "q-value: " + formatP(row.q_value),
      "Incidence rate, " + escapeHtml(row.group_one_label || "primary group") + ": " + formatRate(row.incidence_rate_group_one),
      "Incidence rate, " + escapeHtml(row.group_two_label || "comparison group") + ": " + formatRate(row.incidence_rate_group_two)
    ].join("<br>");
  }

  function shortCancerName(name) {
    return String(name || "")
      .replace("Population-based analysis (female)", "Population, female")
      .replace("Population-based analysis (male)", "Population, male")
      .replace("Population-based analysis", "Population-based")
      .replace("Sibling-based analysis (female)", "Sibling, female")
      .replace("Sibling-based analysis (male)", "Sibling, male")
      .replace("Sibling-based analysis", "Sibling-based");
  }

  function cancerLabelLines(name) {
    const label = shortCancerName(name);
    if (label.includes(", ")) {
      const parts = label.split(", ");
      return [parts[0], parts.slice(1).join(", ")];
    }
    return [label];
  }

  function setOptions(select, options, firstLabel) {
    select.innerHTML = "";
    if (firstLabel) {
      const allOption = document.createElement("option");
      allOption.value = "all";
      allOption.textContent = firstLabel;
      select.appendChild(allOption);
    }
    options.forEach((option) => {
      const element = document.createElement("option");
      element.value = option.value;
      element.textContent = option.label;
      select.appendChild(element);
    });
  }

  function presentSystems(rows) {
    const ids = new Set(rows.map((row) => row.disease_system));
    return SYSTEMS.filter((system) => ids.has(system.id)).map((system) => ({
      value: system.id,
      label: system.name
    }));
  }

  function initScatter() {
    const cancerSelect = byId("cancerFilter");
    const systemSelect = byId("systemFilter");
    const searchInput = byId("scatterSearch");
    const plot = byId("phewasPlot");
    let lastScatterYRange = null;
    let lastScatterXRange = null;

    if (!plot || !cancerSelect || !systemSelect || !searchInput) {
      return;
    }

    setOptions(cancerSelect, CANCERS.map((cancer) => ({ value: cancer.id, label: cancer.name })), null);
    const requestedCancer = new URLSearchParams(window.location.search).get("analysis");
    if (requestedCancer && CANCERS.some((cancer) => cancer.id === requestedCancer)) {
      cancerSelect.value = requestedCancer;
    }

    function currentCancer() {
      return CANCERS.find((item) => item.id === cancerSelect.value) || CANCERS[0];
    }

    function currentCancerRows() {
      const cancer = currentCancer();
      return ROWS.filter((row) => row.cancer_id === cancer.id).sort(compareRows);
    }

    function refreshSystemOptions() {
      const selectedSystem = systemSelect.value || "all";
      const pageRows = currentCancerRows();
      setOptions(systemSelect, presentSystems(pageRows), "All physiological systems");
      if (selectedSystem && Array.from(systemSelect.options).some((option) => option.value === selectedSystem)) {
        systemSelect.value = selectedSystem;
      }
    }

    function highlightLabel(row) {
      const name = String(row.disease_name || "");
      return name.length > 30 ? name.slice(0, 29) + "..." : name;
    }

    function scatterFilename() {
      const cancer = currentCancer();
      return "phewas_scatter_" + String(cancer.id || "cancer").toLowerCase();
    }

    function downloadScatter(format) {
      Plotly.downloadImage(plot, {
        format: format,
        filename: scatterFilename(),
        width: 1400,
        height: 860,
        scale: format === "png" ? 2 : 1
      });
    }

    function zoomScatter(multiplier) {
      const currentRange = plot.layout && plot.layout.yaxis && plot.layout.yaxis.range
        ? plot.layout.yaxis.range.slice()
        : (lastScatterYRange ? lastScatterYRange.slice() : null);
      if (!currentRange) {
        return;
      }
      const center = (Number(currentRange[0]) + Number(currentRange[1])) / 2;
      const halfWidth = (Number(currentRange[1]) - Number(currentRange[0])) * multiplier / 2;
      Plotly.relayout(plot, {
        "yaxis.range": [center - halfWidth, center + halfWidth]
      });
    }

    function resetScatterZoom() {
      if (!lastScatterYRange || !lastScatterXRange) {
        return;
      }
      Plotly.relayout(plot, {
        "xaxis.range": lastScatterXRange.slice(),
        "yaxis.range": lastScatterYRange.slice()
      });
    }

    function redraw() {
      const cancer = currentCancer();
      const pageRows = currentCancerRows();
      const system = systemSelect.value;
      const search = searchInput.value.trim().toLowerCase();
      const filtered = pageRows.filter((row) => {
        const systemMatch = system === "all" || row.disease_system === system;
        return systemMatch;
      }).sort(compareRows);
      filtered.forEach((row) => {
        const searchTarget = (row.disease_name + " " + row.phecode).toLowerCase();
        row._highlight = Boolean(search && searchTarget.includes(search));
      });
      const highlightedCount = filtered.filter((row) => row._highlight).length;
      const hasEffects = Boolean(DATA.metadata && DATA.metadata.has_effect_estimates);

      if (byId("cancerName")) {
        byId("cancerName").textContent = cancer ? cancer.name : "";
      }
      if (byId("retainedCount")) {
        byId("retainedCount").textContent = pageRows.length.toLocaleString();
      }
      if (byId("significantCount")) {
        byId("significantCount").textContent = pageRows.filter((row) => row.significant).length.toLocaleString();
      }

      if (byId("visibleCount")) {
        byId("visibleCount").textContent = filtered.length.toLocaleString();
      }
      if (byId("highlightedCount")) {
        byId("highlightedCount").textContent = highlightedCount.toLocaleString();
      }

      if (!filtered.length) {
        Plotly.purge(plot);
        plot.innerHTML = '<div class="empty-state">No results are available for the selected analysis and physiological system.</div>';
        return;
      }

      filtered.forEach((row, index) => {
        row._scatterX = index + 1;
      });

      const groupRanges = [];
      let currentSystem = null;
      let groupStart = 1;
      filtered.forEach((row, index) => {
        if (row.disease_system !== currentSystem) {
          if (currentSystem !== null) {
            groupRanges.push({ system: currentSystem, start: groupStart, end: index });
          }
          currentSystem = row.disease_system;
          groupStart = index + 1;
        }
      });
      groupRanges.push({ system: currentSystem, start: groupStart, end: filtered.length });

      const traces = SYSTEMS.map((systemInfo) => {
        const groupRows = filtered.filter((row) => row.disease_system === systemInfo.id);
        if (!groupRows.length) {
          return null;
        }
        return {
          type: "scatter",
          mode: "markers+text",
          name: systemInfo.name,
          x: groupRows.map((row) => row._scatterX),
          y: groupRows.map((row) => hasEffects ? row.hr_plot : row.count),
          text: groupRows.map((row) => row._highlight ? highlightLabel(row) : ""),
          textposition: "top center",
          textfont: {
            size: 11,
            color: "#111827"
          },
          customdata: groupRows.map((row) => [
            row.disease_name,
            row.phecode,
            row.disease_system_display,
            row.cancer_name,
            formatHR(row.hr_original),
            formatCI(row),
            formatP(row.p_value),
            formatP(row.q_value),
            formatRate(row.incidence_rate_group_one),
            formatRate(row.incidence_rate_group_two),
            row.group_two_label || "comparison group",
            row.effect_measure || EFFECT_LABEL,
            row.group_one_label || "primary group",
            formatNumber(row.count, 0)
          ]),
          hovertemplate: (hasEffects ? [
            "Medical condition: %{customdata[0]}",
            "Phecode: %{customdata[1]}",
            "Physiological system: %{customdata[2]}",
            "Analysis: %{customdata[3]}",
            "%{customdata[11]}: %{customdata[4]} (%{customdata[5]})",
            "p-value: %{customdata[6]}",
            "q-value: %{customdata[7]}",
            "Incidence rate, %{customdata[12]}: %{customdata[8]}",
            "Incidence rate, %{customdata[10]}: %{customdata[9]}<extra></extra>"
          ] : [
            "Medical condition: %{customdata[0]}",
            "Phecode: %{customdata[1]}",
            "Physiological system: %{customdata[2]}",
            "Individuals with condition: %{customdata[13]}<extra></extra>"
          ]).join("<br>"),
          marker: {
            color: groupRows.map((row) => hexToRgba(row.system_color, row.significant ? 0.92 : 0.34)),
            size: groupRows.map((row) => row._highlight ? 13 : (row.significant ? 9 : 7)),
            symbol: groupRows.map((row) => hasEffects ? (row.direction === "protective" ? "triangle-down" : "triangle-up") : "circle"),
            line: {
              color: groupRows.map((row) => row._highlight ? "#111827" : "#20242a"),
              width: groupRows.map((row) => row._highlight ? 1.8 : (row.significant ? 0.8 : 0))
            }
          }
        };
      }).filter(Boolean);

      const tickVals = groupRanges.map((range) => (range.start + range.end) / 2);
      const tickText = groupRanges.map((range) => {
        const systemInfo = SYSTEMS.find((item) => item.id === range.system);
        return systemInfo ? systemInfo.name : range.system;
      });
      lastScatterXRange = [0, filtered.length + 1];
      const plottedValues = filtered.map((row) => hasEffects ? row.hr_plot : Number(row.count || 0));
      const minEffect = Math.max(0.25, Math.min.apply(null, plottedValues) * 0.88);
      lastScatterYRange = hasEffects
        ? [Math.log10(minEffect), Math.log10(SCATTER_HR_CAP * 1.12)]
        : [0, Math.max.apply(null, plottedValues) * 1.12 || 1];
      const shapes = hasEffects ? [{
        type: "line",
        xref: "paper",
        x0: 0,
        x1: 1,
        yref: "y",
        y0: 1,
        y1: 1,
        line: { color: "#6b7280", width: 1.4, dash: "dash" }
      }] : [];
      groupRanges.slice(0, -1).forEach((range) => {
        shapes.push({
          type: "line",
          xref: "x",
          x0: range.end + 0.5,
          x1: range.end + 0.5,
          yref: "paper",
          y0: 0,
          y1: 1,
          line: { color: "#eeeeee", width: 1 }
        });
      });

      const layout = {
        paper_bgcolor: "rgba(0,0,0,0)",
        plot_bgcolor: "#fbfbfa",
        height: 740,
        margin: { l: 72, r: 210, t: 92, b: 150 },
        title: {
          text: cancer.name || DATA.metadata.subject_name || "DiNetxify analysis",
          x: 0.02,
          xanchor: "left",
          font: { size: 18, color: "#24282d" }
        },
        hovermode: "closest",
        showlegend: true,
        legend: {
          title: { text: "Physiological system" },
          x: 1.02,
          y: 1,
          xanchor: "left",
          yanchor: "top",
          bgcolor: "rgba(255,255,255,0.86)"
        },
        xaxis: {
          title: { text: "Physiological systems ordered by Phecode" },
          tickmode: "array",
          tickvals: tickVals,
          ticktext: tickText,
          tickangle: -35,
          showgrid: false,
          zeroline: false,
          range: lastScatterXRange
        },
        yaxis: {
          title: { text: hasEffects ? EFFECT_LABEL : (DATA.metadata.count_label || "PheWAS count") },
          type: hasEffects ? "log" : "linear",
          tickmode: hasEffects ? "array" : "auto",
          tickvals: hasEffects ? [0.25, 0.5, 1, 2, 5, 10, 20, 30] : undefined,
          ticktext: hasEffects ? ["0.25", "0.5", "1", "2", "5", "10", "20", "30"] : undefined,
          range: lastScatterYRange,
          gridcolor: "#e8eaed",
          zeroline: false
        },
        shapes: shapes,
        annotations: [{
          xref: "paper",
          yref: "paper",
          x: 0,
          y: 1.05,
          showarrow: false,
          align: "left",
          text: hasEffects ? "Triangle up = estimate > 1; triangle down = estimate < 1. Darker points indicate statistically significant results. Search terms label matching Phecodes without filtering the plot." : "Vertical position shows the PheWAS count. Search terms label matching Phecodes without filtering the plot.",
          font: { size: 12, color: "#667085" }
        }]
      };

      const config = {
        responsive: true,
        displayModeBar: false,
        displaylogo: false,
        modeBarButtonsToRemove: ["lasso2d", "select2d", "autoScale2d"]
      };
      Plotly.react(plot, traces, layout, config);
    }

    cancerSelect.addEventListener("change", () => {
      refreshSystemOptions();
      redraw();
      const url = new URL(window.location.href);
      url.searchParams.set("analysis", cancerSelect.value);
      history.replaceState(null, "", url);
    });
    systemSelect.addEventListener("change", redraw);
    searchInput.addEventListener("input", redraw);
    byId("scatterDownloadPng").addEventListener("click", () => downloadScatter("png"));
    byId("scatterDownloadSvg").addEventListener("click", () => downloadScatter("svg"));
    byId("scatterZoomIn").addEventListener("click", () => zoomScatter(0.75));
    byId("scatterZoomOut").addEventListener("click", () => zoomScatter(1.25));
    byId("scatterResetZoom").addEventListener("click", resetScatterZoom);
    refreshSystemOptions();
    redraw();
  }

  function diseaseIndex() {
    const map = new Map();
    ROWS.forEach((row) => {
      const key = row.phecode;
      if (!map.has(key)) {
        map.set(key, {
          key: key,
          phecode: row.phecode,
          disease_name: row.disease_name,
          disease_system: row.disease_system,
          disease_system_display: row.disease_system_display,
          system_color: row.system_color,
          cancerCount: 0,
          significantCount: 0
        });
      }
    });
    map.forEach((disease, key) => {
      const diseaseRows = ROWS.filter((row) => row.phecode === key);
      disease.cancerCount = new Set(diseaseRows.map((row) => row.cancer_id)).size;
      disease.significantCount = diseaseRows.filter((row) => row.significant).length;
    });
    return Array.from(map.values()).sort((a, b) => {
      if (b.significantCount !== a.significantCount) {
        return b.significantCount - a.significantCount;
      }
      if (b.cancerCount !== a.cancerCount) {
        return b.cancerCount - a.cancerCount;
      }
      return Number(a.phecode) - Number(b.phecode);
    });
  }

  function initForest() {
    const systemFilter = byId("systemFilter");
    const diseaseSearch = byId("diseaseSearch");
    const diseaseSelect = byId("diseaseSelect");
    const addButton = byId("addDiseaseButton");
    const clearDiseaseButton = byId("clearDiseases");
    const cancerBox = byId("cancerCheckboxes");
    const plotContainer = byId("plotContainer");
    const tooltip = byId("tooltip");
    if (!systemFilter || !diseaseSearch || !diseaseSelect || !addButton || !clearDiseaseButton || !cancerBox || !plotContainer) {
      return;
    }

    const ADD_ALL_VISIBLE_OPTION = "__add_all_visible__";
    const diseases = diseaseIndex();
    const diseaseMap = new Map(diseases.map((disease) => [disease.key, disease]));
    const rowByDiseaseCancer = new Map(ROWS.map((row) => [row.phecode + "|" + row.cancer_id, row]));
    const fallbackDefault = diseases.slice(0, 6).map((disease) => disease.key);
    const manuscriptDefaults = ["300.4", "966", "327.4", "338.2", "250.2", "626.8"]
      .filter((key) => diseaseMap.has(key));
    let selectedDiseases = (manuscriptDefaults.length ? manuscriptDefaults : fallbackDefault).slice(0, 6);
    let selectedCancers = CANCERS.map((cancer) => cancer.id);
    let forestZoom = 1;
    let currentDiseaseOptions = [];

    setOptions(systemFilter, SYSTEMS.map((system) => ({ value: system.id, label: system.name })), "All physiological systems");

    function updateDiseaseActionState() {
      const atDiseaseLimit = selectedDiseases.length >= MAX_DISEASES;
      const selectedValue = diseaseSelect.value;
      const hasUnselectedVisibleDisease = currentDiseaseOptions.some((disease) => !selectedDiseases.includes(disease.key));
      const canAddAllVisible = selectedValue === ADD_ALL_VISIBLE_OPTION && hasUnselectedVisibleDisease && !atDiseaseLimit;
      const canAddSingle = selectedValue && selectedValue !== ADD_ALL_VISIBLE_OPTION && !selectedDiseases.includes(selectedValue) && !atDiseaseLimit;
      addButton.disabled = !(canAddAllVisible || canAddSingle);
      clearDiseaseButton.disabled = selectedDiseases.length === 0;
    }

    function updateDiseaseOptions() {
      const system = systemFilter.value;
      const query = diseaseSearch.value.trim().toLowerCase();
      const filtered = diseases.filter((disease) => {
        const systemMatch = system === "all" || disease.disease_system === system;
        const searchTarget = (disease.disease_name + " " + disease.phecode).toLowerCase();
        return systemMatch && (!query || searchTarget.includes(query));
      }).slice(0, 250);
      currentDiseaseOptions = filtered;
      diseaseSelect.innerHTML = "";
      const addAllOption = document.createElement("option");
      addAllOption.value = ADD_ALL_VISIBLE_OPTION;
      addAllOption.textContent = "Add all conditions in current list (" + filtered.length + ")";
      diseaseSelect.appendChild(addAllOption);
      filtered.forEach((disease) => {
        const option = document.createElement("option");
        option.value = disease.key;
        option.textContent = disease.disease_name + " [" + disease.phecode + "]";
        diseaseSelect.appendChild(option);
      });
      updateDiseaseActionState();
    }

    function renderCancerSelections() {
      cancerBox.innerHTML = "";
      CANCERS.forEach((cancer) => {
        const label = document.createElement("label");
        label.className = "checkbox-item";
        const checkbox = document.createElement("input");
        checkbox.type = "checkbox";
        checkbox.value = cancer.id;
        checkbox.checked = selectedCancers.includes(cancer.id);
        checkbox.addEventListener("change", () => {
          if (checkbox.checked && !selectedCancers.includes(cancer.id)) {
            selectedCancers.push(cancer.id);
          }
          if (!checkbox.checked) {
            selectedCancers = selectedCancers.filter((id) => id !== cancer.id);
          }
          renderAll();
        });
        label.appendChild(checkbox);
        label.appendChild(document.createTextNode(cancer.name));
        cancerBox.appendChild(label);
      });
    }

    function renderChips() {
      const diseaseChips = byId("selectedDiseaseChips");
      const cancerChips = byId("selectedCancerChips");
      diseaseChips.innerHTML = "";
      selectedDiseases.forEach((key) => {
        const disease = diseaseMap.get(key);
        if (!disease) {
          return;
        }
        const chip = document.createElement("span");
        chip.className = "chip";
        chip.innerHTML = "<span>" + escapeHtml(disease.disease_name) + " [" + escapeHtml(disease.phecode) + "]</span>";
        const remove = document.createElement("button");
        remove.className = "remove";
        remove.type = "button";
        remove.textContent = "x";
        remove.addEventListener("click", () => {
          selectedDiseases = selectedDiseases.filter((item) => item !== key);
          renderAll();
        });
        chip.appendChild(remove);
        diseaseChips.appendChild(chip);
      });
      cancerChips.innerHTML = "";
      selectedCancers.forEach((id) => {
        const cancer = CANCERS.find((item) => item.id === id);
        if (!cancer) {
          return;
        }
        const chip = document.createElement("span");
        chip.className = "chip";
        chip.innerHTML = "<span>" + escapeHtml(cancer.name) + "</span>";
        const remove = document.createElement("button");
        remove.className = "remove";
        remove.type = "button";
        remove.textContent = "x";
        remove.addEventListener("click", () => {
          selectedCancers = selectedCancers.filter((item) => item !== id);
          renderCancerSelections();
          renderAll();
        });
        chip.appendChild(remove);
        cancerChips.appendChild(chip);
      });
      if (byId("diseaseCount")) {
        byId("diseaseCount").textContent = selectedDiseases.length;
      }
      if (byId("cancerCount")) {
        byId("cancerCount").textContent = selectedCancers.length;
      }
      updateDiseaseActionState();
    }

    function clamp(value, min, max) {
      return Math.max(min, Math.min(max, value));
    }

    function logScale(value, x0, width) {
      const clamped = clamp(Number(value), FOREST_MIN, FOREST_MAX);
      const minLog = Math.log(FOREST_MIN);
      const maxLog = Math.log(FOREST_MAX);
      return x0 + ((Math.log(clamped) - minLog) / (maxLog - minLog)) * width;
    }

    function svgEl(name, attrs, text) {
      const element = document.createElementNS("http://www.w3.org/2000/svg", name);
      Object.entries(attrs || {}).forEach(([key, value]) => {
        element.setAttribute(key, value);
      });
      if (text != null) {
        element.textContent = text;
      }
      return element;
    }

    function showTooltip(html, event) {
      if (!tooltip) {
        return;
      }
      tooltip.innerHTML = html;
      tooltip.style.display = "block";
      tooltip.style.left = Math.min(event.clientX + 16, window.innerWidth - 360) + "px";
      tooltip.style.top = Math.min(event.clientY + 16, window.innerHeight - 220) + "px";
    }

    function hideTooltip() {
      if (tooltip) {
        tooltip.style.display = "none";
      }
    }

    function currentForestSvg() {
      return plotContainer.querySelector("svg");
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

    function downloadForestSvg() {
      const svg = currentForestSvg();
      if (!svg) {
        return;
      }
      const clone = svg.cloneNode(true);
      clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
      const style = document.createElementNS("http://www.w3.org/2000/svg", "style");
      style.textContent = ".axis-label{fill:#344054;font:700 13px Arial}.grid-line{stroke:#e5e7eb}.disease-label{fill:#1f2937;font:700 13px Arial}.phecode-label,.muted-label{fill:#667085;font:11px Arial}.reference-line{stroke:#98a2b3;stroke-dasharray:3 3}.ci-line{stroke-width:3.4;stroke-linecap:round}.point{stroke:#20242a;stroke-width:1.3}.column-label{fill:#344054;font:700 12px Arial}";
      clone.insertBefore(style, clone.firstChild);
      const text = new XMLSerializer().serializeToString(clone);
      downloadBlob(new Blob([text], { type: "image/svg+xml;charset=utf-8" }), "phewas_forest_plot.svg");
    }

    function downloadForestPng() {
      const svg = currentForestSvg();
      if (!svg) {
        return;
      }
      const clone = svg.cloneNode(true);
      clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
      const style = document.createElementNS("http://www.w3.org/2000/svg", "style");
      style.textContent = ".axis-label{fill:#344054;font:700 13px Arial}.grid-line{stroke:#e5e7eb}.disease-label{fill:#1f2937;font:700 13px Arial}.phecode-label,.muted-label{fill:#667085;font:11px Arial}.reference-line{stroke:#98a2b3;stroke-dasharray:3 3}.ci-line{stroke-width:3.4;stroke-linecap:round}.point{stroke:#20242a;stroke-width:1.3}.column-label{fill:#344054;font:700 12px Arial}";
      clone.insertBefore(style, clone.firstChild);
      const text = new XMLSerializer().serializeToString(clone);
      const blob = new Blob([text], { type: "image/svg+xml;charset=utf-8" });
      const url = URL.createObjectURL(blob);
      const image = new Image();
      const width = Number(svg.getAttribute("width")) || 1200;
      const height = Number(svg.getAttribute("height")) || 800;
      image.onload = () => {
        const canvas = document.createElement("canvas");
        canvas.width = width * 2;
        canvas.height = height * 2;
        const context = canvas.getContext("2d");
        context.fillStyle = "#ffffff";
        context.fillRect(0, 0, canvas.width, canvas.height);
        context.scale(2, 2);
        context.drawImage(image, 0, 0);
        URL.revokeObjectURL(url);
        canvas.toBlob((pngBlob) => {
          if (pngBlob) {
            downloadBlob(pngBlob, "phewas_forest_plot.png");
          }
        }, "image/png");
      };
      image.src = url;
    }

    function drawForest() {
      plotContainer.innerHTML = "";
      if (!selectedDiseases.length || !selectedCancers.length) {
        plotContainer.innerHTML = '<div class="empty-state">Select at least one disease and one analysis to draw the forest plot.</div>';
        return;
      }

      const labelWidth = 320;
      const columnWidth = Math.round(190 * forestZoom);
      const top = 250;
      const rowHeight = 54;
      const bottom = 70;
      const width = labelWidth + selectedCancers.length * columnWidth + 180;
      const height = top + selectedDiseases.length * rowHeight + bottom;
      const svg = svgEl("svg", {
        width: width,
        height: height,
        viewBox: "0 0 " + width + " " + height,
        role: "img",
        "aria-label": "Interactive disease-by-analysis forest plot"
      });

      svg.appendChild(svgEl("text", { x: 16, y: 34, class: "axis-label" }, "Medical conditions and Phecodes"));
      selectedCancers.forEach((id, index) => {
        const x = labelWidth + index * columnWidth;
        svg.appendChild(svgEl("line", {
          x1: x,
          x2: x,
          y1: top - 12,
          y2: height - bottom + 8,
          class: "grid-line"
        }));
      });

      selectedDiseases.forEach((key, rowIndex) => {
        const disease = diseaseMap.get(key);
        if (!disease) {
          return;
        }
        const rowTop = top + rowIndex * rowHeight;
        const centerY = rowTop + rowHeight / 2;
        if (rowIndex % 2 === 0) {
          svg.appendChild(svgEl("rect", {
            x: 0,
            y: rowTop,
            width: width,
            height: rowHeight,
            fill: "#fbfbfa"
          }));
        }
        svg.appendChild(svgEl("text", {
          x: 16,
          y: centerY - 5,
          class: "disease-label"
        }, disease.disease_name));
        svg.appendChild(svgEl("text", {
          x: 16,
          y: centerY + 14,
          class: "phecode-label"
        }, "Phecode " + disease.phecode + " | " + disease.disease_system_display));

        selectedCancers.forEach((cancerId, columnIndex) => {
          const cellX = labelWidth + columnIndex * columnWidth;
          const plotX = cellX + 24;
          const plotWidth = columnWidth - 48;
          const refX = logScale(1, plotX, plotWidth);
          svg.appendChild(svgEl("line", {
            x1: refX,
            x2: refX,
            y1: rowTop + 7,
            y2: rowTop + rowHeight - 7,
            class: "reference-line"
          }));
          const row = rowByDiseaseCancer.get(key + "|" + cancerId);
          if (!row) {
            svg.appendChild(svgEl("text", {
              x: cellX + columnWidth / 2,
              y: centerY + 4,
              class: "muted-label",
              "text-anchor": "middle"
            }, "NA"));
            return;
          }

          const color = row.system_color || disease.system_color || "#555555";
          const ciLower = row.ci_lower == null ? row.hr_original : row.ci_lower;
          const ciUpper = row.ci_upper == null ? row.hr_original : row.ci_upper;
          const lineX1 = logScale(ciLower, plotX, plotWidth);
          const lineX2 = logScale(ciUpper, plotX, plotWidth);
          const pointOpacity = row.significant ? 0.98 : 0.68;
          svg.appendChild(svgEl("line", {
            x1: lineX1,
            x2: lineX2,
            y1: centerY,
            y2: centerY,
            class: "ci-line",
            stroke: color,
            opacity: pointOpacity
          }));
          if (ciLower >= FOREST_MIN) {
            svg.appendChild(svgEl("line", {
              x1: lineX1,
              x2: lineX1,
              y1: centerY - 5,
              y2: centerY + 5,
              stroke: color,
              "stroke-width": 2.4,
              opacity: pointOpacity
            }));
          }
          if (ciUpper <= FOREST_MAX) {
            svg.appendChild(svgEl("line", {
              x1: lineX2,
              x2: lineX2,
              y1: centerY - 5,
              y2: centerY + 5,
              stroke: color,
              "stroke-width": 2.4,
              opacity: pointOpacity
            }));
          }
          if (ciLower < FOREST_MIN) {
            svg.appendChild(svgEl("text", {
              x: plotX - 9,
              y: centerY + 5,
              fill: color,
              "font-weight": "700",
              "text-anchor": "middle"
            }, "<"));
          }
          if (ciUpper > FOREST_MAX || row.hr_original > FOREST_MAX) {
            svg.appendChild(svgEl("text", {
              x: plotX + plotWidth + 9,
              y: centerY + 5,
              fill: color,
              "font-weight": "700",
              "text-anchor": "middle"
            }, ">"));
          }
          const pointX = logScale(row.hr_original, plotX, plotWidth);
          svg.appendChild(svgEl("circle", {
            cx: pointX,
            cy: centerY,
            r: row.significant ? 5.2 : 4.2,
            fill: color,
            opacity: pointOpacity,
            class: "point"
          }));
          const hover = svgEl("rect", {
            x: cellX,
            y: rowTop,
            width: columnWidth,
            height: rowHeight,
            class: "cell-hover"
          });
          hover.addEventListener("mousemove", (event) => showTooltip(rowTooltip(row), event));
          hover.addEventListener("mouseleave", hideTooltip);
          svg.appendChild(hover);
        });
      });

      const axisY = height - 25;
      selectedCancers.forEach((id, index) => {
        const cellX = labelWidth + index * columnWidth;
        const plotX = cellX + 24;
        const plotWidth = columnWidth - 48;
        FOREST_TICKS.forEach((tick) => {
          const x = logScale(tick, plotX, plotWidth);
          svg.appendChild(svgEl("line", {
            x1: x,
            x2: x,
            y1: axisY - 7,
            y2: axisY - 2,
            stroke: "#6b7280",
            "stroke-width": 1
          }));
          svg.appendChild(svgEl("text", {
            x: x,
            y: axisY + 12,
            class: "muted-label",
            "text-anchor": "middle"
          }, String(tick)));
        });
      });
      svg.appendChild(svgEl("text", {
        x: labelWidth,
        y: axisY + 12,
        class: "muted-label",
        "text-anchor": "end"
      }, "Log " + EFFECT_LABEL + " scale: 0.7 to 30"));

      selectedCancers.forEach((id, index) => {
        const cancer = CANCERS.find((item) => item.id === id);
        const cellX = labelWidth + index * columnWidth;
        const plotX = cellX + 24;
        const plotWidth = columnWidth - 48;
        const labelX = logScale(1, plotX, plotWidth);
        const labelY = top - 26;
        const label = svgEl("text", {
          x: labelX,
          y: labelY,
          class: "column-label",
          "text-anchor": "start",
          transform: "rotate(-52 " + labelX + " " + labelY + ")"
        });
        cancerLabelLines(cancer ? cancer.name : id).forEach((line, lineIndex) => {
          label.appendChild(svgEl("tspan", {
            x: labelX,
            dy: lineIndex === 0 ? 0 : "1.05em"
          }, line));
        });
        svg.appendChild(label);
      });

      plotContainer.appendChild(svg);
    }

    function renderAll() {
      renderChips();
      hideTooltip();
      drawForest();
    }

    systemFilter.addEventListener("change", updateDiseaseOptions);
    diseaseSearch.addEventListener("input", updateDiseaseOptions);
    diseaseSelect.addEventListener("change", updateDiseaseActionState);
    addButton.addEventListener("click", () => {
      const key = diseaseSelect.value;
      if (!key || selectedDiseases.length >= MAX_DISEASES) {
        return;
      }
      if (key === ADD_ALL_VISIBLE_OPTION) {
        const remainingSlots = MAX_DISEASES - selectedDiseases.length;
        const keysToAdd = currentDiseaseOptions
          .map((disease) => disease.key)
          .filter((diseaseKey) => !selectedDiseases.includes(diseaseKey))
          .slice(0, remainingSlots);
        if (!keysToAdd.length) {
          return;
        }
        selectedDiseases = selectedDiseases.concat(keysToAdd);
        renderAll();
        return;
      }
      if (selectedDiseases.includes(key)) {
        return;
      }
      selectedDiseases.push(key);
      renderAll();
    });
    clearDiseaseButton.addEventListener("click", () => {
      selectedDiseases = [];
      renderAll();
    });
    byId("selectAllCancers").addEventListener("click", () => {
      selectedCancers = CANCERS.map((cancer) => cancer.id);
      renderCancerSelections();
      renderAll();
    });
    byId("clearCancers").addEventListener("click", () => {
      selectedCancers = [];
      renderCancerSelections();
      renderAll();
    });
    byId("forestDownloadSvg").addEventListener("click", downloadForestSvg);
    byId("forestDownloadPng").addEventListener("click", downloadForestPng);
    byId("forestZoomIn").addEventListener("click", () => {
      forestZoom = Math.min(1.8, forestZoom + 0.15);
      renderAll();
    });
    byId("forestZoomOut").addEventListener("click", () => {
      forestZoom = Math.max(0.75, forestZoom - 0.15);
      renderAll();
    });
    byId("forestResetZoom").addEventListener("click", () => {
      forestZoom = 1;
      renderAll();
    });

    updateDiseaseOptions();
    renderCancerSelections();
    renderAll();
  }

  window.PhewasApp = {
    initScatter: initScatter,
    initForest: initForest
  };
})();
