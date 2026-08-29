(function () {
  "use strict";

  function escapeHtml(value) {
    return String(value == null ? "" : value).replace(/[&<>"']/g, (char) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#039;" }[char]));
  }

  function fmt(value, digits) {
    return value == null || !Number.isFinite(Number(value)) ? "Not available" : Number(value).toFixed(digits == null ? 2 : digits);
  }

  function fmtP(value) {
    return value == null || !Number.isFinite(Number(value)) ? "Not available" : Number(value) < .001 ? Number(value).toExponential(2) : Number(value).toFixed(3);
  }

  function download(blob, filename) {
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url; link.download = filename; document.body.append(link); link.click(); link.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }

  function serializedSvg(svg) {
    const clone = svg.cloneNode(true);
    clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
    const style = document.createElementNS("http://www.w3.org/2000/svg", "style");
    style.textContent = ".module-boundary{fill-opacity:.07;stroke-opacity:.48;stroke-width:1.5;stroke-dasharray:6 4}.module-title,.node-label{fill:#344054;stroke:#fff;stroke-width:4;paint-order:stroke;font:700 12px Arial;text-anchor:middle}.network-edge,.trajectory-edge{stroke:#667085;stroke-opacity:.3}.network-node,.trajectory-node{stroke:#fff;stroke-width:1.5}.endpoint-node{fill:#111827;stroke:#fff;stroke-width:2}.node-label{display:none}.node-label.visible{display:block}";
    clone.prepend(style);
    return new XMLSerializer().serializeToString(clone);
  }

  function exportSvg(svg, filename) {
    download(new Blob([serializedSvg(svg)], { type: "image/svg+xml;charset=utf-8" }), filename);
  }

  function exportPng(svg, filename) {
    const value = serializedSvg(svg);
    const blob = new Blob([value], { type: "image/svg+xml;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const image = new Image();
    image.onload = function () {
      const rect = svg.getBoundingClientRect();
      const scale = 3;
      const canvas = document.createElement("canvas");
      canvas.width = Math.max(1, Math.round(rect.width * scale));
      canvas.height = Math.max(1, Math.round(rect.height * scale));
      const context = canvas.getContext("2d");
      context.fillStyle = "#ffffff"; context.fillRect(0, 0, canvas.width, canvas.height);
      context.drawImage(image, 0, 0, canvas.width, canvas.height);
      canvas.toBlob((result) => { if (result) download(result, filename); }, "image/png");
      URL.revokeObjectURL(url);
    };
    image.src = url;
  }

  function panZoom(svg, viewport) {
    const state = { x: 0, y: 0, scale: 1, dragging: false, startX: 0, startY: 0 };
    function apply() { viewport.setAttribute("transform", "translate(" + state.x + " " + state.y + ") scale(" + state.scale + ")"); }
    function reset() { state.x = 0; state.y = 0; state.scale = 1; apply(); }
    function zoomBy(factor) { state.scale = Math.min(5, Math.max(.25, state.scale * factor)); apply(); }
    svg.addEventListener("wheel", (event) => { event.preventDefault(); zoomBy(event.deltaY < 0 ? 1.12 : .89); }, { passive: false });
    svg.addEventListener("pointerdown", (event) => { state.dragging = true; state.startX = event.clientX; state.startY = event.clientY; svg.classList.add("dragging"); svg.setPointerCapture(event.pointerId); });
    svg.addEventListener("pointermove", (event) => {
      if (!state.dragging) return;
      const box = svg.viewBox.baseVal;
      const rect = svg.getBoundingClientRect();
      state.x += (event.clientX - state.startX) * box.width / rect.width / state.scale;
      state.y += (event.clientY - state.startY) * box.height / rect.height / state.scale;
      state.startX = event.clientX; state.startY = event.clientY; apply();
    });
    function stop(event) { state.dragging = false; svg.classList.remove("dragging"); if (event.pointerId != null && svg.hasPointerCapture(event.pointerId)) svg.releasePointerCapture(event.pointerId); }
    svg.addEventListener("pointerup", stop); svg.addEventListener("pointercancel", stop);
    return { reset: reset, zoomBy: zoomBy };
  }

  window.DiNetxifyInteractiveCommon = { escapeHtml: escapeHtml, fmt: fmt, fmtP: fmtP, exportSvg: exportSvg, exportPng: exportPng, panZoom: panZoom };
}());

