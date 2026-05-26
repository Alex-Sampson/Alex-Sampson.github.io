/* ============================================================
   site.js — Jekyll build
   The topbar is rendered server-side from _data/nav.yml; this
   script now only handles the live Tweaks panel (accent swap).
   ============================================================ */

const DEFAULTS = /*EDITMODE-BEGIN*/{
  "accent": "garnet"
}/*EDITMODE-END*/;

let TWEAKS = { ...DEFAULTS };

function applyTweaks() {
  document.documentElement.setAttribute("data-accent", TWEAKS.accent);
}

function renderTweaks() {
  const t = document.createElement("div");
  t.className = "tweaks";
  t.id = "tweaks";
  t.innerHTML = `
    <div class="row">
      <h4>Accent</h4>
      <div class="swatches">
        <div class="sw ${TWEAKS.accent==='garnet'?'active':''}"   data-c="garnet"   title="Garnet"></div>
        <div class="sw ${TWEAKS.accent==='ochre'?'active':''}"    data-c="ochre"    title="Ochre"></div>
        <div class="sw ${TWEAKS.accent==='ink-blue'?'active':''}" data-c="ink-blue" title="Ink Blue"></div>
        <div class="sw ${TWEAKS.accent==='sage'?'active':''}"     data-c="sage"     title="Sage"></div>
      </div>
    </div>
    <div class="row" style="display:flex;justify-content:space-between;align-items:center;font-family:var(--sans);font-size:0.75rem;color:var(--ink-faint);">
      <span>Tweaks · live preview</span>
      <button id="tweaks-close" style="background:none;border:0;cursor:pointer;color:var(--ink-faint);font-family:var(--sans);font-size:0.95rem;">×</button>
    </div>
  `;
  document.body.appendChild(t);

  t.querySelectorAll(".sw").forEach(sw => {
    sw.addEventListener("click", () => {
      TWEAKS.accent = sw.dataset.c;
      applyTweaks();
      t.querySelectorAll(".sw").forEach(s => s.classList.toggle("active", s.dataset.c === TWEAKS.accent));
    });
  });

  t.querySelector("#tweaks-close").addEventListener("click", () => {
    t.classList.remove("open");
  });
}

document.addEventListener("DOMContentLoaded", () => {
  renderTweaks();
  applyTweaks();
});
