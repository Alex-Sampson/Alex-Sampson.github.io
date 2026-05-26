---
title: "Markdown patterns I use in notes"
description: "Footnotes, definitions, KaTeX equations, callouts. A small style guide for the way I write here — partly for future me, partly so the rendering stays consistent across posts."
date: 2026-05-23
katex: true
ai_placeholder: true
---

<section>
  <p>
    Most of these notes are short Markdown files. I'd like them to render
    consistently, both visually and semantically, so a returning reader (often me,
    three months later) finds the same conventions on every page. What follows is a
    minimal style guide.
  </p>
</section>

<section>
  <h3 class="section-label">Headings</h3>
  <p>
    Each post opens with a short lede in the page header — written as the
    <code>description</code> in front-matter, not inside the body. Inside the body,
    sections are introduced by the <code>section-label</code> pattern you see above
    (◆&nbsp;HEADING with the hairline rule), not by raw <code>##</code> in Markdown.
    A small Liquid include or custom kramdown post-processor can map <code>##</code>
    to the right span structure so I don't have to think about it while writing.
  </p>
</section>

<section>
  <h3 class="section-label">Inline math &amp; display equations</h3>
  <p>
    Math is rendered with <code>KaTeX</code> at build time. Inline expressions are
    wrapped in single dollar signs in Markdown; display equations are wrapped in
    double dollar signs. For example, the gradient of a smooth function on a
    Riemannian manifold satisfies <span class="math-inline" id="mp-grad"></span>, and
    the Riemannian inner product on tangent vectors is written
    <span class="math-inline" id="mp-inner"></span>.
  </p>
  <p>
    A display equation:
  </p>
  <div class="metric-eq" role="math">
    <p class="eq-label">Sample (eq. demo)</p>
    <div id="mp-display"></div>
  </div>
</section>

<section>
  <h3 class="section-label">Callouts &amp; quotes</h3>
  <p>
    Blockquotes get the warm left-bar treatment for short pulls — usually a sentence
    or two from a paper I'm reading, or a line from a book worth keeping at the top
    of the page.
  </p>
  <blockquote>
    A good notation has a subtlety and suggestiveness which at times make it almost
    seem like a live teacher.
    <cite>— Bertrand Russell</cite>
  </blockquote>
  <p>
    Longer prose-paragraph callouts (like the placeholder notice above) live in a
    tinted panel with hairline rules, top and bottom. I'd reserve those for asides
    that aren't part of the running argument — caveats, definitions, the occasional
    "if you haven't seen this, do this first."
  </p>
</section>

<section>
  <h3 class="section-label">Code &amp; mono</h3>
  <p>
    Short identifiers live in <code>inline code</code> with the warm accent
    background. Longer code blocks are uncommon here — when they show up, they get
    the page's <code>monospace</code> face on a slightly sunk paper background, no
    syntax highlighting unless the code genuinely benefits from it.
  </p>
</section>

<section>
  <h3 class="section-label">Links &amp; footnotes</h3>
  <p>
    Links are underlined in the accent color; the underline thickens on hover. I try
    to keep external links to the things you actually need to follow, and to use
    footnotes for asides that would otherwise interrupt the sentence. Footnote
    rendering is provided by the Jekyll <code>kramdown</code> processor.
  </p>
</section>

<script type="application/json" id="katex-data">
[
  ["mp-grad",    "\\nabla f \\in T_pM"],
  ["mp-inner",   "\\langle u, v \\rangle_g"],
  ["mp-display", "\\int_M \\langle \\nabla f, \\nabla g \\rangle_g \\, \\mathrm{vol} = -\\int_M f \\, \\Delta g \\, \\mathrm{vol}", true]
]
</script>
