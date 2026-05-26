# Jekyll port — deployment guide

This folder is a complete Jekyll site that produces the same HTML as the
`site/` static mockup. Drop the contents into a branch of your existing
GitHub Pages repo, push, and you're live.

---

## What's here

```
jekyll-site/
├── _config.yml                 site config + analytics + CV metadata
├── Gemfile                     gem dependencies (github-pages bundle)
├── _layouts/
│   ├── default.html            <html>…<body> shell + topbar + scripts
│   └── post.html               full post structure (back link, lede, sidebar)
├── _includes/
│   ├── head.html               <head> tags + meta + OG + Umami + favicon
│   ├── topbar.html             sticky nav, driven by _data/nav.yml
│   ├── jsonld-person.html      Person schema (homepage only)
│   └── katex.html              KaTeX render bootstrap
├── _data/
│   ├── nav.yml                 the five nav items
│   ├── projects.yml            project rows (FCM I & II)
│   └── teaching.yml            teaching record
├── _posts/
│   ├── 2026-05-25-welcome.md
│   ├── 2026-05-24-jekyll.md
│   └── 2026-05-23-markdown-patterns.md
├── assets/
│   ├── css/site.css            full stylesheet
│   ├── js/site.js              live Tweaks panel (accent swap)
│   ├── images/                 banner, portrait, surface figure
│   └── pdf-documents/cv-updates/AlexSampson-CV.pdf
├── index.html                  homepage (layout: default, includes math)
├── projects.html               loops over _data/projects.yml
├── teaching.html               loops over _data/teaching.yml
├── notes.html                  loops over site.posts
└── cv.html                     CV summary + PDF download
```

The CSS and JS are byte-identical to `site/assets/`. The rendered HTML
will match the static mockup essentially exactly (same DOM, same classes,
same content) — Jekyll just *generates* it from the Markdown / YAML
sources instead of you hand-writing each page.

---

## Deploying to your GitHub repo

You said your existing repo is `Alex-Sampson/Alex-Sampson.github.io`.
Here's the cleanest path:

### 1. Create a new branch on your local clone

```bash
cd /path/to/Alex-Sampson.github.io
git checkout master
git pull
git checkout -b redesign
```

### 2. Remove the old Jekyll content from this branch

```bash
git rm -rf _layouts _includes _sass _posts assets _Private-Notes \
          index.md cv.md projects.md teaching.md blog.html archive.html 404.html
```

(Leave `.github/`, `.gitignore`, `Gemfile`, `_config.yml`,
`google77d6ea25d4901736.html`, `robots.txt`, `README.md`, and
`UNLICENSE.txt` in place — we'll overwrite the config and gemfile.)

### 3. Copy this folder's contents into the repo root

Everything inside `jekyll-site/` (not the folder itself) goes at the repo
root:

```bash
cp -r /path/to/jekyll-site/. .
```

After this, `index.html`, `_config.yml`, `_layouts/`, `_posts/`, `assets/`,
etc. should all be at the repo root.

### 4. Test locally

```bash
bundle install
bundle exec jekyll serve
```

Visit `http://127.0.0.1:4000/`. You should see the same site you've been
reviewing.

### 5. Push and switch GitHub Pages source

```bash
git add -A
git commit -m "Redesign: warm scholarly site"
git push -u origin redesign
```

In your repo on github.com → **Settings → Pages → Build and deployment**:
- **Source**: Deploy from a branch
- **Branch**: `redesign` / `(root)`
- Save

GitHub will rebuild and serve the new site at `alex-sampson.github.io`
within a minute or two.

When you're happy with it, you can either keep deploying from `redesign`
or merge `redesign` into `master` and switch the Pages source back.

---

## Day-to-day updates

| You want to … | Edit |
|---|---|
| Add a new project | append an item to `_data/projects.yml` |
| Add a new course or teaching material | append an item to `_data/teaching.yml` |
| Add a new note / blog post | drop a Markdown file in `_posts/` named `YYYY-MM-DD-slug.md` |
| Edit homepage prose | edit `index.html` (the body, not the layout) |
| Edit CV summary | edit `cv.html` |
| Add a new top-level page | create `<name>.md` with `layout: default` and a permalink, then add a row to `_data/nav.yml` |
| Update the PDF CV | replace `assets/pdf-documents/cv-updates/AlexSampson-CV.pdf` and bump `cv.updated` in `_config.yml` |
| Change site title / description / Umami ID | edit `_config.yml` |
| Change visual design (colors, type, spacing) | edit `assets/css/site.css` |

---

## Notes on the port

- **Topbar** is now rendered by Liquid (`_includes/topbar.html`) from
  `_data/nav.yml`, not by JS. The "current page" highlight uses the
  `section:` value from each page's front-matter. Posts inherit
  `section: notes` from the `defaults:` block in `_config.yml`, so the
  Notes tab stays active on every post page.

- **Meta tags + Umami** live in one place (`_includes/head.html`) and
  pull title/description from each page's front-matter. To change the
  global OG image, the Umami site ID, or the JSON-LD schema, edit there
  once.

- **KaTeX** is opt-in per page via `katex: true` front-matter. The
  homepage and the markdown-patterns post are the two pages that
  currently use it.

- **Permalinks**: pages use the `permalink:` field in their front-matter
  (e.g. `/projects/`). Posts use `/notes/<slug>/` from
  `_config.yml`'s site-wide `permalink:` setting.

- **Search engines**: `jekyll-sitemap` generates `sitemap.xml` at build
  time. `jekyll-feed` generates `feed.xml` at `/feed.xml`. Both are
  linked from the `<head>` via `{% feed_meta %}`.

---

## If anything breaks

- **Liquid syntax error** while building: most often a `{% raw %}` block
  in a post that wraps something Jekyll already wouldn't try to parse.
  Run `bundle exec jekyll build --trace` to see the file + line.
- **Layout not applied**: check that the page's front-matter `layout:`
  value matches a file in `_layouts/`.
- **Topbar not highlighting**: check that the page's `section:`
  front-matter matches a `section:` value in `_data/nav.yml`.
- **PDF link 404**: filename in `_config.yml` (`cv.pdf:`) must match the
  actual file in `assets/pdf-documents/cv-updates/`.
