# Author: Xu Chen, University of Washington
"""Write index.html (the landing page) and an HTML copy of every lab README.md and homework.md.

usage:  python tools/make_index.py          (needs pandoc for the HTML copies; the .md files stay the source)

Everything is read from the lab folders: titles and overviews from README.md, the problems and points from
homework.md, the demo titles from the first markdown heading of each demo_*.py, the lecture map from the table
in the top-level README.md, and the figures from labs/*/figures/scene.png and out/ (made by run_all.py).
"""
import glob, html, os, re, shutil, subprocess

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LECTURES = range(2, 27)
REPO = 'https://github.com/xchen-me/control-in-practice'
CSS = """
:root { --bg:#faf8f2; --card:#fff; --ink:#1d1d1f; --muted:#62605a; --line:#ddd6c8; --accent:#1f5fa8; --hit:#c9dcf2; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --bg:#17181a; --card:#202124; --ink:#e8e6e1; --muted:#a19e96; --line:#38362f; --accent:#7fb0ea; --hit:#2c4b70; } }
:root[data-theme="dark"] { --bg:#17181a; --card:#202124; --ink:#e8e6e1; --muted:#a19e96; --line:#38362f; --accent:#7fb0ea; --hit:#2c4b70; }
* { box-sizing:border-box; }
body { background:var(--bg); color:var(--ink); margin:0; font:16px/1.55 "Charter","Iowan Old Style",Georgia,serif; }
main { max-width:1080px; margin:0 auto; padding:32px 16px 64px; }
h1,h2,h3 { font-family:-apple-system,"Helvetica Neue",sans-serif; line-height:1.25; }
h1 { font-size:30px; margin:0 0 4px; } h2 { font-size:21px; margin:40px 0 12px; } h3 { font-size:18px; margin:0 0 6px; }
a { color:var(--accent); } .muted { color:var(--muted); } img { max-width:100%; height:auto; }
.byline { color:var(--muted); margin:0 0 20px; font-family:-apple-system,sans-serif; font-size:15px; }
.lab { background:var(--card); border:1px solid var(--line); border-radius:8px; padding:20px; margin:0 0 20px; }
.lab-top { display:grid; grid-template-columns:minmax(0,1fr) minmax(0,340px); gap:20px; align-items:start; }
.scene { background:#fff; border-radius:6px; padding:6px; border:1px solid var(--line); }
.links { font-family:-apple-system,sans-serif; font-size:14px; margin:8px 0 0; }
.links a { margin-right:14px; }
details { margin-top:14px; } summary { cursor:pointer; font-family:-apple-system,sans-serif; font-size:15px; font-weight:600; }
.demos { display:grid; grid-template-columns:repeat(auto-fill,minmax(230px,1fr)); gap:12px; margin-top:10px; }
.demo { font-size:14px; } .demo img { background:#fff; border:1px solid var(--line); border-radius:4px; }
.demo code { font-size:12px; color:var(--muted); }
ol.hw { margin:8px 0 0; padding-left:22px; font-size:15px; } ol.hw .pts { color:var(--muted); white-space:nowrap; }
.map { overflow-x:auto; } table.map-t { border-collapse:collapse; font:13px -apple-system,sans-serif; }
.map-t th, .map-t td { border:1px solid var(--line); padding:3px 0; text-align:center; min-width:24px; }
.map-t th.row { text-align:left; padding:3px 8px; white-space:nowrap; font-weight:500; }
.map-t td.on { background:var(--hit); }
pre { background:var(--card); border:1px solid var(--line); border-radius:6px; padding:10px 12px; overflow-x:auto; font-size:13px; }
@media (max-width:720px) { .lab-top { grid-template-columns:1fr; } h1 { font-size:25px; } }
"""


def read(path):
    with open(path, encoding='utf-8') as f:
        return f.read()


def first_paragraph(md):
    """The overview: the first plain paragraph after the title (skips images, headings and tables)."""
    for block in re.split(r'\n\s*\n', md):
        b = block.strip()
        if b and not b.startswith(('#', '!', '|', '**The sample')):
            return ' '.join(b.split())
    return ''


def md_inline(text):
    """Escape, then the little markdown the overviews use: **bold**, `code`, $math$ shown as plain text."""
    t = html.escape(text)
    t = re.sub(r'\*\*(.+?)\*\*', r'<strong>\1</strong>', t)
    t = re.sub(r'`(.+?)`', r'<code>\1</code>', t)
    return t.replace('$', '')


def lecture_set(spec):
    out = set()
    for part in re.findall(r'\d+(?:-\d+)?', spec):
        a, _, b = part.partition('-')
        out.update(range(int(a), int(b or a) + 1))
    return out


def problems(md):
    """(title, points, lectures) for each problem; the three homework formats in the labs are all handled."""
    pats = [r'^## Problem \d+\. (?P<t>.+?) \((?P<p>\d+) points?(?:; (?P<l>[^)]*))?\)',
            r'^## \d+\. (?P<t>.+?) \((?P<p>\d+) pts\)(?: — (?P<l>.+))?$',
            r'^\d+\. \*\*(?P<t>.+?) \((?P<p>\d+) pts(?:, (?P<l>[^)]*))?\)\.?\*\*']
    out = []
    for line in md.splitlines():
        for pat in pats:
            m = re.match(pat, line)
            if m:
                out.append((m['t'], int(m['p']), (m['l'] or '').strip()))
                break
    return out


def demo_title(path):
    for line in read(path).splitlines():
        if line.startswith('# # '):
            return re.sub(r'^Lab \d+, demo \d+: ', '', line[4:]).strip()
    return os.path.basename(path)


def to_html(md_path, title):
    """README.html / homework.html next to the .md, so its relative figure links keep working."""
    if not shutil.which('pandoc'):
        return None
    out = md_path[:-3] + '.html'
    with open(os.path.join(HERE, 'docs', 'page.css'), 'w') as f:
        f.write(CSS + 'body{max-width:860px;margin:0 auto;padding:24px 16px 64px}'
                ' table{border-collapse:collapse} td,th{border:1px solid var(--line);padding:4px 8px}'
                ' header#title-block-header .title{font-size:26px}\n')
    # -V pagetitle keeps the title in the tab without a second heading (the .md already starts with one)
    subprocess.run(['pandoc', md_path, '-f', 'markdown-implicit_figures', '-s', '--math-method=mathjax', '-V', 'pagetitle=' + title,
                    '-c', '../../docs/page.css', '-o', out], check=True)
    return os.path.relpath(out, HERE)


# the lecture map lives in the table of the top-level README
top = read(os.path.join(HERE, 'README.md'))
lect = {m[1]: lecture_set(m[2]) for m in re.finditer(r'^\| \[(\d\d) [^]]+\]\([^)]*\) \|.*\| ([\d\-, ]+) \|$', top, re.M)}

cards, rows = [], []
for d in sorted(glob.glob(os.path.join(HERE, 'labs', '[0-9][0-9]_*'))):
    lab = os.path.basename(d)
    num = lab[:2]
    readme, hwmd = os.path.join(d, 'README.md'), os.path.join(d, 'homework.md')
    if not (os.path.exists(readme) and os.path.exists(hwmd)):     # a lab still being written
        continue
    rmd, hw = read(readme), read(hwmd)
    title = re.search(r'^# (.+)$', rmd, re.M)[1]
    short = re.sub(r'^Lab \d+: ', '', title)
    hw_title = re.search(r'^# (.+)$', hw, re.M)[1]
    hw_short = re.sub(r'^(Homework, )?Lab \d+( homework)?: ', '', hw_title)
    readme_link = to_html(readme, title) or os.path.relpath(readme, HERE)
    hw_link = to_html(hwmd, hw_title) or os.path.relpath(hwmd, HERE)

    demos = []
    for p in sorted(glob.glob(os.path.join(d, 'demo_*.py'))):
        stem = os.path.basename(p)[:-3]
        fig = os.path.join('out', '%s_%s_1.png' % (lab, stem))
        img = ('<img loading="lazy" src="%s" alt="first figure of %s">' % (fig, stem)
               if os.path.exists(os.path.join(HERE, fig)) else '')
        demos.append('<div class="demo">%s<div>%s</div><code>%s</code></div>'
                     % (img, md_inline(demo_title(p)), html.escape(os.path.relpath(p, HERE))))

    probs = ''.join('<li>%s <span class="pts">(%d pts%s)</span></li>'
                    % (md_inline(t), pts, ', ' + html.escape(l) if l else '') for t, pts, l in problems(hw))
    ls = lect.get(num, set())
    lec_text = ', '.join(str(n) for n in sorted(ls))
    cards.append(f"""
<section class="lab" id="lab{num}">
  <div class="lab-top">
    <div>
      <h3>Lab {num}: {md_inline(short)}</h3>
      <p>{md_inline(first_paragraph(rmd.split(chr(10), 1)[1]))}</p>
      <p class="muted" style="font-size:14px">Lectures {lec_text}</p>
      <p class="links"><a href="{readme_link}">Lab guide</a><a href="{hw_link}">Homework</a>
        <a href="{os.path.relpath(os.path.join(d, 'hw_starter.py'), HERE)}">hw_starter.py</a>
        <a href="notebooks/{lab}/">Notebooks</a><a href="{REPO}/tree/main/labs/{lab}">On GitHub</a></p>
    </div>
    <img class="scene" src="labs/{lab}/figures/scene.png" alt="Scene figure of lab {num}">
  </div>
  <details><summary>Demos for class ({len(demos)})</summary><div class="demos">{''.join(demos)}</div></details>
  <details><summary>Homework: {md_inline(hw_short)}</summary><ol class="hw">{probs}</ol></details>
</section>""")
    rows.append('<tr><th class="row"><a href="#lab%s">%s %s</a></th>%s</tr>'
                % (num, num, html.escape(short.split(':')[0].split(' (')[0]),
                   ''.join('<td class="on" title="L%d"></td>' % n if n in ls else '<td></td>' for n in LECTURES)))

head = ''.join('<th>%d</th>' % n for n in LECTURES)
page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Control in Practice</title><style>{CSS}</style></head><body><main>
<h1>Control in Practice</h1>
<p class="byline">Xu Chen · Introduction to Control Theory · practical demos and homework</p>
<p>Six labs take the course tools off the whiteboard. Each one starts from a logged experiment on a real kind of
machine and carries it through the whole toolchain: a physical model, a transfer function identified from data,
a root-locus or Bode design, a PID, lead or lag controller that has to live with saturation, noise, quantization
and delay, and a Nyquist or margin check of how robust the result is. Every lab has worked demos for class, a
100-point homework on a different unit or test than the demos, and a starter script.</p>
<p class="muted">Lab 04 uses a measurement you make with your own phone. The other logs are simulated from
physical models with realistic sensor noise, resolution, delays and actuator limits; each data file says so.</p>
<p>All lab materials (demos, homework sheets, data, starter scripts, notebooks and the shared tools) are in the
GitHub repository <a href="{REPO}">github.com/xchen-me/control-in-practice</a>; clone it or
<a href="{REPO}/archive/refs/heads/main.zip">download it as a zip file</a>.</p>

<h2>Which lectures each lab uses</h2>
<div class="map"><table class="map-t"><tr><th class="row">Lab \\ lecture</th>{head}</tr>{''.join(rows)}</table></div>

<h2>The labs</h2>
{''.join(cards)}

<h2>Getting started</h2>
<p>Open this page and the lab pages it links to directly in your web browser (double-click
<code>index.html</code> in Finder or File Explorer). Inside JupyterLab the links and figures of HTML pages do not
load; there, open the <code>README.md</code> of a lab instead, which JupyterLab shows as a formatted page.</p>
<p>Get the materials with <code>git clone {REPO}.git</code> or as the
<a href="{REPO}/archive/refs/heads/main.zip">zip file</a>, which unpacks as <code>control-in-practice-main</code>.
With Docker (recommended, the same environment for everyone): put the folder on a local disk, not inside
Google Drive, OneDrive or iCloud, then run</p>
<pre>cd control-in-practice*          # the folder from the clone or the zip
docker compose up --build        # then open http://localhost:8888 and see a lab README for the steps
docker compose run --rm check    # runs every demo and starter, figures in out/</pre>
<p>Without Docker (Python 3.11 or newer):</p>
<pre>python -m venv .venv &amp;&amp; source .venv/bin/activate
python -m pip install -r requirements.txt
export PYTHONPATH=$PWD
python labs/01_thermal/demo_1_identify.py</pre>
<p>The shared tools in <code>ctrlprac/</code> cover identification (<code>fit_fopdt</code>,
<code>fit_second_order_step</code>, <code>log_decrement</code>, <code>etfe</code>, <code>fit_tf</code>),
simulation of a sampled-data loop with saturation, noise, quantization and delay (<code>PID</code>,
<code>LTIController</code>, <code>simulate_loop</code>), and loading the logs. Full instructions:
<a href="README.md">README.md</a>.</p>
</main></body></html>
"""
with open(os.path.join(HERE, 'index.html'), 'w', encoding='utf-8') as f:
    f.write(page)
print(os.path.join(HERE, 'index.html'))
