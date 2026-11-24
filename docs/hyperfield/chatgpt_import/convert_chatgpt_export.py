#!/usr/bin/env python3
"""Convert an official ChatGPT data export into one Markdown file per conversation plus INDEX.md.

Input:  the export .zip, the extracted folder, or a conversations.json file.
        Handles conversations.json and numbered files (conversations-000.json, ...), which OpenAI uses for large exports.
Output: <out>/YYYY-MM-DD_slug.md for each kept conversation, and <out>/INDEX.md.

Usage:
  python convert_chatgpt_export.py EXPORT.zip -o docs/chatgpt
  python convert_chatgpt_export.py EXPORT_FOLDER -o docs/chatgpt --filter            # default physics keywords
  python convert_chatgpt_export.py EXPORT.zip -o out --filter gravity E8 tesseract    # your own keywords
Options: --include-tools (show tool/code-execution messages), --all-branches (append alternate branches),
         --copy-assets (copy referenced images from the export next to the markdown).
Standard library only (Python 3.8+).
"""
import argparse, datetime as dt, io, json, os, re, shutil, sys, unicodedata, zipfile

DEFAULT_KEYWORDS = ["physics", "field", "gravity", "quantum", "dimension", "dimensional", "hypersphere", "hyperfield",
                    "E8", "unified", "tesseract", "brane", "lattice", "relativity", "Einstein", "particle", "proton",
                    "neutron", "electron", "photon", "gluon", "boson", "fermion", "spacetime", "string theory",
                    "Planck", "Kaluza", "Klein", "Higgs", "dark matter", "cosmolog", "wavefunction", "orbital",
                    "Schrodinger", "Schrödinger", "spin", "torus", "hyperspher", "4D", "13D", "11D", "M-theory"]

# ---------------------------------------------------------------- loading
def load_conversations(src):
    """Return (list_of_conversations, asset_lookup) where asset_lookup(name)->bytes or None."""
    convs = []
    if os.path.isdir(src):
        names = sorted(f for f in os.listdir(src) if re.fullmatch(r"conversations(-\d+)?\.json", f))
        if not names:  # search one level down (zip extracted into a subfolder)
            for d in sorted(os.listdir(src)):
                p = os.path.join(src, d)
                if os.path.isdir(p) and any(re.fullmatch(r"conversations(-\d+)?\.json", f) for f in os.listdir(p)):
                    return load_conversations(p)
            raise SystemExit(f"No conversations.json found in {src}")
        for n in names:
            with open(os.path.join(src, n), encoding="utf-8") as f:
                convs += _as_list(json.load(f))
        files = {}
        for root, _, fs in os.walk(src):
            for f in fs: files.setdefault(f, os.path.join(root, f))
        def asset(name):
            p = _find_asset(name, files)
            return open(p, "rb").read() if p else None
        return convs, asset
    if zipfile.is_zipfile(src):
        z = zipfile.ZipFile(src)
        names = sorted(n for n in z.namelist() if re.search(r"(^|/)conversations(-\d+)?\.json$", n))
        if not names: raise SystemExit(f"No conversations.json inside {src}")
        for n in names:
            convs += _as_list(json.loads(z.read(n).decode("utf-8")))
        files = {}
        for n in z.namelist(): files.setdefault(os.path.basename(n), n)
        def asset(name):
            p = _find_asset(name, files)
            return z.read(p) if p else None
        return convs, asset
    with open(src, encoding="utf-8") as f:
        return _as_list(json.load(f)), (lambda name: None)

def _as_list(obj):
    if isinstance(obj, list): return obj
    if isinstance(obj, dict) and "conversations" in obj: return obj["conversations"]
    if isinstance(obj, dict): return [obj]
    return []

def _find_asset(name, files):
    if not name: return None
    key = name.split("://")[-1].split("/")[-1]
    if key in files: return files[key]
    for f, p in files.items():          # export files are often named "<file-id>-<original name>"
        if f.startswith(key): return p
    return None

# ---------------------------------------------------------------- tree walking
def main_branch(conv):
    """Node ids from root to current_node (the branch the user last saw). Fallback: follow last child from the root."""
    mapping = conv.get("mapping") or {}
    cur = conv.get("current_node")
    path = []
    if cur in mapping:
        seen = set()
        while cur and cur in mapping and cur not in seen:
            seen.add(cur); path.append(cur); cur = mapping[cur].get("parent")
        return path[::-1]
    roots = [k for k, v in mapping.items() if not v.get("parent") or v.get("parent") not in mapping]
    if not roots: return []
    node = roots[0]
    while node:
        path.append(node)
        ch = [c for c in (mapping[node].get("children") or []) if c in mapping]
        node = ch[-1] if ch else None
    return path

def alternate_branches(conv, main):
    """Leaf-to-root paths of branches not on the main path (edited/regenerated turns)."""
    mapping = conv.get("mapping") or {}; mainset = set(main); out = []
    for k, v in mapping.items():
        if not (v.get("children")) and k not in mainset:
            p = []; cur = k
            while cur and cur in mapping and cur not in mainset:
                p.append(cur); cur = mapping[cur].get("parent")
            out.append((cur, p[::-1]))  # (fork point on main path, nodes)
    return out

# ---------------------------------------------------------------- rendering
ROLE = {"user": "User", "assistant": "ChatGPT", "tool": "Tool", "system": "System"}

def render_part(p, asset_names):
    if isinstance(p, str): return p
    if isinstance(p, dict):
        ct = p.get("content_type", "")
        if ct == "image_asset_pointer" or "asset_pointer" in p:
            ptr = p.get("asset_pointer", ""); fid = ptr.split("://")[-1]
            asset_names.append(fid)
            w, h = p.get("width"), p.get("height")
            size = f" {w}x{h}" if w and h else ""
            return f"![image{size}](assets/{fid})"
        if ct in ("audio_transcription",) or "text" in p and isinstance(p.get("text"), str):
            return p.get("text", "")
        if ct == "audio_asset_pointer": return "[audio clip]"
        return f"[{ct or 'attachment'}]"
    return str(p)

def message_text(msg, asset_names):
    c = msg.get("content") or {}
    ct = c.get("content_type", "text")
    if ct in ("text", "multimodal_text"):
        return "\n\n".join(t for t in (render_part(p, asset_names) for p in c.get("parts") or []) if t and t.strip())
    if ct == "code":
        lang = c.get("language") or ""; lang = "" if lang == "unknown" else lang
        return f"```{lang}\n{c.get('text', '')}\n```"
    if ct == "execution_output":
        return f"```\n{c.get('text', '')}\n```"
    if ct == "tether_quote":
        return f"> {c.get('title', '')} {c.get('url', '')}\n> " + (c.get("text", "") or "").replace("\n", "\n> ")
    if ct in ("thoughts", "reasoning_recap", "user_editable_context", "model_editable_context", "system_error"):
        return ""      # internal / settings content: skipped
    if "parts" in c:
        return "\n\n".join(t for t in (render_part(p, asset_names) for p in c["parts"]) if t)
    return c.get("text", "") or ""

def is_hidden(msg):
    md = msg.get("metadata") or {}
    return bool(md.get("is_visually_hidden_from_conversation") or md.get("is_user_system_message"))

def turns(conv, ids, include_tools=False):
    mapping = conv.get("mapping") or {}; out = []
    for i in ids:
        msg = (mapping.get(i) or {}).get("message")
        if not msg or is_hidden(msg): continue
        role = (msg.get("author") or {}).get("role", "unknown")
        if role == "system": continue
        if role == "tool" and not include_tools: continue
        names = []
        txt = message_text(msg, names).strip()
        if not txt: continue
        out.append(dict(role=role, text=txt, time=msg.get("create_time"), assets=names))
    return out

def conv_time(conv, tlist):
    for t in (conv.get("create_time"), conv.get("update_time")) + tuple(x["time"] for x in tlist):
        if isinstance(t, (int, float)) and t > 0:
            return dt.datetime.fromtimestamp(t, dt.timezone.utc)
    return None

def slugify(s, n=60):
    s = unicodedata.normalize("NFKD", s or "untitled").encode("ascii", "ignore").decode()
    s = re.sub(r"[^A-Za-z0-9]+", "-", s).strip("-").lower()
    return (s[:n].rstrip("-")) or "untitled"

def keyword_hits(text, keywords):
    hits = []
    for k in keywords:
        pat = r"(?<![A-Za-z0-9])" + re.escape(k)
        if len(k) <= 3: pat += r"(?![A-Za-z0-9])"
        if re.search(pat, text, flags=0 if len(k) <= 3 else re.IGNORECASE): hits.append(k)
    return hits

def fmt_time(t):
    return dt.datetime.fromtimestamp(t, dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC") if isinstance(t, (int, float)) and t > 0 else ""

def to_markdown(conv, tl, when, hits, alts=None):
    title = conv.get("title") or "Untitled"
    cid = conv.get("conversation_id") or conv.get("id") or ""
    L = [f"# {title}", "", f"- Date: {when.strftime('%Y-%m-%d %H:%M UTC') if when else 'unknown'}",
         f"- Conversation id: {cid}" if cid else "", f"- Turns: {len(tl)}",
         f"- Matched keywords: {', '.join(hits)}" if hits else "", "", "---", ""]
    for t in tl:
        stamp = fmt_time(t["time"])
        L += [f"### {ROLE.get(t['role'], t['role'].title())}" + (f"  <sub>{stamp}</sub>" if stamp else ""), "", t["text"], ""]
    for fork, nodes in (alts or []):
        L += ["---", "", f"<details><summary>Alternate branch (edited or regenerated)</summary>", ""]
        for t in nodes:
            L += [f"**{ROLE.get(t['role'], t['role'])}:**", "", t["text"], ""]
        L += ["</details>", ""]
    return "\n".join(x for x in L if x is not None) + "\n"

# ---------------------------------------------------------------- main
def convert(src, out_dir, keywords=None, include_tools=False, all_branches=False, copy_assets=False, quiet=False):
    convs, asset = load_conversations(src)
    os.makedirs(out_dir, exist_ok=True)
    rows, used, skipped = [], set(), 0
    for conv in convs:
        ids = main_branch(conv); tl = turns(conv, ids, include_tools)
        if not tl: skipped += 1; continue
        title = conv.get("title") or "Untitled"
        hits = keyword_hits(title + "\n" + "\n".join(t["text"] for t in tl), keywords) if keywords else []
        if keywords and not hits: skipped += 1; continue
        when = conv_time(conv, tl)
        base = f"{when.strftime('%Y-%m-%d') if when else 'undated'}_{slugify(title)}"
        name = base; k = 2
        while name in used: name = f"{base}-{k}"; k += 1
        used.add(name)
        alts = []
        if all_branches:
            for fork, nodes in alternate_branches(conv, ids):
                at = turns(conv, nodes, include_tools)
                if at: alts.append((fork, at))
        with open(os.path.join(out_dir, name + ".md"), "w", encoding="utf-8") as f:
            f.write(to_markdown(conv, tl, when, hits, alts))
        if copy_assets:
            for t in tl:
                for a in t["assets"]:
                    data = asset(a)
                    if data:
                        os.makedirs(os.path.join(out_dir, "assets"), exist_ok=True)
                        open(os.path.join(out_dir, "assets", a), "wb").write(data)
        rows.append((when, title, name, len(tl), hits))
    rows.sort(key=lambda r: (r[0] is None, r[0] or dt.datetime.min.replace(tzinfo=dt.timezone.utc)))
    with open(os.path.join(out_dir, "INDEX.md"), "w", encoding="utf-8") as f:
        f.write("# ChatGPT conversations\n\n")
        f.write(f"Source: `{os.path.basename(src)}` — {len(rows)} kept, {skipped} skipped"
                + (f" (keyword filter: {', '.join(keywords)})" if keywords else "") + ".\n\n")
        f.write("| Date | Title | Turns | Keywords |\n|---|---|---|---|\n")
        for when, title, name, n, hits in rows:
            t = title.replace("|", "\\|")
            f.write(f"| {when.strftime('%Y-%m-%d') if when else 'undated'} | [{t}]({name}.md) | {n} | {', '.join(hits[:6])} |\n")
    if not quiet:
        print(f"conversations in export: {len(convs)}; written: {len(rows)}; skipped: {skipped}; output: {out_dir}")
    return rows

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src", help="export .zip, extracted folder, or conversations.json")
    ap.add_argument("-o", "--out", default="docs/chatgpt", help="output folder (default docs/chatgpt)")
    ap.add_argument("--filter", nargs="*", default=None, metavar="KEYWORD",
                    help="keep only chats mentioning a keyword; with no words uses the default physics list")
    ap.add_argument("--include-tools", action="store_true"); ap.add_argument("--all-branches", action="store_true")
    ap.add_argument("--copy-assets", action="store_true")
    a = ap.parse_args(argv)
    kw = None if a.filter is None else (a.filter or DEFAULT_KEYWORDS)
    convert(a.src, a.out, kw, a.include_tools, a.all_branches, a.copy_assets)

if __name__ == "__main__":
    main()
