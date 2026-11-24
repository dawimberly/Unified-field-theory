"""Builds a small synthetic ChatGPT export (zip with numbered conversation files + an image) and checks the converter."""
import json, os, shutil, sys, tempfile, zipfile
sys.path.insert(0, os.path.dirname(__file__))
import convert_chatgpt_export as C

def node(i, parent, children, role=None, parts=None, t=None, ct="text", extra=None, hidden=False):
    msg = None
    if role:
        content = {"content_type": ct}
        if parts is not None: content["parts"] = parts
        if extra: content.update(extra)
        msg = {"id": i, "author": {"role": role}, "create_time": t, "content": content,
               "metadata": {"is_visually_hidden_from_conversation": True} if hidden else {}}
    return i, {"id": i, "message": msg, "parent": parent, "children": children}

def conv_physics():
    m = dict([
        node("root", None, ["sys"]),
        node("sys", "root", ["u1"], "system", [""], hidden=True),
        node("u1", "sys", ["a1"], "user", ["Is the hydrogen orbital a slice of a 4D hypersphere?"], 1700000000),
        node("a1", "u1", ["u2_old", "u2_new"], "assistant", ["Fock showed hydrogen has SO(4) symmetry."], 1700000010),
        node("u2_old", "a1", ["a2_old"], "user", ["OLD question that was edited"], 1700000020),
        node("a2_old", "u2_old", [], "assistant", ["OLD answer"], 1700000030),
        node("u2_new", "a1", ["th"], "user", [{"content_type": "image_asset_pointer", "asset_pointer": "file-service://file-ABC123", "width": 640, "height": 480},
                                              "Here is my sketch of the E8 lattice."], 1700000040, ct="multimodal_text"),
        node("th", "u2_new", ["code"], "assistant", None, 1700000045, ct="thoughts", extra={"thoughts": [{"summary": "hidden"}]}),
        node("code", "th", ["out"], "assistant", None, 1700000050, ct="code", extra={"language": "python", "text": "print(240*4/3)"}),
        node("out", "code", ["a2"], "tool", None, 1700000055, ct="execution_output", extra={"text": "320.0"}),
        node("a2", "out", [], "assistant", ["The E8 root system has 240 roots; alpha = 320."], 1700000060),
    ])
    return {"title": "Hypersphere | E8 & gravity", "create_time": 1700000000.5, "update_time": 1700000060, "mapping": m,
            "current_node": "a2", "conversation_id": "c-phys"}

def conv_cooking():   # no conversation timestamps: falls back to message time
    m = dict([node("r", None, ["u"]), node("u", "r", ["a"], "user", ["Best lasagna recipe?"], 1690000000),
              node("a", "u", [], "assistant", ["Layer pasta, sauce, cheese."], 1690000005)])
    return {"title": "Lasagna", "create_time": None, "update_time": None, "mapping": m, "current_node": "a"}

def conv_undated():   # no timestamps anywhere, no current_node (fallback: last child), unicode title
    m = dict([node("r", None, ["u"]), node("u", "r", ["a"], "user", ["Does quantum gravity need 13 dimensions?"]),
              node("a", "u", [], "assistant", ["Bars' 2T physics uses 11+2."])])
    return {"title": "Quántum grávity — 13D?", "mapping": m}

def conv_dup():       # same title and date as physics chat -> slug collision
    c = conv_physics(); c["conversation_id"] = "c-dup"; return c

def conv_empty():
    return {"title": "empty", "create_time": 1700000000, "mapping": {}, "current_node": None}

def build(tmp):
    d = os.path.join(tmp, "export"); os.makedirs(d)
    json.dump([conv_physics(), conv_cooking()], open(os.path.join(d, "conversations-000.json"), "w"))
    json.dump([conv_undated(), conv_dup(), conv_empty()], open(os.path.join(d, "conversations-001.json"), "w"))
    open(os.path.join(d, "file-ABC123-sketch.png"), "wb").write(b"\x89PNG fake")
    zp = os.path.join(tmp, "chatgpt-export.zip")
    with zipfile.ZipFile(zp, "w") as z:
        for f in os.listdir(d): z.write(os.path.join(d, f), "export-2026-11-24/" + f)
    return d, zp

def run():
    tmp = tempfile.mkdtemp()
    d, zp = build(tmp)
    ok = True
    def check(cond, msg):
        nonlocal ok; print(("PASS " if cond else "FAIL ") + msg); ok &= bool(cond)
    # 1. full conversion from zip, all branches, assets
    out1 = os.path.join(tmp, "all")
    rows = C.convert(zp, out1, None, include_tools=True, all_branches=True, copy_assets=True)
    files = sorted(os.listdir(out1)); print("  files:", files)
    check(len(rows) == 4, "4 non-empty conversations written (empty one skipped)")
    phys = open(os.path.join(out1, "2023-11-14_hypersphere-e8-gravity.md")).read()
    check("OLD question" not in phys.split("Alternate branch")[0], "main branch follows current_node (edited turn excluded)")
    check("Alternate branch" in phys and "OLD answer" in phys, "--all-branches appends the edited branch")
    check("![image 640x480](assets/file-ABC123)" in phys, "multimodal image part rendered")
    check(os.path.exists(os.path.join(out1, "assets", "file-ABC123")), "--copy-assets copied the image from the zip")
    check("```python\nprint(240*4/3)\n```" in phys and "320.0" in phys, "code + execution output rendered")
    check("hidden" not in phys, "thoughts / hidden system message skipped")
    check("2023-11-14_hypersphere-e8-gravity-2.md" in files, "slug collision resolved")
    check("2023-07-22_lasagna.md" in files, "missing conversation timestamps fall back to message time")
    check("undated_quantum-gravity-13d.md" in files, "no timestamps at all -> undated_, unicode slugified")
    idx = open(os.path.join(out1, "INDEX.md")).read()
    check("Hypersphere \\| E8 & gravity" in idx, "INDEX escapes '|' in titles")
    # 2. physics filter from extracted folder
    out2 = os.path.join(tmp, "phys")
    rows2 = C.convert(d, out2, C.DEFAULT_KEYWORDS)
    check(sorted(r[2] for r in rows2) == ["2023-11-14_hypersphere-e8-gravity", "2023-11-14_hypersphere-e8-gravity-2", "undated_quantum-gravity-13d"],
          "default physics filter keeps 3, drops lasagna")
    # 3. custom filter via CLI on a single json file
    out3 = os.path.join(tmp, "cli")
    C.main([os.path.join(d, "conversations-001.json"), "-o", out3, "--filter", "13D", "Bars"])
    check(os.listdir(out3) and "undated_quantum-gravity-13d.md" in os.listdir(out3), "CLI single-file + custom keywords")
    print("\n--- sample output (physics chat, main branch) ---\n" + phys.split("---\n\n<details>")[0][:1200])
    print("--- INDEX.md ---\n" + idx)
    shutil.copytree(out1, os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_output"), dirs_exist_ok=True)
    print("ALL PASS" if ok else "SOME FAILED")
    return ok

if __name__ == "__main__":
    sys.exit(0 if run() else 1)
