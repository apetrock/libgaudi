"""Cached text for PDFs in this directory.

`pages(path)` extracts once and stores the result in
`reference/.cache/<name>.json`, keyed by size and mtime. Later calls load
that file. `search(path, keys)` prints a short window around each hit.

`chunk_toc(path)` writes one markdown file per table-of-contents entry
under `reference/.cache/<stem>/`, plus `INDEX.md`. Open the index, then
the one section file. The cache is gitignored.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

CACHE_DIR = Path(__file__).resolve().parent / ".cache"


def _cache_path(pdf: Path) -> Path:
    return CACHE_DIR / f"{pdf.name}.json"


def pages(pdf: str | Path) -> list[str]:
    path = Path(pdf).resolve()
    st = path.stat()
    stamp = {"size": st.st_size, "mtime_ns": st.st_mtime_ns}
    cache = _cache_path(path)
    if cache.is_file():
        data = json.loads(cache.read_text())
        if data.get("size") == stamp["size"] and data.get("mtime_ns") == stamp["mtime_ns"]:
            return list(data["pages"])

    from pypdf import PdfReader

    reader = PdfReader(str(path))
    text = [(page.extract_text() or "") for page in reader.pages]
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps({**stamp, "pages": text}))
    return text


def search(pdf: str | Path, keys: tuple[str, ...] | list[str], window: int = 200) -> None:
    for i, page in enumerate(pages(pdf), start=1):
        low = page.lower()
        if not any(k.lower() in low for k in keys):
            continue
        print(f"--- page {i} hits ---")
        for key in keys:
            j = low.find(key.lower())
            if j < 0:
                continue
            snippet = page[max(0, j - 80) : j + window].replace("\n", " ")
            print(snippet)
            print()


def _slug(title: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", title.lower()).strip("-")
    return slug[:64] or "section"


def _unwrap(text: str) -> str:
    """Join line-break hyphenation. Blank lines stay paragraph breaks."""
    out: list[str] = []
    buf = ""
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            if buf:
                out.append(buf)
                buf = ""
            if not out or out[-1] != "":
                out.append("")
            continue
        if buf.endswith("-") and line[:1].islower():
            buf = buf[:-1] + line
        elif buf:
            buf = buf + " " + line
        else:
            buf = line
    if buf:
        out.append(buf)
    return "\n".join(out).strip()


def chunk_toc(pdf: str | Path) -> Path:
    """Write one markdown file per TOC section. Returns the index path."""
    import fitz

    path = Path(pdf).resolve()
    doc = fitz.open(path)
    toc = doc.get_toc()
    if not toc:
        raise RuntimeError(f"no table of contents: {path.name}")
    n_pages = doc.page_count
    # Next entry owns the shared start page, so a heading with no pages of
    # its own is skipped.
    spans: list[tuple[str, int, int]] = []
    for i, (_level, title, start) in enumerate(toc):
        end = toc[i + 1][2] - 1 if i + 1 < len(toc) else n_pages
        if end < start:
            continue
        spans.append((title.strip(), start, end))

    out_dir = CACHE_DIR / path.stem
    out_dir.mkdir(parents=True, exist_ok=True)
    for old in out_dir.glob("*.md"):
        old.unlink()

    rows: list[str] = []
    for n, (title, start, end) in enumerate(spans):
        parts = [_unwrap(doc[p - 1].get_text() or "") for p in range(start, end + 1)]
        body = "\n\n".join(part for part in parts if part)
        if not body:
            continue
        name = f"{n:02d}-{_slug(title)}.md"
        page = f"{start}" if start == end else f"{start}–{end}"
        (out_dir / name).write_text(
            f"# {title}\n\n"
            f"- source: `{path.name}`\n"
            f"- pages: {page}\n\n"
            f"{body}\n"
        )
        rows.append(f"| {page} | {title} | [{name}]({name}) |")

    index = out_dir / "INDEX.md"
    index.write_text(
        f"# {path.stem}\n\n"
        f"Section text from `{path.name}`. Open one file.\n\n"
        f"| pages | section | file |\n"
        f"|------:|---------|------|\n"
        + "\n".join(rows)
        + "\n"
    )
    return index


def main(argv: list[str]) -> None:
    if len(argv) == 3 and argv[1] == "chunk":
        print(chunk_toc(argv[2]))
        return
    raise SystemExit("usage: pdf_text.py chunk <pdf>")


if __name__ == "__main__":
    main(sys.argv)
