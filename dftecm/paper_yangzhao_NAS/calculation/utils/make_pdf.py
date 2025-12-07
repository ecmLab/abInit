#!/usr/bin/env python3
"""Minimal dependency-free text-to-PDF renderer."""

import math
import textwrap
from pathlib import Path
import sys

PAGE_WIDTH = 612.0  # Letter size (points)
PAGE_HEIGHT = 792.0
MARGIN = 54.0
LINE_HEIGHT = 14.0
FONT_SIZE = 11
CHARS_PER_LINE = 90


def sanitize(text: str) -> str:
    return text.replace('\\', r'\\').replace('(', r'\(').replace(')', r'\)')


def paginate(lines):
    usable_height = PAGE_HEIGHT - 2 * MARGIN
    lines_per_page = max(1, int(usable_height // LINE_HEIGHT))
    wrapped_lines = []
    wrapper = textwrap.TextWrapper(width=CHARS_PER_LINE, replace_whitespace=False, drop_whitespace=False)
    for ln in lines:
        ln = ln.rstrip('\n')
        if not ln:
            wrapped_lines.append('')
            continue
        wrapped_lines.extend(wrapper.wrap(ln) or [''])
    pages = [wrapped_lines[i:i + lines_per_page] for i in range(0, len(wrapped_lines), lines_per_page)]
    if not pages:
        pages = [[]]
    return pages


def make_content_stream(page_lines):
    cmds = ['BT', f'/F1 {FONT_SIZE} Tf', f'1 0 0 1 {MARGIN:.2f} {PAGE_HEIGHT - MARGIN:.2f} Tm']
    for idx, line in enumerate(page_lines):
        cmds.append(f'({sanitize(line)}) Tj')
        if idx != len(page_lines) - 1:
            cmds.append(f'0 {-LINE_HEIGHT:.2f} Td')
    cmds.append('ET')
    stream = '\n'.join(cmds).encode('latin-1', errors='ignore')
    header = f"<< /Length {len(stream)} >>\nstream\n".encode('ascii')
    return header + stream + b"\nendstream\n"


class PdfBuilder:
    def __init__(self):
        self.objects = []  # list of bytes

    def add_object(self, body: bytes) -> int:
        self.objects.append(body)
        return len(self.objects)

    def set_object(self, obj_id: int, body: bytes) -> None:
        self.objects[obj_id - 1] = body

    def write(self, path: Path, root_obj_id: int):
        with open(path, 'wb') as fh:
            fh.write(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
            offsets = []
            for idx, body in enumerate(self.objects, start=1):
                offsets.append(fh.tell())
                fh.write(f"{idx} 0 obj\n".encode('ascii'))
                fh.write(body)
                if not body.endswith(b"\n"):
                    fh.write(b"\n")
                fh.write(b"endobj\n")
            xref_pos = fh.tell()
            fh.write(f"xref\n0 {len(self.objects)+1}\n0000000000 65535 f \n".encode('ascii'))
            for off in offsets:
                fh.write(f"{off:010d} 00000 n \n".encode('ascii'))
            fh.write(b"trailer\n")
            fh.write(f"<< /Size {len(self.objects)+1} /Root {root_obj_id} 0 R >>\n".encode('ascii'))
            fh.write(f"startxref\n{xref_pos}\n%%EOF".encode('ascii'))


def render_pdf(lines, out_path: Path):
    pages = paginate(lines)
    pdf = PdfBuilder()
    font_obj = pdf.add_object(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>")
    pages_obj_id = pdf.add_object(b"")  # placeholder
    page_obj_ids = []
    content_obj_ids = []
    for page_lines in pages:
        content = make_content_stream(page_lines)
        content_id = pdf.add_object(content)
        page_dict = f"<< /Type /Page /Parent {pages_obj_id} 0 R /MediaBox [0 0 {PAGE_WIDTH:.2f} {PAGE_HEIGHT:.2f}] /Contents {content_id} 0 R /Resources << /Font << /F1 {font_obj} 0 R >> >> >>".encode('ascii')
        page_id = pdf.add_object(page_dict)
        page_obj_ids.append(page_id)
        content_obj_ids.append(content_id)
    kids_str = ' '.join(f"{pid} 0 R" for pid in page_obj_ids)
    pages_dict = f"<< /Type /Pages /Kids [{kids_str}] /Count {len(page_obj_ids)} >>".encode('ascii')
    pdf.set_object(pages_obj_id, pages_dict)
    catalog_obj = pdf.add_object(f"<< /Type /Catalog /Pages {pages_obj_id} 0 R >>".encode('ascii'))
    pdf.write(out_path, catalog_obj)


def main():
    if len(sys.argv) != 3:
        print("Usage: make_pdf.py input.txt output.pdf")
        sys.exit(1)
    inp = Path(sys.argv[1])
    out = Path(sys.argv[2])
    lines = inp.read_text(encoding='utf-8', errors='ignore').splitlines()
    render_pdf(lines, out)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()

