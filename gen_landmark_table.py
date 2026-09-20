"""Regenerate the Observed-vs-Fitted landmark table (supplementary.tex, tab:s15).

Every fitted value is read through model.survival_at_month(), so array index ==
month for all three endpoints and both arms.

Usage:
    python gen_landmark_table.py            # print the new rows (dry run)
    python gen_landmark_table.py --write    # rewrite the table in supplementary.tex
"""
import re
import sys

from model import ModelParams, survival_at_month
from rerun_primary import CEAModelV2

TEX = "supplementary.tex"
CURVE = {"OS": "os", "DMFS": "dmfs", "RFS": "rfs"}
TABLE_START = "\\label{tab:s15}"


def table_span(tex: str) -> tuple[int, int]:
    start = tex.index(TABLE_START)
    return start, tex.index("\\end{table}", start)


def parse_rows(chunk: str) -> list[tuple]:
    body = chunk.split("\\midrule\n", 1)[1].split("\\midrule\n", 1)[0]
    rows = []
    for line in body.strip().split("\\\\\n"):
        cells = [c.strip() for c in line.split("&")]
        if len(cells) == 8:
            rows.append((cells[0], cells[1], int(cells[2]), float(cells[3].rstrip("\\%"))))
    return rows


def main() -> int:
    with open(TEX) as fh:
        tex = fh.read()
    start, end = table_span(tex)
    chunk = tex[start:end]
    rows = parse_rows(chunk)

    sur = CEAModelV2(ModelParams(constraint_general_pop=True))._survival()
    lines, residuals = [], []
    for endpoint, arm, month, observed in rows:
        fitted = 100 * survival_at_month(sur[f"{CURVE[endpoint]}_{arm.lower()}"], month)
        res = fitted - observed
        residuals.append(res)
        lines.append(f"{endpoint} & {arm} & {month} & {observed:.1f}\\% & {fitted:.1f}\\% & "
                     f"{res:+.1f}\\% & {abs(res):.1f}\\% & {abs(res) / observed * 100:.1f}\\% \\\\")

    mae = sum(abs(r) for r in residuals) / len(residuals)
    rmse = (sum(r * r for r in residuals) / len(residuals)) ** 0.5
    worst_i = max(range(len(residuals)), key=lambda i: abs(residuals[i]))
    worst = abs(residuals[worst_i])
    print("\n".join(lines))
    print(f"\nMAE {mae:.2f} | RMSE {rmse:.2f} | max |res| {worst:.1f} "
          f"({rows[worst_i][0]} {rows[worst_i][1]} {rows[worst_i][2]}m)")

    if "--write" not in sys.argv:
        return 0

    old_body = chunk.split("\\midrule\n", 1)[1].split("\\midrule\n", 1)[0]
    new_chunk = chunk.replace(old_body, "\n".join(lines) + "\n")
    new_chunk = re.sub(r"(\\textbf\{MAE\}\} & \\multicolumn\{4\}\{l\}\{\\textbf\{)[\d.]+", rf"\g<1>{mae:.2f}", new_chunk)
    new_chunk = re.sub(r"(\\textbf\{RMSE\}\} & \\multicolumn\{4\}\{l\}\{\\textbf\{)[\d.]+", rf"\g<1>{rmse:.2f}", new_chunk)
    new_chunk = re.sub(r"(\\textbf\{Max \$\\\|Residual\\\|\$\}\} & \\multicolumn\{4\}\{l\}\{\\textbf\{)[\d.]+", rf"\g<1>{worst:.1f}", new_chunk)
    with open(TEX, "w") as fh:
        fh.write(tex[:start] + new_chunk + tex[end:])
    print("wrote supplementary.tex")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
