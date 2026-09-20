"""Build data/life_table_surv.json from the NVSS US Life Tables, 2019.

Source: Arias E, Xu J. United States Life Tables, 2019. National Vital Statistics
Reports; vol 70 no 19. Hyattsville, MD: NCHS; 2022. Tables 1-3 (total, male, female),
shipped as data/nvss/Table0{1,2,3}.xlsx.

Convention used by model.py / rerun_primary.py:
    GP survival conditional on being alive at baseline age 61,
    sex-weighted 64% male / 36% female,
    annual values at every year 0..42, log-linearly interpolated to months downstream.
Ages beyond the table's last row (100) hold the terminal value constant.

Output: {"survival_by_year": {"0": 1.0, "1": ..., "42": ...}}

Usage:  python build_life_table.py [out_path]
"""
import json
import os
import sys

import numpy as np
import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
W_MALE = 0.64
BASE_AGE = 61
N_YEARS = 43
TABLES = ("01", "02", "03")  # total, male, female


def lx(table: str) -> dict[int, float]:
    """l_x by age from a NVSR 70-19 life-table spreadsheet."""
    path = os.path.join(HERE, "data", "nvss", f"Table{table}.xlsx")
    ws = openpyxl.load_workbook(path, data_only=True).active
    out = {}
    for row in ws.iter_rows(min_row=4, max_row=104, values_only=True):
        age, l = str(row[0]).strip(), row[2]
        if l is None or age.startswith("SOURCE"):
            continue
        if age == "100 and over":
            out[100] = float(l)
        else:
            out[int(age.split("–")[0])] = float(l)
    return out


def main() -> int:
    male, female = lx("02"), lx("03")

    def weighted(age: int) -> float:
        return (W_MALE * male.get(age, male[100])
                + (1 - W_MALE) * female.get(age, female[100]))

    base = weighted(BASE_AGE)
    surv = {str(y): round(weighted(BASE_AGE + y) / base, 6) for y in range(N_YEARS)}

    out_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "data", "life_table_surv.json")
    with open(out_path, "w") as fh:
        json.dump({"survival_by_year": surv}, fh)
    print(f"wrote {out_path}")

    # self-check: reproduce the GP column of Supplementary Table S5 (mid-interval read)
    years = np.arange(N_YEARS)
    ann = np.array([surv[str(y)] for y in years])
    published = {5: 94.1, 10: 86.4, 15: 76.1, 20: 61.9, 30: 23.1, 40: 2.0}
    worst = 0.0
    for t, want in published.items():
        got = 100 * float(np.exp(np.interp(t + 0.5 / 12, years, np.log(ann))))
        worst = max(worst, abs(got - want))
        print(f"  year {t:>2}: {got:5.2f}  (published {want:5.1f})")
    assert worst < 0.12, f"GP column mismatch: max {worst:.2f} pp"
    print(f"  self-check OK (max {worst:.2f} pp)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
