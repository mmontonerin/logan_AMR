from pathlib import Path
import shutil

src = Path("../data/SRA_metadata_before20231211_logan.csv")
tmp = src.with_suffix(".tmp")    # temporary output
bak = src.with_suffix(".bak")    # keep original as backup

with src.open("r", encoding="utf-8") as fin, tmp.open("w", encoding="utf-8", newline="") as fout:
    for line_num, line in enumerate(fin, start=1):
        # 1) fix trailing comma
        if line.endswith(",\n"):
            line = line[:-2] + "\n"
        # 2) fix unmatched quote
        if line.count('"') % 2:          # odd number → unbalanced
            line = line.replace('"', "") # drop ALL quotes in that row
        fout.write(line)

# atomically replace original (keep a backup just in case)
shutil.move(src, bak)
shutil.move(tmp, src)
print("✔ Clean copy written. Original saved as", bak.name)

