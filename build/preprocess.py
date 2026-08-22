"""
Preprocess manuscript.md for pandoc: embed tables AND figures at first citation.
"""
import re, shutil

with open("manuscript.md") as f:
    text = f.read()

# 1. Remove title line
text = re.sub(r"^# Cost-Effectiveness[^\n]*\n", "", text)

# 2. Extract all figures (convert to pandoc ![] syntax), store by number
figures = {}
def extract_figure(m):
    num = int(m.group(1))
    caption = m.group(2).strip()
    path = m.group(3)
    figures[num] = f"![]({path})"  # no caption text in the image syntax
    # Return the caption as text (will be placed before the image)
    return f"**Figure {num}:** {caption}"

text = re.sub(
    r'-\s+\*\*Figure\s+(\d+)[:\*]?\s*(.*?)\s*→\s*`([^`]+)`',
    extract_figure, text
)

# Remove ## Figures heading and any remaining empty lines
text = re.sub(r"\n## Figures\n", "\n", text)
text = re.sub(r'^-\s*\n', '', text, flags=re.MULTILINE)

# 3. Extract all tables, store by number
tables = {}
# Save tables first
table_pattern = re.compile(r'### (Table \d+\..*?)\n(.*?)(?=\n### |\n---\n|\Z)', re.DOTALL)
for m in table_pattern.finditer(text):
    title = m.group(1)
    content = m.group(2)
    tnum = int(re.search(r'Table (\d+)', title).group(1))
    tables[tnum] = f"**{title}**\n\n{content.strip()}"

# Remove all table definitions from the text
for tnum in sorted(tables.keys(), reverse=True):
    text = re.sub(rf'### Table {tnum}\..*?(?=\n### |\n---\n|\Z)', '', text, count=1, flags=re.DOTALL)

# Remove ## Tables heading
text = re.sub(r"\n## Tables\n", "\n", text)

# 4. Insert tables at first citation
for tnum in sorted(tables.keys()):
    m = re.search(rf'\bTable\s+{tnum}\b', text)
    if m:
        rest = text[m.end():]
        eop = re.search(r'\n\n', rest)
        insert_at = m.end() + (eop.end() if eop else 0)
        text = text[:insert_at] + "\n\n" + tables[tnum] + "\n\n" + text[insert_at:]

# 5. Insert figures at first citation (process in reverse order to preserve positions)
for fnum in sorted(figures.keys(), reverse=True):
    m = re.search(rf'\bFigure\s+{fnum}\b', text)
    if m:
        rest = text[m.end():]
        eop = re.search(r'\n\n', rest)
        insert_at = m.end() + (eop.end() if eop else 0)
        img = figures[fnum]
        text = text[:insert_at] + "\n\n" + img + "\n\n" + text[insert_at:]

# 6. Clean up
text = re.sub(r'\n---\n', '\n\n', text)
text = re.sub(r'\n{3,}', '\n\n', text)

# 7. Add page break markers
text = re.sub(r'\n\n## Abstract\n', '\n\n\\newpage\n\n## Abstract\n', text, count=1)
# After Keywords (end of abstract), before 1. Introduction
text = re.sub(r'\nKeywords:[^\n]+\n\n', '\n\\newpage\n\n', text, count=1)

with open("build/manuscript_pandoc.md", "w") as f:
    f.write(text)

# Verify
fig_count = len(re.findall(r'!\[\]\(output/fig', text))
tab_count = len(re.findall(r'\*\*Table \d+\.\*\*', text))
print(f"Done. Figures embedded: {fig_count}/{len(figures)}, Tables: {tab_count}/{len(tables)}")
# Show where figures are
for fnum in sorted(figures.keys()):
    pos = text.find(figures[fnum])
    context = text[max(0,pos-50):pos+80]
    print(f"  Figure {fnum}: pos={pos} -> ...{context.strip()[:60]}...")