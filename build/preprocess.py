"""
Preprocess manuscript.md for pandoc conversion.
"""
import re

with open("manuscript.md") as f:
    text = f.read()

# 1. Remove title line
text = re.sub(r"^# Cost-Effectiveness[^\n]*\n", "", text)

# 2. Convert figure references to pandoc image syntax
def replace_figure(m):
    num = m.group(1)
    caption = m.group(2).strip()
    path = m.group(3)
    return f"![{caption}]({path})"

text = re.sub(
    r'-\s+\*\*Figure\s+(\d+)[:\*]?\s*(.*?)\s*→\s*`([^`]+)`',
    replace_figure, text
)

# 3. Remove ## Figures and ## Tables headings
text = re.sub(r"\n## Figures\n", "\n", text)
text = re.sub(r"\n## Tables\n", "\n", text)

# 4. Extract tables and embed at first citation
tables = {}
table_pattern = re.compile(r'### (Table \d+\..*?)\n(.*?)(?=\n### |\n---\n|\Z)', re.DOTALL)
for m in table_pattern.finditer(text):
    title = m.group(1)
    content = m.group(2)
    tnum = int(re.search(r'Table (\d+)', title).group(1))
    tables[tnum] = f"**{title}**\n\n{content.strip()}"

# Remove all table definitions
for tnum in sorted(tables.keys(), reverse=True):
    text = re.sub(rf'### Table {tnum}\..*?(?=\n### |\n---\n|\Z)', '', text, count=1, flags=re.DOTALL)

# Insert tables at first citation
for tnum in sorted(tables.keys()):
    m = re.search(rf'\bTable\s+{tnum}\b', text)
    if m:
        rest = text[m.end():]
        eop = re.search(r'\n\n', rest)
        insert_at = m.end() + (eop.end() if eop else 0)
        text = text[:insert_at] + "\n\n" + tables[tnum] + "\n\n" + text[insert_at:]

# 5. Clean up separators and blank lines
text = re.sub(r'\n---\n', '\n\n', text)
text = re.sub(r'\n{3,}', '\n\n', text)

# 6. Add page breaks with raw_tex
# After title page info (before Abstract)
# Insert a blank line + \newpage between the title page info and Abstract
text = re.sub(r'\n\n## Abstract\n', '\n\n\\newpage\n\n## Abstract\n', text, count=1)

# After Abstract (before 1. Introduction)
text = re.sub(r'\nKeywords:[^\n]+\n\n', '\n\\newpage\n\n', text, count=1)

with open("build/manuscript_pandoc.md", "w") as f:
    f.write(text)

print("Done")