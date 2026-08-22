"""Add page breaks to pandoc-generated docx."""
from docx import Document
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

def add_page_break(paragraph):
    """Insert a page break before a paragraph."""
    run = paragraph.insert_paragraph_before().add_run()
    br = OxmlElement('w:br')
    br.set(qn('w:type'), 'page')
    run._r.append(br)

doc = Document('/Users/cary/cea-intismeran/output/manuscript.docx')

# Find target paragraphs
abstract_idx = None
intro_idx = None

for i, p in enumerate(doc.paragraphs):
    text = p.text.strip()
    # Title page break: before Abstract heading
    if p.style.name == 'Heading 2' and text == 'Abstract':
        abstract_idx = i
    # Abstract to body: before 1. Introduction
    if p.style.name == 'Heading 2' and text.startswith('1. Introduction'):
        intro_idx = i

# Insert page breaks
# Remove duplicate "Yichao Jin, Ph.D." from YAML title (paragraph 1)
# Actually, keep it as is - it's fine

# Add page break before Abstract
if abstract_idx and abstract_idx > 0:
    add_page_break(doc.paragraphs[abstract_idx])
    print(f"Page break added before paragraph {abstract_idx} (Abstract)")

# Add page break before 1. Introduction
if intro_idx:
    add_page_break(doc.paragraphs[intro_idx])
    print(f"Page break added before paragraph {intro_idx} (1. Introduction)")

doc.save('/Users/cary/cea-intismeran/output/manuscript_pandoc.docx')
print("Done - page breaks added")