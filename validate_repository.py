"""Validate public file integrity and syntax without importing scientific modules."""
from pathlib import Path
import ast,hashlib,json,csv,sys
P=Path(__file__).resolve().parent
m=json.loads((P/'FILE_MANIFEST.json').read_text())
for r in m['files']:
 p=P/r['file'];assert p.is_file(),r['file'];b=p.read_bytes();assert len(b)==r['bytes'] and hashlib.sha256(b).hexdigest()==r['sha256'],r['file']
for p in P.rglob('*.py'):ast.parse(p.read_text(),filename=str(p))
for p in P.rglob('*'):
 if '.git' in p.parts:continue
 assert p.suffix.lower() not in {'.png','.jpg','.jpeg','.gif','.tif','.tiff','.webp','.pdf','.svg','.ai','.eps','.docx','.pptx'},p
for p in (P/'results').glob('*.tsv'):
 with p.open() as f:assert len(next(csv.reader(f,delimiter='\t')))>1,p
print('PASS: file hashes, Python syntax, table headers, no image/manuscript files. No scientific analysis executed.')
