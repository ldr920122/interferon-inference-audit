"""Verify the archived release contents against MANIFEST_sha256.txt (MIT)."""
from pathlib import Path
import hashlib
root=Path(__file__).resolve().parent
count=0
for line in (root/'MANIFEST_sha256.txt').read_text().splitlines():
    expected, name = line.split('  ',1)
    p=root/name
    assert p.is_file(), f'Missing: {name}'
    assert hashlib.sha256(p.read_bytes()).hexdigest()==expected, f'Checksum mismatch: {name}'
    count+=1
print(f'PASS: {count} release files match SHA-256 manifest')
