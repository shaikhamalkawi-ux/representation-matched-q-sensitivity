"""Validate the exact public-release member set and SHA-256 manifest."""
from pathlib import Path, PurePosixPath
import hashlib

ROOT=Path(__file__).resolve().parent
def verify():
    entries={}
    for line in (ROOT/'SHA256SUMS.txt').read_text(encoding='utf-8').splitlines():
        expected,name=line.split('  ',1)
        p=PurePosixPath(name)
        if p.is_absolute() or '..' in p.parts or '\\' in name or name in entries:
            raise ValueError('Unsafe or duplicate manifest entry')
        path=ROOT/name
        if not path.is_file() or path.is_symlink() or hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
            raise ValueError('Missing or changed release member: '+name)
        entries[name]=expected
    actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file()}
    if actual!=set(entries)|{'SHA256SUMS.txt'}:
        raise ValueError('Unexpected or missing files in release')
    print('PASS:',len(entries),'manifest entries')
if __name__=='__main__': verify()
