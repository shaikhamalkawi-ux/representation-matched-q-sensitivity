"""Verify the exact released payload; standard library only, no network/writes."""
import hashlib
import json
import os
import re
from pathlib import Path


def unique_pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate manifest key: {key}")
        result[key] = value
    return result


def main():
    root = Path(__file__).resolve().parent
    if os.name == 'nt' and not str(root).startswith('\\\\?\\'):
        root = Path('\\\\?\\' + str(root))
    entries = json.loads((root / "MANIFEST_SHA256.json").read_text(encoding="utf-8"),
                         object_pairs_hook=unique_pairs)
    if not isinstance(entries, dict) or not entries:
        raise SystemExit("Expected a nonempty path-to-SHA256 manifest")
    failures = []
    for name, expected in entries.items():
        if not isinstance(name, str) or not isinstance(expected, str) or not re.fullmatch(r"[0-9a-f]{64}", expected):
            failures.append(f"invalid manifest entry: {name}")
            continue
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts or relative.as_posix() != name or ':' in name or '\\' in name:
            failures.append(f"unsafe manifest path: {name}")
            continue
        path = root / relative
        if path.is_symlink() or not path.resolve().is_relative_to(root):
            failures.append(f"symlink or escaped root: {name}")
        elif not path.is_file():
            failures.append(f"missing: {name}")
        elif hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            failures.append(f"SHA256 mismatch: {name}")
    actual = {p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file()}
    extra = actual - set(entries) - {'MANIFEST_SHA256.json'}
    failures.extend(f"unlisted file: {name}" for name in sorted(extra))
    if failures:
        raise SystemExit("\n".join(failures))
    print(f"PASS: {len(entries)} SHA-256 payload checks; no unlisted files (manifest excludes itself)")


if __name__ == "__main__":
    main()
