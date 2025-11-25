#!/usr/bin/env python3
"""
Move files in fornax not overlapping with video to junk folder for deletion.
After running, check that the correct files have bee moved with fornax_footprints.py

Created: Tuesday 25th November 2025.
"""
#!/usr/bin/env python3
from pathlib import Path
import shutil
import pickle
import re
import argparse
import sys

def load_allowed(pkl_path, keys):
    with open(pkl_path, "rb") as f:
        d = pickle.load(f)
    # Normalize everything to strings (pickle may contain ints or strings)
    allowed = set()
    for k in keys:
        for v in d.get(k, []):
            allowed.add(str(v))
    return allowed

def main():
    p = argparse.ArgumentParser(description="Move CDFS fits not in allowed tile lists to ./junk")
    p.add_argument("--base", default=str(Path.home() / "euclid" / "CDFS"), help="Base directory with .fits files")
    p.add_argument("--pkl", default="euclid_DR1_within_video.pkl", help="Pickle file with CDFS lists")
    p.add_argument("--keys", nargs="+", default=["CDFS1","CDFS2","CDFS3"], help="Keys to gather allowed tiles from")
    p.add_argument("--chunk", type=int, default=50, help="Progress print interval")
    p.add_argument("--commit", action="store_true", help="Actually move files. Without this, script does a dry-run.")
    args = p.parse_args()

    base = Path(args.base).expanduser()
    if not base.exists():
        print("Base directory does not exist:", base, file=sys.stderr)
        sys.exit(1)

    junk = base / "junk"
    if args.commit:
        junk.mkdir(exist_ok=True)

    allowed = load_allowed(args.pkl, args.keys)

    tile_re = re.compile(r"TILE(\d+)")
    files = sorted(base.glob("*.fits"))
    total = len(files)
    moved = 0
    to_move = []  # store paths for commit

    for i, fpath in enumerate(files, 1):
        m = tile_re.search(fpath.name)
        if not m:
            reason = "no-TILE"
            should_move = True
        else:
            tile_str = m.group(1)
            should_move = tile_str not in allowed
            reason = f"tile={tile_str}"

        if should_move:
            to_move.append((fpath, reason))

        if i % args.chunk == 0 or i == total:
            print(f"{i}/{total} scanned, {len(to_move)} flagged for move")

    if not to_move:
        print("Nothing to move.")
        return

    print()
    print("Sample flagged files (up to 10):")
    for pth, reason in to_move[:10]:
        print(f"  {pth.name}  [{reason}]")
    print()
    if not args.commit:
        print("Dry-run (no files moved). Rerun with --commit to perform moves.")
        return

    # Commit: move and log
    log_path = base / "moved_files.log"
    with open(log_path, "a") as logf:
        for i, (fpath, reason) in enumerate(to_move, 1):
            target = junk / fpath.name
            try:
                shutil.move(str(fpath), str(target))
                moved += 1
                logf.write(f"{fpath.name}\t{reason}\n")
            except Exception as e:
                print(f"Failed to move {fpath.name}: {e}", file=sys.stderr)

            if i % args.chunk == 0 or i == len(to_move):
                print(f"{i}/{len(to_move)} moved, {moved} successful")

    print(f"Done. {moved} files moved. Log: {log_path}")

if __name__ == "__main__":
    main()
