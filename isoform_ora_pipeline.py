#!/usr/bin/env python3

import subprocess
import argparse

def run(cmd):
    print("\nRunning:", " ".join(cmd))
    subprocess.run(cmd, check=True)

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--dtu", required=True)
    parser.add_argument("--prefix", default="results")
    parser.add_argument(
        "--go_file",
        required=True,
        help="Path to existing GO mapping file (e.g. IPS_go_map.txt)"
    )

    args = parser.parse_args()

    go_file = args.go_file

    # Step 1
    try:
        open(go_file)
        print("\nUsing provided GO mapping file:", go_file)
    except FileNotFoundError:
        raise FileNotFoundError(
            f"GO mapping file not found: {go_file}. "
            "Please provide a valid --go_file path."
        )

    # Step 2
    dtu_file = f"{args.prefix}_dtu.tsv"
    run(["python", "dtu_classify.py",
         "--dtu", args.dtu,
         "--out", dtu_file])

    # Step 3
    run(["python", "enrich.py",
         "--go", go_file,
         "--dtu", dtu_file,
         "--out_prefix", args.prefix])

    print("\nDone.")

if __name__ == "__main__":
    main()