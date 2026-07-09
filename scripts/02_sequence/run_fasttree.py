#!/usr/bin/env python3

from pathlib import Path
import argparse
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser(
        description="Run FastTree on aligned FASTA files."
    )
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    fasttree = shutil.which("FastTree") or shutil.which("fasttree")

    if fasttree is None:
        raise SystemExit("FastTree was not found in PATH.")

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    alignment_files = sorted(input_dir.glob("*_aligned.fasta"))

    if not alignment_files:
        raise SystemExit(f"No aligned FASTA files found in {input_dir}")

    for alignment in alignment_files:
        output = output_dir / f"{alignment.stem}.nwk"

        command = [
            fasttree,
            "-wag",
            str(alignment),
        ]

        with open(output, "w", encoding="utf-8") as handle:
            subprocess.run(
                command,
                stdout=handle,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )

        print(f"Tree {alignment} -> {output}")


if __name__ == "__main__":
    main()
