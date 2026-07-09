#!/usr/bin/env python3

from pathlib import Path
import argparse
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser(
        description="Run MAFFT on all FASTA files in a directory."
    )
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--threads", default="1")
    args = parser.parse_args()

    mafft = shutil.which("mafft")

    if mafft is None:
        raise SystemExit("MAFFT was not found in PATH.")

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_files = sorted(input_dir.glob("*.fasta"))

    if not fasta_files:
        raise SystemExit(f"No FASTA files found in {input_dir}")

    for fasta in fasta_files:
        output = output_dir / f"{fasta.stem}_aligned.fasta"

        command = [
            mafft,
            "--auto",
            "--thread",
            str(args.threads),
            str(fasta),
        ]

        with open(output, "w", encoding="utf-8") as handle:
            subprocess.run(
                command,
                stdout=handle,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )

        print(f"Aligned {fasta} -> {output}")


if __name__ == "__main__":
    main()
