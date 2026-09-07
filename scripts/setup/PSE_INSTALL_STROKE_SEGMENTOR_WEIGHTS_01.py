#!/usr/bin/env python3
"""
PSE_INSTALL_STROKE_SEGMENTOR_WEIGHTS_01.py

Workaround for the current stroke_segmentor Zenodo archive-download failure.

What it does
------------
1. Reads the Zenodo record configured by stroke_segmentor.
2. Uses the direct file download endpoint instead of the failing files-archive endpoint.
3. Downloads stroke_segmentor_weights.zip with resume support.
4. Verifies the downloaded size and MD5 checksum when available.
5. Extracts the archive.
6. Finds model0.ts ... model14.ts recursively.
7. Places them in the exact local cache folder expected by stroke_segmentor.

Run inside the conda environment:
    conda activate pse-stroke
    python PSE_INSTALL_STROKE_SEGMENTOR_WEIGHTS_01.py
"""

from __future__ import annotations

import hashlib
import os
import shutil
import sys
import time
import zipfile
from pathlib import Path

import requests
import stroke_segmentor.zenodo as z


N_MODELS = 15
CHUNK_SIZE = 8 * 1024 * 1024  # 8 MiB
MIN_FREE_GB_WARN = 12.0


def human_bytes(n: int) -> str:
    units = ["B", "KiB", "MiB", "GiB", "TiB"]
    x = float(n)
    for unit in units:
        if x < 1024 or unit == units[-1]:
            return f"{x:.2f} {unit}"
        x /= 1024
    return f"{n} B"


def md5sum(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        while True:
            block = f.read(16 * 1024 * 1024)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def main() -> int:
    weights_root = Path(z.WEIGHTS_FOLDER)
    weights_root.mkdir(parents=True, exist_ok=True)

    free = shutil.disk_usage(weights_root).free
    print("=" * 68)
    print("STROKE_SEGMENTOR WEIGHTS INSTALLER")
    print("=" * 68)
    print(f"Configured record ID : {z.WEIGHTS_RECORD_ID}")
    print(f"Weights root         : {weights_root}")
    print(f"Free disk space      : {human_bytes(free)}")
    if free < MIN_FREE_GB_WARN * (1024**3):
        print(
            f"WARNING: less than {MIN_FREE_GB_WARN:.0f} GiB free. "
            "The archive is ~5 GB and extraction needs additional space."
        )
        answer = input("Continue anyway? [y/N] ").strip().lower()
        if answer not in {"y", "yes"}:
            return 1

    # Use the package's own Zenodo base URL and configured concept record ID.
    metadata_url = f"{z.ZenodoRecord.BASE_URL}/{z.WEIGHTS_RECORD_ID}"
    print("\nFetching Zenodo metadata...")
    r = requests.get(metadata_url, timeout=60)
    r.raise_for_status()
    data = r.json()

    metadata = data.get("metadata", {})
    version = str(metadata.get("version", "")).strip()
    if not version:
        raise RuntimeError("Zenodo metadata has no 'version' field.")

    files = data.get("files", [])
    if not files:
        raise RuntimeError("Zenodo record contains no files.")

    # Prefer the expected weights zip; otherwise use the only/first zip.
    file_info = None
    for item in files:
        if item.get("key") == "stroke_segmentor_weights.zip":
            file_info = item
            break
    if file_info is None:
        zip_candidates = [x for x in files if str(x.get("key", "")).endswith(".zip")]
        if len(zip_candidates) == 1:
            file_info = zip_candidates[0]
        else:
            raise RuntimeError(
                "Could not uniquely identify stroke_segmentor_weights.zip."
            )

    file_name = str(file_info["key"])
    expected_size = int(file_info.get("size", 0))
    checksum = str(file_info.get("checksum", "")).strip()

    links = file_info.get("links", {})
    direct_url = links.get("self") or links.get("content")
    if not direct_url:
        raise RuntimeError("No direct file download link present in Zenodo metadata.")

    # This mirrors the cache-folder naming used by stroke_segmentor.
    target_dir = weights_root / f"{z.WEIGHTS_RECORD_ID}_v{version}"

    expected_models = [target_dir / f"model{i}.ts" for i in range(N_MODELS)]
    if all(p.is_file() for p in expected_models):
        print("\nAll model weights are already installed:")
        print(f"  {target_dir}")
        print("Nothing to do.")
        return 0

    download_dir = weights_root / "_direct_download"
    extract_dir = weights_root / "_direct_extract"
    download_dir.mkdir(parents=True, exist_ok=True)

    final_zip = download_dir / file_name
    part_zip = download_dir / (file_name + ".part")

    print(f"\nZenodo version       : {version}")
    print(f"Remote file          : {file_name}")
    if expected_size:
        print(f"Expected size        : {human_bytes(expected_size)}")
    if checksum:
        print(f"Expected checksum    : {checksum}")
    print(f"Target cache folder  : {target_dir}")

    # If a previous completed archive exists, reuse it.
    if final_zip.is_file() and expected_size and final_zip.stat().st_size == expected_size:
        print("\nComplete archive already present; skipping download.")
    else:
        if final_zip.exists():
            final_zip.unlink()

        existing = part_zip.stat().st_size if part_zip.exists() else 0
        headers = {}
        mode = "wb"

        if existing > 0:
            headers["Range"] = f"bytes={existing}-"
            print(
                f"\nResuming partial download at {human_bytes(existing)} "
                f"({existing / expected_size * 100:.1f}% if remote size unchanged)."
                if expected_size
                else f"\nResuming partial download at {human_bytes(existing)}."
            )
        else:
            print("\nStarting direct Zenodo file download...")

        with requests.get(
            direct_url,
            stream=True,
            headers=headers,
            timeout=(30, 600),
        ) as resp:
            if existing > 0 and resp.status_code == 206:
                mode = "ab"
                downloaded = existing
            elif resp.status_code == 200:
                # Server ignored Range or this is a fresh download.
                mode = "wb"
                downloaded = 0
                if existing:
                    print("Server did not honor Range; restarting download.")
            else:
                resp.raise_for_status()
                downloaded = existing

            resp.raise_for_status()

            last_print = time.time()
            with part_zip.open(mode) as f:
                for chunk in resp.iter_content(chunk_size=CHUNK_SIZE):
                    if not chunk:
                        continue
                    f.write(chunk)
                    downloaded += len(chunk)

                    now = time.time()
                    if now - last_print >= 2:
                        if expected_size:
                            pct = min(100.0, downloaded / expected_size * 100)
                            print(
                                f"\r  {human_bytes(downloaded)} / "
                                f"{human_bytes(expected_size)} ({pct:5.1f}%)",
                                end="",
                                flush=True,
                            )
                        else:
                            print(
                                f"\r  {human_bytes(downloaded)} downloaded",
                                end="",
                                flush=True,
                            )
                        last_print = now

        print()

        if expected_size and part_zip.stat().st_size != expected_size:
            raise RuntimeError(
                f"Downloaded size mismatch: got {part_zip.stat().st_size} bytes, "
                f"expected {expected_size} bytes. "
                "The .part file was kept so the download can be resumed."
            )

        part_zip.replace(final_zip)

    # Verify checksum if Zenodo supplied md5:<hex>.
    if checksum.lower().startswith("md5:"):
        expected_md5 = checksum.split(":", 1)[1].strip().lower()
        print("\nVerifying MD5 checksum (this may take a little while)...")
        actual_md5 = md5sum(final_zip).lower()
        print(f"  Actual MD5: {actual_md5}")
        if actual_md5 != expected_md5:
            raise RuntimeError(
                "MD5 checksum mismatch. Archive kept for inspection."
            )
        print("  Checksum OK.")

    # Clean a prior temporary extraction, but never silently delete the target cache.
    if extract_dir.exists():
        shutil.rmtree(extract_dir)
    extract_dir.mkdir(parents=True)

    print("\nExtracting archive...")
    with zipfile.ZipFile(final_zip, "r", allowZip64=True) as zf:
        bad = zf.testzip()
        if bad is not None:
            raise RuntimeError(f"Corrupt ZIP member detected: {bad}")
        zf.extractall(extract_dir)

    found = {}
    for i in range(N_MODELS):
        matches = list(extract_dir.rglob(f"model{i}.ts"))
        if len(matches) != 1:
            raise RuntimeError(
                f"Expected exactly one model{i}.ts, found {len(matches)}."
            )
        found[i] = matches[0]

    print(f"Found all {N_MODELS} model checkpoints.")

    if target_dir.exists():
        # Keep safety strict: only remove an incomplete cache folder that matches
        # this exact record/version naming convention.
        existing_models = list(target_dir.glob("model*.ts"))
        print(
            f"Removing incomplete target cache folder "
            f"({len(existing_models)} model files found)."
        )
        shutil.rmtree(target_dir)

    target_dir.mkdir(parents=True)

    for i in range(N_MODELS):
        src = found[i]
        dst = target_dir / f"model{i}.ts"
        shutil.move(str(src), str(dst))

    if not all(p.is_file() for p in expected_models):
        raise RuntimeError("Final verification failed: some model files are missing.")

    print("\nFinal model files:")
    total = 0
    for p in expected_models:
        size = p.stat().st_size
        total += size
        print(f"  {p.name:10s} {human_bytes(size)}")

    print(f"Total model size: {human_bytes(total)}")

    # Remove temporary extraction and downloaded ZIP after successful installation
    # to recover disk space.
    shutil.rmtree(extract_dir, ignore_errors=True)
    try:
        final_zip.unlink()
        download_dir.rmdir()
        print("\nTemporary ZIP removed after successful installation.")
    except OSError:
        pass

    print("\n" + "=" * 68)
    print("INSTALLATION COMPLETE")
    print("=" * 68)
    print(f"Weights installed in:\n  {target_dir}")
    print("\nNow test:")
    print("  python -c \"from stroke_segmentor.inferer import Inferer; Inferer(force_cpu=True); print('Inferer OK')\"")
    print("=" * 68)

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\nInterrupted. Partial download is kept and can be resumed.")
        raise SystemExit(130)
    except Exception as exc:
        print(f"\nERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise SystemExit(1)
