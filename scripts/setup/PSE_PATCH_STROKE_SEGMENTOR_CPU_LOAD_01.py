#!/usr/bin/env python3
"""
PSE_PATCH_STROKE_SEGMENTOR_CPU_LOAD_01.py

Patch stroke_segmentor so TorchScript checkpoints saved with CUDA tensors can
be loaded on a CPU-only machine (e.g. Apple Silicon Mac).

Problem
-------
stroke_segmentor currently does:
    torch.jit.load(checkpoint).to(self.device)

The .to(self.device) happens too late: torch.jit.load() may already try to
restore CUDA tensors and fail before the model can be moved to CPU.

Patch
-----
Replace with:
    torch.jit.load(checkpoint, map_location=self.device).to(self.device)

The script:
- locates the installed stroke_segmentor/model_handler.py
- creates a .bak backup once
- applies the minimal one-line patch
- verifies the patched source

Run inside the pse-stroke environment:
    conda activate pse-stroke
    python PSE_PATCH_STROKE_SEGMENTOR_CPU_LOAD_01.py
"""

from __future__ import annotations

import inspect
import shutil
from pathlib import Path

import stroke_segmentor.model_handler as mh


def main() -> int:
    model_handler_file = Path(inspect.getfile(mh)).resolve()
    backup_file = model_handler_file.with_suffix(model_handler_file.suffix + ".bak")

    old = "torch.jit.load(checkpoint).to(self.device)"
    new = "torch.jit.load(checkpoint, map_location=self.device).to(self.device)"

    print("=" * 68)
    print("PATCH STROKE_SEGMENTOR TORCHSCRIPT CPU LOADING")
    print("=" * 68)
    print(f"File   : {model_handler_file}")
    print(f"Backup : {backup_file}")
    print()

    text = model_handler_file.read_text(encoding="utf-8")

    if new in text:
        print("Patch already present. Nothing to change.")
    elif old in text:
        if not backup_file.exists():
            shutil.copy2(model_handler_file, backup_file)
            print("Backup created.")
        else:
            print("Backup already exists; leaving it untouched.")

        text = text.replace(old, new, 1)
        model_handler_file.write_text(text, encoding="utf-8")
        print("Patch applied.")
    else:
        raise RuntimeError(
            "Expected source line was not found. "
            "Package version may have changed; no file was modified."
        )

    verify = model_handler_file.read_text(encoding="utf-8")
    if new not in verify:
        raise RuntimeError("Verification failed.")

    print()
    print("Verified patched line:")
    print(f"  {new}")
    print()
    print("=" * 68)
    print("PATCH COMPLETE")
    print("=" * 68)
    print("Now test:")
    print(
        '  python -c "from stroke_segmentor.inferer import Inferer; '
        "Inferer(force_cpu=True); print('Inferer CPU OK')\""
    )
    print()
    print("Then rerun:")
    print("  python PSE_DEEPISLES_INFER_03.py P002")
    print("=" * 68)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
