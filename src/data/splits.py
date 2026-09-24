#!/usr/bin/env python
"""Shared split-manifest utilities: frozen train/val membership for SPLIT-A/B/C.

A split manifest records which variant IDs go to train vs. validation, plus a
content hash over that membership, so a split can't silently drift between runs
(e.g. if the source parquet is rebuilt in a different row order). SPLIT-A's split
was never frozen this way -- it's recomputed inline by ``train_test_split(...,
random_state=42, stratify=df.label)`` in run_experiment.py/train_transformer_local.py
every run. That default path is left untouched; ``split_manifest_a.json`` here is a
documentation/audit freeze of what that call currently produces, not an enforced
input.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SPLITS_DIR = ROOT / "data" / "splits"


def variant_ids(df: pd.DataFrame) -> pd.Series:
    """Stable per-row identifier: chrom_pos_ref_alt (unique post drop_duplicates in build_dataset.py)."""
    return (
        df["chrom"].astype(str) + "_" + df["pos"].astype(str) + "_" + df["ref"].astype(str) + "_" + df["alt"].astype(str)
    )


def content_hash(train_ids: list[str], val_ids: list[str]) -> str:
    h = hashlib.sha256()
    h.update("\n".join(sorted(train_ids)).encode())
    h.update(b"\x00---VAL---\x00")
    h.update("\n".join(sorted(val_ids)).encode())
    return "sha256:" + h.hexdigest()


def build_manifest(split_id: str, name: str, description: str, source_data: str,
                    train_ids: list[str], val_ids: list[str], seed: int, extra: dict) -> dict:
    manifest = {
        "split_id": split_id,
        "name": name,
        "description": description,
        "source_data": source_data,
        "variant_key": ["chrom", "pos", "ref", "alt"],
        "seed": seed,
        "n_train": len(train_ids),
        "n_val": len(val_ids),
        "train_ids": sorted(train_ids),
        "val_ids": sorted(val_ids),
        "content_hash": content_hash(train_ids, val_ids),
    }
    manifest.update(extra)
    return manifest


def save_manifest(manifest: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2))


def load_manifest(path: Path) -> dict:
    manifest = json.loads(Path(path).read_text())
    recomputed = content_hash(manifest["train_ids"], manifest["val_ids"])
    if recomputed != manifest["content_hash"]:
        raise ValueError(
            f"{path}: content_hash mismatch (file may have been hand-edited or corrupted); "
            f"stored={manifest['content_hash']} recomputed={recomputed}"
        )
    return manifest


def apply_split(df: pd.DataFrame, manifest: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split df into (train_df, val_df) per a loaded manifest's frozen membership."""
    ids = variant_ids(df)
    train_set = set(manifest["train_ids"])
    val_set = set(manifest["val_ids"])
    train_mask = ids.isin(train_set)
    val_mask = ids.isin(val_set)
    missing = len(df) - int(train_mask.sum()) - int(val_mask.sum())
    if missing:
        raise ValueError(
            f"{missing} rows in {manifest.get('source_data')} are not present in manifest "
            f"'{manifest['split_id']}' (source data has drifted from what the split was built on)"
        )
    return df[train_mask].copy(), df[val_mask].copy()
