#!/usr/bin/env python3
"""Compare observed and radial-corrected DR1 positions near one TARGETID."""
from __future__ import annotations

import argparse
import gzip
import json
import io
import struct
from html import escape
from pathlib import Path

import numpy as np

MPC_TO_MLY = 3.261563777


def read_chunk(path: Path, count: int, scales: tuple[float, float]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    raw = gzip.decompress(path.read_bytes())
    if len(raw) < 64 or raw[:4] != b"DSC3":
        raise ValueError(f"{path}: invalid V3 header")
    version, n, type_off, row_off, pos_off, z_off, attr_off, end = struct.unpack_from("<8I", raw, 4)
    if version != 1 or n != count or end != len(raw) or not (64 <= type_off < row_off < pos_off < z_off < attr_off < end):
        raise ValueError(f"{path}: invalid V3 offsets or count")
    types = np.frombuffer(raw, dtype="u1", count=n, offset=type_off)
    rows = np.frombuffer(raw, dtype="<u4", count=n, offset=row_off)
    axes = np.frombuffer(raw, dtype="<i2", count=3*n, offset=pos_off).reshape(3, n).astype(np.int64)
    native = np.cumsum(axes, axis=1, dtype=np.int64).T
    scale = np.where((types & 64) != 0, scales[1], scales[0])
    return native / scale[:, None], types, rows


def read_ids(path: Path, count: int) -> np.ndarray:
    with gzip.open(path, "rb") as stream:
        ids = np.load(io.BytesIO(stream.read()), allow_pickle=False)
    if ids.dtype != np.dtype("<u8") or len(ids) != count:
        raise ValueError(f"{path}: invalid TARGETID index")
    return ids


def selected_ids(types: np.ndarray, rows: np.ndarray, galaxy: np.ndarray, qso: np.ndarray) -> np.ndarray:
    q = (types & 64) != 0
    if (np.any(rows[~q] >= len(galaxy)) or np.any(rows[q] >= len(qso))):
        raise ValueError("V3 original-row index out of bounds")
    ids = np.empty(len(rows), dtype=np.uint64)
    ids[~q] = galaxy[rows[~q]]
    ids[q] = qso[rows[q]]
    return ids


def read_correction(path: Path, count: int) -> np.ndarray:
    values = np.frombuffer(gzip.decompress(path.read_bytes()), dtype="<f4")
    if values.size != count or not np.isfinite(values).all():
        raise ValueError(f"{path}: invalid correction count or values")
    return values


def basis(center: np.ndarray) -> np.ndarray:
    forward = center / np.linalg.norm(center)
    reference = np.array([0., 0., 1.])
    if abs(np.dot(forward, reference)) > .95:
        reference = np.array([0., 1., 0.])
    right = np.cross(reference, forward)
    right /= np.linalg.norm(right)
    up = np.cross(forward, right)
    return np.stack((right, up, forward), axis=1)


def select_region(observed: np.ndarray, center: np.ndarray, frame: np.ndarray,
                  dimensions: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    local = (observed - center) @ frame * MPC_TO_MLY
    inside = np.all(np.abs(local) <= dimensions / 2, axis=1)
    return inside, local


def corrected_positions(observed: np.ndarray, delta_mpc: np.ndarray) -> np.ndarray:
    radius = np.linalg.norm(observed, axis=1)
    direction = np.divide(observed, radius[:, None], out=np.zeros_like(observed),
                          where=radius[:, None] > 0)
    return observed + direction * delta_mpc[:, None]


def point_color(fog_mpc: float, rsd_mpc: float, scale_mpc: float, component: str) -> str:
    if component in ('rsd', 'both') and rsd_mpc != 0:
        fraction = min(abs(rsd_mpc) / scale_mpc, 1.0)
        if rsd_mpc > 0:
            return f'#ff{int(172 - 108 * fraction):02x}50'
        return f'#50{int(170 - 120 * fraction):02x}{int(210 + 45 * fraction):02x}'
    return '#f59e0b' if component in ('fog', 'both') and fog_mpc != 0 else '#73869b'


def write_svg(path: Path, targetid: int, dimensions: np.ndarray, before: np.ndarray,
              after: np.ndarray, fog: np.ndarray, rsd: np.ndarray, target_after: list[float],
              shown: np.ndarray, component: str) -> None:
    panel_width, panel_height = 600, 440
    margin_x, margin_y = 72, 68
    plot_w, plot_h = panel_width - 115, panel_height - 120
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 940">']
    parts += ['<rect width="1200" height="940" fill="#07111c"/>',
              '<style>text{font-family:Arial,sans-serif;fill:#dbeafe} .title{font-size:17px;font-weight:bold}.axis{font-size:13px;fill:#9bb4cd}</style>',
              f'<text x="600" y="30" text-anchor="middle" class="title">TARGETID {targetid} · '
              f'{dimensions[0]:g} × {dimensions[1]:g} × {dimensions[2]:g} Mly · {len(before):,} objects</text>']
    scale_mpc = max(float(np.percentile(np.abs(rsd[rsd != 0]), 95)) if np.any(rsd != 0) else 1., 1.)
    colors = [point_color(float(fog[index]), float(rsd[index]), scale_mpc, component) for index in shown]
    for row, (x, y, section, xlabel, ylabel) in enumerate(((0, 1, 'Sky plane', 'Right (Mly)', 'Up (Mly)'),
                                                              (0, 2, 'Line of sight', 'Right (Mly)', 'Depth (Mly)'))):
        for col, (points, label) in enumerate(((before, 'Observed'), (after, f'Corrected ({component})'))):
            ox, oy = col * panel_width + margin_x, row * panel_height + margin_y + 38
            parts.append(f'<rect x="{ox}" y="{oy}" width="{plot_w}" height="{plot_h}" fill="#0d1b2a" stroke="#456079"/>')
            parts.append(f'<text x="{ox}" y="{oy-14}" class="title">{escape(label)} — {section}</text>')
            sx = ox + (points[shown, x] / dimensions[x] + .5) * plot_w
            sy = oy + (.5 - points[shown, y] / dimensions[y]) * plot_h
            for color, px, py in zip(colors, sx, sy):
                parts.append(f'<circle cx="{px:.2f}" cy="{py:.2f}" r="1.2" fill="{color}" fill-opacity=".55"/>')
            tx, ty = (0., 0.) if col == 0 else (target_after[x], target_after[y])
            px = ox + (tx / dimensions[x] + .5) * plot_w
            py = oy + (.5 - ty / dimensions[y]) * plot_h
            parts.append(f'<path d="M{px-8:.2f},{py:.2f}h16 M{px:.2f},{py-8:.2f}v16" stroke="#ff5c66" stroke-width="2"/>')
            parts.append(f'<text x="{ox+plot_w/2:.1f}" y="{oy+plot_h+26}" text-anchor="middle" class="axis">{xlabel}</text>')
            parts.append(f'<text transform="translate({ox-48},{oy+plot_h/2}) rotate(-90)" text-anchor="middle" class="axis">{ylabel}</text>')
            parts.append(f'<text x="{ox}" y="{oy+plot_h+44}" class="axis">−{dimensions[x]/2:g}</text>')
            parts.append(f'<text x="{ox+plot_w}" y="{oy+plot_h+44}" text-anchor="end" class="axis">+{dimensions[x]/2:g}</text>')
    legend = (f'Applied RSD: blue = inward, orange/red = outward; p95 = {scale_mpc:.2f} Mpc · '
              if component in ('rsd', 'both') else '')
    legend += 'amber = applied FoG without RSD · gray = unchanged · red cross = TARGETID'
    parts.append(f'<text x="600" y="917" text-anchor="middle" class="axis">{legend}</text>')
    parts.append('</svg>')
    path.write_text('\n'.join(parts), encoding='utf-8')


def histogram(values: np.ndarray, bins: int, signed: bool) -> tuple[np.ndarray, np.ndarray, int]:
    finite = values[np.isfinite(values) & (values != 0)]
    if signed:
        span = float(np.percentile(np.abs(finite), 99)) if len(finite) else 1.
        edges = np.linspace(-max(span, 1e-6), max(span, 1e-6), bins + 1)
    else:
        ceiling = float(np.percentile(finite, 99)) if len(finite) else 1.
        edges = np.linspace(0, max(ceiling, 1e-6), bins + 1)
    counts = np.histogram(finite, bins=edges)[0]
    overflow = int(np.count_nonzero((finite < edges[0]) | (finite > edges[-1])))
    return counts, edges, overflow


def write_histograms(path: Path, rsd: np.ndarray, fog: np.ndarray, component: str,
                     bins: int) -> dict:
    selected = (rsd if component in ('rsd', 'both') else 0) + (fog if component in ('fog', 'both') else 0)
    series = [('Signed RSD displacement', rsd, True, '#55aaff'),
              ('Absolute RSD displacement', np.abs(rsd), False, '#f07a52'),
              (f'Absolute applied displacement ({component})', np.abs(selected), False, '#b993ff')]
    panels, report = [], {}
    for title, values, signed, color in series:
        counts, edges, overflow = histogram(values, bins, signed)
        report[title] = {'binEdgesMpc': edges.tolist(), 'counts': counts.tolist(),
                         'zeroCount': int(np.count_nonzero(values == 0)),
                         'outsidePlottedP99Range': overflow}
        panels.append((title, counts, edges, color, overflow))
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1020 1010">',
             '<rect width="1020" height="1010" fill="#07111c"/>',
             '<style>text{font-family:Arial,sans-serif;fill:#dbeafe}.sub{font-size:13px;fill:#9bb4cd}</style>',
             '<text x="510" y="35" text-anchor="middle" font-size="20">Radial displacement histograms · all selected objects</text>']
    for index, (title, counts, edges, color, overflow) in enumerate(panels):
        x0, y0, width, height = 90, 95 + index * 305, 850, 205
        maximum = max(int(counts.max()), 1)
        parts.append(f'<text x="{x0}" y="{y0-18}" font-size="17">{escape(title)}</text>')
        parts.append(f'<rect x="{x0}" y="{y0}" width="{width}" height="{height}" fill="#0d1b2a" stroke="#456079"/>')
        for j, count in enumerate(counts):
            bar_h = count / maximum * (height - 8)
            bar_w = width / len(counts)
            parts.append(f'<rect x="{x0+j*bar_w:.2f}" y="{y0+height-bar_h:.2f}" width="{max(bar_w-1,.2):.2f}" height="{bar_h:.2f}" fill="{color}"/>')
        parts.append(f'<text x="{x0}" y="{y0+height+23}" class="sub">{edges[0]:.2f} Mpc</text>')
        parts.append(f'<text x="{x0+width}" y="{y0+height+23}" text-anchor="end" class="sub">{edges[-1]:.2f} Mpc</text>')
        parts.append(f'<text x="{x0+width/2}" y="{y0+height+23}" text-anchor="middle" class="sub">nonzero values · p99 display range · {overflow:,} outliers outside</text>')
    parts.append('</svg>')
    path.write_text('\n'.join(parts), encoding='utf-8')
    return report


def run(args: argparse.Namespace) -> dict:
    manifest = json.loads((args.catalogs / "desiV3-consolidated-manifest.json").read_text())
    chunks = manifest["chunks"]
    if manifest.get("encoding") != "desi-consolidated-1" or sum(c["count"] for c in chunks) != manifest["count"]:
        raise ValueError("invalid V3 manifest")
    galaxy = read_ids(args.identity / "DR1_GALAXY_targetid.npy.gz", manifest["catalogs"][0]["count"])
    qso = read_ids(args.identity / "DR1_QSO_targetid.npy.gz", manifest["catalogs"][1]["count"])
    scales = (manifest["catalogs"][0]["unitsPerMpc"], manifest["catalogs"][1]["unitsPerMpc"])
    if not np.isfinite(args.dimensions).all() or np.any(args.dimensions <= 0):
        raise ValueError("width, height and depth must be positive finite Mly")
    if args.targetid <= 0 or args.targetid >= 2**64:
        raise ValueError("TARGETID must be a positive unsigned 64-bit integer")

    target = np.uint64(args.targetid)
    center = None
    center_chunk = None
    for number, spec in enumerate(chunks):
        positions, types, rows = read_chunk(args.catalogs / spec["file"], spec["count"], scales)
        ids = selected_ids(types, rows, galaxy, qso)
        hits = np.flatnonzero(ids == target)
        if hits.size:
            if center is not None or hits.size != 1:
                raise ValueError("TARGETID is not unique")
            center = positions[hits[0]]
            center_chunk = number
    if center is None or np.linalg.norm(center) == 0:
        raise ValueError("TARGETID missing or has zero-distance position")

    frame = basis(center)
    before_parts, after_parts, id_parts, fog_parts, rsd_parts = [], [], [], [], []
    total = 0
    for number, spec in enumerate(chunks):
        count = spec["count"]
        positions, types, rows = read_chunk(args.catalogs / spec["file"], count, scales)
        inside, local_before = select_region(positions, center, frame, args.dimensions)
        if not inside.any():
            continue
        ids = selected_ids(types[inside], rows[inside], galaxy, qso)
        fog = read_correction(args.corrections / f"fog_{number:03d}.bin", count)[inside]
        rsd = read_correction(args.corrections / f"rsd_{number:03d}.bin", count)[inside]
        displacement = (fog if args.component in ("fog", "both") else 0) + (rsd if args.component in ("rsd", "both") else 0)
        after = corrected_positions(positions[inside], displacement)
        before_parts.append(local_before[inside])
        after_parts.append((after - center) @ frame * MPC_TO_MLY)
        id_parts.append(ids)
        fog_parts.append(fog)
        rsd_parts.append(rsd)
        total += int(inside.sum())

    before = np.concatenate(before_parts) if before_parts else np.empty((0, 3))
    after = np.concatenate(after_parts) if after_parts else np.empty((0, 3))
    ids = np.concatenate(id_parts) if id_parts else np.empty(0, dtype=np.uint64)
    fog = np.concatenate(fog_parts) if fog_parts else np.empty(0)
    rsd = np.concatenate(rsd_parts) if rsd_parts else np.empty(0)
    if not np.isfinite(after).all():
        raise ValueError("nonfinite corrected position")
    retained = np.all(np.abs(after) <= args.dimensions / 2, axis=1)
    # Keep the central coordinate fixed so every panel uses the same reference point.
    center_row = np.flatnonzero(ids == target)
    if center_row.size != 1:
        raise ValueError("central TARGETID was not recovered exactly once")
    target_after = after[center_row[0]].tolist()

    args.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    shown = np.arange(total) if total <= args.max_points else np.sort(rng.choice(total, args.max_points, replace=False))
    write_svg(args.out / "comparison.svg", args.targetid, args.dimensions,
              before, after, fog, rsd, target_after, shown, args.component)
    histograms = write_histograms(args.out / 'histograms.svg', rsd, fog, args.component, args.histogram_bins)

    applied = (rsd if args.component in ('rsd', 'both') else 0) + (fog if args.component in ('fog', 'both') else 0)
    nonzero = np.abs(applied[applied != 0])
    largest = np.argsort(np.abs(applied))[-min(25, total):][::-1]
    (args.out / 'largest_displacements.csv').write_text(
        'targetid,fog_mpc,rsd_mpc,applied_mpc\n' +
        ''.join(f'{int(ids[i])},{fog[i]:.7g},{rsd[i]:.7g},{applied[i]:.7g}\n' for i in largest))

    summary = {
        "targetid": str(args.targetid), "centerChunk": center_chunk,
        "centerObservedMpc": center.tolist(), "dimensionsMly": args.dimensions.tolist(),
        "component": args.component, "selection": "observed-position axis-aligned box in Earth line-of-sight frame",
        "objectsInObservedBox": total, "correctedPositionsStillInBox": int(retained.sum()),
        "correctedPositionsOutsideBox": int((~retained).sum()),
        "objectsWithNonzeroFoG": int(np.count_nonzero(fog)),
        "objectsWithNonzeroRSD": int(np.count_nonzero(rsd)),
        "objectsWithAppliedDisplacement": int(np.count_nonzero(applied)),
        "absoluteAppliedDisplacementMpc": {
            "median": float(np.median(nonzero)) if len(nonzero) else 0.,
            "p95": float(np.percentile(nonzero, 95)) if len(nonzero) else 0.,
            "maximum": float(np.max(nonzero)) if len(nonzero) else 0.,
        },
        "histogramData": "histograms.json",
        "samplePlotted": int(len(shown)), "sampleSeed": args.seed,
        "redCross": "same TARGETID in both panels; corrected cross may move radially",
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (args.out / 'histograms.json').write_text(json.dumps(histograms, indent=2) + '\n')
    print(json.dumps(summary, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("targetid", type=int)
    parser.add_argument("--width-mly", type=float, required=True)
    parser.add_argument("--height-mly", type=float, required=True)
    parser.add_argument("--depth-mly", type=float, required=True)
    parser.add_argument("--component", choices=("rsd", "fog", "both"), default="both")
    parser.add_argument("--catalogs", type=Path, default=Path("desiV3/catalogs"))
    parser.add_argument("--identity", type=Path, default=Path("review_candidate/identity"))
    parser.add_argument("--corrections", type=Path, default=Path("review_candidate/payloads"))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--max-points", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--histogram-bins", type=int, default=40)
    args = parser.parse_args()
    args.dimensions = np.array([args.width_mly, args.height_mly, args.depth_mly], dtype=float)
    if args.max_points < 1 or args.histogram_bins < 2:
        parser.error("--max-points must be positive and --histogram-bins must be at least 2")
    run(args)


if __name__ == "__main__":
    main()
