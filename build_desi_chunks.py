#!/usr/bin/env python3
"""
build_desi_chunks.py — Empaqueta el catalogo DESI como trozos de PUNTOS para carga
progresiva, en vez de un unico .sog de 50 MB que hay que esperar entero.

El objetivo es el TIEMPO HASTA EL PRIMER FOTOGRAMA, no el tamano total. El primer trozo
pesa un par de MB y entra en menos de un segundo; el resto llega por detras sin bloquear.

Como se reparten los puntos entre trozos
----------------------------------------
Ordenar por curva de Morton y cortar en bloques contiguos seria lo mejor para comprimir,
pero el primer trozo seria una ESQUINA del universo: verias un pedazo denso y el resto
vacio. Aqui se ordena por Morton y luego se INTERCALA (`i % n_chunks`), de modo que cada
trozo es una muestra repartida por todo el volumen. Desde el primer trozo se ve el survey
entero, y los siguientes solo le suben la densidad.

El precio es compresion: dentro de un trozo los vecinos ya no son adyacentes sino que van
de k en k sobre la curva, asi que los deltas son mayores. Sigue comprimiendo mucho mejor
que el orden original, porque la curva mantiene la localidad espacial. El script imprime
las dos cifras para que se vea cuanto cuesta.

Formato
-------
manifest.json     { version, count, unitsPerMpc, chunks: [ {file, count, bytes} ] }
chunk_NNN.bin     gzip de: int16[count] dx | int16[count] dy | int16[count] dz
                  separado por eje (comprime mejor que intercalado) y en deltas respecto
                  al punto anterior del mismo trozo.

El gzip se deshace en el navegador con DecompressionStream, asi no depende de que el
servidor este configurado con Content-Encoding.

El color NO va aqui: se calcula del radio comovil en el shader (las capas de la lamina de
DESI). Meter g-r y r-z subiria el paquete de 41 a 58 MB, mas que el .sog actual. Si algun
dia se quiere el color fotometrico, que vaya en un archivo aparte que solo se descargue al
activarlo.

Uso:
  python build_desi_chunks.py                  # 12 trozos a ply/chunks/
  python build_desi_chunks.py --chunks 16
"""

import argparse
import gzip
import json
import struct
import sys
from pathlib import Path

import numpy as np

ASSETS = Path(__file__).parent / 'assets'
OUT_DIR = Path(__file__).parent / 'ply' / 'chunks'
LEGACY_TRACERS = ['BGS_BRIGHT', 'LRG', 'ELG', 'QSO']
DR1_TRACERS = ['DR1_GALAXY', 'DR1_QSO']


def read_desi(path: Path):
    raw = path.read_bytes()
    if raw[:4] != b'DESI':
        raise ValueError(f'{path.name}: no empieza por DESI')
    version, count, units_per_mpc, z_min, z_max = struct.unpack('<IIfff', raw[4:24])
    pos = np.frombuffer(raw, np.int16, count * 3, 48).reshape(count, 3)
    base = 48 + count * 8
    gr = np.frombuffer(raw, np.uint8, count, base + count * 2)
    rz = np.frombuffer(raw, np.uint8, count, base + count * 3)
    if version < 3:
        raise ValueError(f'{path.name}: se requiere .desi v3 con TARGETID para el indice FoG')
    target_id = np.frombuffer(raw, np.uint64, count, base + count * 16)
    return pos, float(units_per_mpc), float(z_max), gr, rz, target_id


def morton_order(pos: np.ndarray) -> np.ndarray:
    """Indices que ordenan por curva de Morton sobre una rejilla de 10 bits por eje.

    10 bits bastan: el objetivo es la localidad para el delta, no reproducir la posicion.
    La posicion se guarda entera, en int16.
    """
    q = ((pos.astype(np.int32) + 32768).astype(np.uint32) >> 6)

    def spread(x):
        x = x.astype(np.uint64) & 0x3ff
        x = (x | (x << 16)) & 0x30000ff
        x = (x | (x << 8)) & 0x300f00f
        x = (x | (x << 4)) & 0x30c30c3
        x = (x | (x << 2)) & 0x9249249
        return x

    key = spread(q[:, 0]) | (spread(q[:, 1]) << 1) | (spread(q[:, 2]) << 2)
    return np.argsort(key, kind='stable')


def encode_chunk(pos: np.ndarray) -> bytes:
    """Deltas int16 separados por eje, comprimidos con gzip.

    `prepend` va a CERO, no a pos[0]: el primer valor tiene que ser la posicion absoluta
    para que el decodificador arranque el acumulador en el sitio. Con prepend=pos[:1] el
    primer delta sale 0 y el trozo entero queda desplazado por su propio origen — cada
    trozo por un sitio distinto, y la nube se deforma.
    """
    d = np.diff(pos.astype(np.int32), axis=0, prepend=0).astype(np.int16)
    return gzip.compress(np.ascontiguousarray(d.T).tobytes(), 9)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--chunks', type=int, default=12, help='numero de trozos (def. 12)')
    ap.add_argument('--out', default=str(OUT_DIR))
    ap.add_argument('--tag', default='',
                    help='sufijo de version para invalidar cache, ej. 20260819')
    ap.add_argument('--fog-index-dir', default='',
                    help='directorio SOLO offline para TARGETID por fila de chunk. No se publica.')
    args = ap.parse_args()
    tag = f'_{args.tag}' if args.tag else ''

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    dr1_present = [(ASSETS / f'{name}.desi').exists() for name in DR1_TRACERS]
    if any(dr1_present) and not all(dr1_present):
        raise FileNotFoundError('DR1 requiere ambos archivos: DR1_GALAXY.desi y DR1_QSO.desi')
    tracers = DR1_TRACERS if all(dr1_present) else LEGACY_TRACERS
    tracer_codes = {name: code for code, name in enumerate(tracers)}

    # Cada .desi tiene su PROPIO unitsPerMpc: cada trazador se cuantizo a su propio alcance
    # para exprimir el int16 (BGS llega mucho menos lejos que QSO, asi que su escala es casi
    # el doble). Para juntarlos hay que pasar por Mpc y recuantizar con una escala unica.
    parts, types, target_parts, flag_parts = [], [], [], []
    max_redshift = 0.0
    for tracer in tracers:
        src = ASSETS / f'{tracer}.desi'
        if not src.exists():
            print(f'  [skip]  falta {src.name}')
            continue
        pos_q, u, z_max, gr, rz, target_id = read_desi(src)
        max_redshift = max(max_redshift, z_max)
        mpc = pos_q.astype(np.float32) / u
        parts.append(mpc)
        types.append(np.full(len(mpc), tracer_codes[tracer], dtype=np.uint8))
        target_parts.append(target_id)
        if tracers == DR1_TRACERS:
            flags_path = ASSETS / f'{tracer}.target-flags.bin'
            if not flags_path.exists() or flags_path.stat().st_size != len(mpc):
                raise ValueError(f'{flags_path.name}: missing or wrong row count')
            flag_parts.append(np.memmap(flags_path, dtype='u1', mode='r'))
        print(f'  [read]  {tracer}: {len(mpc):,}  (unitsPerMpc propio {u:.3f}, '
              f'radio max {np.linalg.norm(mpc, axis=1).max():.0f} Mpc)')

    if not parts:
        print('sin datos')
        return 1

    mpc = np.concatenate(parts)
    type_all = np.concatenate(types)
    target_flags = np.concatenate(flag_parts) if flag_parts else None
    target_all = np.concatenate(target_parts)
    n = len(mpc)

    # Escala global: la coordenada mas extrema manda, con un 2 % de margen para que el
    # redondeo no se salga del int16.
    extent = float(np.abs(mpc).max()) * 1.02
    units = 32767.0 / extent
    pos = np.rint(mpc * units).astype(np.int16)
    err = float(np.abs(pos.astype(np.float32) / units - mpc).max())
    print(f'\nTotal: {n:,} objetos')
    print(f'Escala global: {units:.4f} u/Mpc  (extensión {extent:.0f} Mpc, '
          f'error máx {err * 1000:.1f} kpc)')

    order = morton_order(pos)
    pos_m = pos[order]

    # Referencia: cuanto pesaria en bloques contiguos (mejor compresion, peor experiencia).
    contiguous = sum(len(encode_chunk(c)) for c in np.array_split(pos_m, args.chunks))

    manifest = {
        'version': 3,
        'count': int(n),
        'zMax': max_redshift,
        'unitsPerMpc': units,
        'partition': 'position-hash-v1',
        'tracers': tracers,
        'tracerCounts': {name: int((type_all == code).sum())
                         for name, code in tracer_codes.items()},
        'chunks': []
    }
    if target_flags is not None:
        manifest['typeEncoding'] = 'dr1-target-bits-v1'
        manifest['targetCounts'] = {name: int(np.count_nonzero(target_flags & (1 << bit)))
                                    for bit, name in enumerate(
                                        ('BGS', 'BGS_BRIGHT', 'LRG', 'ELG', 'QSO_TARGET', 'OTHER'))}
        # Bits 0-5 retain overlapping targeting labels; bit 6 is spectral QSO.
        type_all = target_flags | (type_all << 6)
    fog_index_dir = Path(args.fog_index_dir) if args.fog_index_dir else None
    if fog_index_dir:
        fog_index_dir.mkdir(parents=True, exist_ok=True)
    total = 0
    ttotal = 0
    position_hash = np.abs(pos[:, 0].astype(np.int64) * 17 +
                           pos[:, 2].astype(np.int64) * 59 +
                           pos[:, 1].astype(np.int64) * 101)
    for i in range(args.chunks):
        # Partición espacial determinista. Cualquier cliente puede calcular, solo desde la
        # posición, qué chunk contiene un punto; así el picker móvil nunca selecciona filas de
        # los chunks que decidió no cargar. Dentro de cada bucket mantenemos Morton para gzip.
        chosen = np.flatnonzero(position_hash % args.chunks == i)
        local_order = morton_order(pos[chosen])
        chosen = chosen[local_order]
        chunk = pos[chosen]
        blob = encode_chunk(chunk)
        name = f'chunk_{i:03d}{tag}.bin'
        (out_dir / name).write_bytes(blob)

        # Un byte por objeto conserva el tipo original tras mezclar los catálogos.
        # Se descarga después de las posiciones para no retrasar el primer fotograma.
        tname = f'type_{i:03d}{tag}.bin'
        tblob = gzip.compress(type_all[chosen].tobytes(), 9)
        (out_dir / tname).write_bytes(tblob)

        manifest['chunks'].append({'file': name, 'count': int(len(chunk)),
                                   'bytes': len(blob),
                                   'type': tname, 'typeBytes': len(tblob)})
        if fog_index_dir:
            # Este archivo vive solo en el entorno de procesamiento: permite generar la
            # corrección FoG con el orden EXACTO de la geometría publicada, sin enviar
            # TARGETID ni un hash de millones de filas al navegador.
            iname = f'targetid_{i:03d}{tag}.npy'
            np.save(fog_index_dir / iname, target_all[chosen])
            manifest['chunks'][-1]['fogTargetIndex'] = iname
        total += len(blob)
        ttotal += len(tblob)
        print(f'  [write] {name}: {len(chunk):,} puntos, {len(blob) / 1e6:.2f} MB'
              f'  + tipo {len(tblob) / 1e6:.2f} MB')

    (out_dir / 'manifest.json').write_text(json.dumps(manifest, indent=1))

    first = manifest['chunks'][0]
    print(f'\nTotal   : {total / 1e6:.1f} MB en {args.chunks} trozos')
    print(f'Tipos   : {ttotal / 1e6:.1f} MB aparte (asíncrono, no bloquea el arranque)')
    print(f'Contiguo: {contiguous / 1e6:.1f} MB '
          f'(+{(total - contiguous) / 1e6:.1f} MB es lo que cuesta el intercalado)')
    print(f'Primero : {first["bytes"] / 1e6:.2f} MB, {first["count"]:,} puntos')
    print(f'Salida  : {out_dir}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
