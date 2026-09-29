# Auditoría de máscara RSD DR1 v1.5 — 2026-09-24

## Decisión para los assets

**Los assets de extensión están activos en la escena PlayCanvas para evaluación visual, sin aprobación científica.** El usuario pidió implementar los cambios antes de la validación posterior y luego mantener solo los nuevos archivos en el proyecto. La carpeta `ply/reconstruction_dr1/staged_extensions_v1` contiene 12 payloads RSD por chunk, 12 FoG equivalentes al baseline en bytes descomprimidos y dos pares LOD FoG/RSD. `playcanvas-review-deployment.json` documenta los 15 assets nuevos (incluido overlay), los 15 anteriores retirados y el Launch de comprobación. El baseline retirado está respaldado localmente con hashes coincidentes con PlayCanvas. `asset-staging-manifest.json` conserva `readyForUpload: false` como gate científico, aunque los archivos ya están subidos para inspección. La extensión tiene reproducibilidad de campo y soporte 3-D proxy, pero carece de validación independiente de selección angular/radial. `materialize_dr1_corrections.py --stage-unvalidated-rsd` habilita salida local de revisión; el modo de exportación aprobado sigue exigiendo `--mask-approval-report` ligado por hash al manifest exacto.

## Fuentes y trazabilidad

- Selección: DESI DR1 LSS `iron/LSScats/v1.5`, cinco tracers en NGC/SGC, random índice 0 asociado a cada clustering catalog, cortes z semiabiertos de `rsd_config.json` y peso positivo finito `WEIGHT × WEIGHT_FKP`.
- Archivos locales y tamaños: `D:\desi\data_external\RSD\source-manifest.json`. Las diez fuentes random y diez fuentes data usadas por la reconstrucción ya estaban disponibles. Tras liberar espacio, se descargó la realización random índice 1 para los diez campos en `tools/rsd/validation_sources`: **11.208.458.880 bytes** en total. `source-manifest-index1.json` registra URL, tamaño, SHA256, filas y `DESIDR=dr1` de cada FITS; columnas y estructura verificadas. Esta segunda realización contrasta el muestreo del índice 0, sin reemplazar una máscara oficial independiente. Los archivos `full`/`noveto`/`HPmapcut` requieren todavía una elección científica precisa de productos antes de descargar variantes de varios GB.
- DESI recomienda los productos `clustering` para análisis de clustering y aplicar a los randoms los mismos cortes que a los datos: https://data.desi.lbl.gov/doc/releases/dr1/ y https://desidatamodel.readthedocs.io/en/latest/DESI_ROOT/vac/RELEASE/lss/VERSION/LSScats/clustering/ . El catálogo `full` distingue explícitamente versiones antes y después de vetos angulares: https://desidatamodel.readthedocs.io/en/latest/DESI_ROOT/survey/catalogs/RELEASE/LSS/SPECPROD/LSScats/VERSION/index.html .
- El directorio oficial v1.5 contiene catálogos `full`, `full_noveto`, `full_HPmapcut` y mapas de propiedades HEALPix; esos mapas de propiedades no equivalen, por sí solos, a una máscara de aceptación por tracer, hemisferio y z. Antes de cualquier descarga masiva, un equipo DESI debe fijar qué versión de veto y función de selección reproduce exactamente cada `clustering` de esta corrida.

## Procedimiento ejecutado

1. Confirmar `extensions/extension-manifest.json` con diez campos completos, hashes de cada `.npz` y los gates de origen/malla/caché. `audit_rsd_extensions.py` produjo `final-extension-audit.json` con estado `complete`; esto solo comprueba integridad y reproducción del mismo campo.
2. Usar `audit_rsd_support_by_field.py` y `rsd-field-support-audit.json` para soporte tridimensional: randoms oficiales índice 0 completos, split par/impar, radio p99 de vecino más cercano y límites de la malla. El p99 es una heurística, no el borde DESI.
3. Ejecutar `audit_extension_angular_mask.py`: muestra fija de 10.000 TARGETID por campo (semilla 20260924); reconstruye su posición desde la geometría `.desi` v2 y el índice exacto por chunk; elige hasta 500.000 randoms válidos por campo; divide train/holdout par/impar; calcula en la esfera los percentiles 99 de distancia al primer y octavo vecino del train para el holdout. Consulta la muestra de extensión y registra ambos criterios. La distancia al octavo vecino detecta baja densidad angular; no clasifica inequívocamente un agujero de veto.
4. La auditoría escrita en `extensions/angular-mask-audit.json` revisó 100.000 posiciones. **5.798** superaron al menos uno de los dos límites p99. Por campo, rechazos/10.000: BGS NGC 445, SGC 505; LRG NGC 545, SGC 1.058; LRG+ELG NGC 334, SGC 584; ELG NGC 390, SGC 541; QSO NGC 564, SGC 832. Son tasas de la muestra, no conteos exactos de las extensiones completas.
5. Se repitió la misma muestra con random índice 1: `extensions/angular-mask-audit-index1.json`. Señaló **5.738/100.000**. `compare_angular_random_realizations.py` cruzó TARGETID exactos de los archivos privados de flags: **3.447** señalados por ambos índices, **8.089** por al menos uno. LRG SGC concentra 761 alertas persistentes de 10.000. Los 3.447 casos persistentes merecen prioridad de inspección, pero tampoco se clasifican automáticamente como fuera de máscara oficial. Informe: `extensions/angular-random-realization-comparison.json`.

## Reproducción local

En Windows, desde `D:\desiV2` con las fuentes en las rutas anteriores:

```powershell
$env:PYTHONPATH='D:\desi\.python-packages'
& 'C:\Program Files\Blender Foundation\Blender 4.5\4.5\python\bin\python.exe' tools/rsd/audit_extension_angular_mask.py
```

Para repetir el audit global desde el entorno DESI WSL cuando esté disponible:

```bash
cd /mnt/d/desiV2
/home/khyro/.venvs/desi-rsd/bin/python tools/rsd/audit_rsd_extensions.py
```

La ejecución Windows de la auditoría angular terminó sin errores; el script usa NumPy 2.4.6, SciPy 1.17.1 y Astropy 8.0.1. El entorno WSL no estuvo disponible en esta sesión (`E_ACCESSDENIED`), por lo que no se repitió aquí el audit global.

Para descargar y verificar la segunda realización random cuando se repita el trabajo:

```powershell
$env:PYTHONPATH='D:\desi\.python-packages'
& 'C:\Program Files\Blender Foundation\Blender 4.5\4.5\python\bin\python.exe' tools/rsd/download_validation_randoms_parallel.py
& 'C:\Program Files\Blender Foundation\Blender 4.5\4.5\python\bin\python.exe' tools/rsd/verify_validation_sources.py
& 'C:\Program Files\Blender Foundation\Blender 4.5\4.5\python\bin\python.exe' tools/rsd/audit_extension_angular_mask.py --random-index 1
& 'C:\Program Files\Blender Foundation\Blender 4.5\4.5\python\bin\python.exe' tools/rsd/compare_angular_random_realizations.py
```

El descargador compara `Content-Length` oficial con tamaño local, verifica la respuesta `Content-Range` de cada parte, reutiliza archivos completos y ensambla las partes comprobando el tamaño final. El descargador PowerShell secuencial `download_validation_randoms.ps1` es una alternativa simple. El verificador revisa columnas FITS, filas y SHA256, y escribe `validation_sources/source-manifest-index1.json`. Los dos audits usan la misma semilla y muestra de TARGETID; sus `.npz` de flags contienen identificadores privados y no son assets para subir.

## Revisión requerida por el equipo astronómico

1. Fijar la función de selección angular y radial oficial aplicable a los diez `clustering` v1.5, incluidos vetos, `HPmapcut`, completitud y tratamientos dependientes de tracer/hemisferio. Registrar versión, código, parámetros, hashes y procedencia de cada insumo. Si se requieren archivos adicionales, descargar los `full`/`full_noveto`/`full_HPmapcut` y auxiliares concretos de la versión aprobada, no asumir que un mapa de propiedades es la máscara.
2. Clasificar los 5.798 casos señalados por la muestra y una muestra de control no señalada: borde externo, hueco interno/veto, densidad aleatoria baja, posible cuantización de posición, o dentro de selección. Medir distancia angular al borde oficial, no solo al random más cercano.
3. Aplicar la máscara oficial a los **8.298.917** TARGETID, desglosando aceptación/rechazo por campo, hemisferio, z y distancia al borde. Comprobar en holdout espacial y bins de z los deltas reconstruidos contra una referencia independiente; comparar sensibilidad a random stride, pesos y malla.
4. Revisar los assets ya preparados conforme a `STAGED_ASSETS_2026-09-24.md`. Tras aprobación científica, escribir un informe JSON con `status: "approved"`, `extensionManifestSha256` del manifest exacto, autoría, fecha, versión de máscara, resultados por campo y referencia a la revisión. Regenerar con `--rsd-extensions`, `--mask-approval-report` y `--out` vacío, comparar hashes descomprimidos y conteos con el staging, revisar visualmente y después crear el overlay con los IDs reales de los assets subidos. La precedencia por tracer y z ya está codificada; los 4.166 cruces con cachés extranjeros conservan el delta de la extensión asignada.

## Límites

Los randoms `clustering` son una representación de la selección, no una prueba independiente de máscara; ambas realizaciones comparten la misma construcción y los mismos vetos. El muestreo de 500.000 randoms por campo ensancha las distancias de vecinos respecto al catálogo completo. Los 3.447 casos señalados de forma persistente no son automáticamente objetos inválidos, y los demás no quedan certificados. La posición reconstruida desde `.desi` está cuantizada. La coincidencia del campo con caché no valida velocidades verdaderas ni permite publicación científica.
