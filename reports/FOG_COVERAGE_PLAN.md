# FoG DR1: auditoría de cobertura y plan de ampliación

Fecha: 2026-09-23  
Alcance: cobertura FoG del catálogo PlayCanvas de 15.786.217 TARGETID. Este informe usa los FITS Gfinder DR1 v1.0 ya descargados en `D:\desi\data_external\desi\dr1\gfinder\v1.0` y las salidas locales de `D:\desi\data_external\fog`. No cambia binarios publicados ni DELTAS activos.

## Resultado actual

El materializador DR1 reporta 595.821 TARGETID FoG presentes en el catálogo exacto de 15.786.217 filas (3,77 %); 195.200 reciben un delta radial no nulo (1,24 % del catálogo). El análisis de origen produjo 595.823 TARGETID; al cruzarlos con las 12 listas exactas del catálogo DR1, dos no están presentes y quedan 595.821. La cobertura FoG actual, por tanto, describe el alcance del archivo de deltas ya calculado, no el alcance máximo de los FITS Gfinder.

El cruce por clase da 589.608 filas `DR1_GALAXY` (194.006 no nulas) y 6.213 filas `DR1_QSO` (1.194 no nulas). Conviene reportar ambas clases por separado y validar los QSO: Gfinder es un catálogo de galaxias y la calidad/interpretación de su redshift no necesariamente es equivalente para un QSO.

El delta no nulo es pequeño para la mayoría de los objetos. En valor absoluto, sobre los 195.200 valores no nulos: mediana 3,06 Mpc; p90 9,27; p95 11,84; p99 17,94; p99,9 29,09; máximo 175,81 Mpc. El 48,6 % de los 28.080 grupos corregidos queda en el piso `alpha=0,08` (`13.645` grupos). El piso no produce un límite de 0,08 Mpc: define la razón entre dispersión transversal y radial usada en la contracción.

El corte actual `maxDeltaMpc=300` no es el motivo de la baja cobertura ni limita los valores generados: el máximo observado es 175,81 Mpc; hay 0 deltas por encima de 200 o 300 Mpc. Sí hay 117 objetos por encima de 35 Mpc, 52 por encima de 50, 11 por encima de 100 y dos por encima de 150 Mpc. El máximo es TARGETID `39628423219380830`, ZSPEC `0,21545`, asociado por Gfinder al grupo `IGRP=3129` (RICH=71, zgrupo `0,1679`); en el análisis ese grupo tiene 18 miembros espectroscópicos, `sigma_parallel=7,07 Mpc`, `sigma_perp=0,377 Mpc`, elongación 18,76 y `alpha=0,08`. Su delta equivale aproximadamente a 10.900 km/s usando `v ≈ H(z) |delta_chi|/(1+z)`. Es un caso extremo para auditar como posible miembro interlopante/error de redshift; no se debe presentar como un desplazamiento peculiar confirmado.

## Por qué se queda en 595.821

La cadena actual tiene dos restricciones separadas:

1. `tools/fog/match_gfinder.py` se ejecutó contra `viewer-target-index.npz` del visualizador LSS anterior: 9.751.955 filas. Escaneó 134.731.580 filas de Gfinder y obtuvo 5.507.091 TARGETID Gfinder únicos, equivalentes a 5.549.432 filas viewer (56,91 %, por filas con duplicados en viewer). Esos 5.507.091 IDs sí se cruzaron casi todos con el catálogo actual: 5.388.525 `DR1_GALAXY` + 117.430 `DR1_QSO` = 5.505.955 filas. Así que el salto a 15.786.217 filas no explica por sí solo los 595.821; una reindexación completa del mismo conjunto de coincidencias sólo recuperaría 1.136 filas más.
2. El análisis FoG produjo deltas sólo para miembros de grupos resueltos a partir de esas coincidencias. `collect_members.py` dejó fuera grupos con RICH<5. Entre los 5.137.633 grupos relevantes, 3.197.763 tienen RICH=1: al no tener otro miembro, no se puede medir dispersión interna para corregir FoG. 352.351 grupos tienen RICH>=5 y se recogieron 3.593.784 miembros; el subconjunto de TARGETID cruzados con la selección del visualizador que cae en esos grupos produjo 595.823 filas en `targetid-delta.npy` (dos no están en el catálogo DR1 actual). Luego `analyze_fog.py` exigió al menos cinco miembros con redshift espectroscópico (`ZSRC>0`) y elongación radial/transversal `>=1,5`; pasaron 28.080 grupos y 195.200 filas recibieron delta no nulo. El archivo contiene además filas de grupos de riqueza elegible que reciben delta cero.

Esto explica también por qué no basta con bajar `elongation`: en una sensibilidad ejecutada sobre los mismos miembros, cambiar `1,5` a `1,0` con `min-spec=5` subió de 28.080 a 28.217 grupos y de 195.200 a 195.769 valores no nulos. En cambio, cambiar `min-spec=5` a `3`, manteniendo elongación `1,5`, dio 65.117 grupos y 283.000 deltas no nulos; con elongación `1,0`, 65.907 grupos y 285.051 deltas no nulos. En las tres pruebas el cruce exacto al catálogo siguió en 595.821 TARGETID: amplían los objetos con delta no nulo dentro del conjunto ya procesado, pero no añaden nuevas coincidencias. Son sensibilidades exploratorias, no parámetros adoptados.

La tabla de grupos local muestra además 1.017.771 grupos de RICH=2, 390.517 de RICH=3 y 179.231 de RICH=4. Sus miembros no están materializados en `members/`, por el corte RICH>=5. Ampliar el match completo a los 5.505.955 objetos que ya tienen ID Gfinder es viable y amplía el número de coincidencias con Gfinder, pero muchos de esos objetos están en grupos de riqueza 1–4 y no obtienen por ello una medición de dispersión fiable. La extensión a RICH>=3 es una oportunidad de cobertura que debe validarse con los ZSPEC de cada grupo; RICH=1/2 no se debe convertir automáticamente en corrección.

## Caso de control: Abell 2162

Usando como consulta la posición de catálogo Abell A2162, aproximadamente RA=243,125°, Dec=+29,533° y z≈0,0310 ([referencia de la lista de supercúmulos](https://www.icc.dur.ac.uk/~tt/Lectures/Galaxies/LocalGroup/Back/nearsc.html)), el grupo Gfinder cercano en redshift es `IGRP=4641737`, RICH=2, z=0,0320, a 0,045°; queda descartado directamente por `minimum-richness=5`. También aparece `IGRP=52780607`, RICH=1, z=0,0295, a 0,016°. El grupo RICH=5 cercano `IGRP=964226` está a z=0,0549 y sólo tiene 1 miembro ZSPEC de 5 en los miembros recogidos, por lo que falla `min-spec=5`; no es un grupo FoG validado para Abell 2162. Esto explica los cortes locales, pero no demuestra que alguno de esos grupos sea físicamente parte de A2162. Para cerrar el caso hay que cruzar membresías Gfinder individuales, ZSPEC, máscara y una referencia externa del cúmulo antes de asignar correcciones.

## Auditoría de límites físicos

Un límite uniforme de 35 Mpc no sale directamente de la bibliografía de FoG revisada. Como chequeo de escala, si se adopta hipotéticamente un máximo de velocidad peculiar `v_max`, la conversión de primer orden para desplazamiento comóvil radial es:

```text
|delta_chi|max(z) ≈ (1 + z) v_max / H_Planck18(z)
```

Con Planck18, `v_max=2.000–2.500 km/s` equivale a 29,7–37,1 Mpc en z=0,008; 31,0–38,7 Mpc en z=0,1; y 33,2–41,4 Mpc en z=1. Por tanto, 35 Mpc aproxima un corte de velocidad de alrededor de 2.300 km/s, pero no representa exactamente el mismo máximo a todo z. Además, una velocidad peculiar individual de 2.000–2.500 km/s sería un límite duro que requiere evidencia para este catálogo; no se adoptó. El exceso del caso máximo sobre esa escala debe investigarse en membresía, `ZERR`/calidad espectral y dispersión del grupo.

Tegmark et al. describen una compresión FoG por grupos: identifican grupos con un algoritmo FoF, comparan dispersión radial y transversal, comprimen radialmente cuando la primera supera a la segunda y prueban la sensibilidad al umbral de densidad. No proponen un límite universal de 300 Mpc ni un tope universal de 35 Mpc. El algoritmo local usa membresías Gfinder y MAD robusta en vez del mismo FoF/desvío estándar, por lo que es una aproximación visual exploratoria, no una reproducción exacta de Tegmark.

## Plan de implementación

### 1. Medir el techo de coincidencias Gfinder en el catálogo actual

- Construir un índice ordenado y deduplicado de los TARGETID de las 15.786.217 filas DR1; conservar el tipo/trazador y redshift de cada fila.
- Repetir el escaneo de las 134.731.580 filas Gfinder usando el índice actual, codificando el TARGETID a partir de `RELEASE`, `BRICKID` y `OBJID`; no usar tolerancia angular ni emparejamiento por redondeo.
- Registrar filas Gfinder escaneadas, coincidencias, TARGETID únicos, duplicados, coincidencias por trazador, por redshift y por `ZSRC`. Éste da la cobertura potencial de objetos que Gfinder realmente identifica dentro del catálogo actual; no es todavía la cobertura con delta FoG.
- Comprobación de regresión: con el índice LSS anterior, reproducir 5.507.091 filas Gfinder / 5.549.432 filas viewer antes de aceptar el índice nuevo. Validar unicidad del índice nuevo y contar explícitamente duplicados antes de producir deltas.

### 2. Recuperar membresías oficiales completas para las coincidencias actuales

- Unir Gfinder `IGAL -> IGRP` por `iDESIDR9.y1.v1_1.fits` y las propiedades por `IGRP` en `DESIDR9.y1.v1_group.fits`. Sólo asignar grupo cuando existe esa relación oficial.
- Recolectar las membresías de grupos RICH>=3 para una sensibilidad basada en el rango de uso del método publicado; mantener aparte RICH=1/2 y miembros sin ZSPEC. No inferir grupos con vecinos espaciales.
- Usar sólo miembros con `ZSRC>0` al estimar la dimensión radial. `ZSRC=0` identifica redshift fotométrico en la documentación oficial Gfinder y no tiene la precisión radial requerida para medir dispersión FoG.

### 3. Calibrar el estimador y los cortes, sin convertir sensibilidad en resultado

- Calcular distancias comóviles con Planck18 y dispersión radial/perpendicular respecto del centro y zgrupo Gfinder. Comparar MAD escalada (actual) contra desviación estándar robusta/truncada, y registrar efectos de `ZERR`/calidad espectral e interlopantes.
- Publicar una matriz exploratoria de RICH>=3/5, `min-spec`=3/5 y elongación=1,0/1,5/2,0. Los cortes actuales 5, 1,5 y `alpha-min=0,08` son decisiones locales sin calibración publicada específica para DR1; no cambiar el perfil activo hasta contrastarlo con mocks, cúmulos conocidos y estabilidad al remuestreo.
- Comparar corrección por grupo con dispersión radial observada, dispersión residual tras corrección, miembros retenidos y fracción de grupos que cambia al retirar uno de sus espectros. Grupos con muestra insuficiente o dispersión no positiva reciben delta cero.
- Tratar el límite de desplazamiento por separado. Reportar percentiles y número de casos arriba de 35/50/100 Mpc. Si se estudia un cap de velocidad, convertirlo por redshift mediante Planck18, comparar sensibilidad a 2.000 y 2.500 km/s y confirmar contra distribución espectroscópica del grupo/mocks antes de fijarlo. El cap actual 300 Mpc no actúa en los valores existentes; cambiarlo por un valor visual sería injustificado.

### 4. Criterio de validación antes de sustituir FoG

- No se crea ningún grupo: cada delta debe trazar a un `IGAL`, relación `IGAL/IGRP` y fila de propiedades Gfinder verificables.
- Reportar denominadores separados: TARGETID actual; match con Gfinder; miembro de grupo; miembro ZSPEC; grupo que supera geometría; miembro con delta no nulo. Incluir totales y desglose por tracer/redshift.
- El archivo de salida debe indicar el perfil y parámetros, conservar delta cero para los no elegibles, tener todos los valores finitos y respetar el orden TARGETID/chunk exacto. Repetir con al menos un caso rico conocido, incluyendo Abell 2162, y mostrar sus membresías/cortes fallidos antes y después.
- Mostrar cobertura y distribución de deltas junto a la comparación geométrica. Aceptar un perfil sólo si pasa la validación de membresía y dispersión; si no, mantener la corrección como parcial y exploratoria.

## Fundamento y límites de las fuentes

- [Documentación oficial Gfinder DR1](https://data.desi.lbl.gov/doc/releases/dr1/vac/gfinder/): define los tres archivos, campos `IGAL`, `IGRP`, `RICH`, `GRP_Z`, y que `ZSRC=0` es `PHOTZ` y `ZSRC>0` es `SPECZ`. Describe Gfinder DR1 como actualización del catálogo halo-based extendido con espectros Iron; la selección se limita a MAG_Z<21.
- [Yang et al., método de grupo halo-based extendido](https://arxiv.org/abs/2012.14998): evaluación sobre mocks Legacy DR8 (z<=1, MAG_Z<=21); informa completitud de miembros y pureza condicionadas a masa y grupos con al menos 3 miembros. Sirve como antecedente del umbral RICH>=3 para una prueba de sensibilidad, no como calibración directa del catálogo DR1.
- [Tegmark et al. 2004, sección 3.3 y figura 7](https://web.physics.rutgers.edu/grad/690/Tegmark-etal-2004.pdf): define compresión radial de grupos FoG según dispersión radial y transversal y explora umbrales. Apoya la forma del estimador geométrico y la obligación de estudiar sensibilidad, no los cortes RICH>=5, cinco espectros, elongación>=1,5, alpha-min=0,08 o el límite 300 Mpc usados por esta implementación local.

## Artefactos consultados

- `D:\desiV2\ply\reconstruction_dr1\payloads\coverage-report.json`
- `D:\desiV2\tools\rsd\materialize_dr1_corrections.py`
- `D:\desi\data_external\fog\matches\summary.json`
- `D:\desi\data_external\fog\groups\summary.json`
- `D:\desi\data_external\fog\members\summary.json`
- `D:\desi\data_external\fog\analysis\summary.json`, `targetid-delta.npy`, `group-stats.npy`
- Sensibilidades temporales ejecutadas con `analyze_fog.py`: (5 espectros, elongación 1,0), (3, 1,5), (3, 1,0). No modificaron los deltas activos ni los artefactos publicados.

Comandos PowerShell para reproducir las tres sensibilidades (requieren el Python Blender y paquetes aislados indicados en el handoff FoG):

```powershell
$py = 'C:\Program Files\Blender Foundation\Blender 4.5\4.5\python\bin\python.exe'
$env:PYTHONPATH = 'D:\desi\.python-packages'
$script = 'D:\desiV2\tools\fog\analyze_fog.py'
$common = @('--members','D:\desi\data_external\fog\members', '--groups','D:\desi\data_external\fog\groups', '--matches','D:\desi\data_external\fog\matches')
& $py $script @common --out "$env:TEMP\fog_sensitivity_s5_e1" --min-spec 5 --elongation 1.0
& $py $script @common --out "$env:TEMP\fog_sensitivity_s3_e1p5" --min-spec 3 --elongation 1.5
& $py $script @common --out "$env:TEMP\fog_sensitivity_s3_e1" --min-spec 3 --elongation 1.0
```
