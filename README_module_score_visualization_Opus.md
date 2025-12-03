# Module Score Visualization (Opus)

Este README cubre el script `module_score_visualization_Opus.py` — una implementación adaptada de la lógica Seurat `AddModuleScore` para datos de Bulk RNA-seq, con visualizaciones estilo figura 3G/3H de artículos (boxplot/violin) y gráficos complementarios (heatmap, barplot, dotplot).

**Resumen**
- Método: Seurat-like AddModuleScore adaptado a Bulk RNA-seq.
- Cálculo: Score = mean(genes del set) - mean(genes control pareados por expresión).
- Adaptación para Bulk: control:target = 1:1 por defecto (en Seurat scRNA-seq es 100:1).
- Visualizaciones: `boxplot`, `violin`, `heatmap`, `barplot`, `dotplot`.
- Atención: el script incluye advertencias y consejos cuando hay pocos replicados (p. ej. n=2 por grupo).

Requiere: `pandas`, `numpy`, `seaborn`, `matplotlib`, `scipy`.
Instalación rápida:

```powershell
pip install pandas numpy seaborn matplotlib scipy
```

---

## Propósito

Generar scores de módulos (module scores) para conjuntos de genes y producir figuras listas para publicación que permitan comparar condiciones en datos de Bulk RNA-seq, siguiendo la filosofía de `AddModuleScore` (pareo por nivel de expresión para seleccionar genes control).

## Archivos de entrada

- Matriz VST/normalizada (`TSV`): filas = genes, columnas = muestras. Debe contener una columna con nombres de genes (por defecto `external_gene_name`).
- Listas de genes: uno por archivo, un gen por línea.
- Archivo de grupos (opcional, TSV): columnas `sample_id` y `group`.

Formato de ejemplo `groups.txt` (tab-sep):

```
sample_id	group
sample1	shLuci
sample2	shLuci
sample3	shNsun2
sample4	shNsun2
```

Si no proporcionas `--sample-groups`, puedes usar `--infer-groups` para que el script intente inferir grupos por nombre de muestra.

---

## Uso (ejemplos)

Comando completo (generar todas las visualizaciones):

```powershell
python module_score_visualization_Opus.py \
  --vst condition_21mo_shNsun2_vs_21mo_shLuci.txt \
  --gene-lists 2c_DOWN.txt 2c_UP.txt 7c_UP.txt 7c_DOWN.txt \
  --gene-column external_gene_name \
  --sample-groups groups.txt \
  --output-prefix results/module_scores \
  --plots all
```

Inferir grupos desde nombres de muestras:

```powershell
python module_score_visualization_Opus.py \
  --vst expression_vst.txt \
  --gene-lists genelist1.txt genelist2.txt \
  --infer-groups \
  --output-prefix results/figure3_style
```

Generar sólo boxplot o violin (estilo Fig 3G/3H):

```powershell
python module_score_visualization_Opus.py \
  --vst data.txt \
  --gene-lists genelist.txt \
  --sample-groups groups.txt \
  --plots boxplot \
  --output-prefix fig3G_style

python module_score_visualization_Opus.py \
  --vst data.txt \
  --gene-lists genelist.txt \
  --sample-groups groups.txt \
  --plots violin \
  --output-prefix fig3H_style
```

Opciones clave:
- `--method`: `seurat-like` (por defecto) o `mean`.
- `--bins`: número de bins para parear expresión (por defecto 24).
- `--seed`: semilla aleatoria (por defecto 42).
- `--format`: `png`, `pdf` o `svg`.

---

## Salidas generadas

Con `--output-prefix results/module_scores` se crean:
- `results/module_scores.tsv` — tabla de module scores (gene sets × samples)
- `results/module_scores_boxplot.png` — boxplots por condición
- `results/module_scores_violin.png` — violin plots por condición
- `results/module_scores_heatmap.png` — heatmap con clustering
- `results/module_scores_barplot.png` — barras con error bars (mean ± SEM)
- `results/module_scores_dotplot.png` — dotplot tipo Seurat

---

## Interpretación y recomendaciones (especialmente para n pequeño)

- Score > 0: el set de genes está relativamente sobre-expresado.
- Score < 0: el set de genes está relativamente sub-expresado.
- Score ≈ 0: sin desviación del background pareado.

Advertencias cuando hay pocos replicados:
- El script imprime una advertencia si algún grupo tiene n ≤ 2.
- Con n=2 por grupo los p-valores tienen baja potencia; interpreta los resultados como exploratorios.
- Recomendaciones para n=2:
  - Prioriza la magnitud del efecto (|mean(group1) - mean(group2)|) y la dirección consistente entre sets.
  - Evita basar conclusiones firmes únicamente en p-valores; valida con más muestras cuando sea posible.
  - Reporta claramente en métodos que los resultados son preliminares y basados en n pequeño.

---

## Lógica técnica (resumen)

1. Calcula la media de expresión por gen y crea bins por cuartiles/quantiles (hasta `--bins`).
2. Para cada gen del set objetivo selecciona genes control del mismo bin (excluyendo genes del set).
3. Control por defecto: 1 gen control por cada gen del set (parámetro `ctrl_genes_per_gene` en la función interna).
4. Score por muestra: media(target genes) - media(control genes).

Para mayor robustez en scRNA-seq se puede aumentar el ratio control:target (p. ej. 100:1), pero en Bulk esto puede sobremuestrear y no es recomendable.

---

## Dependencias

```powershell
pip install pandas numpy seaborn matplotlib scipy
```

Opcional (si usas análisis de enriquecimiento alternativo):

```powershell
pip install gseapy
```

---

## Referencias

- AddModuleScore (Seurat): Tirosh et al., Science (2016)
- GSVA: Hänzelmann et al., BMC Bioinformatics (2013)
- Ejemplo de paper con figuras tipo 3G/3H: Lu et al., Cell (2025)

---

Si quieres, actualizo también la tabla de scripts en el README principal para referenciar este archivo y su uso. ¿Lo añado?