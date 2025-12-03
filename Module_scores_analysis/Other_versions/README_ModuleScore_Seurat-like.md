# Module Score Analysis for Bulk RNA-seq

Este repositorio contiene scripts para calcular **module scores** en datos de Bulk RNA-seq, similar a la función `AddModuleScore` de Seurat para scRNA-seq.

## Contexto

El método `AddModuleScore` de Seurat calcula un score para cada célula basado en la expresión media de un set de genes, restando la expresión media de genes control pareados por nivel de expresión. Esto corrige el sesgo de que genes altamente expresados contribuyan más al score.

Para **Bulk RNA-seq**, adaptamos este enfoque:
- **Score = mean(genes del set) - mean(genes control pareados por expresión)**

## Scripts disponibles

| Script | Método | Descripción |
|--------|--------|-------------|
| `module_score.py` | Mean / Seurat-like | Cálculo básico con heatmap |
| `module_score_gsva.py` | GSVA | Gene Set Variation Analysis |
| `module_score_visualization.py` | Seurat-like | **Visualizaciones tipo paper (boxplot, violin, etc.)** |

## Uso rápido

### 1. Preparar archivos

**Matriz de expresión VST** (`condition_21mo_shNsun2_vs_21mo_shLuci.txt`):
```
external_gene_name    sample1    sample2    sample3    ...
Actb                  12.5       11.8       12.2       ...
Gapdh                 14.2       13.9       14.1       ...
...
```

**Listas de genes** (un gen por línea):
```
# 2c_DOWN.txt
Gene1
Gene2
Gene3
```

**Archivo de grupos** (opcional, `groups.txt`):
```
sample_id	group
sample1	shLuci
sample2	shLuci
sample3	shNsun2
sample4	shNsun2
```

### 2. Ejecutar análisis

```bash
# Visualizaciones completas (similar a Figura 3G/3H)
python module_score_visualization.py \
    --vst condition_21mo_shNsun2_vs_21mo_shLuci.txt \
    --gene-lists 2c_DOWN.txt 2c_UP.txt 7c_UP.txt 7c_DOWN.txt \
    --gene-column external_gene_name \
    --sample-groups groups.txt \
    --output-prefix results/module_scores \
    --plots all

# O inferir grupos desde nombres de muestras
python module_score_visualization.py \
    --vst condition_21mo_shNsun2_vs_21mo_shLuci.txt \
    --gene-lists 2c_DOWN.txt 2c_UP.txt 7c_UP.txt 7c_DOWN.txt \
    --gene-column external_gene_name \
    --infer-groups \
    --output-prefix results/module_scores
```

### 3. Outputs generados

- `module_scores.tsv` - Tabla de scores (gene sets × samples)
- `module_scores_boxplot.png` - **Boxplots por condición** (estilo Fig 3G/3H)
- `module_scores_violin.png` - Violin plots por condición
- `module_scores_heatmap.png` - Heatmap con clustering
- `module_scores_barplot.png` - Barras agrupadas con error bars
- `module_scores_dotplot.png` - Dot plot estilo Seurat

## Visualizaciones

### Boxplot (recomendado para comparaciones)
Similar a las figuras 3G/3H del paper, muestra la distribución de module scores por condición con test estadístico (t-test):

```bash
python module_score_visualization.py \
    --vst data.txt \
    --gene-lists genelist.txt \
    --sample-groups groups.txt \
    --plots boxplot \
    --output-prefix fig3G_style
```

### Violin plot
Muestra la distribución completa de los scores:

```bash
python module_score_visualization.py \
    --vst data.txt \
    --gene-lists genelist.txt \
    --sample-groups groups.txt \
    --plots violin \
    --output-prefix fig3H_style
```

## Métodos de scoring

### Seurat-like (default)
```
Score = mean(genes en set) - mean(genes control)
```
Los genes control se seleccionan de bins de expresión similar para corregir sesgos.

### Mean simple
```
Score = mean(genes en set) - mean global por muestra
```
Más simple pero puede tener sesgos si los genes del set tienen expresión atípica.

## Parámetros importantes

| Parámetro | Default | Descripción |
|-----------|---------|-------------|
| `--method` | seurat-like | Método de scoring |
| `--bins` | 24 | Número de bins de expresión para controles |
| `--seed` | 42 | Semilla para reproducibilidad |
| `--format` | png | Formato de salida (png/pdf/svg) |

## Requisitos

```bash
pip install pandas numpy seaborn matplotlib scipy
# Para GSVA (opcional):
pip install gseapy
```

## Interpretación

- **Score > 0**: El set de genes está sobre-expresado en esa muestra/condición
- **Score < 0**: El set de genes está sub-expresado
- **Score ≈ 0**: Expresión similar al background

Los p-valores en los boxplots indican si la diferencia entre grupos es estadísticamente significativa.

## Referencia

- Método AddModuleScore: Tirosh et al., Science (2016)
- GSVA: Hänzelmann et al., BMC Bioinformatics (2013)
- Paper de ejemplo: Lu et al., Cell (2025) - "Prevalent mesenchymal drift in aging and disease is reversed by partial reprogramming"
