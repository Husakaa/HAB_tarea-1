# Análisis de Enriquecimiento Funcional de Genes

Pipeline completo Python + R para análisis de enriquecimiento funcional utilizando g:Profiler.

## Descripción

Este proyecto realiza análisis de enriquecimiento funcional de genes consultando tres bases de datos principales:

- **GO:BP** (Gene Ontology - Biological Process)
- **KEGG** (Kyoto Encyclopedia of Genes and Genomes)
- **REAC** (Reactome Pathway Database)

El análisis se realiza en dos fases: Python para la obtención de datos y R para la visualización y análisis estadístico.

## Instalación (recomendado con entorno virtual)

1. `python -m venv .venv`
2. `source .venv/bin/activate`  (Linux/Mac) ó 
   `.venv\Scripts\activate`     (Windows)
3. `pip install -r requirements.txt`

### Requisitos

- Python 3.7+
- R 4.0+
- RStudio (recomendado)
 


### Dependencias Python
```bash
pip install -r requirements.txt
```

### Dependencias R
```r
install.packages(c("knitr", "readr", "tidyr", "dplyr", 
                   "kableExtra", "ggplot2", "scales", "RColorBrewer"))
```

## Estructura del Proyecto
```
proyecto/
├── data/
│   └── genes_input.txt
├── scripts/
│   └── enrichment_analysis.Rmd    
│   └── functional_analysis.py
├── results/
│   ├── validacion_genes.csv
│   ├── resultados_completos.csv
│   ├── GO_BP_resultados.csv
│   ├── KEGG_resultados.csv
│   ├── REAC_resultados.csv
│   └── R/                         
│       ├── *.png
│       └── *.csv
├── requirements.txt
├── README.md
```

## Uso

### Paso 1: Análisis funcional (Python)
```bash
python functional_analysis.py -i genes_input.txt
```

**Opciones disponibles:**
```bash
# Con archivo de entrada
python scripts/functional_analysis.py -i data/genes_input.txt

# Con genes directos
python scripts/functional_analysis.py -g COX4I2 ND1 ATP6

# Especificar directorio de salida
python scripts/functional_analysis.py -i data/genes_input.txt -o mi_directorio

# Cambiar organismo (default: hsapiens)
python scripts/functional_analysis.py -i data/genes_input.txt -org mmusculus
```

### Paso 2: Análisis estadístico y visualización (R)
```bash
cd scripts
# Abrir enrichment_analysis.Rmd en RStudio y hacer clic en "Knit"
```

O desde terminal:
```bash
cd scripts
Rscript -e "rmarkdown::render('enrichment_analysis.Rmd')"
```

## 📊 Resultados

### Archivos generados por Python (`results/`)

- `validacion_genes.csv` - Validación y conversión de genes
- `resultados_completos.csv` - Todos los términos enriquecidos
- `GO_BP_resultados.csv` - Términos Gene Ontology
- `KEGG_resultados.csv` - Vías metabólicas KEGG
- `REAC_resultados.csv` - Vías Reactome
- `resultados_completos.xlsx` - Excel con todas las hojas

### Archivos generados por R (`results/R/`)

**Visualizaciones (PNG, 300 DPI):**
- `distribucion_pvalores.png`
- `top20_terminos.png`
- `go_bp_top15.png`
- `kegg_pathways.png`
- `reactome_pathways.png`
- `comparacion_bases_datos.png`
- `dotplot_integrado.png`

**Tablas:**
- `resumen_estadisticas_bd.csv`
- `top50_terminos.csv`
- `terminos_mitocondriales.csv`

**Informe:**
- `enrichment_analysis.html` - Informe completo interactivo

## Formato de Entrada

El archivo de genes puede tener dos formatos:

**Opción 1: Genes separados por comas**
```
COX4I2, ND1, ATP6
```

**Opción 2: Un gen por línea**
```
COX4I2
ND1
ATP6
```

## Genes Mitocondriales

El script busca automáticamente genes mitocondriales con el prefijo `MT-`:

- `ND1` → se busca como `MT-ND1`
- `ATP6` → se busca como `MT-ATP6`
- `COX1` → se busca como `MT-CO1`

## Ejemplo Completo
```bash
# 1. Crear archivo con genes
echo "COX4I2, ND1, ATP6" > genes_input.txt

# 2. Ejecutar análisis funcional
python scripts/functional_analysis.py -i data/genes_input.txt

# 3. Generar informe R
cd scripts
# Abrir enrichment_analysis.Rmd en RStudio → Knit

# 4. Visualizar resultados
open ../results/enrichment_analysis.html  # macOS
xdg-open ../results/enrichment_analysis.html  # Linux
```

## Solución de Problemas

**Error: "No se encuentra validacion_genes.csv"**
- Ejecutar primero el script Python

**Error: "Paquete no encontrado"**
- Instalar paquetes R faltantes: `install.packages("nombre_paquete")`

**No se encontraron términos enriquecidos**
- Normal si los genes no están relacionados funcionalmente
- Verificar ortografía de los genes

## 📚 Referencias

- [g:Profiler](https://biit.cs.ut.ee/gprofiler/)
- [Gene Ontology](http://geneontology.org/)
- [KEGG](https://www.genome.jp/kegg/)
- [Reactome](https://reactome.org/)

## 📄 Licencia

Este proyecto está bajo la Licencia MIT.

## 👥 Autor

Hugo Salas Calderón - [hugosalascalderon@gmail.com](mailto:hugosalascalderon@gmail.com)

