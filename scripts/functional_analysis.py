#!/usr/bin/env python3
"""
Script de Análisis Funcional de Genes con Validación de IDs
===========================================================

Realiza análisis de enriquecimiento funcional con validación y conversión
automática de identificadores de genes.

Flujo:
    - Acepta múltiples tipos de IDs (símbolos, Entrez, Ensembl, UniProt)
    - Control de errores para genes mal escritos o no encontrados
    - Valida y normaliza identificadores usando MyGene.info
    - Análisis de enriquecimiento con g:Profiler

Bases de datos utilizadas:
    - Gene Ontology (GO:BP): Procesos biológicos
    - KEGG: Vías metabólicas
    - Reactome: Vías de reacción biológicas

Uso:
    python analisis_funcional.py -i data/genes_input.txt
    python analisis_funcional.py -g COX4I2 ND1 ATP6
    python analisis_funcional.py -g 1327 4535 4508  # Entrez IDs
"""

import argparse
import pandas as pd
import sys
import os

try:
    from gprofiler import GProfiler
    import mygene
except ImportError as e:
    print("Error: Faltan dependencias. Instala con:")
    print("  pip install -r requirements.txt")
    sys.exit(1)


def leer_genes(archivo):
    """
    Lee genes desde un archivo de texto.
    
    Soporta múltiples formatos:
    1. Un gen por línea:
       COX4I2
       ND1
       ATP6
    
    2. Genes separados por comas (en una o varias líneas):
       COX4I2, ND1, ATP6
    
    3. Formato mixto:
       COX4I2, ND1
       ATP6
    
    Args:
        archivo: Ruta al archivo con genes
        
    Returns:
        Lista de genes
    """
    try:
        with open(archivo, 'r') as f:
            contenido = f.read()
        
        genes = []
        
        # Detectar si hay comas en el archivo
        if ',' in contenido:
            print(f"[INFO] Formato detectado: genes separados por comas")
            # Reemplazar saltos de línea por espacios, dividir por comas
            genes = [gene.strip() for gene in contenido.replace('\n', ',').split(',') if gene.strip()]
        else:
            print(f"[INFO] Formato detectado: un gen por línea")
            # Un gen por línea
            genes = [line.strip() for line in contenido.split('\n') if line.strip()]
        
        print(f"[INFO] Leídos {len(genes)} genes desde {archivo}")
        return genes
        
    except FileNotFoundError:
        print(f"[ERROR] Archivo no encontrado: {archivo}")
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] Error al leer archivo: {e}")
        sys.exit(1)


def validar_y_convertir_genes(genes, organismo='human'):
    """
    Valida y convierte identificadores de genes a símbolos oficiales.
    
    Usa MyGene.info para:
    - Validar que los genes existan
    - Convertir entre diferentes tipos de IDs (Entrez, Ensembl, símbolos)
    - Normalizar a símbolos oficiales (HGNC)
    - Detectar genes mal escritos o no encontrados
    - Buscar automáticamente variantes con prefijo MT- para genes mitocondriales
      (e.g., si "ND1" no se encuentra, intenta con "MT-ND1")
    
    Args:
        genes: Lista de identificadores de genes (pueden ser símbolos, Entrez IDs, etc.)
        organismo: Especie ('human', 'mouse', etc.)
        
    Returns:
        dict con:
            - 'validos': Lista de símbolos validados
            - 'mapping': Diccionario de conversión {input -> símbolo}
            - 'no_encontrados': Lista de genes no encontrados
            - 'advertencias': Lista de mensajes de advertencia
    """
    print(f"\n{'='*70}")
    print("VALIDACIÓN Y CONVERSIÓN DE IDENTIFICADORES")
    print(f"{'='*70}")
    print(f"[INFO] Validando {len(genes)} genes con MyGene.info...")
    
    # Inicializar cliente de MyGene
    mg = mygene.MyGeneInfo()
    
    # Determinar species code para MyGene
    species_map = {
        'human': 'human',
        'mouse': 'mouse',
        'rat': 'rat',
        'hsapiens': 'human',
        'mmusculus': 'mouse',
        'rnorvegicus': 'rat'
    }
    species = species_map.get(organismo.lower(), 'human')
    
    # Consultar MyGene.info
    # scopes: buscar en múltiples tipos de IDs
    try:
        results = mg.querymany(
            genes, 
            scopes='symbol,entrezgene,ensembl.gene,uniprot.Swiss-Prot',
            fields='symbol,entrezgene,name',
            species=species,
            returnall=True
        )
    except Exception as e:
        print(f"[ERROR] Error al conectar con MyGene.info: {e}")
        print("[INFO] Continuando sin validación...")
        return {
            'validos': genes,
            'mapping': {g: g for g in genes},
            'no_encontrados': [],
            'advertencias': ['No se pudo validar genes (sin conexión)']
        }
    
    # Procesar resultados
    genes_validos = []
    mapping = {}
    no_encontrados = []
    advertencias = []
    genes_a_reintentar = []  # Para genes no encontrados que intentaremos con MT-
    
    for gene_input, result in zip(genes, results['out']):
        if 'notfound' in result and result['notfound']:
            # Gen no encontrado - guardar para reintentar con MT-
            genes_a_reintentar.append(gene_input)
            
        elif 'symbol' in result:
            # Gen encontrado - usar símbolo oficial
            simbolo = result['symbol']
            genes_validos.append(simbolo)
            mapping[gene_input] = simbolo
        else:
            # Resultado ambiguo o incompleto
            genes_a_reintentar.append(gene_input)
    
    # Reintentar genes no encontrados con prefijo MT- (genes mitocondriales)
    if genes_a_reintentar:
        print(f"[INFO] Buscando variantes mitocondriales (MT-) para {len(genes_a_reintentar)} genes...")
        
        # Crear variantes con MT-
        variantes_mt = [f"MT-{gene}" for gene in genes_a_reintentar]
        
        try:
            results_mt = mg.querymany(
                variantes_mt,
                scopes='symbol,entrezgene,ensembl.gene,uniprot.Swiss-Prot',
                fields='symbol,entrezgene,name',
                species=species,
                returnall=True
            )
            
            for gene_original, gene_mt, result_mt in zip(genes_a_reintentar, variantes_mt, results_mt['out']):
                if 'notfound' not in result_mt and 'symbol' in result_mt:
                    # Encontrado con prefijo MT-
                    simbolo = result_mt['symbol']
                    genes_validos.append(simbolo)
                    mapping[gene_original] = simbolo
                else:
                    # Definitivamente no encontrado
                    no_encontrados.append(gene_original)
                    advertencias.append(f"Gen no encontrado: '{gene_original}'")
                    
        except Exception as e:
            # Si falla el reintento, agregar todos a no encontrados
            for gene_original in genes_a_reintentar:
                no_encontrados.append(gene_original)
                advertencias.append(f"Gen no encontrado: '{gene_original}'")
    
    # Resumen
    print(f"Validación completada: {len(genes_validos)} genes válidos")
    if no_encontrados:
        print(f"Advertencia: {len(no_encontrados)} genes no encontrados: {', '.join(no_encontrados)}")
    
    # Si no hay genes válidos, detener
    if not genes_validos:
        print("[ERROR] No se encontró ningún gen válido. Verifica los nombres/IDs.")
        sys.exit(1)
    
    print(f"{'='*70}\n")
    
    return {
        'validos': genes_validos,
        'mapping': mapping,
        'no_encontrados': no_encontrados,
        'advertencias': advertencias
    }


def analizar_genes(genes, organismo='hsapiens'):
    """
    Realiza análisis de enriquecimiento funcional usando g:Profiler.
    
    g:Profiler consulta múltiples bases de datos para identificar procesos
    biológicos y vías metabólicas enriquecidas.
    
    Bases de datos consultadas:
        - GO:BP: Gene Ontology - Biological Process
        - KEGG: Kyoto Encyclopedia of Genes and Genomes (vías metabólicas)
        - REAC: Reactome Pathway Database (vías de reacción biológicas)
    
    Args:
        genes: Lista de símbolos de genes validados
        organismo: Código del organismo (default: 'hsapiens')
        
    Returns:
        DataFrame con resultados del análisis
    """
    print(f"\n{'='*70}")
    print("ANÁLISIS DE ENRIQUECIMIENTO FUNCIONAL")
    print(f"{'='*70}")
    print(f"[INFO] Analizando {len(genes)} genes validados")
    print("[INFO] Bases de datos: GO:BP, KEGG, REAC")
    
    # Inicializar cliente de g:Profiler
    gp = GProfiler(return_dataframe=True)
    
    # Bases de datos a consultar
    # GO:BP - Gene Ontology Biological Process
    # KEGG - Kyoto Encyclopedia of Genes and Genomes (vías metabólicas)
    # REAC - Reactome Pathway Database (vías de reacción)
    sources = ['GO:BP', 'KEGG', 'REAC']
    
    try:
        # Realizar análisis de enriquecimiento
        # - user_threshold: nivel de significancia (p-valor ajustado < 0.05)
        # - significance_threshold_method: corrección FDR (False Discovery Rate)
        resultados = gp.profile(
            organism=organismo,
            query=genes,
            sources=sources,
            user_threshold=0.05,
            significance_threshold_method='fdr'
        )
        
        if resultados is not None and not resultados.empty:
            print(f"Análisis completado: {len(resultados)} términos enriquecidos encontrados")
        else:
            print("No se encontraron términos enriquecidos significativos (p < 0.05)")
            resultados = pd.DataFrame()
            
    except Exception as e:
        print(f"[ERROR] Error en el análisis: {str(e)}")
        resultados = pd.DataFrame()
    
    print(f"{'='*70}\n")
    return resultados


def exportar_resultados(df, validacion_info, output_dir='results'):
    """
    Exporta los resultados a archivos separados por base de datos para análisis en R.
    
    Estructura de salida:
        results/
        ├── validacion_genes.csv         # Info de validación
        ├── resultados_completos.csv     # Todos los resultados
        ├── resultados_completos.xlsx    # Todos los resultados (Excel)
        ├── GO_BP_resultados.csv         # Solo GO:BP
        ├── KEGG_resultados.csv          # Solo KEGG
        └── REAC_resultados.csv          # Solo Reactome
    
    Args:
        df: DataFrame con resultados del análisis
        validacion_info: Información de validación de genes
        output_dir: Directorio de salida
    """
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"\n{'='*70}")
    print("EXPORTACIÓN DE RESULTADOS")
    print(f"{'='*70}")
    
    # 1. Exportar información de validación en formato CSV
    val_data = {
        'gen_input': [],
        'gen_validado': [],
        'estado': [],
        'notas': []
    }
    
    for orig, val in validacion_info['mapping'].items():
        val_data['gen_input'].append(orig)
        val_data['gen_validado'].append(val)
        val_data['estado'].append('valido')
        val_data['notas'].append(f'Convertido de {orig}' if orig != val else 'Validado')
    
    for gene in validacion_info['no_encontrados']:
        val_data['gen_input'].append(gene)
        val_data['gen_validado'].append('NA')
        val_data['estado'].append('no_encontrado')
        val_data['notas'].append('Gen no encontrado en MyGene.info')
    
    df_validacion = pd.DataFrame(val_data)
    val_file = os.path.join(output_dir, 'validacion_genes.csv')
    df_validacion.to_csv(val_file, index=False)
    print("Validación: validacion_genes.csv")
    
    # 2. Exportar resultados del análisis funcional
    if df.empty:
        print("No hay resultados de enriquecimiento para exportar")
        print(f"{'='*70}\n")
        return
    
    # CSV completo
    csv_file = os.path.join(output_dir, 'resultados_completos.csv')
    df.to_csv(csv_file, index=False)
    print(f"CSV completo: resultados_completos.csv ({len(df)} términos)")
    
    # Excel completo
    excel_file = os.path.join(output_dir, 'resultados_completos.xlsx')
    df.to_excel(excel_file, index=False, engine='openpyxl')
    print("Excel completo: resultados_completos.xlsx")
    
    # 3. Archivos separados por base de datos
    bases_datos = {
        'GO:BP': 'GO_BP_resultados.csv',
        'KEGG': 'KEGG_resultados.csv',
        'REAC': 'REAC_resultados.csv'
    }
    
    for source, filename in bases_datos.items():
        df_source = df[df['source'] == source]
        if not df_source.empty:
            source_file = os.path.join(output_dir, filename)
            df_source.to_csv(source_file, index=False)
            print(f"{source}: {filename} ({len(df_source)} términos)")
        else:
            print(f"{source}: Sin resultados")
    
    print(f"\nTodos los archivos exportados a: {output_dir}/")
    print(f"{'='*70}\n")



def main():
    """Función principal con interfaz de línea de comandos."""
    
    print("\n" + "#"*70)
    print("ANÁLISIS FUNCIONAL DE GENES")
    print("#"*70)
    
    # Configurar argumentos de línea de comandos
    parser = argparse.ArgumentParser(
        description='Análisis Funcional de Genes con Validación de IDs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:
    # Con archivo de entrada
    python analisis_funcional.py -i data/genes_input.txt
    
    # Con genes por línea de comandos (símbolos)
    python analisis_funcional.py -g COX4I2 ND1 ATP6
    
    # Con Entrez Gene IDs
    python analisis_funcional.py -g 1327 4535 4508
    
    # Con organismo diferente
    python analisis_funcional.py -i genes.txt -org mmusculus
        """
    )
    
    # Entrada de genes (archivo o lista)
    grupo = parser.add_mutually_exclusive_group(required=True)
    grupo.add_argument('-i', '--input', help='Archivo con genes (uno por línea)')
    grupo.add_argument('-g', '--genes', nargs='+', help='Lista de genes o IDs')
    
    # Opciones adicionales
    parser.add_argument('-o', '--output', default='results', 
                       help='Directorio de salida (default: results)')
    parser.add_argument('-org', '--organism', default='hsapiens',
                       help='Organismo para g:Profiler (default: hsapiens)')
    
    args = parser.parse_args()
    
    # Obtener lista de genes
    if args.input:
        genes_input = leer_genes(args.input)
    else:
        genes_input = args.genes
    
    # Validar y convertir genes
    validacion = validar_y_convertir_genes(genes_input, args.organism)
    genes_validos = validacion['validos']
    
    # Realizar análisis funcional
    resultados = analizar_genes(genes_validos, args.organism)
    
    # Exportar resultados
    exportar_resultados(resultados, validacion, args.output)
    
    print("\n" + "#"*70)
    print("ANÁLISIS COMPLETADO EXITOSAMENTE")
    print("#"*70 + "\n")


if __name__ == '__main__':
    main()