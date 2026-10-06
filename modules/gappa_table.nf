process GAPPA_TABLE {
    tag "Merge_Abundance_Taxonomy"
    label 'process_low'

    // Usamos um container básico com Python 3 instalado
    container 'quay.io/biocontainers/python:3.9--1'
    
    publishDir "${params.outdir}/final_tables", mode: 'copy'

    input:
    path per_query_tsv // O output global.per_query.tsv do GAPPA_ASSIGN
    path shared_file   // O output final_otu_table.shared do VSEARCH_MAP_OTUS

    output:
    path "Final_OTU_Table_with_Taxonomy.csv", emit: final_table

    script:
    """
    python3 -c "
import csv
from collections import defaultdict

# 1. Lê a taxonomia gerada pelo GAPPA
# per_query.tsv tem colunas como: name, taxpath, taxpath_sn, etc.
taxonomy_dict = {}
with open('${per_query_tsv}', 'r') as f_tax:
    reader = csv.DictReader(f_tax, delimiter='\\t')
    for row in reader:
        # A coluna 'name' tem o nome da OTU. A coluna 'taxpath' tem a linhagem.
        otu_name = row['name'].strip()
        taxonomy_dict[otu_name] = row.get('taxpath', 'Unclassified')

# 2. Lê a tabela .shared do VSEARCH e a transpõe (linhas viram colunas)
# Formato .shared: label, Group(Sample), numOtus, OTU1, OTU2...
samples = []
abundances = defaultdict(dict) # abund[otu][sample] = count

with open('${shared_file}', 'r') as f_shared:
    reader = csv.reader(f_shared, delimiter='\\t')
    header = next(reader)
    otu_list = header[3:] # A partir da 4ª coluna são as OTUs
    
    for row in reader:
        sample_name = row[1]
        samples.append(sample_name)
        counts = row[3:]
        for otu, count in zip(otu_list, counts):
            abundances[otu][sample_name] = count

# 3. Escreve a Tabela Final Mestra (OTUs nas linhas, Amostras nas colunas + Taxonomia)
with open('Final_OTU_Table_with_Taxonomy.csv', 'w') as f_out:
    # Cabeçalho
    header_out = ['OTU_ID'] + samples + ['Taxonomy']
    f_out.write(','.join(header_out) + '\\n')
    
    # Preenche os dados
    for otu in otu_list:
        row_data = [otu]
        for sample in samples:
            row_data.append(abundances[otu].get(sample, '0'))
        
        # Puxa a taxonomia (se não achou no GAPPA, coloca Unclassified)
        tax = taxonomy_dict.get(otu, 'Unclassified')
        row_data.append(tax)
        
        f_out.write(','.join(row_data) + '\\n')
"
    """
}