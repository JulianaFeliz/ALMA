#!/usr/bin/env python3
import glob, re, os
from collections import defaultdict

def get_taxon_name(taxopath):
    parts = [p.strip() for p in taxopath.split(';') 
             if p.strip() and p.strip().lower() not in ('fungi', 'taxopath')]
    if not parts:
        return None
    last = parts[-1]
    is_species = last[0].islower() or ' ' in last
    if is_species:
        return f"{parts[-2]}_{last.replace(' ', '_')}" if len(parts) > 1 else last
    else:
        return f"{last}_sp"

def run():
    # Processa por amostra, mantendo OTU IDs separados
    sample_taxa_counts = {}  # {sample: {taxon: count}}
    all_taxa = set()

    # Descobre amostras pelos shared files
    for sf in glob.glob("*.shared"):
        sample = os.path.basename(sf).split('.')[0]
        
        # Lê shared
        with open(sf) as f:
            lines = f.readlines()
        if len(lines) < 2:
            continue
        headers = lines[0].strip().split('\t')
        otu_counts = {}
        for line in lines[1:]:
            vals = line.strip().split('\t')
            if len(vals) < 4:
                continue
            for h, v in zip(headers[3:], vals[3:]):
                otu_id = re.sub(r'Otu0+', 'Otu', h)
                otu_counts[otu_id] = otu_counts.get(otu_id, 0) + int(v)

        # Lê per_query da mesma amostra
        qf = f"{sample}.per_query.tsv"
        if not os.path.exists(qf):
            continue
        
        otu_to_taxon = {}
        with open(qf) as f:
            for line in f.readlines()[1:]:
                if not line.strip():
                    continue
                cols = line.strip().split('\t')
                if len(cols) < 6:
                    continue
                q_name, lwr = cols[0], float(cols[1])
                match = re.search(r'(Otu\d+)', q_name)
                if not match:
                    continue
                otu_id = re.sub(r'Otu0+', 'Otu', match.group(1))
                taxon = get_taxon_name(cols[-1])
                if not taxon:
                    continue
                # Melhor LWR por OTU dentro da amostra
                if otu_id not in otu_to_taxon or lwr > otu_to_taxon[otu_id][0]:
                    otu_to_taxon[otu_id] = (lwr, taxon)

        # Soma reads por táxon para essa amostra
        taxa_counts = defaultdict(int)
        for otu_id, count in otu_counts.items():
            if otu_id in otu_to_taxon:
                taxon = otu_to_taxon[otu_id][1]
                taxa_counts[taxon] += count

        sample_taxa_counts[sample] = taxa_counts
        all_taxa.update(taxa_counts.keys())

    all_taxa = sorted(all_taxa)
    print(f"[INFO] Amostras processadas: {len(sample_taxa_counts)}")
    print(f"[INFO] Táxons únicos: {len(all_taxa)}")

    with open("ALMA_OTU_abundance.tsv", "w") as out:
        out.write('\t'.join(['label', 'Group', 'numOtus'] + all_taxa) + '\n')
        for sample in sorted(sample_taxa_counts.keys()):
            counts = sample_taxa_counts[sample]
            out.write('\t'.join(
                ['0.02', sample, str(len(all_taxa))] +
                [str(counts.get(t, 0)) for t in all_taxa]
            ) + '\n')

    print("[INFO] Tabela salva em ALMA_OTU_abundance.tsv")

if __name__ == '__main__':
    run()
