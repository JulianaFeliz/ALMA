#!/usr/bin/env python3
import glob
import re
import os

def run():
    # =========================================================
    # 1. Carregar os dados de abundância do Mothur (.shared)
    # =========================================================
    shared_data = {}
    for sf in glob.glob("results/**/*.shared", recursive=True):
        sample = os.path.basename(sf).split('.')[0]
        with open(sf, 'r') as f:
            lines = f.readlines()
            if len(lines) < 2:
                continue
            headers = lines[0].strip().split('\t')
            values  = lines[1].strip().split('\t')
            clean_headers = [re.sub(r'Otu0+', 'Otu', h) for h in headers[3:]]
            shared_data[sample] = dict(zip(clean_headers, [int(x) for x in values[3:]]))

    # =========================================================
    # 2. Carregar o melhor hit taxonômico por OTU (.per_query.tsv)
    # =========================================================
    otu_best_hit = {}
    for qf in glob.glob("results/**/*.per_query.tsv", recursive=True):
        with open(qf, 'r') as f:
            lines = f.readlines()
        if len(lines) < 2:
            continue

        for line in lines[1:]:
            if not line.strip():
                continue
            cols = line.strip().split('\t')
            if len(cols) < 3:
                continue

            q_name = cols[0]
            lwr    = float(cols[1])

            # Extrai o OTU ID diretamente do nome da sequência
            # O nome vem do MAFFT como: UFRN_D7_..._Otu14_5...
            match = re.search(r'(Otu\d+)', q_name)
            if not match:
                continue

            otu_id = re.sub(r'Otu0+', 'Otu', match.group(1))

           # Guarda a taxonomia limpando os espaços
            taxopath = cols[-1]
            parts = [p.strip() for p in taxopath.split(';') 
                     if p.strip() and p.strip().lower() != 'fungi']
            
            if parts:
                current_depth = len(parts) # Mede quão profundo chegou (ex: Família=4, Espécie=6)
                
                # Pega os dados que já salvamos dessa OTU (ou zera se for a primeira vez)
                stored_depth = otu_best_hit[otu_id]['depth'] if otu_id in otu_best_hit else 0
                stored_lwr = otu_best_hit[otu_id]['lwr'] if otu_id in otu_best_hit else -1

                # LÓGICA "TIP-FIRST":
                # Atualiza SE chegou mais fundo na árvore (maior depth)
                # OU SE empatou na profundidade, mas tem uma certeza maior (LWR >)
                if current_depth > stored_depth or (current_depth == stored_depth and lwr > stored_lwr):
                    otu_best_hit[otu_id] = {
                        'lwr'      : lwr,
                        'depth'    : current_depth, # Precisamos salvar a profundidade para comparar!
                        'taxon'    : parts[-1],
                        'genus'    : parts[-2] if len(parts) > 1 else 'Unknown',
                        'full_path': taxopath
                    }

    # =========================================================
    # 3. Construir nomes finais para cada OTU
    # =========================================================
    final_otu_names    = {}

    for otu_id, info in otu_best_hit.items():
        last_rank = info['taxon']
        # Heurística: se começa com minúscula ou tem espaço → é epíteto específico
        is_species = last_rank[0].islower() or ' ' in last_rank

        if is_species:
            final_otu_names[otu_id] = f"{info['genus']}_{last_rank.replace(' ', '_')}"
        else:
            # Se não é espécie, usa o Nível + OTU ID (Ex: Dominikiaceae_Otu14)
            final_otu_names[otu_id] = f"{last_rank}_{otu_id}"

    # =========================================================
    # O resto do script (4. Escrever a tabela) continua EXATAMENTE igual!
    # =========================================================
    # 4. Escrever a tabela final
    # =========================================================
    all_names = sorted(set(final_otu_names.values()))

    # Diagnóstico rápido
    matched   = sum(1 for otu in shared_data.get(next(iter(shared_data), ''), {})
                    if otu in final_otu_names)
    unmatched = sum(1 for otu in shared_data.get(next(iter(shared_data), ''), {})
                    if otu not in final_otu_names)
    print(f"[INFO] OTUs com hit taxonômico : {len(final_otu_names)}")
    print(f"[INFO] Táxons únicos na tabela : {len(all_names)}")
    print(f"[INFO] OTUs sem hit (ignorados) : {unmatched}")

    with open('ALMA_OTU_abundance.tsv', 'w') as out:
        out.write('\t'.join(['label', 'Group', 'numOtus'] + all_names) + '\n')
        for sample in sorted(shared_data.keys()):
            sums = {name: 0 for name in all_names}
            for otu_id, count in shared_data[sample].items():
                if otu_id in final_otu_names:
                    sums[final_otu_names[otu_id]] += count
            out.write('\t'.join(
                ['0.02', sample, str(len(all_names))] +
                [str(sums[n]) for n in all_names]
            ) + '\n')

    print("[INFO] Tabela salva em ALMA_OTU_abundance.tsv")

if __name__ == '__main__':
    run()
