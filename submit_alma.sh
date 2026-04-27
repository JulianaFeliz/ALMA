#!/bin/bash
#SBATCH --job-name=ALMA_Analysis        # Nome do trabalho no cluster
#SBATCH --output=ALMAout_%j.log        # Arquivo de log (saída)
#SBATCH --error=ALMAerr_%j.log         # Arquivo de erro
#SBATCH --nodes=1                       # Usar 1 nó do supercomputador
#SBATCH --cpus-per-task=8              # Quantidade de CPUs (Nextflow vai dividir isso)
#SBATCH --mem=32G                       # Memória RAM total (ajuste se necessário)
#SBATCH --time=19:00:00                 # Tempo máximo de execução (24 horas)
#SBATCH --partition=intel-128               # Nome da fila (confirme se é 'batch' no seu cluster)

# --- Notificações por E-mail ---
#SBATCH --mail-user=dearfelixx@gmail.com
#SBATCH --mail-type=ALL               # Receber e-mail para todos os eventos (BEGIN, END, FAIL)

# --- 1. Carregar os módulos necessários ---
# No Ubuntu/WSL você não usa isso, mas no Supercomputador é essencial:
module load compilers/java/23.0.1
module load singularity/3.7.1  # O Singularity é vital para o ALMA rodar os containers
export PATH=$HOME:$PATH

#---- 2. Caminhos (Ajuste aqui)
# Definimos onde está o executável e onde está o pipeline
NEXTFLOW_BIN=$HOME/nextflow
PIPELINE_DIR=$HOME/ALMA


# --- 3. Comando de Execução do ALMA ---

$NEXTFLOW_BIN run main.nf \
    --input_dir data/ \
    --outdir results/ \
    -profile singularity \
    -resume

# --- 4. Dica de Pós-Processamento ---
echo "Pipeline completed. Your results are in the results/ folder :)"

