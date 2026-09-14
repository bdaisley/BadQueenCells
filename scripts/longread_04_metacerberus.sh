#!/bin/bash

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 1: Setup directories
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

working_dir=~/CBQC_project/longread
medaka_dir=${working_dir}/04_medaka_consensus
prodigal_dir=${working_dir}/05_prodigal
metacerberus_dir=${working_dir}/06_metacerberus

mkdir -p $prodigal_dir $metacerberus_dir $metacerberus_dir/fna_symlinks

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 2: Make fna
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
cd $medaka_dir
for i in ${medaka_dir}/*/*.fna; do
  fna_base=$(basename "$i" | awk -F'\\.fna' '{print $1}') ;
  echo $fna_base
  ln -s ${medaka_dir}/${fna_base}.fna $metacerberus_dir/fna_symlinks/${fna_base}.fna
done

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 3: Run metacerberus
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
conda activate metacerberus

#--------------------FIX FOR multiple metacerbs being run
export PYTHONPATH="$HOME/tmp/hydrampp005"
python scripts/longread_04_metacerberus_hydra.py
#-----------------------------------------------
hydra_port=$(
python - <<'PY'
import socket

with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as sock:
    sock.bind(("127.0.0.1", 0))
    print(sock.getsockname()[1])
PY
)

echo "Using HydraMPP port: $hydra_port"

metacerberus.py --prodigal $metacerberus_dir/fna_symlinks \
  --hmm ALL \
  --cpus 24  \
  --address local \
  --port "$hydra_port" \
  --dir_out ${metacerberus_dir}
