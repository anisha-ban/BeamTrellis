#!/bin/bash

pip install --upgrade pip

# install modified fast-ctc-decode
cd ~/ont-bonito/fast-ctc-decode # UPDATE THIS!!!!!!!!!!!!!!!!!!
pip install "maturin>=0.14,<0.15"
python -m maturin build --release --features python
pip install target/wheels/*.whl --force-reinstall

# install modified bonito (v0.1.2)
cd ~/nanopore_dna_storage_bonito/bonito # UPDATE THIS!!!!!!!!!!!!!!!!!!
pip install --no-build-isolation -r requirements.txt
pip install -e .
cd ..

pip install six

# install fast5_research
cd fast5_research
pip install --no-build-isolation .
cd ..

pip install pycparser cffi

# ========= old installation commands from Chandak et al.'s repo ===========
pip3 install crc8 \
            distance \
            h5py \
            numpy \
            scipy \
            scrappie \
            pysam

pip3 install fast5_research


cd RSCode_schifra
make clean
make
cd ../
g++ viterbi/viterbi_convolutional_code_approach_1.cpp -std=c++11 -o viterbi/viterbi_nanopore.out -Wall  -Iviterbi -fopenmp -O3 -march=native
g++ util/read_length_distribution.cpp -std=c++11 -o util/read_length_distribution.out -Wall -O3 -march=native
