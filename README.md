# BeamTrellis

- Project forked from [nanopore_dna_storage](https://github.com/shubhamchandak94/nanopore_dna_storage)
- The `extension-to-other-datasets` branch of this Git repo extends the `bonito` branch of the original project by allowing the use of BeamTrellis [1, 2] over other FAST5 datasets - those that do not necessarily relate to DNA strands corresponding to convolutional codeswords.
- Based on "SynDe: Syndrome-guided Decoding of Raw Nanopore Reads" (link here)
- Keywords: DNA data storage, nanopore sequencing, convolutional codes, basecaller-decoder

## Prerequisites

## Download and installation

It is advisable to run the upcoming commands in a virtual environment:
```bash
python3 -m venv venv3
source venv3/bin/activate
```

Following this, run the following to download and install:
```bash
git clone --recursive -b extension-to-other-datasets https://github.com/anisha-ban/BeamTrellis/
cd BeamTrellis/
./install.sh
```

## Usage of scripts

- Old scripts:
    - The instructions for these scripts can be found in the old README, named `BeamTrellis-init-readme.md`.
- New scripts:
    - `synde_generate_sim_targets.py`: generates simulation targets (JSON files with each entry containing the ID of a raw read, the leading primer, starting position in raw read, etc.) for Chandak et al.'s raw data (HDF5 files)
    - `generate_decoded_lists_for_chandak.py`: (from Chandak et al.'s original repo)  performs decoding on each raw read in the given HDF5 file by Chandak et al.'s procedure (see [1,2])
    - `compute_error_rate_from_decoded_lists_for_stanford_dataset.py`: (from Chandak et al.'s original repo) Used to compute the number of errors (due to no CRC match being found and due to incorrect CRC match) from the decoded lists.
    - `compute_error_rate_from_decoded_lists_for_other_datasets.py`: adapts `compute_error_rate_from_decoded_lists_for_stanford_dataset.py` to work on other datasets.
    - `generate_decoded_lists_other_datasets.py`: adapts `generate_decoded_lists_for_chandak.py` to work on other datasets. The modifications applied to accomplish this are explained in the next section.
    - `add_stanford_primer_positions.py` (is this used or integrated into synde_generate_sim_targets?)
    - `my_utils.py`: contains helper functions for the aforementioned new scripts

Usage of these scripts can be found in xxxxxxxxxxx

## How new FAST5 datasets are accommodated

To simulate our decoding algorithm Synde on other FAST5 datasets, say for an error correction code $C \in `{0,1,2,3}`^m$, we adopt the following approach.
(We assume the presence of primer sequences, each of 25 nucleotides, on either side of the payload)

1.  A given raw read is first basecalled, and the resulting sequence is matched to the causal reference sequence by utilizing alignment tools such as SAMtools [4, 5]. If the read is reverse-completeded, or has a secondary or a chimeric alignment, it is ignored (see [4] for details).
2.   Following this, the part of the reference that is expressed in the raw read, say $y \in `{A,C,G,T}`^z$ is extracted. If $z<2 \cdot 25 + m$, then the process is aborted. We do so since such a read cannot accommodate two primer sequences and a DNA codeword of $130$ symbols.
3.  Next, we choose an index $i \in `{0,1,...,z-50-m}`$ uniformly at random, and set the two primer sequences as $v^{(1)}=y_{i+1}^{i+25}$ and $v^{(2)}=y_{i+m+26}^{i+m+50}$, while the payload is assumed to be $x=y_{i+26}^{i+m+25}$.
4.  A random codeword $u \in C$ is drawn from the chosen error-correction scheme, and the offset with respect to the true payload, say $o$, is computed as $o = f^{-1}(x) - u \pmod{4}$.
5.  Our decoding pipeline is then provided the two primers $v^{(1)}, v^{(2)}$ and the computed offset $o$: PrimerSeeker searches for the leading primer $v^{(1)}$ and the estimate of the starting position of $v^{(1)}$ in the raw read that PrimerSeeker produces, is leveraged by \syndec~to initiate decoding on the corresponding syndrome trellis, while incorporating the offset $o$ appropriately.
6.  If the decoder fails to produce a codeword sequence of the requisite length $m$, an erasure is declared. Otherwise, if the decoded codeword, say $x' \in C$, does exactly matches the initial codeword $u$, the decoder is said to have worked correctly. Else, a frame error is declared.



# References
[1] S. Chandak et al., “Overcoming High Nanopore Basecaller Error Rates for DNA Storage via Basecaller-Decoder Integration and Convolutional Codes,” in ICASSP 2020 - 2020 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), May 2020, pp. 8822–8826. doi: 10.1109/ICASSP40776.2020.9053441.

[2] B. Lau et al., “Magnetic DNA random access memory with nanopore readouts and exponentially-scaled combinatorial addressing,” Sci Rep, vol. 13, no. 8514, Art. no. 1, May 2023, doi: 10.1038/s41598-023-29575-z.

[3] A. Banerjee et al., ............

[4] H. Li et al., “The Sequence Alignment/Map format and SAMtools,” Bioinformatics, vol. 25, no. 16, pp. 2078–2079, Aug. 2009, doi: 10.1093/bioinformatics/btp352.

[5] P. Danecek et al., “Twelve years of SAMtools and BCFtools,” Gigascience, vol. 10, no. 2, p. giab008, Feb. 2021, doi: 10.1093/gigascience/giab008.
