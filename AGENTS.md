## Cursor Cloud specific instructions

This is a collection of standalone bioinformatics CLI scripts (Python 3 + Perl 5). There are no services, databases, or build steps.

### Dependencies

- **Python 3.x** with: `numpy`, `pandas`, `scikit-learn`, `scipy`, `biopython`, `matplotlib`, `seaborn`
- **Perl 5.x** (standard library only, no CPAN modules needed)

Install Python deps: `pip install numpy pandas scikit-learn scipy biopython matplotlib seaborn`

### Running scripts

All scripts are in the repo root. See `README.md` for usage of each script. Key notes:

- `calculate_position_importance.py` requires a FASTA alignment file (`test.aln.fas`) and a tab-separated labels file (`labels.txt`) in the working directory. Use `matplotlib.use('Agg')` before importing pyplot when running headless (no display server).
- `convert_position.py` reads from stdin for position input; pipe values or use `echo "1,2,3" | python3 convert_position.py <fasta> <gene>`.
- `miniprot_GFF_2_EVM_GFF3.py` writes to stdout; redirect with `> output.gff3`.
- Perl scripts expect specific input files (GFF, FASTA, protein files) in the current working directory; see each script's header for expected filenames.

### Testing

There are no automated tests or linting configured in this repository. Verification is done by running individual scripts against appropriate input data.
