# This the motif-mark assignment for BGMP@UO

## Objective: visualize motifs across genes or gene-like sequences

### Input
- a fasta file (-f, --file; required)
- a motifs text file (-m, --motifs; required)
- an output name (-o, --out; optional)

### Output
- one .svg image with:
  - labels
  - one or more 'genes'
    - accommodates multiple, potentially ambiguous motifs

### Requires
- Python 3
- PyCairo
