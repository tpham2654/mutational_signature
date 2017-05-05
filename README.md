# Mutational Signature Library

Library for computing enrichment of mutational signature (like APOBEC signature) :alien:

# Usage

(Awesome & easy-to-follow documentation is coming - very soon)

## Input

- Mutations: VCF or MAF
- Reference genome: FASTA or FA

## Output

- Results (n_features-x-n_samples): TXT (tab-separated)
- Ranking & distribution of the enrichment): PDF

# Comparison with Firehose's P-MACD

[Comparison results](https://github.com/KwatME/mutational-signature/blob/master/media/compare_with_p-macd.pdf)

# References

Algorithm was implemented based on:

- <http://firebrowse.org/> (P-MACD tool)
- <http://www.nature.com/nature/journal/v500/n7463/full/nature12477.html>
- <http://www.nature.com/ng/journal/v45/n9/abs/ng.2702.html>
