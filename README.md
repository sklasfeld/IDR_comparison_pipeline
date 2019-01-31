# IDR_comparison_pipeline

This script generates self-pseudo replicates and pooled-pseudo replicates for a ChIP-seq
experiment. It then run MACS2 (Zhang et al 2008) to call ChIP peaks on the true replicates,
pooled replicates, self-pseudo replicates, and pooled-pseudo replicates.
Next, IDR (Li et al 2011) is run on all pairs of true replicates, self-pseudoreplicates,
ad pooled pseudo replicates. The results file gives the number of peaks
that pass IDR for each pairwise comparison. Using the maximum number of
peaks that pass IDR given two true replicates, a conservative set of peaks
is output from the pooled MACS2 results.

This pipeline is based on the ENCODE paper (Landt et al 2012)

Dependencies:
* python 2.7
* python libraries:
 * pandas (0.19.2)
 * numpy (1.12.0)
* bedtools
* IDR (https://github.com/nboley/idr)
* MACS2 (https://github.com/taoliu/MACS)
* samtools

Citations:
* Landt, S. G., Marinov, G. K., Kundaje, A., Kheradpour, P., Pauli, F., Batzoglou, S., ... & Chen, Y. (2012). ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia. Genome research, 22(9), 1813-1831.
* Li, Q., Brown, J. B., Huang, H., & Bickel, P. J. (2011). Measuring reproducibility of high-throughput experiments. The annals of applied statistics, 5(3), 1752-1779.
* Zhang, Y., Liu, T., Meyer, C. A., Eeckhoute, J., Johnson, D. S., Bernstein, B. E., ... & Liu, X. S. (2008). Model-based analysis of ChIP-Seq (MACS). Genome biology, 9(9), R137.
