Lux
===================

Overview
-------------
This version supports BS-seq, TAB-seq and fCAB-seq data with fixed experimental parameters. For more details, please see the main README document.

Features
-------------
Please see the main README document.

Installation
-------------
Please see the main README document.

### Generating input files
A python script **parse.py** is supplied for generating input files for the use with **Lux**. 
Basically, it transforms the user-supplied data files into the data format read by CmdStan. In addition, it initializes the model parameters with values sampled from the corresponding priors

    usage: parse.py [-h] -d DATA -p PRIOR -b BSEFF [BSEFF ...] -i BSBEFF [BSBEFF ...] -o OXEFF [OXEFF ...] -l LABEFF [LABEFF ...] -pro PROEFF [PROEFF ...] -s SEQERR [SEQERR ...] -pr PREFIX [-v]

    Generates data and init files in the dump format for Lux

    optional arguments:
      -h, --help                                             show this help message and exit
      -d DATA, --data DATA                                   noncontrol cytosine data
      -p PRIOR, --prior PRIOR                                prior of the noncontrol cytosines
      -b BSEFF [BSEFF ...], --bseff BSEFF [BSEFF ...]        bisulfite conversion efficiencies for each replicate
      -i BSBEFF [BSBEFF ...], --bsbeff BSBEFF [BSBEFF ...]   inaccurate bisulfite conversion efficiencies for each replicate
      -o OXEFF [OXEFF ...], --oxeff OXEFF [OXEFF ...]        oxidation efficiencies for each replicate
      -l LABEFF [LABEFF ...], --labeff LABEFF [LABEFF ...]   labeling efficiencies for each replicate
      -pro PROEFF [PROEFF ...], --proeff PROEFF [PROEFF ...] protection efficiencies for each replicate
      -s SEQERR [SEQERR ...], --seqerr SEQERR [SEQERR ...]   sequencies errors for each replicate
      -pr PREFIX, --prefix PREFIX                            prefix of the output files
      -v, --version                                          show program's version number and exit

For instance, the script **parse.py** is called in the case of one sample as 

    ./parse.py --data data.tsv --prior prior.tsv -b 0.99 -i 0.001 -o 0.8 -l 0.9 -pro 0.9 -s 0.001 --prefix sample

This command will generate the two files, **sample_data.R** and **sample_init.R**, which are ready to be used in **Lux**.

As you notice **parse.py** takes two files as input, so next, we will go through their content and format in detail.

#### Noncontrol cytosines
The files **data.tsv** and **prior.tsv** have the count data and prior information on the noncontrol cytosines of interest, respectively.

Each of the cytosines has its own line in the files, and thus the files are required to have the same number of lines, moreover, the order of the cytosines is assumed to be the same in the files. 

Each line in the file **data.tsv** is composed of one or more replicate-specific tab-separated blocks of six tab-separated nonnegative integers
>N<sub>BS</sub><sup>C</sup>\tN<sub>BS</sub>\tN<sub>TAB</sub><sup>C</sup>\tN<sub>TAB</sub>\tN<sub>fCAB</sub><sup>C</sup>\tN<sub>fCAB</sub>

That is, the number of Cs and the total BS-seq read-outs are listed first, which are followed by the number of Cs and total TAB-seq read-outs and the number of Cs and total fCAB-seq read-outs.

Moreover, on each line there should be exactly 6×N<sub>replicates</sub> values separated with the tabs.

The replicate-level in the hierarchical model formulation explicitly assumes that the replicates share the same prior. Thus, the prior of μ is defined by giving the parameter α for all the noncontrol cytosines in terms of four pseudo-counts α<sub>1</sub>, α<sub>2</sub>, α<sub>3</sub> and α<sub>4</sub> corresponding to p(C), p(5mC), p(5hmC) and p(5fC), respectively.

That is, all the lines in the file **prior.tsv** are required to have the following format
>α<sub>1</sub>\tα<sub>2</sub>\tα<sub>3</sub>\tα<sub>4</sub>

### References
[1] M. J. Betancourt, “Generalizing the No-U-Turn Sampler to Riemannian Manifolds,” arXiv, vol. 1304, no. 1920, 2013. 

[2] K. D. Hansen, B. Langmead and R. A. Irizarry, “BSmooth: from whole genome bisulfite sequencing reads to differentially methylated regions.,” Genome Biol, vol. 13, no. 10, pp. R83, Oct 2012. 

[3] M. D. Hoffman and A. Gelman, “The No-U-Turn Sampler: Adaptively Setting Path Lengths in Hamiltonian Monte Carlo,” Journal of Machine Learning Research, vol. in press, 2013. 

[4] F. Krueger and S. R. Andrews, “Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications.,” Bioinformatics, vol. 27, no. 11, pp. 1571-1572, Jun 2011. 

[5] Stan Development Team, “Stan Modeling Language Users Guide and Reference Manual, Version 2.2,” 2014. 

[6] Stan Development Team, “Stan: A C++ Library for Probability and Sampling, Version 2.2,” 2014. 

[7] D. Sun, Y. Xi, B. Rodriguez, H. J. Park, P. Tong, M. Meong, M. A. Goodell and W. Li, “MOABS: model based analysis of bisulfite sequencing data.,” Genome Biol, vol. 15, no. 2, pp. R38, Feb 2014. 
