Lux
===================

Overview
-------------
This locus based version supports BS-seq and oxBS-seq data. For more details, please see the main README document.

Features
-------------
Please see the main README document.

Installation
-------------
Please see the main README document.

### Generating input files
A python script **parse.py** is supplied for generating input files for the use with **Lux**. 
Basically, it transforms the user-supplied data files into the data format read by CmdStan. In addition, it initializes the model parameters with values sampled from the corresponding priors

    usage: parse.py [-h] -d DATA -p PRIOR -cd CONTROL_DATA -cp CONTROL_PRIOR -pr PREFIX [-v]
    
    Generates data and init files in the dump format for Lux
    
     optional arguments:
     -h, --help                                         show this help message and exit
     -d DATA, --data DATA                               noncontrol cytosine data
     -p PRIOR, --prior PRIOR                            prior of the noncontrol cytosines
      -cd CONTROL_DATA, --control-data CONTROL_DATA     control cytosine data
      -cp CONTROL_PRIOR, --control-prior CONTROL_PRIOR  priors of the control cytosines
      -pr PREFIX, --prefix PREFIX                       prefix of the output files
     -v, --version                                      show program's version number and exit

For instance, the script **parse.py** is called as

    ./parse.py --data data.tsv --prior prior.tsv --control-data control_data.tsv --control-prior control_prior.tsv --prefix sample

This command will generate the two files, **sample_data.R** and **sample_init.R**, which are ready to be used in **Lux**.

As you notice **parse.py** takes four files as input, so next, we will go through their content and format in detail.

#### Noncontrol cytosines
The files **data.tsv** and **prior.tsv** have the count data and prior information on the noncontrol cytosines of interest, respectively.

This version requires that loci are analyzed separately. Each of the cytosines has its own line in **data.txt**. 

Each line in the file **data.tsv** is composed of one or more replicate-specific tab-separated blocks of four tab-separated nonnegative integers
>N<sub>BS</sub><sup>C</sup>\tN<sub>BS</sub>\tN<sub>oxBS</sub><sup>C</sup>\tN<sub>oxBS</sub>

That is, the number of Cs and the total BS-seq read-outs are listed first, which are followed by the number of Cs and total oxBS-seq read-outs.

Moreover, on each line there should be exactly 4×N<sub>replicates</sub> values separated with the tabs.

The replicate-level in the hierarchical model formulation explicitly assumes that the cytosines in the locus and replicates share the same prior. The prior is defined in terms of three pseudo-counts α<sub>1</sub>, α<sub>2</sub> and α<sub>3</sub> corresponding to p(C), p(5mC) and p(5hmC), respectively.

That is, the file **prior.tsv** is required to have one line in the following format
>α<sub>1</sub>\tα<sub>2</sub>\tα<sub>3</sub>

#### Control cytosines
The files **control_data.tsv** and **control_prior.tsv** have the count data and prior information on the control cytosines of interest, respectively.

Luckily, the data for the control cytosines is supplied in the same format as noncontrol data (see *above*). That is, each line in the file **control_data.tsv** is composed of one or more replicate-specific tab-separated blocks of four tab-separated nonnegative integers.

As in **data.tsv**, on each line there should be exactly 4×N<sub>replicates</sub> values separated with the tabs. Each of the replicate specified in **data.csv** should have its own control data. Moreover, the order of the replicate-specific blocks between the noncontrol and control data is assumed to be the same.

The prior knowledge on the control cytosines is supplied in the file **control_prior.tsv**. Although the hierarchical model allows that the control cytosines would have different priors between replicates, this is not implemented in the current version. Therefore, each line in **control_prior.tsv** should have exactly three tab-separated values and the order of the rows, i.e., control cytosines, should be the same in **control_data.tsv** and **control_prior.tsv**.

### References
[1] M. J. Betancourt, “Generalizing the No-U-Turn Sampler to Riemannian Manifolds,” arXiv, vol. 1304, no. 1920, 2013. 

[2] K. D. Hansen, B. Langmead and R. A. Irizarry, “BSmooth: from whole genome bisulfite sequencing reads to differentially methylated regions.,” Genome Biol, vol. 13, no. 10, pp. R83, Oct 2012. 

[3] M. D. Hoffman and A. Gelman, “The No-U-Turn Sampler: Adaptively Setting Path Lengths in Hamiltonian Monte Carlo,” Journal of Machine Learning Research, vol. in press, 2013. 

[4] F. Krueger and S. R. Andrews, “Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications.,” Bioinformatics, vol. 27, no. 11, pp. 1571-1572, Jun 2011. 

[5] Stan Development Team, “Stan Modeling Language Users Guide and Reference Manual, Version 2.2,” 2014. 

[6] Stan Development Team, “Stan: A C++ Library for Probability and Sampling, Version 2.2,” 2014. 

[7] D. Sun, Y. Xi, B. Rodriguez, H. J. Park, P. Tong, M. Meong, M. A. Goodell and W. Li, “MOABS: model based analysis of bisulfite sequencing data.,” Genome Biol, vol. 15, no. 2, pp. R38, Feb 2014. 

