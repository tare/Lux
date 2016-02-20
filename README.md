Lux
===================

Overview
-------------
Lux is an integrative hierarchical Bayesian model for analyzing BS-seq and oxBS-seq data enabling accurate and unbiased quantification of different cytosine modifications and differential methylation at individual cytosines or loci, with or without replicates, while taking imperfect experimental parameters into account.

Navigation
-------------
[TOC]


Features
-------------
- Model-based integration and analysis of BS-seq and oxBS-seq data from whole genome, reduced representation or targeted experiments
- Consideration of nonideal experimental parameters through modeling
- Model-based integration of biological replicates
- Detection of differential methylation at individual cytosines or regions using Bayes factors (DMRs)
- Full Bayesian inference with Hamiltonian Monte Carlo (HMC) using No-U-Turn sampler (NUTS) implemented in **Stan**

Downloads
-------------
#### Example data

Here we provide the input data sets and the obtained HMC chains used in the publication.

Input data
[v65_t2kd_data.tar.bz2](https://dl.dropboxusercontent.com/u/2392459/v65_t2kd_data.tar.bz2)

Output chains
[v65_t2kd_chains.tar.bz2](https://dl.dropboxusercontent.com/u/2392459/v65_t2kd_chains.tar.bz2)

Quick introduction
-------------
An usual Lux pipeline has the following steps
1. Alignment of BS-seq and oxBS-seq data (e.g., [Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/) or [BSmooth](http://rafalab.jhsph.edu/bsmooth/)
2.  Extraction of converted and unconverted counts (e.g., [Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/) or [BSmooth](http://rafalab.jhsph.edu/bsmooth/))
3. Generation of input files for **Lux** (**parse.py**)
4. Integrative methylation analysis (**Lux**)
5. Calculation of Bayes factors (**bf.py**)

This documentation focus on the the points three, four, and five.

Installation
-------------
Below the environmental parameters **\$STAN_HOME** and **\$LUX_HOME** refers to the source directories of **CmdStan** and **Lux**.

- **Lux** requires **CmdStan** (tested on version 2.2.0)
- **parse.py** and **bf.py** requires **Python** (tested on version 2.7.5), **NumPy** (tested on version 1.7.1) and **Scipy** (tested on version 0.13.0)

#### Installing CmdStan
To install **CmdStan**, please first download the source code of **CmdStan** from http://mc-stan.org/interfaces/cmdstan.html and follow the included installation instructions.

To cut the long story short, after extracting the source code running the following commands should be sufficient for our purposes

    cd $STAN_HOME
    make bin/print
    make bin/stanc

The binaries stanc and print are used for translating **Stan** codes into executables and summarizing HMC chains, respectively.

For more details on installing **CmdStan**, please see **\$STAN_HOME/README.txt** and the documentation of **CmdStan** available at [http://mc-stan.org/manual.html].

#### Compiling Lux
The Lux source code is supplied in the file **\$LUX_HOME/lux.stan**.
The use of Lux requires that it is first compiled for the platform in use, which can be done using the installed **CmdStan**

    cd $STAN_HOME
    make $LUX_HOME/lux

After the compilation is done the directory **$LUX_HOME** should have a binary with the name lux.

For more details on the compilation, such as code optimization, please see the documentation of **CmdStan** available at http://mc-stan.org/documentation/.

Using Lux
-------------
After completing the previous steps **Lux** is ready to be used. As an example, one can run **Lux** with the supplied cytosine-level data sets used in the manuscript as follows

    cd $LUX_HOME
    ./lux method=sample algorithm=hmc engine=nuts max_depth=8 stepsize=0.02 data file=data/v65_data.R init=data/v65_init.R output file=v65_output.csv
    ./lux method=sample algorithm=hmc engine=nuts max_depth=8 stepsize=0.02 data file=data/t2kd_data.R init=data/t2kd_init.R output file=t2kd_output.csv

**Stan** uses the dump data format, and thus the data and init files have to be in that format. A procedure for generating these data files is described below.

For more details, please see the documentation of **CmdStan** available at http://mc-stan.org/documentation/.

### Generating input files
A python script **parse.py** is supplied for generating input files for the use with **Lux**. Basically it transforms the user-supplied data files into the data format read by **CmdStan**, in addition, it initializes the model parameters with values sampled from the corresponding priors

    usage: parse.py [-h] -d DATA -p PRIOR -cd CONTROL_DATA -cp CONTROL_PRIOR -pr PREFIX [-v]
    
    Generates data and init files in the dump format for Lux
    
     optional arguments:
     -h, --help    show this help message and exit
     -d DATA, --data DATA    noncontrol cytosine data
     -p PRIOR, --prior PRIOR    prior of the noncontrol cytosines
      -cd CONTROL_DATA, --control-data CONTROL_DATA    control cytosine data
      -cp CONTROL_PRIOR, --control-prior CONTROL_PRIOR    priors of the control cytosines
      -pr PREFIX, --prefix PREFIX    prefix of the output files
     -v, --version    show program's version number and exit

For instance, the script **parse.py** is called as

    ./parse.py --data data.tsv --prior prior.tsv --control-data control_data.tsv --control-prior control_prior.tsv --prefix sample

This command will generate the two files, **sample_data.R** and **sample_init.R**, which are ready to be used in **Lux**.

As you notice **parse.py** takes four files as input, so next, we will go through their content and format in detail.

#### Noncontrol cytosines
The files **data.tsv** and **prior.tsv** have the count data and prior information on the noncontrol cytosines of interest, respectively.

Each of the cytosines has its own line in the files, and thus the files are required to have the same number of lines, moreover, the order of the cytosines is assumed to be the same in the files. 

Each line in the file **data.tsv** is composed of one or more replicate-specific tab-separated blocks of four tab-separated nonnegative integers
>$N_\text{BS}^\text{C}$\t$N_\text{BS}$\t$N_\text{oxBS}^\text{C}$\t$N_\text{oxBS}$

That is, the number of Cs and the total BS-seq read-outs are listed first, which are followed by the number of Cs and total oxBS-seq read-outs.

Moreover, on each line there should be exactly $4 \times N_\text{replicates}$ values separated with the tabs.

The replicate-level in the hierarchical model formulation explicitly assumes that the replicates share the same prior. Thus, the prior of $\mu$ is defined by giving the parameter $\alpha$ for all the noncontrol cytosines in terms of three pseudo-counts $\alpha_1$, $\alpha_2$ and $\alpha_3$ corresponding to $p(\text{C})$, $p(\text{5mC})$ and $p(\text{5hmC})$, respectively.

That is, all the lines in the file **prior.tsv** are required to have the following format
>$\alpha_1$\t$\alpha_2$\t$\alpha_3$

In the original publication these were $\alpha_1=\alpha_2=\alpha_3=0.8$ for all the noncontrol cytosines.

#### Control cytosines
The files **control_data.tsv** and **control_prior.tsv** have the count data and prior information on the control cytosines of interest, respectively.

Luckily, the data for the control cytosines is supplied in the same format as noncontrol data (see *above*). That is, each line in the file **control_data.tsv** is composed of one or more replicate-specific tab-separated blocks of four tab-separated nonnegative integers.

As in **data.tsv**, on each line there should be exactly $4 \times N_{\text{replicates}}$ values separated with the tabs. Each of the replicate specified in **data.csv** should have its own control data. Moreover, the order of the replicate-specific blocks between the noncontrol and control data is assumed to be the same.

The prior knowledge on the control cytosines is supplied in the file **control_prior.tsv**. Although, the hierarchical model allows that the control cytosines would have different priors between replicates but this is not implemented in the current version. Therefore, each line in **control_prior.tsv** should have exactly three tab-separated values and the order of the rows, i.e., control cytosines, should be the same in **control_data.tsv** and **control_prior.tsv**.

### Summarizing HMC chains
After the runs are completed their output can be summarized using **print** as follows

    $STAN_HOME/bin/print v65_output.csv
    $STAN_HOME/bin/print t2kd_output.csv

For instance, the summary includes posterior means, standard deviations, quantiles and convergence diagnostic statistics for each of the parameters.

For more details on the output of **print**, please see the documentation of **CmdStan** available at http://mc-stan.org/documentation/.

#### Interpreting the parameters in the implementation
The replicate-specific experimental parameters $\text{BS}_\text{eff}$, $\text{ox}_\text{eff}$, $\text{BS}^*_\text{eff}$ and $\text{seq}_\text{err}$ are stored in the parameter vectors **bsEff**, **oxEff**, **bsBEff** and **seqErr**, respectively.

The methylation proportion parameters of noncontrol cytosines $\mu$ and $\theta$ have the following variable names in the implementation **mu** and **theta**, respectively.
In the implementation **mu** is a two-dimensional matrix ($N_\text{cytosines} \times N_\text{replicates}$), where the elements are 3-simplices. Whereas, **mu** is a vector ($N_\text{cytosines}$) in which the elements are 3-simplices.

The methylation proportion parameters of control cytosines $\theta_\text{control}$ has the following variable names in the implementation **theta_control**. The parameters **theta_control** is a vector ($N_\text{control cytosines}$) in which the elements are 3-simplices.

For more details on this, please see the parameter declarations in the source code of **Lux** in the file **\$LUX_HOME/lux.stan**.

### Detecting differential methylation
A python script **bf.py** is supplied for calculating Bayes factors

    usage: bf.py [-h] -c1 CHAIN1 -c2 CHAIN2 [-v]
    
    Calculates Bayes factor between two conditions based on the Stan HMC output chains
    
    optional arguments:
      -h, --help    show this help message and exit
      -c1 CHAIN1, --chain-1 CHAIN1    output chain of the first condition
      -c2 CHAIN2, --chain-2 CHAIN2    output chain of the second condition
      -v, --version    show program's version number and exit

For instance, to calculate Bayes factor from the supplied HMC chains used in the publication, type

    ./bf.py --chain-1 v65_output.csv --chain-2 t2kd_output.csv

Then the first four output lines are
> mu\[1]\t0.000087\n
> mu\[2]\t0.000092\n
> mu\[3]\t0.000876\n
> mu\[4]\t0.000130

The script calculates Bayes factor for every **mu** variable.

Notably, the chains are required to have the same number of **mu** variables, i.e., the same cytosines, and they should be supplied in the same order. However, calculation of Bayes factors between conditions with different number of replicates is supported. Additionally, **bf.py** assumes that $\alpha_1=\alpha_2=\alpha_3=0.8$.

For more details on the the calculation of Bayes factor, please see the original publication in which the procedure is explained in detail.

### Advanced uses
Some users might want to consider of using **Lux** through the **RStan** or **PyStan** interfaces (see http://mc-stan.org/interfaces/).

### References
[1 ] M. J. Betancourt, “Generalizing the No-U-Turn Sampler to Riemannian Manifolds,” arXiv, vol. 1304, no. 1920, 2013. 
[2] K. D. Hansen, B. Langmead and R. A. Irizarry, “BSmooth: from whole genome bisulfite sequencing reads to differentially methylated regions.,” Genome Biol, vol. 13, no. 10, pp. R83, Oct 2012. 
[3] M. D. Hoffman and A. Gelman, “The No-U-Turn Sampler: Adaptively Setting Path Lengths in Hamiltonian Monte Carlo,” Journal of Machine Learning Research, vol. in press, 2013. 
[4] F. Krueger and S. R. Andrews, “Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications.,” Bioinformatics, vol. 27, no. 11, pp. 1571-1572, Jun 2011. 
[5] Stan Development Team, “Stan Modeling Language Users Guide and Reference Manual, Version 2.2,” 2014. 
[6] Stan Development Team, “Stan: A C++ Library for Probability and Sampling, Version 2.2,” 2014. 
[7] D. Sun, Y. Xi, B. Rodriguez, H. J. Park, P. Tong, M. Meong, M. A. Goodell and W. Li, “MOABS: model based analysis of bisulfite sequencing data.,” Genome Biol, vol. 15, no. 2, pp. R38, Feb 2014. 

