
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SynSigRun

Run Mutational Signature Analysis Software Packages Using Mutational
Spectra generated by
[SynSigGen](https://github.com/steverozen/SynSigGen) package.

## Purpose

An easy-to-use package for non-experts which runs software packages
reproducibly with synthetic tumors generated by
[SynSigGen](https://github.com/steverozen/SynSigGen). SynSigRun gives
necessary information to mutational-signature analysis programs. These
programs used catalogs of synthetic mutational spectra created by
package [SynSigGen](https://github.com/steverozen/SynSigGen), and
results were assessed by
[SynSigEval](https://github.com/WuyangFF95/SynSigEval).

## Installation

Install the development version of SynSigRun from
[GitHub](https://github.com/WuyangFF95/SynSigRun) with the R command
line:

``` r
install.packages("devtools")
devtools::install_github("WuyangFF95/SynSigRun")
```

## Usage

-   `1_data_generation` - Calls `SynSigGen` generation script to
    generate 20 SBS1-SBS5 datasets at `data-raw/` or other repositories.
-   `2_running_approaches` - running computational approaches using
    `SynSigRun` wrapper functions, and output results as a 5-level
    folder structure:

Level 1: Datasets (e.g. `S.0.1.Rsq.0.1`);  
Layer 2: Folder `sp.sp`:

-   Meaning of 1st “sp”: synthetic tumors were generated using signature
    profiles extracted by SigProfiler in [Alexandrov et al.,
    2020](https://www.nature.com/articles/s41586-020-1943-3). The
    signature profiles are available
    [here](https://www.synapse.org/#!Synapse:syn12025148).
-   Meaning of 2nd “sp”: according to [Alexandrov et al.,
    2015](https://www.nature.com/articles/ng.3441), many tumor types
    have a mean SBS1 exposure of \~300, inferred by SigProfiler.
    Therefore, we chose the `log10(exposure of SBS1)` in 20 synthetic
    datasets to be `2.5`.

Level 3: De-novo extraction (`ExtrAttr`), or extraction with
ground-truth signatures known (`ExtrAttrExact`);  
Level 4: Results of computational approaches (e.g. `hdp.results`);  
Level 5: Results of runs with seeds (e.g. `seed.1`, `run.1`).

-   `3_evaluation` - evaluating performance of signature extraction by
    calling evaluation functions in `SynSigEval`.

## Reference manual

<https://github.com/WuyangFF95/SynSigRun/blob/master/data-raw/SynSigRun_0.1.2.pdf>
