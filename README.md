
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SynSigRun

Run mutational signature analysis software packages Packages and
benchmarking the performance of these packages.

## Purpose

A package to 1. wrap R-based signature analysis packages in functions
handy for non-expert users, by wrapping default argument values and all
necessary steps in the function bodies. 2. reproduce benchmarking
analysis of signature analysis packages in papers by [Rozen
Lab](https://github.com/Rozen-Lab).

Typically, a benchmarking analysis to evaluation accuracy of signature
extraction and/or exposure inference involves the following steps:

-   Generation of synthetic tumor spectra based on signatures and
    synthetic tumor exposures using wrapper functions in
    [`SynSigGen`](https://github.com/steverozen/SynSigGen). Usually,
    -   signatures are real signatures downloaded from
        [COSMIC](https://cancer.sanger.ac.uk/signatures/),
    -   synthetic tumor exposures are drawn from a distribution which
        mimics the distribution of a real tumor type.
-   Run of computational approaches (can be an R/Python/Julia/C++
    package) on generated data sets. It involves two steps:
    -   Number of signatures (*K*) to be extracted is estimated
        -   heuristically - starting from an user-provided or default
            number, without providing the range of possible *K* s.
        -   semi-automatically - selecting from a range of possible
            *K* s.
        -   manually - requiring users to choose the best *K* based on
            the diagnostic plot or table generated by the software.
    -   Extract a specific number of signatures, AND
    -   Infer the exposures of extracted signatures For computational
        approaches based on R and can do
    -   Signature extraction which heuristically or semi-automatically
        selects *K* AND/OR
    -   exposure inference (attribution), we wrote wrapper functions in
        **`R/`** folder of this package for non-expert users to run
        these approaches with a simple function call.
-   Evaluation of accuracy on signature extraction AND/OR exposure
    inference. Many of the evaluation functions are in package
    [`SynSigEval`](https://github.com/steverozen/SynSigEval).

The full

## Installation

Install the development version of `SynSigRun` from
[GitHub](https://github.com/WuyangFF95/SynSigRun) with the R command
line:

``` r
install.packages("devtools")
devtools::install_github("WuyangFF95/SynSigRun")
```

## Usage

### Benchmarking analysis in Alexandrov et al. (2020)

*Nature* paper “The repertoire of mutational signatures in human cancer”
([link](https://doi.org/10.1038/s41586-020-1943-3)) involves
benchmarking analysis compared to
[`SigProfiler`](https://www.mathworks.com/matlabcentral/fileexchange/38724-sigprofiler)
(the ancestor of
[`SigProfilerExtractor`](https://pypi.org/project/SigProfilerExtractor/))
and `SignatureAnalyzer`.

It used some functions and top-level codes in this package. Some of the
codes are in `data-raw/Alexandrov_2020`.

### Re-produce benchmarking analysis in Wu et al. (2022)

Manuscript under revision in *Scientific Reports*.

In order to reproduce benchmarking of signature extraction accruacy on
synthetic spectra with correlated exposures to SBS1 and SBS5 signatures,
users can go to `data-raw/Wu_2022/1_scripts.for.SBS1SBS5` to generate
the main figure and the full data of this analysis. The sub-folders hold
scripts for:

-   `1_data_generation` - Calls `SynSigGen` generation script to
    generate 20 SBS1-SBS5 datasets at `data-raw/` or other repositories.

-   `2_running_approaches` - running computational approaches directly
    or using `SynSigRun` wrapper functions. The results are generated as
    a 5-level folder structure:

Level 1: Datasets (e.g. `S.0.1.Rsq.0.1`);  
Level 2: De-novo extraction without specifying `K = 2` (`ExtrAttr`), or
extraction with number of ground-truth signature `K = 2` provided to
computational approaches (`ExtrAttrExact`);  
Level 3: Results of computational approaches (e.g. `hdp.results`);  
Level 4: Results of runs with seeds (e.g. `seed.1`, `run.1`).

-   `3_evaluation` - evaluating performance of signature extraction by
    calling evaluation functions in `SynSigEval`.

### Re-produce benchmarking analysis in Liu et al. (2022)

The paper for new computational approach `mSigHdp`, “mSigHdp:
hierarchical Dirichlet processes in mutational signature extraction”,
Liu et al. (2022) (Manuscript in revision) includes a benchmarking study
on real-tumor-based synthetic spectra with SBS or indel mutations.

The benchmarking code of this study calls the wrapper function in
`SynSigRun` to run computational approaches
[`signeR`](https://bioconductor.org/packages/release/bioc/html/signeR.html)
and [`SignatureAnalyzer`]().

## Reference manual

<https://github.com/WuyangFF95/SynSigRun/blob/master/data-raw/SynSigRun_1.0.0.pdf>
