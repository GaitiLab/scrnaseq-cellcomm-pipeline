# GaitiLab/scrnaseq-cellcomm-pipeline

<!-- [![GitHub Actions CI Status](https://github.com/GaitiLab/scrnaseq-cellcomm-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/GaitiLab/scrnaseq-cellcomm-pipeline/actions/workflows/ci.yml) -->
<!-- [![GitHub Actions Linting Status](https://github.com/GaitiLab/scrnaseq-cellcomm-pipeline/actions/workflows/linting.yml/badge.svg)](https://github.com/GaitiLab/scrnaseq-cellcomm-pipeline/actions/workflows/linting.yml) -->
<!-- [![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX) -->
<!-- [![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com) -->
<!-- [![GitHub Actions CI Status](https://github.com/GaitiLab/scrnaseq-cellcomm-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/GaitiLab/scrnaseq-cellcomm-pipeline/actions/workflows/ci.yml) -->
<!-- [![GitHub Actions Linting Status](https://github.com/GaitiLab/scrnaseq-cellcomm-pipeline/actions/workflows/linting.yml/badge.svg)](https://github.com/GaitiLab/scrnaseq-cellcomm-pipeline/actions/workflows/linting.yml) -->
<!-- [![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX) -->
<!-- [![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com) -->

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/) [![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

<!-- TODO More testing needed
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)  

<!-- [![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/GaitiLab/scrnaseq-cellcomm-pipeline) -->

## Introduction

**GaitiLab/scrnaseq-cellcomm-pipeline** is a bioinformatics pipeline that infers cell-cell interactions from scRNAseq data using various publicly available tools.

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.

### Requirements

* Unix-like operating system (Linux, macOS, etc)
* Java 18
* Nextflow 24.10.5

> Disclaimer: pipeline has been only been tested the abovementioned versions.

First, clone this GitHub repository:

```bash
git clone https://github.com/GaitiLab/scrnaseq-cellcomm-pipeline.git
```

If you run the pipeline offline, then please install the required plugin.

```bash
nextflow plugin install nf-schema@2.3.0
```

Then specify the parameters in `params.yml`, which contains the minimal parameters that need to be set:

* `input_file`, a Seurat object containing multiple samples.
* `annot`, column in Seurat object's metadata containing the annotation labels.
* `sample_var`, column in Seurat object's metadata containing the sample IDs.

> NOTE: ensure that the metadata of your Seurat object, does **not** have a column `cell_type` if `annot` is **not** "cell_type".

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run scrnaseq-cellcomm-pipeline \
   -profile <docker/singularity/.../institute> \
   --outdir <OUTDIR> -params-file "params.yml" 
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

<!-- ## Credits

GaitiLab/scrnaseq-cellcomm-pipeline was originally written by Joan Kant.

We thank the following people for their extensive assistance in the development of this pipeline: -->

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

<!-- ## Contributions and Support
If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md). -->

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use GaitiLab/scrnaseq-cellcomm-pipeline for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
