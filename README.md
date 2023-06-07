# NextFlow Pipeline for Inferring Relatedness with BioMark X SNP Panel Data
_Please note! This workflow is under development and subject to frequent near-term revisions._

## Overview
In early 2023, it was revealed that some of the research colony monkeys were illegally harvested from the wild, rather than being the offspring of a documented research animal. The resulting controversy was a major stress-test for the biomedical research community. It also brought organizations who maintain or ship research animals under scrutiny, even when they did things "by the books," with high ethical standards.

To address this situation and serve the research community, the Genomic Services Unit at the Wisconsin National Primate Research Center began developing a panel of SNPs that could, in conjunction with a BioMark X SNP chip instrument, be used to determine whether two animals that are expected to be related are in fact unrelated. This makes it possible to determine whether any animals, unbeknownst to the colony managers, are in fact of dubious origin.

This pipeline, in its final form, will automate the processing of BioMark X intensity data and use animal metadata to generate informative reports. Specifically, it will:
1. Automatically select optimal assay cycles from which to draw data (a CSV spreadsheet of preferred cycles can also be provided).
2. Generate a wide pivot table, complete with helpful color-coding and other formatting/
3. Generate confidence statistics that contextualize the likelihood that two animal are _not_ related given concordance on the SNP panel.
4. Produce an easy-to-read final report and certificate of relatedness for any desired animals.

At this stage, the automated selection of cycles, spreadsheet formatting, and report/certificate generation are still on the horizon. The scaffolding of the pipeline is, however, now established. The workflow is fully containerized and reproducible, and can be run on a variety of POSIX-ct compatible machines ranging from laptops to AWS servers to HPC clusters. We will also produce a jupyter notebook version of the pipeline in the near future, so stay tuned!
