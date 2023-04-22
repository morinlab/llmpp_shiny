# Manual variant review

## Preamble

This repository contains the source code for two interactive Shiny sites that support the variant review process. _Reviewers do not need the files housed here to contribute to this project_. That said, feel free to dig into the source code if you are into! The repository also contains several thousand IGV screenshots from a single large cancer exome data set [Reddy et al, 2017](https://pubmed.ncbi.nlm.nih.gov/28985567/). The screenshots are split between two directories (`false_positives` and `false_negatives`). Respectively, these represent sets of variants that were not identified in a reanalysis of this data set with modern pipelines (`false positives`) or were newly identified by these pipelines (`false_negatives`). This division should be ignored while performing the review process to avoid biasing the opinion of reviewers as they rate variants. 

### Getting Started: Training & Calibration

Use [this page](https://shiny.rcg.sfu.ca/u/rdmorin/calibrate/) to learn the basics of the scoring system.

### Participating in Variant Review

When you are ready, use [this page](https://shiny.rcg.sfu.ca/u/rdmorin/llmpp_shiny/) to review the variants. 
