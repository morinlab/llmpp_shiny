# Manual variant review

## Preamble

This repository contains the source code for two interactive Shiny sites that support the variant review process. _Reviewers do not need the files housed here to contribute to this project_. That said, feel free to dig into the source code if you are into! The repository also contains several thousand IGV screenshots from a single large cancer exome data set [Reddy et al, 2017](https://pubmed.ncbi.nlm.nih.gov/28985567/). The screenshots are split between two directories (`false_positives` and `false_negatives`). Respectively, these represent sets of variants that were not identified in a reanalysis of this data set with modern pipelines (`false positives`) or were newly identified by these pipelines (`false_negatives`). This division should be ignored while performing the review process to avoid biasing the opinion of reviewers as they rate variants. 

### Getting Started: Training & Calibration

The reviewing will follow the guidelines detailed in [this paper](https://www.gimjournal.org/article/S1098-3600(21)00974-6/fulltext#ec0015) along with a streamlined scoring system we devised to simplify the process. Use [this page](https://shiny.rcg.sfu.ca/u/rdmorin/calibrate/) to learn the basics of the scoring system. Essentially, the system uses a 5-point scale to grade variants based on their quality. The maximal score is reserved for variants with an ideal amount of support and no issues that cause reduction in the reviewers confidence in its accuracy. Hence, a variant with good support but one or more confounders should be down-graded to a lower score.

0. Our scoring scale allows for variants to be given a score of 0 for when it has absolutely no evidence in the data according to IGV. This scenario should _theoretically_ never happen.
1. Reseved for variants that appear to have the bare minimum of support. Either a single read supporting the variant or, in the case of short-insert libraries, you will likely see two reads (F and R from the same molecule). Variants should rarely be downgraded to this category otherwise.
2. For variants with slightly above the minimal possible support (i.e. 2-3 molecules from up to 6 total reads, in the case of short insert libraries). You may downgrade a variant with higher support to this score based on the presence of other confounders 
3. For variants with less-than-ideal (i.e. modest) support but exceeding that in the category below. You may downgrade a variant with stronger support to this score based on the presence of other confounders 
4. For variants with almost ideal support or with ideal support and possibly one confounder
5. For the ideal variant with no confounders


### Participating in Variant Review

When you are ready, use [this page](https://shiny.rcg.sfu.ca/u/rdmorin/llmpp_shiny/) to review the variants. 
