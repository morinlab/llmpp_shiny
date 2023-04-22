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

When you are ready, use [this page](https://shiny.rcg.sfu.ca/u/rdmorin/llmpp_shiny/) to review the variants. The interface should be relatively self-explanatory. Be sure to enter a user ID in the box to avoid your submissions being tracked under `anonymous`. We recommend using the first part of the email address associated with your GitHub account (everything before the `@`). If you have done this correctly, your ID should appear on the leaderboard on the bottom of the side pane after you have submitted at least one review. 

Things to note: 

* None of the samples have a matched normal available. The page auto-selects the No count normal (NCN) tag for you to track this. 
* Pretty much all libraries have a short insert size. This means that many of the candidate variants will be supported by overlapping F and R from the same pair. Try to take this into consideration when using the scoring and aim to consider how many distinct molecules support the mutation (rather than reads), when possible.
* Some regions of certain genes have been noted to have evidence of contaminating PCR amplicon. Examples of what this looks like can be found [in the PubPeer post describing them](https://pubpeer.com/publications/E61AC72AE0402C6A62A84E36ED2AEA#1)

