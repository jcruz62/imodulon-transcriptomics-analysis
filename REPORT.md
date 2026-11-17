# Independently Modulated Gene Sets in *Acinetobacter baumannii* Across 1,498 Transcriptomes

**Jade Cruz**
Original Investigation · Microbial Genomics

---

## Key Points

**Question** Which regulatory programs in *Acinetobacter baumannii* turn on and off together, and do those relationships hold up after accounting for the fact that samples come in studies?

**Findings** Independent component analysis of 1,498 RNA-seq samples from 136 studies identified 252 independently modulated gene sets (iModulons) that explained 66.2% of expression variance. Ribosome and ATP synthase activity were strongly correlated (*r* = 0.71; 95% CI, 0.61-0.79), and the correlation held in 86% of individual studies. Arsenic resistance and type IV pili activity were inversely correlated (*r* = −0.32; 95% CI, −0.43 to −0.22), but in only 63% of studies.

**Meaning** Growth-related programs in *A. baumannii* are tightly coordinated. The proposed trade-off between metal detoxification and pili-driven biofilm formation is a hypothesis worth testing experimentally, not an established finding.

---

## Abstract

**Importance** *A. baumannii* is a hospital-acquired pathogen that the World Health Organization lists as a critical priority because many strains resist last-line antibiotics. Its regulatory network is poorly mapped: only a few dozen regulatory interactions have been confirmed experimentally. That makes large expression datasets hard to interpret.

**Objective** To break the *A. baumannii* transcriptome into independently regulated gene sets, and to test which of these programs are coordinated across conditions using methods that respect study structure.

**Design, Setting, and Participants** This was a secondary analysis of 1,498 public RNA-seq samples from 136 projects in the National Center for Biotechnology Information Sequence Read Archive and Gene Expression Omnibus. The samples covered antibiotic exposure, stress, virulence and gene-knockout experiments.

**Exposures** Experimental condition (738 project-condition groups).

**Main Outcomes and Measures** iModulon gene membership and activity. Associations between iModulon activities were measured with Pearson correlation at the condition level, with 95% CIs from a project-level cluster bootstrap (2,000 resamples) and consistency within each project.

**Results** The 252 iModulons covered 3,718 genes and explained 66.2% of expression variance. Of these, 107 (42%) were single-gene components. Ribosome (IM_109) and ATP synthase (IM_139) activities were strongly correlated (*r* = 0.71; 95% CI, 0.61-0.79; median within-project *r* = 0.64). Arsenic resistance (IM_104) and type IV pili (IM_45) activities were inversely correlated (*r* = −0.32; 95% CI, −0.43 to −0.22; median within-project *r* = −0.35). However, the pili correlation had the same sign in only 44 of 70 projects (63%). Principal component analysis and k-means clustering separated a distinct group of 80 conditions, mainly from viable-but-nonculturable and PmrAB studies.

**Conclusions and Relevance** Coordinated regulation of protein synthesis and energy production in *A. baumannii* was robust across studies. The arsenic-pili relationship was moderate and inconsistent across studies, so whether metal stress suppresses biofilm formation will need controlled experiments to confirm.

---

## Introduction

*Acinetobacter baumannii* is a gram-negative bacterium that causes pneumonia, bloodstream infections and wound infections, mostly in patients who are already critically ill. What makes it especially concerning is that it can survive for long periods on dry hospital surfaces and readily picks up resistance genes from other bacteria. Many strains now resist colistin, a last-resort antibiotic.

Despite its clinical importance, little is known about how *A. baumannii* controls its genes. For the ATCC 19606 reference strain, only 29 regulatory interactions and 7 transcription factors have been confirmed experimentally. Without that map, it is hard to tell from expression data alone which genes respond to which signals.

Independent component analysis (ICA) helps fill this gap. ICA is an unsupervised machine learning method that splits a gene expression matrix **X** into two parts, **X ≈ M · A**. Each column of **M** lists gene weights for one component, and the genes with large weights form an *iModulon* (an independently modulated gene set). Each row of **A** gives that iModulon's *activity* in every sample. Principal component analysis (PCA) looks for orthogonal directions of maximum variance. ICA instead looks for components that are statistically independent, which tend to line up with real regulators, because each regulator controls its own group of genes.

This study asked three questions. First, which regulatory programs are active together? Second, do conditions group by their overall regulatory state? Third, do the relationships hold *within* individual studies, or only across them? The third question matters because samples from the same lab, strain and protocol are not independent. Ignoring that can make correlations look stronger and more certain than they are.

## Methods

### Data and Decomposition
RNA-seq samples were collected with BioProject and BioSample identifiers, then processed and quality-controlled with the iModulonMiner pipeline. ICA was run on log-transformed transcripts per million (log-TPM) values. iModulon membership was set by gene-weight thresholds. Metadata were curated for every project and condition, UniProt identifiers were mapped (about 100% coverage), and the dataset was uploaded to iModulonDB (0 failed validation checks). iModulons were annotated by function and regulator, for example ribosome, ATP synthase, arsenic resistance, PmrA and Fur.

### Statistical Analysis
Biological replicates were averaged within each project-condition pair, so each condition counted once (738 conditions). The following were then calculated:

1. **Condition-level correlations** (Pearson and Spearman) between iModulon activities.
2. **Project-level cluster bootstrap.** Whole projects were resampled with replacement 2,000 times to build 95% CIs that reflect between-study variability.
3. **Within-project correlations** for projects with at least 6 samples and real variation in both iModulons. This separates true co-regulation from differences between studies, such as strain or growth medium.
4. **PCA and k-means clustering** on condition-level activity profiles. The number of clusters (*k* = 2-8) was chosen by silhouette score.

Analyses used Python 3 (pandas, NumPy, SciPy, scikit-learn, matplotlib). Code is in `analysis/imodulon_analysis.py`.

## Results

### Decomposition
The 252 iModulons explained 66.2% of total expression variance across 3,718 genes. The largest iModulon contained 101 genes. The 107 single-gene iModulons (42%) each captured variation in a single gene.

### Coordinated Growth Program
Ribosome (IM_109; 38 genes) and ATP synthase (IM_139; 23 genes) activities were strongly and consistently correlated (**Table**). The correlation was almost unchanged after averaging replicates (sample level, *r* = 0.716; condition level, *r* = 0.715). It also held inside 76 of 88 projects (86%).

**Table. Associations Between Key iModulon Activities**

| Pair | Condition-level *r* (95% CI)ᵃ | Spearman ρ | Within-project median *r* | Same sign, No./total (%) |
|---|---|---|---|---|
| Ribosome vs ATP synthase | 0.71 (0.61 to 0.79) | 0.67 | 0.64 | 76/88 (86) |
| Arsenic resistance vs type IV pili | −0.32 (−0.43 to −0.22) | −0.33 | −0.35 | 44/70 (63) |

ᵃ 95% CI from a project-level cluster bootstrap (2,000 resamples); n = 738 conditions.

### Arsenic Resistance and Type IV Pili
Arsenic resistance (IM_104; 6 genes; 1.9% of variance) was inversely correlated with type IV pili (IM_45; 35 genes), and the bootstrap CI excluded zero. Two features call for caution, though. First, IM_104 activity was close to on/off (bimodal), so the correlation mostly reflects a difference between "on" and "off" conditions rather than a smooth relationship. Second, IM_104 was most active in conditions where no arsenic was added, such as viable-but-nonculturable cells and stationary phase. This suggests the iModulon tracks strain background or a general stress state, and not only arsenic exposure.

### Condition Structure
The first 2 principal components explained 17.5% and 11.4% of condition-level variance. The silhouette score was highest at *k* = 2 (0.36). One cluster of 80 conditions, mainly from viable-but-nonculturable and PmrAB studies, separated clearly from the other 658. It was driven by IM_7, IM_32 and the arsenic iModulon.

## Discussion

This study found that, across almost 1,500 transcriptomes, *A. baumannii* tightly coordinates the machinery it needs to grow: ribosomes to make protein and ATP synthase to supply energy. This coupling is expected biologically, which makes it a useful check that the decomposition recovers real regulation. It also stayed strong once replicates and same-study samples were no longer treated as independent.

The arsenic-pili relationship is more interesting but less certain. Type IV pili help *A. baumannii* twitch across surfaces and form biofilms, which are a major reason it persists in hospitals. If metal stress suppresses pili, that would point to a trade-off between detoxification and colonization. The robustness checks show the association is real across the dataset but inconsistent between studies. It is best treated as a **hypothesis**: that arsenic exposure reduces type IV pili expression and biofilm formation. A controlled experiment that adds arsenite to one strain and measures pili gene expression and biofilm mass would test this directly.

### Limitations
First, 42% of iModulons were single-gene components, which suggests the decomposition dimension was set higher than needed and some components capture noise. Second, the samples mixed several strains (ATCC 17978, ATCC 19606, AB5075 and clinical isolates), and differences in gene content between strains can create apparent co-regulation. Third, all findings are correlational, and regulator assignments relied on a small curated regulatory network.

## Conclusions

ICA recovered a robust, biologically sensible growth program in *A. baumannii* and suggested a possible trade-off between metal stress and pili-driven biofilm formation. Checking results at the condition and project level, not just across all samples, changed how strongly each finding can be claimed.
