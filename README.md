# iModulon Transcriptomics Analysis of Acinetobacter baumannii

## Overview
This project applies Independent Component Analysis (ICA) to a large transcriptomic dataset to identify independently regulated gene modules (iModulons) in *Acinetobacter baumannii*. Using 1,498 RNA-seq samples collected from 136 independent studies, this analysis identifies transcriptional patterns that explain how the bacterium responds to different environmental conditions, including antibiotic stress and infection related processes.
## Dataset & ICA Results
- **1,498 RNA-seq samples**
- **136 experimental projects**
- **3,718 genes analyzed**
- **252 iModulons identified**
- **66.23% of total transcriptomic variance explained**

These iModulons represent distinct regulatory groups, including known transcription factors, stress responses, and previously uncharacterized gene groups.

## Project Structure

### Phase 1: iModulonMiner Execution
- Ran ICA to generate M (gene weight) and A (activity) matrices  
- Identified iModulon gene sets and evaluated explained variance  

### Phase 2: Metadata Curation & Database Integration
- Curated metadata across 136 projects and 705 conditions  
- Mapped UniProt IDs with ~100% coverage  
- Uploaded dataset to iModulonDB with successful validation  

### Phase 3: Biological Curation
- Annotated iModulons with functions and regulators  
- Identified known regulons (e.g., BfmR, Fur, PmrA)  
- Characterized novel and partially annotated gene modules  

### Phase 4: Deep Analysis & Interpretation
- Performed correlation analysis, PCA, and clustering  
- Identified condition-specific activity patterns  
- Developed biological hypotheses from regulatory relationships  

## Key Findings

- **Co-regulation of growth-related processes**  
  Ribosome (IM_109) and ATP-Synthase (IM_139) iModulons showed strong positive correlation (r = 0.716), suggesting coordinated regulation of protein synthesis and energy production.

- **Stress vs virulence trade-off**  
  A negative correlation (r = -0.337) was observed between Arsenic (IM_104) and Type-IV-Pili (IM_45), indicating a potential trade-off between detoxification and surface motility.

- **Novel biological hypothesis**  
  These results suggest that arsenic exposure may suppress biofilm formation by reducing pili-related gene activity.

## Why This Matters

*Acinetobacter baumannii* is a highly antibiotic-resistant pathogen responsible for serious hospital-acquired infections. Understanding its transcriptional regulatory structure is critical for identifying new therapeutic targets and understanding how it adapts to stress and treatment.

## Methods & Tools
- Independent Component Analysis (ICA)
- iModulonMiner pipeline
- RNA-seq transcriptomic data
- PCA and clustering analysis
- iModulonDB platform

## Author

Jade Cruz
