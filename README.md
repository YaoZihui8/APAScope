# APAScope: A Modular Workflow for Quantitative Analysis of APA Events Based on Bulk RNA-Seq Data

APAScope is a bioinformatics workflow designed and optimized for the quantitative analysis of alternative polyadenylation (APA) events based on bulk RNA-Seq data. The workflow employs a modular design, making it highly extensible and adaptable to more complex analytical requirements.

## Workflow Overview


### **Module 1: High-Quality Poly(A) Site Identification**
- Users can customize filtering criteria based on the specific characteristics of their data.
- Multiple filtering steps ensure the accuracy and reliability of identified poly(A) sites.
- High-quality poly(A) sites are validated through comparison with authoritative databases such as PolyASite and PolyA DBv3.
- The output includes high-quality poly(A) site data annotated with the number of supporting RNA-seq reads as an indicator of site usage intensity.

### **Module 2: APA Reference Library Construction**
- The high-quality poly(A) site data from Module 1 is input into this module.
- APA reference libraries, including 3' UTR and intronic polyadenylation (IPA) sites, are generated using two quantitative tools: QAPA and APAlyzer.

### **Module 3: Reference Library Integration and Update**
- The APA reference libraries generated in Module 2 are integrated with the original reference libraries.
- This step overcomes the limitations of existing annotation libraries, ensuring comprehensive and specific APA regulation analysis for target tissues or diseases.

### **Module 4: Quantitative Analysis of APA Events**
- Enables precise quantification of APA events tailored to specific research objectives.

## Key Advantages
- **Versatility:** The workflow is applicable to any eukaryotic organism. Other eukaryotic organisms can directly use the identified poly(A) sites to construct reference libraries for 3' UTR APA and IPA event analysis.
- **Reliability:** The 3' UTR APA analysis incorporates a dual-platform cross-validation mechanism, enhancing the reliability and comprehensiveness of the results.
- **Scalability:** The modular design allows for easy expansion and customization to accommodate more complex analytical needs.

## Availability
The APAScope source code is openly available on GitHub, accompanied by detailed documentation and environment dependencies. Users can directly download and use the workflow for their studies.

Through this workflow, we aim to uncover more biological insights into APA from large-scale bulk RNA-seq datasets.

## Citation
If you use APAScope in your research, please cite our work.
