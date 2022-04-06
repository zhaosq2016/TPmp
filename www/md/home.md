
## Exploring miRNAs and phasiRNAs of tea plants (*Camellia Sinensis* L.)
<p>&nbsp;</p>
&nbsp;&nbsp;&nbsp;&nbsp;Since tea made from the fresh leaves of tea plant  (<em>Camellia sinensis</em> L.) contains a large number of secondary metabolites that are closely related to human health, tea has become one of the most beneficial non-alcoholic beverages in the world. miRNAs and phasiRNAs are involved in the regulation of various processes, such as plant hormone regulation, morphogenesis, development and stress response. However, there is still a lack of a genome-wide understanding of the miRNA of tea plant. In addition, the phasiRNA in tea plant has not been identified.


&nbsp;&nbsp;&nbsp;&nbsp;In this database, We collected 91 sequencing data, including 81 small RNA sequencing data and 10 degradome sequencing data.


<p>&nbsp;</p>
<p align="center">
<img src="img/miRNA_and_phasiRNA.png" width="800" hegiht="1000">
</p>
<p style="text-align:center">Fig. 1 Biogenesis of miRNA and phasiRNA, <a href="https://www.nature.com/articles/s41438-018-0072-8">[1]</a>. </p>


&nbsp;&nbsp;&nbsp;&nbsp;We used the following pipeline (Fig. 2) to process the data:


<p align="center">
<img src="img/pipeline.png" width="600" hegiht="1000">
</p>
<p style="text-align:center">Fig. 2 The pipeline for data processing procedures </p>

1. All raw single-end reads were trimmed for quality control using fastx-toolkit.
2. The remaining high-quality reads from sRNA sample were used to identify miRNAs, phasiRNAs and hcsiRNAs through sRNAminer software.
3. The remaining high-quality reads from degradome were used to identified miRNA targets by Cleaveland 4.


&nbsp;&nbsp;&nbsp;&nbsp;Finally, 952 MIRNAs, 427 PHAS locis (336 21-nt and 91 24-nt), 898,516 hcsiRNAs were identified (Fig. 3). Then we identified its target gene or trigger. The expression levels of these small RNAs in different samples were also analyzed. The response element or transcription factor binding site on the promoter sequence of miRNAs was identified.


<p align="center">
<img src="img/chromosome1.png" width="800" hegiht="500">
</p>
<p style="text-align:center">Fig. 3 The distribution of miRNAs and phasiRNAs in <em>Camellia sinensis</em> genomic regions </p>

