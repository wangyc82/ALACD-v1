# ALACD-v1
ALACD is to get the drug-lncRNA assoications via bi-level optimization and decoded in R.
It includes SVM standard algorithm for argumenting the genes' co-expressed lncRNAs (enhancering-lncRNA-mRNA-coexpression.R), and the 
optimization procedure for optimal identification of genes aassociated with both lncRNAs and drug to get the drug-lncRNA assoications (optimizating.R).
crossvalSVM.R and getperf.R are for SVM cross-validtion and generating the evluation criteria, respectively.
It also includes the DEGanalysis.R for differential expression analysis, and
survival-analysis.R for surivival analysis of lncRNAs.

The example data for running the ALACD (BRCA-example.RData) is available at https://mail.cstnet.cn/coremail/XT3/nf/doGetFile.jsp?mid=7%3a1tbiBwMCBllo9S9XCwAAs3&sid=SAghVAOOLPeXIIzjkJOOpTGMHWwkWpOf
