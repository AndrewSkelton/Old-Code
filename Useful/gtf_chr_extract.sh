awk '($1 == "15")' merged.gtf > Ensembl_merged_chr15.gtf
awk '($1 == "15" && $4 >= 89304530 && $5 <= 89460729)' merged.gtf > Ensembl_merged_ACAN.gtf
