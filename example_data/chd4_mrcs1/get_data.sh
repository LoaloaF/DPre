
wget -c ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE108nnn/GSE108087/suppl/GSE108087_knockdowns_expr.csv.gz
gunzip GSE108087_knockdowns_expr.csv.gz

python3 extract_data.py
Rscript doDESeq2.R