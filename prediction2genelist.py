import sys
prediction_txt = sys.argv[1]
threshold = sys.argv[2]
gene_go = {}
for line in open(prediction_txt,'r').readlines():
    terms = line.strip().split()
    if terms[0] in ['AUTHOR','MODEL','KEYWORDS','END']: continue
    if float(terms[2]) < float(threshold): continue
    if not terms[0] in gene_go:
        gene_go[terms[0]] = terms[1]
    else:
        gene_go[terms[0]] = gene_go[terms[0]]+' '+terms[1]
gene_list_txt = prediction_txt.replace('.txt','_genelist.txt')
out_fh = open(gene_list_txt,'w')
for key in gene_go:
    out_fh.write(f'{key} {gene_go[key]}\n')
out_fh.close()
print(f'write to {gene_list_txt}')
