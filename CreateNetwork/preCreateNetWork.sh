GTF="/gpfs/gibbs/project/noonan/yj345/mm10_fasta_gtf/gencode.vM10.annotation.gtf"

awk 'BEGIN{OFS="\t"} $3=="gene"{
  if (match($0,/gene_name "([^"]+)"/,n) && match($0,/gene_id "([^"]+)"/,g))
    print n[1],g[1]
}' "$GTF" | sort -u > symbol_to_id.tsv

for f in \
  FL_E10_merged_peaks_annotated_withEnh.bed \
  HL_E10_merged_peaks_annotated_withEnh.bed \
  P1_E10_merged_peaks_annotated_withEnh.bed \
  P2_E10_merged_peaks_annotated_withEnh.bed
do
  gawk -F'\t' -v OFS='\t' '
    NR==FNR {
      sym2id[$1]=$2;
      next
    }

    {
      idlist = $7;
      symlist = $8;

      n1 = split(idlist, ids, /,/);
      n2 = split(symlist, syms, /,/);

      for(i=1;i<=n1;i++){ sub(/^ +| +$/,"",ids[i]) }
      for(i=1;i<=n2;i++){ sub(/^ +| +$/,"",syms[i]) }

      if (n1==1 && n2==1) {
        print;
        next
      }

      if (n1==n2 && n1>0) {
        for(i=1;i<=n1;i++){
          id = (ids[i]=="" || ids[i]=="NA" || ids[i]==".") ? ((syms[i] in sym2id) ? sym2id[syms[i]] : "NA") : ids[i];
          sym = (syms[i]=="" ? "." : syms[i]);

          print $1,$2,$3,$4,$5,$6,id,sym;
        }
        next
      }

      for(i=1;i<=n2;i++){
        sym = (syms[i]=="" ? "." : syms[i]);
        id  = (sym in sym2id) ? sym2id[sym] : "NA";

        print $1,$2,$3,$4,$5,$6,id,sym;
      }
    }
  ' symbol_to_id.tsv "$f" > "${f}.tmp" && mv "${f}.tmp" "$f"

done
