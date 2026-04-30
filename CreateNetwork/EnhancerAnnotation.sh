module load BEDTools 

GTF="/gpfs/gibbs/project/noonan/yj345/mm10_fasta_gtf/gencode.vM10.annotation.gtf"

declare -A PCHIC
PCHIC[PA1]="GSM6505208_5209_merged.txt"
PCHIC[FL]="GSM6505208_5209_merged.txt"
PCHIC[PA2]="GSM6505210_5211_merged.txt"
PCHIC[HL]="GSM6505210_5211_merged.txt"

declare -A PEAKS
PEAKS[PA1]="P1_E10_merged_peaks_annotated.bed"
PEAKS[PA2]="P2_E10_merged_peaks_annotated.bed"
PEAKS[FL]="FL_E10_merged_peaks_annotated.bed"
PEAKS[HL]="HL_E10_merged_peaks_annotated.bed"

declare -A K27AC
K27AC[PA1]="PA1E10_K27ac.narrowPeak"
K27AC[PA2]="PA2E10_K27ac.narrowPeak"
K27AC[FL]="FLE10_K27ac.narrowPeak"
K27AC[HL]="HLE10_K27ac.narrowPeak"

PROMOTERS="promoters_2kb_up_0.5kb_down.bed"

awk 'BEGIN{OFS="\t"}
$3=="transcript"{
  gene=""; gid=""; tid="";
  if (match($0,/gene_name "([^"]+)"/,m)) gene=m[1];
  if (match($0,/gene_id "([^"]+)"/,g))   gid=g[1];
  if (match($0,/transcript_id "([^"]+)"/,t)) tid=t[1];

  chr=$1; strand=$7; start=$4; end=$5;

  tss = (strand=="+") ? start : end;

  up=tss-2000; if (up<0) up=0;
  down=tss+500;

  name = gene";"gid";"tid;

  print chr, up, down, name, 0, strand
}' "$GTF" \
| sort -k1,1 -k2,2n > "$PROMOTERS"

awk 'BEGIN{OFS="\t"} $3=="gene"{
  if (match($0,/gene_name "([^"]+)"/,n) && match($0,/gene_id "([^"]+)"/,g))
    print n[1], g[1]
}' "$GTF" | sort -u \
| awk 'BEGIN{FS=OFS="\t"}{
    m[$1]=(m[$1]?m[$1]","$2:$2)
  }
  END{
    for(s in m) print s,m[s]
  }' \
> symbol_id_pairs.tsv

awk -F'\t' 'BEGIN{OFS="\t"}{
  split($2,a,",");
  print $1,a[1]
}' symbol_id_pairs.tsv > symbol_to_id.tsv

run_per_tissue () {
  tissue="$1"

  pchic="${PCHIC[$tissue]}"
  peaks="${PEAKS[$tissue]}"
  k27ac="${K27AC[$tissue]}"
  bindetect_dir="${BINDETECT[$tissue]}"

  outdir="${tissue}_enhancer_annotation"
  mkdir -p "$outdir"

  echo "Running tissue: $tissue"
  echo "Output dir: $outdir"
 
  awk 'BEGIN{OFS="\t"}{
    split($1,a,":");
    split(a[2],b,"-");
    print a[1],b[1],b[2],"A1_"NR
  }' "$pchic" > "$outdir/anchor1.bed"

  awk 'BEGIN{OFS="\t"}{
    split($2,a,":");
    split(a[2],b,"-");
    print a[1],b[1],b[2],"A2_"NR
  }' "$pchic" > "$outdir/anchor2.bed"

  bedtools intersect -wa -wb \
    -a "$outdir/anchor1.bed" \
    -b "$PROMOTERS" \
  | awk 'BEGIN{OFS="\t"}{print $4,$8}' \
  | awk -F'\t' 'BEGIN{OFS="\t"}{
      a[$1]=a[$1]","$2
    }
    END{
      for(k in a){
        sub(/^,/,"",a[k]);
        print k,a[k]
      }
    }' \
  > "$outdir/anchor1_to_genes.tsv"

  nl -ba "$pchic" \
  | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' \
  > "$outdir/interactions_with_id.tsv"

  awk 'BEGIN{OFS="\t"}
    NR==FNR{
      map[$1]=$2;
      next
    }
    {
      id=$1;
      a1="A1_"id;
      gene=(a1 in map ? map[a1] : "NA");
      print $2,$3,$4,gene
    }' "$outdir/anchor1_to_genes.tsv" "$outdir/interactions_with_id.tsv" \
  > "$outdir/merged_with_promoterGene.tsv"

  bedtools intersect -u \
    -a "$peaks" \
    -b "$k27ac" \
  | cut -f1-3 \
  | sort -k1,1 -k2,2n \
  | bedtools merge \
  > "$outdir/active_enhancers.bed"

  awk 'BEGIN{OFS="\t"}{
    sub(/^A1_/,"",$1);
    print $1,$2
  }' "$outdir/anchor1_to_genes.tsv" \
  > "$outdir/id_to_genes.tsv"

  bedtools intersect -wa -wb \
    -a "$outdir/active_enhancers.bed" \
    -b "$outdir/anchor2.bed" \
  > "$outdir/enh_x_a2.tsv"

  awk 'BEGIN{OFS="\t"}{
    enh=$1 OFS $2 OFS $3;
    id=$7;
    sub(/^A2_/,"",id);
    print enh,id
  }' "$outdir/enh_x_a2.tsv" \
  > "$outdir/enh_id.tsv"

  awk 'BEGIN{OFS="\t"}
    NR==FNR{
      g[$1]=$2;
      next
    }
    {
      id=$4;
      genes=(id in g ? g[id] : "NA");
      print $1,$2,$3,genes
    }' "$outdir/id_to_genes.tsv" "$outdir/enh_id.tsv" \
  > "$outdir/enh_id_genes_raw.bed"

  awk 'BEGIN{OFS="\t"}{
    k=$1 OFS $2 OFS $3;
    if($4!="NA") agg[k]=(k in agg ? agg[k]","$4 : $4)
  }
  END{
    for(k in agg){
      split(k,a,FS);
      print a[1],a[2],a[3],agg[k]
    }
  }' "$outdir/enh_id_genes_raw.bed" \
  > "$outdir/active_enhancers_annotated.bed"

  awk -F'\t' 'BEGIN{OFS="\t"}{
    split($4,sets,",");
    delete seen;
    genes="";
    for(i in sets){
      split(sets[i],p,";");
      g=p[1];
      if(!(g in seen)){
        seen[g]=1;
        genes=(genes?genes","g:g)
      }
    }
    print $1,$2,$3,(genes==""?"NA":genes)
  }' "$outdir/active_enhancers_annotated.bed" \
  > "$outdir/active_enhancers_genes_only.bed"

  bedtools intersect -wa -wb \
    -a "$peaks" \
    -b "$outdir/active_enhancers_genes_only.bed" \
  | awk -F'\t' 'BEGIN{OFS="\t"}{
      peak=$4;
      n=split($12,a,",");
      for(i=1;i<=n;i++) {
        if(a[i]!="") print peak,a[i]
      }
    }' \
  | sort -u \
  | awk -F'\t' 'BEGIN{OFS="\t"}{
      m[$1]=(m[$1]?m[$1]","$2:$2)
    }
    END{
      for(k in m) print k,m[k]
    }' \
  > "$outdir/peak_to_enhancer_genes.tsv"

  gawk -F'\t' -v OFS='\t' '
    ARGIND==1 {
      sym2id[$1]=$2;
      next
    }
    ARGIND==2 {
      pk2sym[$1]=$2;
      next
    }
    ARGIND==3 {
      peak=$4;

      if (peak in pk2sym) {
        new_syms = pk2sym[peak];

        split(new_syms, S, ",");
        delete idseen;
        new_ids="";

        for (i in S) {
          s=S[i];
          if (s in sym2id) {
            n=split(sym2id[s], IDA, ",");
            for (j=1;j<=n;j++) {
              id=IDA[j];
              if (!(id in idseen)) {
                idseen[id]=1;
                new_ids=(new_ids?new_ids","id:id)
              }
            }
          }
        }

        if (new_ids!="") {
          if ($7=="" || $7=="NA" || $7==".") $7=new_ids;
          else                                $7=$7","new_ids;
        }

        if (new_syms!="") {
          if ($8=="" || $8=="NA" || $8==".") $8=new_syms;
          else                                $8=$8","new_syms;
        }
      }

      print
    }' symbol_id_pairs.tsv "$outdir/peak_to_enhancer_genes.tsv" "$peaks" \
  > "$outdir/${tissue}_E10_merged_peaks_annotated_withEnh.tsv"

  awk 'BEGIN{OFS="\t"}{
    print $1,$2,$3,$4,$5,$6,$7,$8
  }' "$outdir/${tissue}_E10_merged_peaks_annotated_withEnh.tsv" \
  | sort -k1,1 -k2,2n \
  > "$outdir/${tissue}_E10_merged_peaks_annotated_withEnh.bed"

  PEAKS_WITH_ENH="$outdir/${tissue}_E10_merged_peaks_annotated_withEnh.bed"

  done

  echo "Finished tissue: $tissue"
}

for tissue in PA1 PA2 FL HL; do
  run_per_tissue "$tissue"
done
