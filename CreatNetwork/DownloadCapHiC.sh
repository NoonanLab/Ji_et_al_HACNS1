wget -O GSE211902_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE211902&format=file"
tar -xvf GSE211902_RAW.tar
gunzip *.gz

awk '{
  split($1,a,","); split($2,b,",");
  key = a[1]":"a[2]"-"a[3]","b[1]":"b[2]"-"b[3];
  sum[key]+=$3; count[key]++;
}
END{
  for(k in sum){
    split(k,s,",");
    print s[1]"\t"s[2]"\t"sum[k]/count[k];
  }
}' GSM6505208_Md_CapHiC_Merged_All-Tr_2ndBioRep_Chicago_output_washU_text.txt \
   GSM6505209_Md_CapHiC_3rdBioRep_Chicago_output_washU_text.txt \
   > GSM6505208_5209_merged.txt

awk '{
  split($1,a,","); split($2,b,",");
  key = a[1]":"a[2]"-"a[3]","b[1]":"b[2]"-"b[3];
  sum[key]+=$3; count[key]++;
}
END{
  for(k in sum){
    split(k,s,",");
    print s[1]"\t"s[2]"\t"sum[k]/count[k];
  }
}' GSM6505210_PA2_CapHiC_2ndBioRep_Chicago_output_washU_text.txt \
   GSM6505211_PA2_E10.5wt_CapHiC_3rdBioRep_Tr1_Chicago_output_washU_text.txt \
   > GSM6505210_5211_merged_avg.txt

awk '{
  split($1,a,","); split($2,b,",");
  key = a[1]":"a[2]"-"a[3]","b[1]":"b[2]"-"b[3];
  sum[key]+=$3; count[key]++;
}
END{
  for(k in sum){
    split(k,s,",");
    print s[1]"\t"s[2]"\t"sum[k]/count[k];
  }
}' GSM6505210_PA2_CapHiC_2ndBioRep_Chicago_output_washU_text.txt \
   GSM6505211_PA2_E10.5wt_CapHiC_3rdBioRep_Tr1_Chicago_output_washU_text.txt \
   > GSM6505210_5211_merged_avg.txt

awk '{
  split($1,a,","); split($2,b,",");
  key = a[1]":"a[2]"-"a[3]","b[1]":"b[2]"-"b[3];
  sum[key]+=$3; count[key]++;
}
END{
  for(k in sum){
    split(k,s,",");
    print s[1]"\t"s[2]"\t"sum[k]/count[k];
  }
}' GSM6505210_PA2_CapHiC_2ndBioRep_Chicago_output_washU_text.txt \
   GSM6505211_PA2_E10.5wt_CapHiC_3rdBioRep_Tr1_Chicago_output_washU_text.txt \
   > GSM6505210_5211_merged.txt
