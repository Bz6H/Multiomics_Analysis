#!/bin/bash
# Usage: ./jaspar2homer.sh input.jaspar output.motif

awk 'BEGIN{OFS="\t"}
  /^>/ { name=$2; next }
  /^A/ { gsub(/[\[\]A]/,""); split($0,a) }
  /^C/ { gsub(/[\[\]C]/,""); split($0,c) }
  /^G/ { gsub(/[\[\]G]/,""); split($0,g) }
  /^T/ { gsub(/[\[\]T]/,""); split($0,t); 
    n=length(a); consensus="";
    for(i=1;i<=n;i++){
      max=a[i]; base="A";
      if(c[i]>max){max=c[i]; base="C"}
      if(g[i]>max){max=g[i]; base="G"}
      if(t[i]>max){max=t[i]; base="T"}
      consensus = consensus base;
    }
    printf ">%s\t%s\t6.0\n", consensus, name;
    for(i=1;i<=n;i++){
      tot=a[i]+c[i]+g[i]+t[i];
      printf "%.3f\t%.3f\t%.3f\t%.3f\n", a[i]/tot, c[i]/tot, g[i]/tot, t[i]/tot;
    }
  }' "$1" > "$2"