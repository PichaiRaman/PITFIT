#!/bin/bash

cd /srv/shiny-server/PITFIT/data
echo "changed directory"
firehose_get -b -tasks segmented_scna_minus_germline_cnv_hg19__seg.Level_3 Mutation_Packager_Oncotated_Calls.Level_3 Merge_Clinical.Level_1 rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3 stddata latest OV PAAD PRAD LIHC READ

#Expand all tar gz
cd */OV/*
for a in `ls -1 *.tar.gz`; do tar -zxvf $a; done
rm *.md5
rm *.gz
cd /srv/shiny-server/PITFIT/data
cd */PAAD/*
for a in `ls -1 *.tar.gz`; do tar -zxvf $a; done
rm *.md5
rm *.gz
cd /srv/shiny-server/PITFIT/data
cd */PRAD/*
for a in `ls -1 *.tar.gz`; do tar -zxvf $a; done
rm *.md5
rm *.gz
cd /srv/shiny-server/PITFIT/data
cd */LIHC/*
for a in `ls -1 *.tar.gz`; do tar -zxvf $a; done
rm *.md5
rm *.gz
cd /srv/shiny-server/PITFIT/data
cd */READ/*
for a in `ls -1 *.tar.gz`; do tar -zxvf $a; done
rm *.md5
rm *.gz




