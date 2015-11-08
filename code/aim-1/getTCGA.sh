#!/bin/bash

cd /srv/shiny-server/PITFIT/data
echo "changed directory"
firehose_get -b -tasks segmented_scna_minus_germline_cnv_hg19__seg.Level_3 Mutation_Packager_Oncotated_Calls.Level_3 Merge_Clinical.Level_1 rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3 stddata latest OV PAAD PRAD LIHC READ

#Expand all tar gz
cd z
cd OV
cd *
for a in `ls -1 *.tar.gz`; do tar -zxvf $a; done
rm *.md5
rm *.gz
cd /srv/shiny-server/PITFIT/data
cd *
cd PAAD
cd *
for a in `ls -1 *.tar.gz`; do tar -zxvf $a; done
rm *.md5
rm *.gz
cd /srv/shiny-server/PITFIT/data
cd *
cd PRAD
cd *
for a in `ls -1 *.tar.gz`; do tar -zxvf $a; done
rm *.md5
rm *.gz
cd /srv/shiny-server/PITFIT/data
cd *
cd LIHC
cd *
for a in `ls -1 *.tar.gz`; do tar -zxvf $a; done
rm *.md5
rm *.gz
cd /srv/shiny-server/PITFIT/data
cd *
cd READ
cd *
for a in `ls -1 *.tar.gz`; do tar -zxvf $a; done
rm *.md5
rm *.gz




