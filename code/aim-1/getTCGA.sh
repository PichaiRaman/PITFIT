#!/bin/bash

cd /srv/shiny-server/PITFIT/data
firehose_get -b -tasks segmented_scna_minus_germline_cnv_hg19__seg.Level_3 Mutation_Packager_Oncotated_Calls.Level_3 Merge_Clinical.Level_1 rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3 stddata latest OV PAAD PRAD LIHC READ

