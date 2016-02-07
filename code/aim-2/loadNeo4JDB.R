###########################################
#Code to load Neo4J
#
#
###########################################

#Install Neo4J first and make sure you have JDK8
#http://tecadmin.net/install-java-8-on-centos-rhel-and-fedora/

#using package RNeo4J
#https://github.com/nicolewhite/RNeo4j
#Please note when installing you may need to install openssldevel first  : yum install -y openssl-devel
#Also make sure to disable authentication in conf/neo4j-server : dbms.security.auth_enabled=false


#Call library
library(RNeo4j)

#Start a graph & if its full get rid of 
graph = startGraph("http://localhost:7474/db/data/")
clear(graph, F)


#First let's start parse genemania data, there are 553 files
setwd("/home/ramanp/pitfit/data/GeneMania/genemania.org/data/current/Homo_sapiens");
allFiles <-list.files();















