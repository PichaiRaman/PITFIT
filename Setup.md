This server was setup on the AWS cloud with the following config [CentOS 6.5 x86_64 SequenceIQ (upgraded) (ami-48e8cc20)].
The server has 7.5 GB of Ram and 20 gigs of space. 


Have to download R and Shiny Server

```
sudo su -c 'rpm -Uvh http://download.fedoraproject.org/pub/epel/6/i386/epel-release-6-8.noarch.rpm'
sudo yum update
sudo yum install R
sudo su - \
-c "R -e \"install.packages('shiny', repos='http://cran.rstudio.com/')\""

wget http://download3.rstudio.org/centos-5.9/x86_64/shiny-server-1.3.0.403-rh5-x86_64.rpm
sudo yum install --nogpgcheck shiny-server-1.3.0.403-rh5-x86_64.rpm
```
After this it also doesn't hurt to install rmarkdown, so go into R and 

```
install.packages("rmarkdown");
```
I chose repo 92 but feel free to change that up, go nuts. 

also should install git

```
sudo yum install git
```
Make sure to configure it, i.e. user.name, user.email, and add ssh key. Also get twitter bootstrap.
```
wget https://github.com/twbs/bootstrap/releases/download/v3.3.4/bootstrap-3.3.4-dist.zip
```

Also needed to download broad firehose to get TCGA data

```
https://confluence.broadinstitute.org/display/GDAC/Download
```





