
#!/bin/bash
# This script is to set up the required softwares for RNA-seq analysis.
# The script is intended as a wrapper to install all required software packages locally without external dependencies.
# It is important to check error messages (if any) to make sure that all packages required to run this pipeline are installed correctly.
 
 
set -e 
set -u
set -o pipefail


module load apps/gcc/4.9.0
module load apps/python/2.7.13
module load apps/java/8u112
module load apps/parallel/20180322 
module load apps/biopython/1.56
module load apps/bioperl/1.6.924
module load apps/samtools
module load apps/cmake
echo "PICARD_PATH=/gne/research/data/dnaseq/analysis/aplle/software/picard-tools-1.126/picard.jar" >> ../paths.txt

mkdir tools
cd tools

## install bam-readcount
if [ -z ${BAM_READCOUNT_PATH} ]
then 
  git clone https://github.com/genome/bam-readcount.git
  cd bam-readcount
  cmake -Wno-dev $RNA_HOME/tools/bam-readcount
  make
  cd ..
  export BAM_READCOUNT=$('pwd')bam-readcount/bin
  sed -i '/BAM_READCOUNT_PATH/d' ../../paths.txt
  echo "BAM_READCOUNT_PATH=${BAM_READCOUNT}" >> ../../paths.txt
fi


## install HiSAT2
if [ -z ${HISAT2_PATH} ]
then
  wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip
  unzip hisat2-2.1.0-Linux_x86_64.zip
  export HISAT2=$('pwd')/hisat2-2.1.0
  sed -i '/HISAT2_PATH/d' ../../paths.txt
  echo "HISAT2_PATH=${HISAT}" >> ../../paths.txt
fi


## install StringTie
if [ -z ${STRINGTIE_PATH} ]
then
  wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.4b.Linux_x86_64.tar.gz
  tar -xzf stringtie-1.3.4b.Linux_x86_64.tar.gz
  rm stringtie-1.3.4b.Linux_x86_64.tar.gz
  export STRINGTIE=$('pwd')/stringtie-1.3.4b.Linux_x86_64
  sed -i '/STRINGTIE_PATH/d' ../../paths.txt
  echo "STRINGTIE_PATH=${STRINGTIE}" >> ../../paths.txt
fi  


## install HTSEQ
if [ -z ${HTSEQ_PATH} ] 
then
  wget https://github.com/simon-anders/htseq/archive/release_0.9.1.tar.gz
  tar -zxf release_0.9.1.tar.gz
  rm release_0.9.1.tar.gz
  cd htseq-release_0.9.1/
  python setup.py install --user
  chmod +x scripts/htseq-count
  cd ..
  export HTSEQ=$('pwd')/scripts
  sed -i '/HTSEQ_PATH/d' ../../paths.txt
  echo "HTSEQ_PATH=${HTSEQ}" >> ../../paths.txt
fi


## install TopHat
if [ -z ${TOPHAT_PATH} ]
then
  wget https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz
  tar -zxf tophat-2.1.1.Linux_x86_64.tar.gz
  rm tophat-2.1.1.Linux_x86_64.tar.gz
  export TOPHAT=$('pwd')/tophat-2.1.1.Linux_x86_64
  sed -i '/TOPHAT_PATH/d' ../../paths.txt
  echo "TOPHAT_PATH=${TOPHAT}" >> ../../paths.txt
fi


## install Kallisto
if [ -z ${KALLISTO_PATH} ]
then
  wget https://github.com/pachterlab/kallisto/releases/download/v0.44.0/kallisto_linux-v0.44.0.tar.gz
  tar -zxf kallisto_linux-v0.44.0.tar.gz
  rm kallisto_linux-v0.44.0.tar.gz
  export KALLISTO=$('pwd')/kallisto_linux-v0.44.0
  sed -i '/KALLISTO_PATH/d' ../../paths.txt
  echo "KALLISTO_PATH=${KALLISTO}" >> ../../paths.txt
fi
