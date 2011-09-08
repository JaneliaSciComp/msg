#!/bin/sh

# Download sample data for MSG
#wget "https://github.com/tinathu/msg/raw/master/example/msg.cfg" --no-check-certificate -O msg.cfg
#wget "https://github.com/tinathu/msg/raw/master/example/barcodes_file" --no-check-certificate -O barcodes_file.txt

FILES="ftp://anonymous@ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX029/SRX029935/SRR071201/SRR071201.sra ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX032/SRX032362/SRR074287/SRR074287.sra ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX032/SRX032363/SRR074288/SRR074288.sra ftp://anonymous@ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX029/SRX029935/SRR071201/SRR071201.sra"

for F in $FILES
do
	wget "$F" -N
done
