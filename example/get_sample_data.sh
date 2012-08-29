#!/bin/sh

# Download sample data for MSG

FILES="ftp://anonymous@ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX029/SRX029935/SRR071201/SRR071201.sra ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX032/SRX032362/SRR074287/SRR074287.sra ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX032/SRX032363/SRR074288/SRR074288.sra ftp://anonymous@ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX029/SRX029935/SRR071201/SRR071201.sra"

for F in $FILES
do
	wget "$F" -N
done
