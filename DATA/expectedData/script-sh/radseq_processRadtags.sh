rm -rf ../radseq/outputRadseq/
mkdir ../radseq/outputRadseq/

/usr/local/stacks-1.29/bin/process_radtags -p ../radseq/radseqInput/lane5/ -o ../radseq/outputRadseq/ -b ../radseq/radseqInput/barcode/barcode_run01_lane5.txt -e apeKI -i fastq

mv ../radseq/outputRadseq/process_radtags.log ../radseq/outputRadseq/process_radtags_lane5.log

/usr/local/stacks-1.29/bin/process_radtags -p ../radseq/radseqInput/lane6/ -o ../radseq/outputRadseq/ -b ../radseq/radseqInput/barcode/barcode_run01_lane6.txt -e apeKI -i fastq

mv ../radseq/outputRadseq/process_radtags.log ../radseq/outputRadseq/process_radtags_lane6.log
