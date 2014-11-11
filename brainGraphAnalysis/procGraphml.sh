case  in
    
#takes a file with names of zip graphml files as input
#unzips the graphml file, runs R script and deletes the unzipped file
#
#$1: file with list of zip file names
#$2: output directory

while read zipFile
do
    graphFile=$2$(echo $zipFile | awk -F "/" '{print $NF}' \
        | awk -F "." '{OFS=FS}{NF=NF-1; print}')
    outputFile=$2$(echo $zipFile | awk -F "/" '{print $NF}' \
        | awk -F "." '{OFS=FS}{NF=NF-2; print}')"_attr.txt"
    unzip -o $zipFile -d $2
    R CMD BATCH "--args $graphFile $outputFile" procGraphml.R
    rm -f $graphFile
done < $1
