INPUT_CSV=$1
HEADER="Chr	Start	End	Ref	Alt	Zygozity	Func.refGene	Gene.refGene	ExonicFunc.refGene	AAChange.refGene	EUR_1000g2012apr	AMR_1000g2012apr	ASN_1000g2012apr	AFR_1000g2012apr	dbSNP_138	omim_phenotype	avsift_score	clinvar_20140211	ljb23_pp2hvar	CADDgt10	vcf_QUAL	vcf_INFO	GT:AD:DP:GQ:PL"

echo "$HEADER"

awk '{
    tail=$0
    do {
        match(tail,/\"[^\"]*\"/)
        if (RSTART <= 0) {
                gsub(/,/,"\t",tail)
                printf "%s\n",tail
        } else {
                head = substr(tail,1,RSTART-1)
                gsub(/,/,"\t",head)
                printf "%s%s",head,substr(tail,RSTART,RLENGTH)
                tail = substr(tail,RSTART+RLENGTH)
        }
    } while (RSTART > 0)
}' $INPUT_CSV | sed 's/"//g' | tail -n +2 | while read line
do
	one2five=`echo "$line" | cut -d"	" -f1-5`
	six2nineteen=`echo "$line" | cut -d"	" -f6-19`
	zygozity=`echo "$line" | cut -d"	" -f20`
	ending=`echo "$line" | cut -d"	" -f28,30,32`
	echo "$one2five	$zygozity	$six2nineteen	$ending"
done

