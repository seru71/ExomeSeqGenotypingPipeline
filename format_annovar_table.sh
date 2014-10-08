INPUT_CSV=$1
HEADER="Gene	Zygozity	Chr	Start	End	Ref	Alt	Func.refGene	Gene.refGene	ExonicFunc.refGene	AAChange.refGene	EUR_1000g2012apr	AMR_1000g2012apr	ASN_1000g2012apr	AFR_1000g2012apr	dbSNP_138	avsift_score	clinvar_20140211	ljb23_pp2hvar	CADDgt10	vcf_QUAL	vcf_INFO	GT:AD:DP:GQ:PL"

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
	one2five=`echo "$line" | cut -f1-5`
	gene=`echo "$line" | cut -f7`
	six2nineteen=`echo "$line" | cut -f6-18`
	zygozity=`echo "$line" | cut -f19`
	ending=`echo "$line" | cut -f27,29,31`
	echo "$gene	$zygozity	$one2five	$six2nineteen	$ending"
done

