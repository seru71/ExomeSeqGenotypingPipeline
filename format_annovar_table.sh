INPUT_CSV=$1
#HEADER="Chr	Start	End	Ref	Alt	Zygozity	Func.refGene	Gene.refGene	ExonicFunc.refGene	AAChange.refGene	EUR_1000g2012apr	AMR_1000g2012apr	ASN_1000g2012apr	AFR_1000g2012apr	dbSNP_138	omim_phenotype	avsift_score	clinvar_20140211	ljb23_pp2hvar	CADDgt10	vcf_QUAL	vcf_INFO	GT:AD:DP:GQ:PL"
HEADER="Chr	Start	End	Ref	Alt	Zygozity	Func.refGene	Gene.refGene	GeneDetail.refGene	ExonicFunc.refGene	AAChange.refGene	EUR_1000g2014oct	AMR_1000g2014oct	EAS_1000g2014oct	SAS_1000g2014oct	AFR_1000g2014oct	ExAC_ALL	ExAC_AFR	ExAC_AMR	ExAC_EAS	ExAC_FIN	ExAC_NFE	ExAC_OTH	ExAC_SAS	snp138	omim_phenotype	clinvar_20150330	SIFT_score	SIFT_pred	Polyphen2_HDIV_score	Polyphen2_HDIV_pred	Polyphen2_HVAR_score	Polyphen2_HVAR_pred	LRT_score	LRT_pred	MutationTaster_score	MutationTaster_pred	MutationAssessor_score	MutationAssessor_pred	FATHMM_score	FATHMM_pred	RadialSVM_score	RadialSVM_pred	LR_score	LR_pred	VEST3_score	CADD_raw	CADD_phred	GERP++_RS	phyloP46way_placental	phyloP100way_vertebrate	SiPhy_29way_logOdds	vcf_QUAL	vcf_FILTER	vcf_INFO	GT:AD:DP:GQ:PL"

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
	variantLocation=`echo "$line" | cut -d"	" -f1-5`
	dbAnnotations1=`echo "$line" | cut -d"	" -f6-14,16-25`
	omim=`echo "$line" | cut -d"	" -f15`
	dbAnnotations2=`echo "$line" | cut -d"	" -f26-51`
	zygozity=`echo "$line" | cut -d"	" -f52`
	callInfo=`echo "$line" | cut -d"	" -f60-62,64`
	echo "$variantLocation	$zygozity	$dbAnnotations1	$omim	$dbAnnotations2	$callInfo"
done

