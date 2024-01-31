#!/bin/sh
#SBATCH -J pGoM.sh
#SBATCH -o /home/%u/MosaiC-All/Log/pGoM.slurm-%j.out
#SBATCH -p skylake,icelake
#SBATCH -N 1
#SBATCH -n 5
#SBATCH --time=12:30:00
#SBATCH --mem=30GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=%u@adelaide.edu.au

## Hard coded paths you may wish to change for your system
LOGDIR="/home/${USER}/MosaiC-All/Log"

usage()
{
echo "
#
#README
#This script is designed for MosaiC-All  pipeline
#It requires family VCF in trios
#10-father
#11-mother
#12-child (proband/sib)
#
#SCRIPT EXECUTION
#
#sbatch $0 -v /path/to/familytrio.vcf -s /path/to/ID.list.txt -o /path/to/output/
#
# Options
#-v 	<file>		REQUIRED: familyVCF
#-s 	<file>      	REQUIRED: A file e.g. sampleID.list (one header row and then tab-delimited columns \$BAMdir,\$ProbandID,\$Gender,\$Mother,\$Father)
#-o 	<directory> 	REQUIRED: output directory
"
}

# Load modules
module purge
module use /apps/Load modules/all
module load BCFtools/1.17-GCC-11.2.0

## Set Variables ##
while [ "$1" != "" ]; do
        case $1 in
                -v )                    shift
                                        VCFDIR=$1
                                        ;;
                -s )                    shift
                                        SAMPLELIST=$1
                                        ;;
                -o )                    shift
                                        OUTDIR=$1
                                        ;;
                * )                     usage
                                        exit 1
        esac
        shift
done


#Check if all information are provided

if [ -z "$VCFDIR" ]; then
    usage
    echo "## ERROR: You need to provide the directory where we can find all the FamilyVCFs listed in the $SAMPLELIST" > $LOGDIR/pGoM.log
    exit 1
fi

if [ -z "$SAMPLELIST" ]; then
    usage
    echo "## ERROR: You need to provide a sample list" 
    echo "#-s REQUIRED sampleID.list (one header row and then tab-delimited columns \$BAMdir,\$ProbandID,\$Gender,\$Mother,\$Father,\$FamilyVCFID)" >> $LOGDIR/pGoM.log
	exit 1
fi

if [ ! -d "${OUTDIR}" ]; then
    mkdir -p ${OUTDIR}
	echo "## INFO: output directory created, you'll find all of the outputs and log files in here: ${OUTDIR}" >> $LOGDIR/pGoM.log
fi

#Array from a list of Samples (ignoring the header of the file)
mapfile -t SAMPLEID < <(tail -n +2 "$SAMPLELIST")

# Iteration starts here
for SAMPLEID in "${SAMPLEID[@]}"; do

    #Defining variables from each row

    	ProbandID=$(awk '{print $2}' <<< "$SAMPLEID ")
   	MotherID=$(awk '{print $4}' <<< "$SAMPLEID ")
   	FatherID=$(awk '{print $5}' <<< "$SAMPLEID ")
    	FamilyVCF=$(awk '{print $6}' <<< "$SAMPLEID ")
     	VCF=$VCFDIR/$FamilyVCF

     	echo "pGoM Pipeline for $ProbandID, $MotherID, $FatherID based on familyvcf:$VCF" >> $OUTDIR/$ProbandID.pipeline.log

     	#Check if the VCF is available in the specified directory or else dont proceed
       	if [ -s "$VCF" ]; then
    	echo "$VCF exists and is not empty." >> $LOGDIR/$ProbandID.pGoM.pipeline.log
	else
    	echo "$VCF DOES NOT exist or is empty." >> $LOGDIR/$ProbandID.pGoM.pipeline.log
	fi

	# Find column numbers for sample IDs and store in variables
	col_P=$(grep -m 1 "^#CHROM" "$VCF" | awk -v proband="$ProbandID" '{for (i=10; i<=NF; i++) {if ($i == proband) print i}}')
	col_M=$(grep -m 1 "^#CHROM" "$VCF" | awk -v mother="$MotherID" '{for (i=10; i<=NF; i++) {if ($i == mother) print i}}')
	col_F=$(grep -m 1 "^#CHROM" "$VCF" | awk -v father="$FatherID" '{for (i=10; i<=NF; i++) {if ($i == father) print i}}')

	#Filtering for father
 	
  	## Filter 1: All inherited variants from the father
 	bcftools view "$VCF" | grep "^#" > $OUTDIR/$ProbandID.fa.all.inherited.vcf
 	bcftools view "$VCF" | awk -v col_P="$col_P" -v col_M="$col_M" -v col_F="$col_F" 'BEGIN {OFS="\t"} col_F ~ /^0\/1/ && col_P ~ /^0\/1/ && col_M ~ /^0\/0/ {print}' >> $OUTDIR/$ProbandID.fa.all.inherited.vcf

	# Filter 2: DP20 of inherited variants for all sample 
	bcftools view -O v -f PASS -i 'FORMAT/DP>=20' "$VCF" > $OUTDIR/$ProbandID.fa.DP20.vcf
			
   	# Filter 3: aaf filter in father
	bcftools view -O v -f PASS -i '(AD[0]/FORMAT/DP[0])<=0.4 | (AD[0]/FORMAT/DP[0])>=0.7' --samples $FatherID $OUTDIR/$ProbandID.fa.DP20.vcf > $OUTDIR/$ProbandID.fa.mosaic.DP20.vcf

	# Filter 4: aaf filter children is within 0.4-0.7
	bcftools view -O v -f PASS -i '(AD[0]/FORMAT/DP[0])>=0.4 && (AD[0]/FORMAT/DP[0])<=0.7' --samples "$ProbandID" $OUTDIR/$ProbandID.fa.mosaic.DP20.vcf > $OUTDIR/$ProbandID.fa.final.vcf

	#COUNTS
	V=$( bcftools view $VCF | grep -v "^#" | wc -l)
	A=$(grep -v "^#" $OUTDIR/$ProbandID.fa.all.inherited.vcf | wc -l)
	B=$(grep -v "^#" $OUTDIR/$ProbandID.fa.DP20.vcf| wc -l)
	C=$(grep -v "^#" $OUTDIR/$ProbandID.fa.mosaic.DP20.vcf | wc -l)
	D=$(grep -v "^#" $OUTDIR/$ProbandID.fa.final.vcf  | wc -l)
					
	echo -e "$ProbandID\tfather\tTrios\t$V\t$A\t$B\t$C\t$D" | tr " " "\t" >> $OUTDIR/mosaic.father.trios.counts.txt

	#Filtering for mother
 	
  	## Filter 1: All inherited variants from the mother
 	bcftools view "$VCF" | grep "^#" > $OUTDIR/$ProbandID.mo.all.inherited.vcf
 	bcftools view "$VCF" | awk -v col_P="$col_P" -v col_M="$col_M" -v col_F="$col_F" 'BEGIN {OFS="\t"} col_M ~ /^0\/1/ && col_P ~ /^0\/1/ && col_F ~ /^0\/0/ {print}' >> $OUTDIR/$ProbandID.mo.all.inherited.vcf

	# Filter 2: DP20 of inherited variants for all sample 
	bcftools view -O v -f PASS -i 'FORMAT/DP>=20' "$VCF" > $OUTDIR/$ProbandID.mo.DP20.vcf
			
   	# Filter 3: aaf filter in father
	bcftools view -O v -f PASS -i '(AD[0]/FORMAT/DP[0])<=0.4 | (AD[0]/FORMAT/DP[0])>=0.7' --samples $MotherID $OUTDIR/$ProbandID.mo.DP20.vcf > $OUTDIR/$ProbandID.mo.mosaic.DP20.vcf

	# Filter 4: aaf filter children is within 0.4-0.7
	bcftools view -O v -f PASS -i '(AD[0]/FORMAT/DP[0])>=0.4 && (AD[0]/FORMAT/DP[0])<=0.7' --samples "$ProbandID" $OUTDIR/$ProbandID.mo.mosaic.DP20.vcf > $OUTDIR/$ProbandID.mo.final.vcf

	#COUNTS
	V=$( bcftools view $VCF | grep -v "^#" | wc -l)
	A=$(grep -v "^#" $OUTDIR/$ProbandID.mo.all.inherited.vcf | wc -l)
	B=$(grep -v "^#" $OUTDIR/$ProbandID.mo.DP20.vcf| wc -l)
	C=$(grep -v "^#" $OUTDIR/$ProbandID.mo.mosaic.DP20.vcf | wc -l)
	D=$(grep -v "^#" $OUTDIR/$ProbandID.mo.final.vcf  | wc -l)
					
	echo -e "$ProbandID\tfather\tTrios\t$V\t$A\t$B\t$C\t$D" | tr " " "\t" >> $OUTDIR/mosaic.mother.trios.counts.txt

 	echo "pGoM Pipeline for $ProbandID is completed" >> $LOGDIR/$ProbandID.pipeline.log

done


						
	


