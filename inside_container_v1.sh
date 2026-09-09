#!/bin/bash

PATH=$PATH:/srv/bin

cd /mnt

./populate_hpc.py

for num in $(seq -w 001 064); do
	sRNAbench input=/mnt/data/31064-${num}.R1.fastq.gz output=/mnt/QNS31064/results/${num} adapter=AACTGTAGGCACCATCAAT umi=3pA12 	graphics=true microRNA=hsa species=GRCh38_p13_mp libs=GRCh38_p13_ncRNA libs=GRCh38_p13_RNAcentral 	tRNA=GRCh38_p13_genomic_tRNA.fa libs=GRCh38_p13_cdna 'libsStringTypes=mature#sense;hairpin#sense|snRNA#sense|	snoRNA#sense|tRNA#sense|rRNA#sense|ncRNA#sense|protein_coding#sense;cdna#sense;cds#sense|	protein_coding#antisense;cdna#antisense;cds#antisense|piRNA#sense|repeat#sense' 'libsStringNames=miRBase(sense)|snRNA|snoRNA|	tRNA|rRNA|ncRNA|cdna(sense)|mRNA(antisense)|piRNA|repeats(sense)' isoMiR=true fullIsoStat=true predict=true

done

echo "concatenate columns seqVariants microRNAname_type_position"

for num in $(seq -w 001 064); do

	awk '{print $1"_"$2"_"$7"\t"$3"\t"$4"\t"$5"\t"$6}' < /scratch/p309221/micro_rna/QNS31064/results/${num}/microRNA_seqVariants.txt > 	/scratch/p309221/micro_rna/QNS31064/results/${num}/microRNA_seqVariants_concat.txt

done

echo "perform differential expression profiling"

echo "mTBI vs HC"

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBI_vs_HC grpString=001:002:003:004:005:006:007:008:009:010:011:012:013:014:015:016:017:018:019:020:021:022:023:024:025:026:027:028:029:030:034:036:037:038:042:046:047:048:049:050:051:052:053:054:055:056:057:058:059:060:061:062:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=mTBI#controls minRCexpr=5 isoSummary=true 

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBI_vs_HC grpString=001:002:003:004:005:006:007:008:009:010:011:012:013:014:015:016:017:018:019:020:021:022:023:024:025:026:027:028:029:030:034:036:037:038:042:046:047:048:049:050:051:052:053:054:055:056:057:058:059:060:061:062:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=mTBI#controls minRCexpr=5 diffExpr=true diffExprFiles=GRCh38_p13_genomic_tRNA_sense.grouped makeSingleAssignDE=false annotName=GRCh38_p13_genomic_tRNA

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBI_vs_HC grpString=001:002:003:004:005:006:007:008:009:010:011:012:013:014:015:016:017:018:019:020:021:022:023:024:025:026:027:028:029:030:034:036:037:038:042:046:047:048:049:050:051:052:053:054:055:056:057:058:059:060:061:062:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=mTBI#controls minRCexpr=5 diffExpr=true diffExprFiles=GRCh38_p13_ncRNA_sense.grouped makeSingleAssignDE=true annotName=GRCh38_p13_ncRNA

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBI_vs_HC grpString=001:002:003:004:005:006:007:008:009:010:011:012:013:014:015:016:017:018:019:020:021:022:023:024:025:026:027:028:029:030:034:036:037:038:042:046:047:048:049:050:051:052:053:054:055:056:057:058:059:060:061:062:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=mTBI#controls minRCexpr=5 diffExpr=true diffExprFiles=GRCh38_p13_RNAcentral_sense.grouped makeSingleAssignDE=true annotName=GRCh38_p13_RNAcentral

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBI_vs_HC grpString=001:002:003:004:005:006:007:008:009:010:011:012:013:014:015:016:017:018:019:020:021:022:023:024:025:026:027:028:029:030:034:036:037:038:042:046:047:048:049:050:051:052:053:054:055:056:057:058:059:060:061:062:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=mTBI#controls minRCexpr=5 diffExpr=true diffExprFiles=tRNA_mature_sense.grouped makeSingleAssignDE=true annotName=GRCh38_p13_genomic_tRNA

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBI_vs_HC grpString=001:002:003:004:005:006:007:008:009:010:011:012:013:014:015:016:017:018:019:020:021:022:023:024:025:026:027:028:029:030:034:036:037:038:042:046:047:048:049:050:051:052:053:054:055:056:057:058:059:060:061:062:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=mTBI#controls minRCexpr=5 diffExpr=false seqStat=false stat=true statFiles=mappingStat.txt minRCdata=0 folderData=stat 

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBI_vs_HC grpString=001:002:003:004:005:006:007:008:009:010:011:012:013:014:015:016:017:018:019:020:021:022:023:024:025:026:027:028:029:030:034:036:037:038:042:046:047:048:049:050:051:052:053:054:055:056:057:058:059:060:061:062:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=mTBI#controls minRCexpr=5 stat=true statFiles=isomiR_NTA.txt colData=3 folderData=stat

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBI_vs_HC grpString=001:002:003:004:005:006:007:008:009:010:011:012:013:014:015:016:017:018:019:020:021:022:023:024:025:026:027:028:029:030:034:036:037:038:042:046:047:048:049:050:051:052:053:054:055:056:057:058:059:060:061:062:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=mTBI#controls minRCexpr=5 diffExpr=true diffExprFiles=mature_sense.grouped makeSingleAssignDE=true annotName=mature

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBI_vs_HC grpString=001:002:003:004:005:006:007:008:009:010:011:012:013:014:015:016:017:018:019:020:021:022:023:024:025:026:027:028:029:030:034:036:037:038:042:046:047:048:049:050:051:052:053:054:055:056:057:058:059:060:061:062:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=mTBI#controls minRCexpr=5 stat=true statFiles=isomiR_otherVariants.txt colData=3 folderData=stat

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBI_vs_HC grpString=001:002:003:004:005:006:007:008:009:010:011:012:013:014:015:016:017:018:019:020:021:022:023:024:025:026:027:028:029:030:034:036:037:038:042:046:047:048:049:050:051:052:053:054:055:056:057:058:059:060:061:062:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=mTBI#controls minRCexpr=5 stat=true statFiles=microRNA_seqVariants_concat.txt colData=4

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBI_vs_HC grpString=001:002:003:004:005:006:007:008:009:010:011:012:013:014:015:016:017:018:019:020:021:022:023:024:025:026:027:028:029:030:034:036:037:038:042:046:047:048:049:050:051:052:053:054:055:056:057:058:059:060:061:062:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=mTBI#controls minRCexpr=5 stat=true statFiles=microRNA_seqVariants_concat.txt colData=1

echo "ctPos vs CTneg"

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_CTneg grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062 grpDesc=CTpos#CTneg minRCexpr=5 isoSummary=true 

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_CTneg grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062 grpDesc=CTpos#CTneg minRCexpr=5 diffExpr=true diffExprFiles=mature_sense.grouped makeSingleAssignDE=true annotName=mature

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_CTneg grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062 grpDesc=CTpos#CTneg minRCexpr=5 diffExpr=true diffExprFiles=GRCh38_p13_genomic_tRNA_sense.grouped makeSingleAssignDE=false annotName=GRCh38_p13_genomic_tRNA

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_CTneg grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062 grpDesc=CTpos#CTneg minRCexpr=5 diffExpr=true diffExprFiles=GRCh38_p13_ncRNA_sense.grouped makeSingleAssignDE=true annotName=GRCh38_p13_ncRNA

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_CTneg grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062 grpDesc=CTpos#CTneg minRCexpr=5 diffExpr=true diffExprFiles=GRCh38_p13_RNAcentral_sense.grouped makeSingleAssignDE=true annotName=GRCh38_p13_RNAcentral

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_CTneg grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062 grpDesc=CTpos#CTneg minRCexpr=5 diffExpr=true diffExprFiles=tRNA_mature_sense.grouped makeSingleAssignDE=true annotName=GRCh38_p13_genomic_tRNA

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_CTneg grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062 grpDesc=CTpos#CTneg minRCexpr=5 diffExpr=false seqStat=false stat=true statFiles=mappingStat.txt minRCdata=0 folderData=stat

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_CTneg grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062 grpDesc=CTpos#CTneg minRCexpr=5 stat=true statFiles=isomiR_NTA.txt colData=3 folderData=stat

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_CTneg grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062 grpDesc=CTpos#CTneg minRCexpr=5 stat=true statFiles=isomiR_otherVariants.txt colData=3 folderData=stat

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_CTneg grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062 grpDesc=CTpos#CTneg minRCexpr=5 stat=true statFiles=microRNA_seqVariants_concat.txt colData=4

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_CTneg grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062 grpDesc=CTpos#CTneg minRCexpr=5 stat=true statFiles=microRNA_seqVariants_concat.txt colData=1

echo "CTpos vs healthy controls"

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_HC grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=CTpos#HC minRCexpr=5 isoSummary=true

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_HC grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=CTpos#HC minRCexpr=5 diffExpr=true diffExprFiles=mature_sense.grouped makeSingleAssignDE=true annotName=mature

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_HC grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=CTpos#HC minRCexpr=5 diffExpr=true diffExprFiles=GRCh38_p13_genomic_tRNA_sense.grouped makeSingleAssignDE=false annotName=GRCh38_p13_genomic_tRNA

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_HC grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=CTpos#HC minRCexpr=5 diffExpr=true diffExprFiles=GRCh38_p13_ncRNA_sense.grouped makeSingleAssignDE=true annotName=GRCh38_p13_ncRNA

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_HC grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=CTpos#HC minRCexpr=5 diffExpr=true diffExprFiles=GRCh38_p13_RNAcentral_sense.grouped makeSingleAssignDE=true annotName=GRCh38_p13_RNAcentral

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_HC grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=CTpos#HC minRCexpr=5 diffExpr=true diffExprFiles=tRNA_mature_sense.grouped makeSingleAssignDE=true annotName=GRCh38_p13_genomic_tRNA

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_HC grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=CTpos#HC minRCexpr=5 diffExpr=false seqStat=false stat=true statFiles=mappingStat.txt minRCdata=0 folderData=stat

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_HC grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=CTpos#HC minRCexpr=5 stat=true statFiles=isomiR_NTA.txt colData=3 folderData=stat

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_HC grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=CTpos#HC minRCexpr=5 stat=true statFiles=isomiR_otherVariants.txt colData=3 folderData=stat

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_HC grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=CTpos#HC minRCexpr=5 stat=true statFiles=microRNA_seqVariants_concat.txt colData=4

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_HC grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=CTpos#HC minRCexpr=5 stat=true statFiles=microRNA_seqVariants_concat.txt colData=1

echo "CTneg vs healthy controls"

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTneg_vs_HC grpString=001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062#031:032:033:035:039:040:041:043:044:045 grpDesc=CTneg#HC minRCexpr=5 isoSummary=true

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTneg_vs_HC grpString=001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062#031:032:033:035:039:040:041:043:044:045 grpDesc=CTneg#HC minRCexpr=5 diffExpr=true diffExprFiles=mature_sense.grouped makeSingleAssignDE=true annotName=mature

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTneg_vs_HC grpString=001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062#031:032:033:035:039:040:041:043:044:045 grpDesc=CTneg#HC minRCexpr=5 diffExpr=true diffExprFiles=GRCh38_p13_genomic_tRNA_sense.grouped makeSingleAssignDE=false annotName=GRCh38_p13_genomic_tRNA

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTneg_vs_HC grpString=001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062#031:032:033:035:039:040:041:043:044:045 grpDesc=CTneg#HC minRCexpr=5 diffExpr=true diffExprFiles=GRCh38_p13_ncRNA_sense.grouped makeSingleAssignDE=true annotName=GRCh38_p13_ncRNA

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTneg_vs_HC grpString=001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062#031:032:033:035:039:040:041:043:044:045 grpDesc=CTneg#HC minRCexpr=5 diffExpr=true diffExprFiles=GRCh38_p13_RNAcentral_sense.grouped makeSingleAssignDE=true annotName=GRCh38_p13_RNAcentral

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTneg_vs_HC grpString=001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062#031:032:033:035:039:040:041:043:044:045 grpDesc=CTneg#HC minRCexpr=5 diffExpr=true diffExprFiles=tRNA_mature_sense.grouped makeSingleAssignDE=true annotName=GRCh38_p13_genomic_tRNA

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTneg_vs_HC grpString=001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062#031:032:033:035:039:040:041:043:044:045 grpDesc=CTneg#HC minRCexpr=5 diffExpr=false seqStat=false stat=true statFiles=mappingStat.txt minRCdata=0 folderData=stat

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTneg_vs_HC grpString=001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062#031:032:033:035:039:040:041:043:044:045 grpDesc=CTneg#HC minRCexpr=5 stat=true statFiles=isomiR_NTA.txt colData=3 folderData=stat

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTneg_vs_HC grpString=001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062#031:032:033:035:039:040:041:043:044:045 grpDesc=CTneg#HC minRCexpr=5 stat=true statFiles=isomiR_otherVariants.txt colData=3 folderData=stat

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTneg_vs_HC grpString=001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062#031:032:033:035:039:040:041:043:044:045 grpDesc=CTneg#HC minRCexpr=5 stat=true statFiles=microRNA_seqVariants_concat.txt colData=4

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTneg_vs_HC grpString=001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062#031:032:033:035:039:040:041:043:044:045 grpDesc=CTneg#HC minRCexpr=5 stat=true statFiles=microRNA_seqVariants_concat.txt colData=1

echo "mTBI males vs females"

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBImale_vs_mTBIfemale grpString=003:007:011:012:014:016:017:019:020:021:022:023:025:026:027:036:037:038:047:048:052:053:057:058:059:062:064#001:002:004:005:006:008:009:010:013:015:018:024:028:029:030:034:042:046:049:050:051:054:055:056:060:061:063 grpDesc=male#female minRCexpr=5 isoSummary=true

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBImale_vs_mTBIfemale grpString=003:007:011:012:014:016:017:019:020:021:022:023:025:026:027:036:037:038:047:048:052:053:057:058:059:062:064#001:002:004:005:006:008:009:010:013:015:018:024:028:029:030:034:042:046:049:050:051:054:055:056:060:061:063 grpDesc=male#female minRCexpr=5 diffExpr=true diffExprFiles=mature_sense.grouped makeSingleAssignDE=true annotName=mature

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBImale_vs_mTBIfemale grpString=003:007:011:012:014:016:017:019:020:021:022:023:025:026:027:036:037:038:047:048:052:053:057:058:059:062:064#001:002:004:005:006:008:009:010:013:015:018:024:028:029:030:034:042:046:049:050:051:054:055:056:060:061:063 grpDesc=male#female minRCexpr=5 diffExpr=true diffExprFiles=GRCh38_p13_genomic_tRNA_sense.grouped makeSingleAssignDE=false annotName=GRCh38_p13_genomic_tRNA

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBImale_vs_mTBIfemale grpString=003:007:011:012:014:016:017:019:020:021:022:023:025:026:027:036:037:038:047:048:052:053:057:058:059:062:064#001:002:004:005:006:008:009:010:013:015:018:024:028:029:030:034:042:046:049:050:051:054:055:056:060:061:063 grpDesc=male#female minRCexpr=5 diffExpr=true diffExprFiles=GRCh38_p13_ncRNA_sense.grouped makeSingleAssignDE=true annotName=GRCh38_p13_ncRNA

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBImale_vs_mTBIfemale grpString=003:007:011:012:014:016:017:019:020:021:022:023:025:026:027:036:037:038:047:048:052:053:057:058:059:062:064#001:002:004:005:006:008:009:010:013:015:018:024:028:029:030:034:042:046:049:050:051:054:055:056:060:061:063 grpDesc=male#female minRCexpr=5 diffExpr=true diffExprFiles=GRCh38_p13_RNAcentral_sense.grouped makeSingleAssignDE=true annotName=GRCh38_p13_RNAcentral

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBImale_vs_mTBIfemale grpString=003:007:011:012:014:016:017:019:020:021:022:023:025:026:027:036:037:038:047:048:052:053:057:058:059:062:064#001:002:004:005:006:008:009:010:013:015:018:024:028:029:030:034:042:046:049:050:051:054:055:056:060:061:063 grpDesc=male#female minRCexpr=5 diffExpr=true diffExprFiles=tRNA_mature_sense.grouped makeSingleAssignDE=true annotName=GRCh38_p13_genomic_tRNA

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBImale_vs_mTBIfemale grpString=003:007:011:012:014:016:017:019:020:021:022:023:025:026:027:036:037:038:047:048:052:053:057:058:059:062:064#001:002:004:005:006:008:009:010:013:015:018:024:028:029:030:034:042:046:049:050:051:054:055:056:060:061:063 grpDesc=male#female minRCexpr=5 diffExpr=false seqStat=false stat=true statFiles=mappingStat.txt minRCdata=0 folderData=stat

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBImale_vs_mTBIfemale grpString=003:007:011:012:014:016:017:019:020:021:022:023:025:026:027:036:037:038:047:048:052:053:057:058:059:062:064#001:002:004:005:006:008:009:010:013:015:018:024:028:029:030:034:042:046:049:050:051:054:055:056:060:061:063 grpDesc=male#female minRCexpr=5 stat=true statFiles=isomiR_NTA.txt colData=3 folderData=stat

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBImale_vs_mTBIfemale grpString=003:007:011:012:014:016:017:019:020:021:022:023:025:026:027:036:037:038:047:048:052:053:057:058:059:062:064#001:002:004:005:006:008:009:010:013:015:018:024:028:029:030:034:042:046:049:050:051:054:055:056:060:061:063 grpDesc=male#female minRCexpr=5 stat=true statFiles=isomiR_otherVariants.txt colData=3 folderData=stat

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBImale_vs_mTBIfemale grpString=003:007:011:012:014:016:017:019:020:021:022:023:025:026:027:036:037:038:047:048:052:053:057:058:059:062:064#001:002:004:005:006:008:009:010:013:015:018:024:028:029:030:034:042:046:049:050:051:054:055:056:060:061:063 grpDesc=male#female minRCexpr=5 stat=true statFiles=microRNA_seqVariants_concat.txt colData=4

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBImale_vs_mTBIfemale grpString=003:007:011:012:014:016:017:019:020:021:022:023:025:026:027:036:037:038:047:048:052:053:057:058:059:062:064#001:002:004:005:006:008:009:010:013:015:018:024:028:029:030:034:042:046:049:050:051:054:055:056:060:061:063 grpDesc=male#female minRCexpr=5 stat=true statFiles=microRNA_seqVariants_concat.txt colData=1

echo "seed sequences

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBI_vs_HC grpString=001:002:003:004:005:006:007:008:009:010:011:012:013:014:015:016:017:018:019:020:021:022:023:024:025:026:027:028:029:030:034:036:037:038:042:046:047:048:049:050:051:052:053:054:055:056:057:058:059:060:061:062:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=mTBI#controls minRCexpr=5 stat=true statFiles=seedExpression.tsv colData=4

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_CTneg grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062 grpDesc=CTpos#CTneg minRCexpr=5 stat=true statFiles=seedExpression.tsv colData=4

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTpos_vs_HC grpString=004:005:008:011:019:022:023:024:027:029:034:036:037:038:042:046:048:049:051:052:054:055:057:059:060:063:064#031:032:033:035:039:040:041:043:044:045 grpDesc=CTpos#HC minRCexpr=5 stat=true statFiles=seedExpression.tsv colData=4

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/CTneg_vs_HC grpString=001:002:003:006:007:009:010:012:013:014:015:016:017:018:020:021:025:026:028:030:047:050:053:056:058:061:062#031:032:033:035:039:040:041:043:044:045 grpDesc=CTneg#HC minRCexpr=5 stat=true statFiles=seedExpression.tsv colData=4

sRNAde input=/mnt/QNS31064/results output=/mnt/QNS31064/results/mTBImale_vs_mTBIfemale grpString=003:007:011:012:014:016:017:019:020:021:022:023:025:026:027:036:037:038:047:048:052:053:057:058:059:062:064#001:002:004:005:006:008:009:010:013:015:018:024:028:029:030:034:042:046:049:050:051:054:055:056:060:061:063 grpDesc=male#female minRCexpr=5 stat=true statFiles=seedExpression.tsv colData=4

exit
