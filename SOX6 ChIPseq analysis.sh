###########################################################
# DOWNLOADING DATA and MAPPING QC
###########################################################

nohup wget -O unfiltered_control.bam https://www.encodeproject.org/files/ENCFF063YKZ/@@download/ENCFF063YKZ.bam &
nohup wget -O unfiltered_rep1.bam https://www.encodeproject.org/files/ENCFF756IFT/@@download/ENCFF756IFT.bam &
nohup wget -O unfiltered_rep2.bam https://www.encodeproject.org/files/ENCFF362APT/@@download/ENCFF362APT.bam &

top

nohup samtools flagstat unfiltered_control.bam > control_flagstat.txt &
nohup samtools flagstat unfiltered_rep1.bam > rep1_flagstat.txt &
nohup samtools flagstat unfiltered_rep2.bam > rep2_flagstat.txt &


nohup samtools view -bq 1 unfiltered_control.bam > filtered_unique_control.bam &
nohup samtools view -bq 1 unfiltered_rep1.bam > filtered_unique_rep1.bam &
nohup samtools view -bq 1 unfiltered_rep2.bam > filtered_unique_rep2.bam &

top

nohup samtools flagstat filtered_unique_control.bam >> control_flagstat.txt &
nohup samtools flagstat filtered_unique_rep1.bam >> rep1_flagstat.txt &
nohup samtools flagstat filtered_unique_rep2.bam >> rep2_flagstat.txt &

###########################################################
# CREATING ALL SCENARIOS FOR PEAK HANDLING
###########################################################

nohup macs2 callpeak -t filtered_unique_rep1.bam -c filtered_unique_control.bam -g hs -n rep1 -q 0.01 &
nohup macs2 callpeak -t filtered_unique_rep2.bam -c filtered_unique_control.bam -g hs -n rep2 -q 0.01 &
nohup macs2 callpeak -t filtered_unique_rep1.bam filtered_unique_rep2.bam -c filtered_unique_control.bam -g hs -n merged -q 0.01 &

nohup R < rep1_model.r --vanilla &
nohup R < rep2_model.r --vanilla &
nohup R < merged_model.r --vanilla &

bedtools intersect -wa -u -a rep1_peaks.narrowPeak -b rep2_peaks.narrowPeak > intersect_peaks.narrowPeak
awk '{OFS="\t"; print $1, $2+$10, $2+$10+1,"encode", $9}' intersect_peaks.narrowPeak > intersect_summits.bed          # maybe I should change the encode label since irrelevant
bedtools closest -a intersect_summits.bed -b rep2_summits.bed -d -t first > intersect_rep2_closest.bed
awk '{ if ($11 < 100 && $11 >=0) { print } }' intersect_rep2_closest.bed > intersect_rep2_closest_filtered100.bed

bedtools intersect -wa -u -a intersect_summits.bed -b intersect_rep2_closest_filtered100.bed> common_summits.bed
bedtools intersect -wa -u -a rep1_peaks.narrowPeak -b common_summits.bed > common_peaks.narrowPeak 



###### applying directly the summit proximity for the common ones without intersecting first    
bedtools closest -a rep1_summits.bed -b rep2_summits.bed -d -t first > rep1_rep2_closest.bed
awk '{ if ($11 < 100 && $11 >=0) { print } }' rep1_rep2_closest.bed > rep1_rep2_closest_filtered100.bed
bedtools intersect -wa -u -a rep1_summits.bed -b rep1_rep2_closest_filtered100.bed> common_summits.bed
bedtools intersect -wa -u -a rep1_peaks.narrowPeak -b common_summits_TEST.bed > common_peaks.narrowPeak
########## 


wc -l rep1_peaks.narrowPeak                       # or	wc -l rep1_summits.bed 
wc -l rep2_peaks.narrowPeak 
wc -l common_peaks.narrowPeak 
wc -l merged_peaks.narrowPeak 

wc -l rep1_ENCODE_closest_filtered100.bed
wc -l rep2_ENCODE_closest_filtered100.bed
wc -l common_ENCODE_closest_filtered100.bed
wc -l merged_ENCODE_closest_filtered100.bed


###########################################################
# SETTING UP THE ENCODE FILES
###########################################################

nohup wget -O ENCODE_blacklist.bed.gz  https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz &
gunzip ENCODE_blacklist.bed.gz

nohup wget -O ENCODE_peaks.narrowPeak.gz  https://www.encodeproject.org/files/ENCFF863JMY/@@download/ENCFF863JMY.bed.gz &
gunzip ENCODE_peaks.narrowPeak.gz

sort -k1,1 -k2,2n ENCODE_peaks.narrowPeak > ENCODE_bashsorted_peaks.narrowPeak                                     ####NO NEED TO STORE ANOTHER FILE  MODIFY NEXT LINE TOO 
awk '{OFS="\t"; print $1, $2+$10, $2+$10+1,"encode", $9}' ENCODE_bashsorted_peaks.narrowPeak > ENCODE_summits.bed
sort -k1,1 -k2,2n ENCODE_summits.bed > ENCODE_bashsorted_summits.bed		          ####NO NEED TO STORE ANOTHER FILE  MODIFY NEXT LINE TOO
 

###########################################################
# SUMMIT PROXIMITY AND JACCARD
###########################################################
 
bedtools intersect -v -a rep1_summits.bed -b ENCODE_blacklist.bed > rep1_blacklistless_summits.bed

bedtools closest -a rep1_blacklistless_summits.bed -b ENCODE_bashsorted_summits.bed -d -t first > rep1_ENCODE_closest.bed         		# summit proximity applied
awk '{ if ($11 < 100 && $11 >=0) { print } }' rep1_ENCODE_closest.bed > rep1_ENCODE_closest_filtered100.bed			                    # summit proximity applied

bedtools intersect -v -a rep1_peaks.narrowPeak -b ENCODE_blacklist.bed > rep1_blacklistless_peaks.narrowPeak
bedtools jaccard -a rep1_blacklistless_peaks.narrowPeak -b ENCODE_bashsorted_peaks.narrowPeak > rep1_jaccard.txt

# repeat the last 5 lines of code for all the other scenarios (rep2, common, merged) 


###########################################################
# BOX PLOTS
###########################################################
# could be done also with the summit file and cut –f 5, but in this way the code is more readable and intuitive


bedtools intersect -wa -u -a merged_peaks.narrowPeak -b merged_ENCODE_closest_filtered100.bed |  cut -f 9 > qval_merged_overlapping_summits.bed
bedtools intersect -v -a merged_peaks.narrowPeak -b merged_ENCODE_closest_filtered100.bed |  cut -f 9 > qval_merged_non_overlapping_summits.bed

bedtools intersect -wa -u -a rep1_peaks.narrowPeak -b rep1_ENCODE_closest_filtered100.bed |  cut -f 9 > qval_rep1_overlapping_summits.bed
bedtools intersect -v -a rep1_peaks.narrowPeak -b rep1_ENCODE_closest_filtered100.bed |  cut -f 9 > qval_rep1_non_overlapping_summits.bed

bedtools intersect -wa -u -a rep2_peaks.narrowPeak -b rep2_ENCODE_closest_filtered100.bed |  cut -f 9 > qval_rep2_overlapping_summits.bed
bedtools intersect -v -a rep2_peaks.narrowPeak -b rep2_ENCODE_closest_filtered100.bed |  cut -f 9 > qval_rep2_non_overlapping_summits.bed

bedtools intersect -wa -u -a common_peaks.narrowPeak -b common_ENCODE_closest_filtered100.bed |  cut -f 9 > qval_common_overlapping_summits.bed
bedtools intersect -v -a common_peaks.narrowPeak -b common_ENCODE_closest_filtered100.bed |  cut -f 9 > qval_common_non_overlapping_summits.bed

# see R script for continuation 


###########################################################
# ChromHMM
###########################################################

bedtools intersect -a HepG2_chromHMM_18states_hg38.bed -b merged_summits.bed > HepG2_chromHMM_merged_summits.bed             #(MAKE SURE TO NOT PUT UNIQUE)

# see R script for continuation 


###########################################################
# GREAT
###########################################################

sort -k5,5nr merged_summits.bed > merged_bashsorted_summits.bed
awk 'BEGIN {OFS=FS="\t"} {NF--} 1' merged_bashsorted_summits.bed> GREAT_bashsorted_merged_summit.bed

head merged_bashsorted_summits.bed -n 10000 > top10k_merged_summit.bed
awk 'BEGIN {OFS=FS="\t"} {NF--} 1' top10k_merged_summit.bed > GREAT_top10k_merged_summit.bed

head merged_bashsorted_summits.bed -n 5000 > top5k_merged_summit.bed
awk 'BEGIN {OFS=FS="\t"} {NF--} 1' top5k_merged_summit.bed > GREAT_top5k_merged_summit.bed


# upload GREAT_top10k_merged_summit.bed on GREAT, then download and rename the tables needed 


while IFS= read -r line; do
    numbers=$(echo "$line" | grep -oP '\(\K[+-]?\d+(?=\))' | sed 's/)//' | tr '\n' '\t')
    echo -e "${numbers}" >> extracted_numbers.txt
done < alphabetically_ordered_GREAT_10K_gene_table.txt


# upload extracted_numbers.txt on R and proceed with the script 



###########################################################
# GENOME BROWSER
###########################################################

bedtools intersect -v -a merged_peaks.narrowPeak -b ENCODE_peaks.narrowPeak| sort -k7,7nr | head
bedtools intersect -v -a ENCODE_peaks.narrowPeak -b merged_peaks.narrowPeak | sort -k7,7nr | head


###########################################################
# PscaChIP
###########################################################

# Upload  merged_bashsorted_summits.bed 
















------------------------------------------------------------------- WRONG ATTEMPTS

#####WRONG#######sed 's/ \+/\t/g' ENCODE_sorted_peaks.narrowPeak > ENCODE_fixed_sorted_peaks.narrowPeak

bedtools closest -a rep1_peaks.narrowPeak -b ENCODE_bashsorted_peaks.narrowPeak -d -t first > rep1_ENC_closest.bed

ps aux | grep  BCG2023_ramarie
USER         PID %CPU %MEM    VSZ   RSS TTY      STAT START   TIME COMMAND
root     1731860  0.0  0.0  14092  9132 ?        Ss   12:41   0:00 sshd: BCG2023_ramarie [priv]
BCG2023+ 1731967  0.0  0.0  14228  6036 ?        S    12:41   0:00 sshd: BCG2023_ramarie@pts/15
BCG2023+ 1732669  0.0  0.0   6436   724 pts/15   S+   12:48   0:00 grep BCG2023_ramarie

517  bedtools closest -a rep1_summits.bed -b ENCODE_bashsorted_summits.bed -d -t first > rep1_ENCODE_closest.bed
  518  sort -k1,1 -k2,2n ENCODE_bashsorted_summits.bed > ENCODE_resorted_summits.bed
  519  bedtools closest -a rep1_summits.bed -b ENCODE_resorted_summits.bed -d -t first > rep1_ENCODE_closest.bed


bedtools closest -a rep1_summits.bed -b rep1_summits.bed -d -t first > rep1_rep2_closest.bed
bedtools intersect –wa -u -a rep1_summits.bed -b rep2_closest.bed> rep1_blacklistless_summits.bed
awk '{ if ($11 < 100 && $11 >=0) { print } }' rep1_rep2_closest.bed > rep1_rep2_closest_filtered100.bed
############################bedtools intersect –a rep1_summits.bed –b ENCODE_blacklist.bed > rep1_blacklistless_summits.bed 

bedtools intersect -wa -u -a merged_summits.bed -b merged_ENCODE_closest_filtered100.bed |  awk '{print $5}'  > qval_merged_summits.bed
bedtools intersect –v -a merged_summits.bed -b merged_ENCODE_closest_filtered100.bed |  awk '{print $5}’ >> qval_merged_summits.bed

file1.txt > merged_file.txt
file2.txt >> merged_file.txt

@@@@@@@@@ to create the files from which to assign the regulation type to GREAT’s target genes@@@@@@@
grep -v “#” significantly_ordered_GREAT_10k_gene_table.txt | cut -f 2 |  grep -v "NONE"| cut -d " " -f 2,4,6,8 | grep -E "+|-[1-9][0-9]{0,2}“ > sign_test.txt
grep -v “#” alphabetically_ordered_GREAT_10K_gene_table.txt | cut -f 2 |  cut -d " " -f 2,4,6,8,10,12,14,16,18,20 | grep -E "+|-[1-9][0-9]{0,2}"> alpha_test.txt
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



@@@@@@@@@@ Extended reg-ex instead of the perl regex @@@@@@
while IFS= read -r line; do
    numbers=$(echo "$line" | grep -oE '\([+-]?[0-9]+(?=\))' | sed 's/)//' | tr '\n' '\t')
    echo -e "${numbers}" >> extracted_numbers.txt
done < alphabetically_ordered_GREAT_10K_gene_table.txt
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



@@@@@@@@@ finding peaks unique to either my experiment or encode --- WRONG SINCE DONE ON SUMMITS rather than peaks@@@@@
bedtools intersect -v -a merged_bashsorted_summits.bed -b ENCODE_bashsorted_summits.bed | head
bedtools intersect -v -a ENCODE_bashsorted_summits.bed -b merged_ENCODE_closest_filtered100.bed | sort -k5,5nr merged_summits.bed | head
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

