import os
import subprocess
import argparse

######### parsing arguments
parser = argparse.ArgumentParser()
parser.add_argument('--rawdata')
parser.add_argument('--flowcell')
parser.add_argument('--work_dir')
args = parser.parse_args()

######### Place where raw reads are set
rawdata = args.rawdata
flowcell = args.flowcell

######### Working directory
work_dir = args.work_dir

######### Processing

def trimming (wd, id):
    TRIMM = "/home/MegaBOLT/anaconda3/bin/trimmomatic"
    subprocess.check_output(
        f"{TRIMM} PE -threads 60 {wd}/{id}_1.fq {wd}/{id}_2.fq {wd}/{id}_1P.fq {wd}/{id}_1U.fq {wd}/{id}_2P.fq {wd}/{id}_2U.fq LEADING:30 TRAILING:30",
        shell=True)
    return

def alignment (wd, id):
    samtools = "/mnt/ssd/MegaBOLT_scheduler/bin/samtools"
    bwa = "/data/SOFT/bwa/bwa"
    subprocess.check_output(
        f'{bwa} mem -t 60 /mnt/ssd/MegaBOLT/reference/hg38.fa {wd}/{id}_1P.fq {wd}/{id}_2P.fq | {samtools} sort -@ 60 -n -O BAM -o {wd}/{id}_sorted.bam -',
        shell=True)
    return

def mark_dup (wd, id): # REPLACE WITH SAMTOOLS AFTER TEST
    samtools = "/mnt/ssd/MegaBOLT_scheduler/bin/samtools"
    subprocess.check_output(f'{samtools} fixmate -@ 60 -m {wd}/{id}_sorted.bam {wd}/{id}_fixmated.bam',
        shell=True)
    subprocess.check_output(f'{samtools} sort -@ 60 -o {wd}/{id}_fix_sort.bam {wd}/{id}_fixmated.bam',
        shell=True)
    subprocess.check_output(f'{samtools} markdup -@ 60 -s {wd}/{id}_fix_sort.bam {wd}/{id}_sorted_dedup.bam',
        shell=True)

    return

def bam_indexing (wd, id):
    samtools = "/mnt/ssd/MegaBOLT_scheduler/bin/samtools"
    subprocess.check_output(f'{samtools} index -@ 60 {wd}/{id}_sorted_dedup.bam',
                            shell=True)
    return

def variant_calling (wd, id):
    for chrom in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18',
                  '19', '20', '21', '22', 'X']:
        print('CALLING ID', id, ' CHROM', chrom)
        subprocess.check_output(
            f'/usr/bin/bcftools mpileup --threads 60 -f /data/COVID/hs38Dh/hs38DH.fa -I -E -a "FORMAT/DP" -T /data/COVID/reference_panel_Glimpse/1000GP.chr{chrom}.sites.tsv.gz -r chr{chrom} {wd}/{id}_sorted_dedup.bam -Ou | bcftools call --threads 60 -Aim -C alleles -T  /data/COVID/reference_panel_Glimpse/1000GP.chr{chrom}.sites.tsv.gz -Oz -o {wd}/vcf_chrs/{id}_chr{chrom}.vcf.gz',
            shell=True)
        subprocess.check_output(f'/usr/bin/bcftools index --threads 60 {wd}/vcf_chrs/{id}_chr{chrom}.vcf.gz',
                                shell=True)
    return

def making_ID_list(wd):
    ID_list = []
    with open(f'{wd}/ID_table_unmodified.csv', 'r') as ID_table_unmod:
        for s in ID_table_unmod:
            s = s.strip().split('\t')
            ID_list.append(s[0])
    return ID_list

def merging_vcfs(wd, id_list):
    for chr in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
              '20', '21', '22', 'X']:
        print('making chr lists')
        with open(f'{wd}/vcf_chrs/chr{chr}_list_paths_sed.txt', 'w+') as T:
            for id in id_list:
                T.write(f'{wd}/vcf_chrs/{id}_chr{chr}.vcf.gz\n')
        print('merging chromosome ', chr)
        subprocess.check_output(
            f'bcftools merge --threads 60 -m none -r chr{chr} -Oz -o  {wd}/vcfs_merged/target_samples_merged.chr{chr}.vcf.gz -l {wd}/vcf_chrs/chr{chr}_list_paths_sed.txt',
            shell=True)
        subprocess.check_output(
            f'bcftools index --threads 60 {wd}/vcfs_merged/target_samples_merged.chr{chr}.vcf.gz', shell=True)
    return

######### Main

with open(f'{work_dir}/ID_table.csv', 'r') as ID_table:
    for s in ID_table:
        s = s.strip().split('\t')
        print(s)
        ID = s[0]
        for i in (1, 2):
            if ',' in s[1]:
                barcode1 = s[1].split(',')[0]
                barcode2 = s[1].split(',')[1]
                subprocess.check_output(f"zcat  {rawdata}/L01/{flowcell}_L01_{barcode1}_{i}.fq.gz \
                                               {rawdata}/L01/{flowcell}_L01_{barcode2}_{i}.fq.gz \
                                               {rawdata}/L02/{flowcell}_L02_{barcode1}_{i}.fq.gz \
                                               {rawdata}/L02/{flowcell}_L02_{barcode2}_{i}.fq.gz \
                                               {rawdata}/L03/{flowcell}_L03_{barcode1}_{i}.fq.gz \
                                               {rawdata}/L03/{flowcell}_L03_{barcode2}_{i}.fq.gz \
                                               {rawdata}/L04/{flowcell}_L04_{barcode1}_{i}.fq.gz \
                                               {rawdata}/L04/{flowcell}_L04_{barcode2}_{i}.fq.gz \
                                                > {work_dir}/{ID}_{i}.fq",
                                                shell=True)
        print(s)
        print('TRIMMING ID  ', ID)
        trimming(work_dir, ID)
        subprocess.check_output(f"rm -rf {work_dir}/{ID}_1.fq  {work_dir}/{ID}_2.fq", shell=True)
        print('ALIGNING ID  ', ID)
        alignment(work_dir, ID)
        subprocess.check_output(f"rm -rf {work_dir}/{ID}_1P.fq  {work_dir}/{ID}_2P.fq {work_dir}/{ID}_1U.fq {work_dir}/{ID}_2U.fq", shell=True)
        print('MARKING DUPLICATES  ', ID)
        mark_dup(work_dir, ID)
        subprocess.check_output(f'rm -rf {work_dir}/{ID}_sorted.bam {work_dir}/{ID}_fixmated.bam {work_dir}/{ID}_fix_sort.bam', shell=True)
        print('INDEXING ID  ', ID)
        bam_indexing(work_dir, ID)
        print('STARTING CALLING ID  ', ID)
        variant_calling(work_dir, ID)
        subprocess.check_output(f'rm -rf {work_dir}/{ID}_sorted_dedup.bam {work_dir}/{ID}_sorted_dedup.bam.bai', shell=True)

merging_vcfs(work_dir, making_ID_list(work_dir))
