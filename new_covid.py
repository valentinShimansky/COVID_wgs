import subprocess
import argparse


def trimming(wd, id_value):
    trimm = "/home/MegaBOLT/anaconda3/bin/trimmomatic"
    subprocess.check_output(
        f"{trimm} PE -threads 60 "
        f"{wd}/{id_value}_1.fq {wd}/{id_value}_2.fq "
        f"{wd}/{id_value}_1P.fq {wd}/{id_value}_1U.fq "
        f"{wd}/{id_value}_2P.fq {wd}/{id_value}_2U.fq "
        f"LEADING:30 TRAILING:30",
        shell=True)
    return


def alignment(wd, id_value):
    samtools = "/mnt/ssd/MegaBOLT_scheduler/bin/samtools"
    bwa = "/data/SOFT/bwa/bwa"
    subprocess.check_output(
        f'{bwa} mem -t 60 /mnt/ssd/MegaBOLT/reference/hg38.fa '
        f'{wd}/{id_value}_1P.fq {wd}/{id_value}_2P.fq | '
        f'{samtools} sort -@ 60 -n -O BAM -o {wd}/{id_value}_sorted.bam -',
        shell=True)
    return


def mark_dup(wd, id_value):  # REPLACE WITH SAMTOOLS AFTER TEST
    samtools = "/mnt/ssd/MegaBOLT_scheduler/bin/samtools"
    subprocess.check_output(f'{samtools} fixmate -@ 60 -m '
                            f'{wd}/{id_value}_sorted.bam '
                            f'{wd}/{id_value}_fixmated.bam',
                            shell=True)
    subprocess.check_output(f'{samtools} sort -@ 60 -o '
                            f'{wd}/{id_value}_fix_sort.bam '
                            f'{wd}/{id_value}_fixmated.bam',
                            shell=True)
    subprocess.check_output(f'{samtools} markdup -@ 60 -s '
                            f'{wd}/{id_value}_fix_sort.bam '
                            f'{wd}/{id_value}_sorted_dedup.bam',
                            shell=True)

    return


def bam_indexing(wd, id_value):
    samtools = "/mnt/ssd/MegaBOLT_scheduler/bin/samtools"
    subprocess.check_output(f'{samtools} index -@ 60 {wd}/{id_value}_sorted_dedup.bam',
                            shell=True)
    return


def variant_calling(wd, id_value):
    list_chr = [str(i) for i in range(1, 23)] + ['X']
    for chrom in list_chr:
        print('CALLING ID', id_value, ' CHROM', chrom)
        subprocess.check_output(
            f'/usr/bin/bcftools mpileup --threads 60 '
            f'-f /data/COVID/hs38Dh/hs38DH.fa -I -E -a "FORMAT/DP" '
            f'-T /data/COVID/reference_panel_Glimpse/1000GP.chr{chrom}.sites.tsv.gz '
            f'-r chr{chrom} {wd}/{id_value}_sorted_dedup.bam -Ou | '
            f'bcftools call --threads 60 -Aim -C alleles '
            f'-T /data/COVID/reference_panel_Glimpse/1000GP.chr{chrom}.sites.tsv.gz '
            f'-Oz -o {wd}/vcf_chrs/{id_value}_chr{chrom}.vcf.gz',
            shell=True)
        subprocess.check_output(f'/usr/bin/bcftools index --threads 60 {wd}/vcf_chrs/{id_value}_chr{chrom}.vcf.gz',
                                shell=True)
    return


def making_id_list(wd):
    id_list = []
    with open(f'{wd}/ID_table_unmodified.csv', 'r') as ID_table_unmod:
        for s in ID_table_unmod:
            s = s.strip().split('\t')
            id_list.append(s[0])
    return id_list


def merging_vcfs(wd, id_list):
    list_chr = [str(i) for i in range(1, 23)] + ['X']
    for chrom in list_chr:
        print('making chr lists')
        with open(f'{wd}/vcf_chrs/chr{chrom}_list_paths_sed.txt', 'w+') as T:
            for id_value in id_list:
                T.write(f'{wd}/vcf_chrs/{id_value}_chr{chrom}.vcf.gz\n')
        print('merging chromosome ', chrom)
        subprocess.check_output(
            f'bcftools merge --threads 60 '
            f'-m none -r chr{chrom} -Oz -o {wd}/vcfs_merged/target_samples_merged.chr{chrom}.vcf.gz '
            f'-l {wd}/vcf_chrs/chr{chrom}_list_paths_sed.txt',
            shell=True)
        subprocess.check_output(
            f'bcftools index --threads 60 {wd}/vcfs_merged/target_samples_merged.chr{chrom}.vcf.gz', shell=True)
    return


def execute(args):
    rawdata = args.rawdata
    flowcell = args.flowcell
    work_dir = args.work_dir

    with open(f'{work_dir}/ID_table.csv', 'r') as ID_table:
        for s in ID_table:
            s = s.strip().split('\t')
            print(s)
            id_value = s[0]
            for i in (1, 2):
                if ',' in s[1]:
                    barcode1 = s[1].split(',')[0]
                    barcode2 = s[1].split(',')[1]
                    subprocess.check_output(f"zcat "
                                            f"{rawdata}/L01/{flowcell}_L01_{barcode1}_{i}.fq.gz "
                                            f"{rawdata}/L01/{flowcell}_L01_{barcode2}_{i}.fq.gz "
                                            f"{rawdata}/L02/{flowcell}_L02_{barcode1}_{i}.fq.gz "
                                            f"{rawdata}/L02/{flowcell}_L02_{barcode2}_{i}.fq.gz "
                                            f"{rawdata}/L03/{flowcell}_L03_{barcode1}_{i}.fq.gz "
                                            f"{rawdata}/L03/{flowcell}_L03_{barcode2}_{i}.fq.gz "
                                            f"{rawdata}/L04/{flowcell}_L04_{barcode1}_{i}.fq.gz "
                                            f"{rawdata}/L04/{flowcell}_L04_{barcode2}_{i}.fq.gz "
                                            f"> {work_dir}/{id_value}_{i}.fq",
                                            shell=True)
            print(s)
            print('TRIMMING ID  ', id_value)
            trimming(work_dir, id_value)
            subprocess.check_output(f"rm -rf {work_dir}/{id_value}_1.fq  {work_dir}/{id_value}_2.fq",
                                    shell=True)
            print('ALIGNING ID  ', id_value)
            alignment(work_dir, id_value)
            subprocess.check_output(
                f"rm -rf {work_dir}/{id_value}_1P.fq  {work_dir}/{id_value}_2P.fq "
                f"{work_dir}/{id_value}_1U.fq {work_dir}/{id_value}_2U.fq",
                shell=True)
            print('MARKING DUPLICATES  ', id_value)
            mark_dup(work_dir, id_value)
            subprocess.check_output(
                f'rm -rf {work_dir}/{id_value}_sorted.bam '
                f'{work_dir}/{id_value}_fixmated.bam '
                f'{work_dir}/{id_value}_fix_sort.bam',
                shell=True)
            print('INDEXING ID  ', id_value)
            bam_indexing(work_dir, id_value)
            print('STARTING CALLING ID  ', id_value)
            variant_calling(work_dir, id_value)
            subprocess.check_output(f'rm -rf {work_dir}/{id_value}_sorted_dedup.bam '
                                    f'{work_dir}/{id_value}_sorted_dedup.bam.bai',
                                    shell=True)

    merging_vcfs(work_dir, making_id_list(work_dir))


def parse():
    parser = argparse.ArgumentParser()
    input_group = parser.add_argument_group("Input parameters")
    input_group.add_argument('--rawdata', type=str, required=True,
                             help='Path to rawdata.')
    input_group.add_argument('--flowcell', type=str, required=True,
                             help='Path to flowcell.')
    input_group.add_argument('--work_dir', type=str, required=True,
                             help='Path to work directory.')
    return parser


def run_alignment_main():
    args = parse().parse_args()
    try:
        execute(args)
    except Exception as error:
        raise error


if __name__ == '__main__':
    run_alignment_main()
