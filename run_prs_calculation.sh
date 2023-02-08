PATH_TO_DATA=""
PATH_TO_OUT=""

fir chr in {1..22}; do
python3 prs_calculation.py --model \
                           --snp-id 1 \
                           --ea 3 \
                           --beta 6 \
                           --bcf ${PATH_TO_DATA}/target_file_chr${chr}_updated_without_extension \
                           --output ${PATH_TO_OUT}

done

awk '{num[FNR]=$1;arr[FNR]+=$2}END{for(i=1;i<=FNR;i+=1){print num[i]"\t"arr[i]}}' ${PATH_TO_OUT}/target_file_chr*_updated_without_extension.sscore \
> ${PATH_TO_OUT}/target_file_all_updated_without_extension.sscore
