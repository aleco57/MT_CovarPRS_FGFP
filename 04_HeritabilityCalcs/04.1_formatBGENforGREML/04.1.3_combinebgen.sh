####################################################################
### This is a script to combine the begn files bychr into one bgen file
####################################################################
data_dir="../../../../data"
geno_prefix="${data_dir}/PhenoPRS/greml/bgen"

cmd=""
for chrom in {01,02,03,04,05,06,07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22}; do

        cmd="${cmd} ${geno_prefix}/qc_chr${chrom}_wchr.bgen"

done

#This is the path to bgenix tool. This is found on my local bp dir. Change if performing this elsewhere
~/bgenix/./build/apps/cat-bgen \
    -g ${cmd} \
    -clobber \
    -og ${geno_prefix}/qcwchr_allchr.bgen

