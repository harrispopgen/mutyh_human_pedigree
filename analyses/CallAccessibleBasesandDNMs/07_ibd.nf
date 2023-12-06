nextflow.enable.dsl=2
​
def getChromosomeNumber( file ) {
    regexp = /chr[0-9,X]+/
    (file =~ regexp)[0]
}
​
process merge_with_onekg {
    publishDir  "$baseDir/nf_output/vcfs/${chr_num}/stats", pattern: "*.stats", mode: 'copy'
    module 'htslib/1.12 : bcftools/1.12'
    input:
        tuple val(chr_num), path(onekg_vcf), path(onekg_vcf_idx), path(pedigree_vcf), path(pedigree_vcf_idx) 
        
    output:
        tuple val(chr_num), path("chr*_merge.vcf.gz"), path("chr*_merge.vcf.gz.tbi"), emit: vcf
        path("*.stats")
        
    shell:
    '''
    vcf_output_filename="!{chr_num}_merge.vcf.gz"
    stats_output_filename="!{chr_num}_merge.stats"
​
    bcftools merge -f PASS -Oz -o "${vcf_output_filename}" !{onekg_vcf} !{pedigree_vcf} 
    tabix -p vcf "${vcf_output_filename}"
    bcftools stats ${vcf_output_filename} >> ${stats_output_filename}
    '''
}
​
// keep only bi-allelic snps, then filter out sites with <10% minor allele frequency
// lastly remove all one kg samples from the vcf
process filter_vcf {
    publishDir  "$baseDir/nf_output/vcfs/${chr_num}/stats", pattern: "*.stats", mode: 'copy'
    module 'htslib/1.12 : bcftools/1.12'
    input:
        each path(pedigree_samples)
        tuple val(chr_num), path(vcf_file), path(vcf_idx_file)
        
    output:
        tuple val(chr_num), path("chr*_filtered.vcf.gz"), path("chr*_filtered.vcf.gz.tbi"), emit: vcf
        path("*.stats")
        
    shell:
    '''
    vcf_output_filename="!{chr_num}_filtered.vcf.gz"
    stats_output_filename="!{chr_num}_filtered.stats"
​
    bcftools view -v snps -m2 -M2 -Ou !{vcf_file} | bcftools view -i "INFO/AN>30 & MAF>=0.1" -Ou | bcftools view -S !{pedigree_samples} -Oz -o ${vcf_output_filename} 
    bcftools stats ${vcf_output_filename} >> ${stats_output_filename}
    tabix -p vcf ${vcf_output_filename}
    '''
}
​
process remove_unphased_and_missing_sites {
    publishDir  "$baseDir/nf_output/vcfs/${chr_num}/stats", pattern: "*.log", mode: 'copy'
    publishDir  "$baseDir/nf_output/vcfs/${chr_num}", pattern: "*filter.vcf.gz*", mode: 'copy'
    module 'htslib/1.12 : bcftools/1.12'
    input:
        tuple val(chr_num), path(vcf_file), path(vcf_idx_file)
    output:
        tuple val(chr_num), path("chr*filter.vcf.gz"), path("chr*filter.vcf.gz.tbi"), emit: vcf
        path("*.log")
        
    shell:
    '''
    vcf_output_filename="!{chr_num}_phm_filter.vcf.gz"
    log_filename="unphased_and_missing_sites.log"
    bcftools view -HG !{vcf_file} | wc -l >> ${log_filename}
    bcftools view -g ^miss -Ou !{vcf_file} | bcftools view --phased -Oz -o ${vcf_output_filename}
    bcftools view -HG ${vcf_output_filename} | wc -l >> ${log_filename}
    tabix -p vcf ${vcf_output_filename}
    '''
}
​
process reformat_map_file {
    publishDir "${output_dir}/map_files/", pattern: "plink.chr*.map", mode: "copy"
    input:
        tuple val(chr_num), path(map_file)
​
    output:
        tuple val(chr_num), path("plink.chr*.map")
​
    shell:
    '''
    # need to reformat the map file by appending "chr" in front of the chromosome number 
    # to match vcf chromosome format
​
    chr_map_file="plink.!{chr_num}.map"
    sed -e 's/^/chr/' !{map_file} > ${chr_map_file}
    '''
}
​
// infers ibd segment using hap-ibd 
// three output files, of suffix *.hbd.gz, *.ibd.gz, *.log
// ibd file (.ibd.gz) contains IBD segments shared between individuals. 
// hbd file (.hbd.gz) contains HBD segments within individuals.
process infer_ibd {
    publishDir "${output_dir}/hap-ibd", pattern: "chr*hap.*", mode: "copy"
    label "hap_ibd"
​
    input:
        each path(hap_ibd_jar)
        tuple val(chr_num), path(vcf_file), path(vcf_idx_file), path(map_file)
        
    output:
        tuple val(chr_num), path("chr*_hap.ibd"), emit: ibd
        path "chr*_hap.*"
    
    shell:
    '''
        output_filename="!{chr_num}_hap"
        java -Xmx6g -jar !{hap_ibd_jar} gt="!{vcf_file}" map="!{map_file}" out="${output_filename}" min-seed=1.0 max-gap=1000 min-extend=0.2 min-output=2 min-markers=100
        gunzip ${output_filename}.ibd.gz
    '''
}
​
process convert_dataformat {
    publishDir "${output_dir}/reformat_ibd/", pattern: "*_reformat.ibd", mode: "copy"
    input:
        tuple val(chr_num), path(hap_ibd_file)
​
    output:
        tuple val(chr_num), path("*_reformat.ibd")
​
    shell:
    '''
    reformat_ibd_output.py !{hap_ibd_file}
    mv reformat.ibd !{chr_num}_reformat.ibd
    '''
}
​
// 15 individual WGS vcf
ped_vcfs = Channel.fromPath(params.phased_pedigree_vcf_dir + "/*.vcf.gz")
                 .map(f -> [getChromosomeNumber(f.name), f, f+".tbi"])
​
one_kg_vcfs = Channel.fromPath(params.onekg_vcf_dir + "/*chr[0-9]*.vcf.gz")
                     .map(f -> [getChromosomeNumber(f.name), f, f+".tbi"])
​
one_kg_vcfs.combine(ped_vcfs, by:0)
             .set { merge_input_ch }
​
map_files_ch = Channel.fromPath(params.map_file_dir + "/*.chr[0-9]*map")
                     .map(f -> [getChromosomeNumber(f.name), f])
​
chr_ch = Channel.from(1..22).map(chr_num -> "chr" + chr_num)
pedigree_sample_file = params.pedigree_sample_file
output_dir = "$baseDir/nf_output"
​
workflow {
    merge_with_onekg( merge_input_ch )
    filter_vcf( pedigree_sample_file, merge_with_onekg.out.vcf )
    remove_unphased_and_missing_sites( filter_vcf.out.vcf )
​
    reformat_map_file( map_files_ch )
    remove_unphased_and_missing_sites.out.vcf.combine(
        reformat_map_file.out, by:0
    ).set { vcf_map_ch }
​
    infer_ibd(params.hap_ibd_jar, vcf_map_ch)
    convert_dataformat( infer_ibd.out.ibd )
}
