#!/usr/bin/env nextflow


def interpolateReferenceFile(Map reference, def chr) {
    assert reference.dir != null : "Missing dir in ${reference}"
    assert reference.prefix != null : "Missing prefix in ${reference}"
    assert reference.postfix != null : "Missing postfix in ${reference}"

    return "${reference.dir}/${reference.prefix}${chr}${reference.postfix}"
}


def getRegionForChr(chr) {
    def chrParamName = "chr${chr}"
    assert params.refs.chrSize[chrParamName] != null : "Could not find size for chr${chr}"
    def chrEnd = params.refs.chr[chrParamName]

    return "${chr}:0-${chrEnd}"
}


process ToConformedBCF {
    tag "${trio_id}.${relationship}"
    publishDir "${params.output_dir}/conformed_vcfs/${trio_id}", mode: "link"

    input:
    tuple val(trio_id), val(relationship), path(vcf), path(vcf_csi)

    output:
    tuple val(trio_id), val(relationship), path("${trio_id}.${relationship}.bcf"), path("${trio_id}.${relationship}.bcf.csi")

    script:
    """
        export TMPDIR=/ifs/temp/
        module load bioinf/bcftools

        echo "${trio_id}.${relationship}" > new_name.txt

        # Filter to autosomal chromosomes only, then remove 'chr' prefix from regions
        bcftools view -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 $vcf | bcftools annotate --rename-chrs $params.chr_map_file -Oz -o "${trio_id}.${relationship}.temp.vcf.gz"
        bcftools reheader -s new_name.txt "${trio_id}.${relationship}.temp.vcf.gz" > "${trio_id}.${relationship}.reheadered.vcf.gz"

        # Convert to BCF extension
        bcftools view -Ob -o "${trio_id}.${relationship}.bcf" "${trio_id}.${relationship}.reheadered.vcf.gz"
        bcftools index "${trio_id}.${relationship}.bcf"
    """
}

process MergeRelatedBCFs {
    tag "${trio_id}"
    publishDir "${params.output_dir}/trio_bcfs/${trio_id}", mode: "link"

    input:
    tuple val(trio_id), path(proband_bcf), path(proband_bcf_csi), path(mother_bcf), path(mother_bcf_csi), path(father_bcf), path(father_bcf_csi)

    output:
    tuple val(trio_id), path("${trio_id.bcf"), path("${trio_id}.bcf.csi")

    script:
    """
        export TMPDIR=/ifs/temp/
        module load bioinf/bcftools

        bcftools merge -Ob -o ${trio_id}.bcf ${proband_vcf} ${mother_vcf} ${father_vcf}
        bcftools index ${trio_id}.bcf
    """
}

process BatchBCFs {
    tag "${chr}"
    publishDir "${params.output_dir}/batched_bcfs/", mode: "link"

    input:
    tuple val(chr), path(bcfs), path(bcf_csis)

    output:
    tuple val(chr), path("chr${chr}.bcf"), path("chr${chr}.bcf.csi")

    script:
    """
        export TMPDIR=/ifs/temp/
        module load bioinf/bcftools

        bcftools merge -Ob -o merged.bcf $bcfs
        bcftools index merged.bcf

        bcftools +fill-tags merged.bcf -Oz -o chr${chr}.bcf -- -t AC,AN,AF
        bcftools index chr${chr}.bcf

        # Cleanup
        rm merged.bcf
    """
}

process SelectRegionFromBCF {
    tag "${trio_id}"
    publishDir "${params.output_dir}/batched_bcfs_chr/${trio_id}/${chr}", mode: "link"

    input:
    tuple val(trio_id), val(chr), path(bcf), path(bcf_csi)

    output:
    tuple val(trio_id), val(chr), path("${trio_id}.chr${chr}.bcf"), path("${trio_id}.chr${chr}.bcf.csi")

    script:
    """
        export TMPDIR=/ifs/temp/
        module load bioinf/bcftools

        bcftools view -r ${chr} -Ob -o ${trio_id}.chr${chr}.bcf $bcf
        bcftools index ${trio_id}.chr${chr}.bcf
    """
}

process PhaseBCF {
    tag "${chr}"
    publishDir "${params.output_dir}/phased_bcfs", mode: "link"

    input:
    tuple val(chr), path(bcf), path(bcf_csi)

    output:
    tuple val(chr), path("chr${chr}.bcf"), path("chr${chr}.bcf.csi")

    script:
    def ref = interpolateReferenceFile(params.interpolated_refs.1000genomes_genotype, chr)
    def map = interpolateReferenceFile(params.interpolated_refs.genetic_maps, chr)
    """
        export TMPDIR=/ifs/temp/
        module load bioinf/bcftools

        ${params.bins.shapeit5} --input ${bcf} --pedigree ${params.ref.ped} --region ${chr} --map ${map} --output chr${chr}.bcf --reference ${ref} --thread ${task.cpus}
        bcftools index chr${chr}.bcf
    """
}

process ImputeBCF {
    tag "${chr}"
    publishDir "${params.output_dir}/imputed_bcfs", mode: "link"

    input:
    tuple val(chr), path(bcf), path(bcf_csi)

    output:
    tuple val(chr), path("chr${chr}.vcf.gz"), path("chr${chr}.vcf.gz.csi")

    script:
    def ref = interpolateReferenceFile(params.interpolated_refs.1000genomes_genotype, chr)
    def map = interpolateReferenceFile(params.interpolated_refs.genetic_maps, chr)
    def region = getRegionForChr(chr)
    """
        export TMPDIR=/ifs/temp/
        module load bioinf/bcftools

        ${params.bins.impute5} --h ${ref} --g ${bcf} --m ${map} --r ${region} --buffer-region ${region} --out-ap-field --o chr${chr}.vcf.gz --threads ${task.cpus}
        bcftools index chr${chr}.vcf.gz
    """
}

process FetchInfoScoreDistribution {
    tag "${chr}"
    publishDir "${params.output_dir}/imputation_scores", mode: "link"

    input:
    tuple val(chr), path(bcf), path(bcf_csi)

    output:
    tuple val(chr), path("chr${chr}.tsv")

    script:
    """
        export TMPDIR=/ifs/temp/
        module load bioinf/bcftools

        bcftools view -i "" $bcf | bcftools query -Oz -o chr${chr}.tsv -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/INFO\\n"
    """
}

process CalculateOffTargetCoverage {
    tag "${trio_id}.${relationship}"
    publishDir "${params.output_dir}/coverage/${trio_id}", mode: "link"

    input:
    tuple val(trio_id), val(relationship), path(cram), path(cram_crai)

    output:
    tuple val(trio_id), val(relationship), path("${trio_id}.${relationship}.tsv")

    script:
    """
        export TMPDIR=/ifs/temp/
        module load bioinf/mosdepth

        mosdepth --thresholds 1 "${trio_id}_${relationship}" "${cram}"
        zcat "${trio_id}_${relationship}.thresholds.bed.gz" | awk 'NR>1 {covered+=$5; total+=($3-$2)} END {print covered"\\t"total"\\t"covered/total"\n"}' > ${trio_id}.${relationship}.tsv
    """
}

workflow {
    trio_dirs = Channel.fromPath("${params.input_dir}/rumc_trio_*/", type: "dir")

    def relationships = Channel.of("proband", "mother", "father")
    def trio_files = trio_dirs.combine(relationships).map { trio_dir, relationship ->
        def trio_id = trio_dir.getName()
        def vcf = file("${trio_dir}/${trio_id}.vcf.gz")
        def vcf_csi = file("${trio_dir}/${trio_id}.vcf.gz.csi")
        def cram = file("${trio_dir}/${trio_id}.cram")
        def cram_crai = file("${trio_dir}/${trio_id}.cram.crai")
        tuple(trio_id, vcf, vcf_csi, cram, cram_crai)
    }

    // Get coverage
    def trio_crams = trio_files.map { (trio_id, relationship, vcf, vcf_csi, cram, cram_crai) -> tuple(trio_id, relationship, cram, cram_crai) }
    trio_crams | CalculateOffTargetCoverage | set { coverage }


    // Get info score distributions
    def trio_vcfs = trio_files.map { (trio_id, relationship, vcf, vcf_csi, cram, cram_crai) -> tuple(trio_id, relationship, vcf, vcf_csi) }

    trio_vcfs | ToConformedBCF | set { conformed_bcfs }
    conformed_bcfs.branch { it ->
        proband: it.relationship == "proband"
        mother: it.relationship == "mother"
        father: it.relationship == "father"
    }.map { (trio_id, relationship, bcf, bcf_csi) -> tuple(trio_id, bcf, csi) }.set { branched_bcfs }

    branched_bcfs.proband
        .join(branched_bcfs.mother)
        .join(branched_bcfs.father)
        .toSortedList(a, b -> a[0] <=> b[0])
        .set { joined_trio_bcfs }

    joined_trio_bcfs | MergeRelatedBCFs | set { trio_bcfs }

    def chromosomes = Channel.from(1..22)
    trio_bcfs.combine(chromosomes).map { (trio_id, bcf, bcf_csi, chr) -> tuple(trio_id, chr, bcf, bcf_csi) }.set { trio_chr_bcfs }
    trio_chr_bcfs.groupTuple(by: 1).map { (chr, trio_id, bcfs, bcf_csis) -> tuple(chr, bcfs, bcf_csis) }.set { chr_bcfs }

    chr_bcfs | BatchBCFs | PhaseBCF | ImputeBCF | FetchInfoScoreDistribution | set { scores }

}

