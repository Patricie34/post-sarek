nextflow.enable.dsl=2

// ============================================
// HELPER FUNCTION: Parse sample info from filename (like Ruby)
// ============================================

def parseSampleInfo(vcf_path) {

    def fullpath = vcf_path.toString()

    // filename only
    def filename = fullpath.tokenize('/')[-1]

    // remove extension + tool suffix
    def basename = filename
        .replace('.freebayes.filtered.norm.sorted_VEP.ann.vcf.gz', '')
        .replace('.vcf.gz', '')

    // split tumor vs normal using reliable delimiter
    def parts = basename.split('_vs_')

    if (parts.size() != 2) {
        error "Expected '_vs_' structure but got: ${filename}"
    }

    def tumor  = parts[0]
    def normal = parts[1]

    // sample_id = stable prefix before first underscore group
    def sample_id = tumor.split('_')[0]

    return [sample_id, tumor, normal]
}

// ============================================
// HELPER FUNCTION: Derive file paths from FreeBayes path (like Ruby .sub())
// ============================================

def deriveFilePaths(freebayes_vcf) {
    def path_str = freebayes_vcf.toString()

    def mutect_vcf = path_str.replace('freebayes/', 'mutect2/').replace('.freebayes.', '.mutect2.')

    def strelka_base = path_str.replace('freebayes/', 'strelka/').replace('.freebayes.', '.strelka.')
    def strelka_snv = strelka_base.replace('.filtered.norm.sorted_VEP.', '.somatic_snvs.norm.sorted_VEP.')
    def strelka_indel = strelka_base.replace('.filtered.norm.sorted_VEP.', '.somatic_indels.norm.sorted_VEP.')

    return [mutect_vcf, strelka_snv, strelka_indel]
}

// ============================================
// PROCESSES
// ============================================

process FREEBAYES_BCBIO {
    tag { sample_id }
    container 'patricie/freebayes-filter:v5'
    publishDir "${params.outdir}/${tumor}_vs_${normal}", mode: 'copy'
    publishDir "${params.outdir2}/annotation/freebayes/${tumor}_vs_${normal}", mode: 'copy'
    label 's_cpu'
    label 'm_mem'
    
    input:
    tuple val(sample_id), val(tumor), val(normal), path(freebayes_vcf)

    output:
    tuple val(sample_id), val(tumor), val(normal), path("*.bcbio.vcf.gz"), path("*.bcbio.vcf.gz.tbi"), emit: freebayes_bcbio

    script:
    """
    python3 /app/freebayes_debug.py \
        -i ${freebayes_vcf} \
        -r ${params.fasta} \
        --tumor ${sample_id}_${tumor} \
        --normal ${sample_id}_${normal} \
        -d
    """
}

process REHEADER_VCFs {
    tag { sample_id }
    container 'staphb/bcftools:1.23.1'
    publishDir "${params.outdir}/${tumor}_vs_${normal}", mode: 'copy'
    publishDir "${params.outdir2}/annotation/freebayes/${tumor}_vs_${normal}", mode: 'copy'
    label 's_cpu'
    label 's_mem'

    input:
    tuple val(sample_id), val(tumor), val(normal), path(freebayes_vcf), path(freebayes_vcf_tbi), path(mutect_vcf), path(strelka_snv), path(strelka_indel)
        
    output:
    tuple val(sample_id), val(tumor), val(normal), \
          path("*.somatic_snvs.rename.ann.vcf.gz"), path("*.somatic_snvs.rename.ann.vcf.gz.tbi"), \
          path("*.somatic_indels.rename.ann.vcf.gz"), path("*.somatic_indels.rename.ann.vcf.gz.tbi"), \
          path("*.mutect2.rename.ann.vcf.gz"), path("*.mutect2.rename.ann.vcf.gz.tbi"), \
          path("*.freebayes.rename.ann.vcf.gz"), path("*.freebayes.rename.ann.vcf.gz.tbi"), \
          emit: reheadered

    script:
    def snv_out   = "${tumor}_vs_${normal}.strelka.somatic_snvs.rename.ann.vcf.gz"
    def indel_out = "${tumor}_vs_${normal}.strelka.somatic_indels.rename.ann.vcf.gz"
    def mutect_out = "${tumor}_vs_${normal}.mutect2.rename.ann.vcf.gz"
    def freebayes_out = "${tumor}_vs_${normal}.freebayes.rename.ann.vcf.gz"
    """
    cat > rename.file << EOF
    TUMOR\t${tumor}
    NORMAL\t${normal}
    EOF

    bcftools reheader -s rename.file ${strelka_snv} -o ${snv_out}
    bcftools index -t ${snv_out}

    bcftools reheader -s rename.file ${strelka_indel} -o ${indel_out}
    bcftools index -t ${indel_out}

    cat > other_samples.txt << EOF
    ${tumor}
    ${normal}
    EOF

    bcftools reheader -s other_samples.txt ${mutect_vcf} -o ${mutect_out}
    bcftools index -t ${mutect_out}

    bcftools reheader -s other_samples.txt ${freebayes_vcf} -o ${freebayes_out}
    bcftools index -t ${freebayes_out}
    """
}

process COMBINE_VARIANTS {
    tag { sample_id }
    container 'broadinstitute/gatk3:3.6-0'
    publishDir "${params.outdir}/${tumor}_vs_${normal}", mode: 'copy'
    publishDir "${params.outdir2}/annotation/freebayes/${tumor}_vs_${normal}", mode: 'copy'
    label 's_cpu'
    label 'xl_mem'

    input:
    tuple val(sample_id), val(tumor), val(normal),
        path(strelka_snv), path(strelka_snv_tbi),
        path(strelka_indel), path(strelka_indel_tbi),
        path(mutect_vcf), path(mutect_tbi),
        path(freebayes_vcf), path(freebayes_tbi)
    
    output:
    tuple val(sample_id), val(tumor), val(normal), path("*.ensemble.vcf.gz"), path("*.ensemble.vcf.gz.tbi"), emit: combined

    script:
    def out_vcf = "${tumor}_vs_${normal}.ensemble.vcf.gz"
    """
    java -jar /usr/GenomeAnalysisTK.jar \
        -T CombineVariants \
        --variant:mutect ${mutect_vcf} \
        --variant:freebayes ${freebayes_vcf} \
        --variant:strelka_snv ${strelka_snv} \
        --variant:strelka_indel ${strelka_indel} \
        -genotypeMergeOptions PRIORITIZE \
        -priority mutect,freebayes,strelka_snv,strelka_indel \
        -R ${params.fasta} \
        -o ${out_vcf}
    """
}

process VCF2TABLE {
    tag { sample_id }
    container 'patricie/vcf2table-image:v1.3'
    publishDir "${params.outdir}/${tumor}_vs_${normal}", mode: 'copy'
    // publishDir "${params.outdir2}/annotation/freebayes/${tumor}_vs_${normal}", mode: 'copy'
    label 'm_cpu'
    label 'xxl_mem'

    input:
    tuple val(sample_id), val(tumor), val(normal), path(filt_vcf), path(filt_tbi)

    output:
    tuple val(sample_id), val(tumor), val(normal), path("*.csv"), emit: table

    script:
    def out_csv = "${tumor}_vs_${normal}.ensemble_norm_filt_VEP_ann.csv"
    """
    python ${params.vcf2table} genome \
        -report_normal \
        --build hg38 \
        --parallel_proc ${task.cpus} \
        --log_time \
        --input ${filt_vcf} \
        --tumor ${tumor} > ${out_csv}
    """
}

// ============================================
// WORKFLOW
// ============================================

workflow {

    // --- 1. Collect files ---
    vcfs_ch = channel.fromPath('/storage2/sarek_ps3/Results_VEP_05_14_2026/annotation/freebayes/*/*_vs_*.filtered.norm.sorted_VEP.ann.vcf.gz')
        .map { fb_vcf ->
            def info = parseSampleInfo(fb_vcf)
            def sample_id = info[0]
            def tumor = info[1]
            def normal = info[2]

            def paths = deriveFilePaths(fb_vcf)
            def mutect_vcf = paths[0]
            def strelka_snv = paths[1]
            def strelka_indel = paths[2]

            if (!file(mutect_vcf).exists()) {
                error "Mutect2 file not found: ${mutect_vcf}"
            }
            if (!file(strelka_snv).exists()) {
                error "Strelka SNV file not found: ${strelka_snv}"
            }
            if (!file(strelka_indel).exists()) {
                error "Strelka Indel file not found: ${strelka_indel}"
            }

            return [sample_id, tumor, normal, fb_vcf, file(mutect_vcf), file(strelka_snv), file(strelka_indel)]
        }
    // vcfs_ch.view { "DISCOVERED: $it" }
    freebayes_filter = vcfs_ch.map { sample_id, tumor, normal, fb_vcf, mt_vcf, sk_snv, sk_ind  -> [sample_id, tumor, normal, fb_vcf] }
    FREEBAYES_BCBIO(freebayes_filter) 

    // // --- 2. Reheader vcfs ---
    vcfs_filtered = vcfs_ch.map { sample_id, tumor, normal, fb_vcf, mt_vcf, sk_snv, sk_ind  -> [sample_id, tumor, normal, mt_vcf, sk_snv, sk_ind] }
    // vcfs_filtered.view()
    // FREEBAYES_BCBIO.out.freebayes_bcbio.view()
    vcfs_reheader = FREEBAYES_BCBIO.out.freebayes_bcbio.join(vcfs_filtered, by: [0,1,2]) // group by sample_id, tumor, normal 
    // vcfs_reheader.view()

    REHEADER_VCFs(vcfs_reheader)
    // REHEADER_VCFs.out.reheadered.view()

    // --- 3. Combine Variants ---

    COMBINE_VARIANTS(REHEADER_VCFs.out.reheadered)

    // --- 4. VCF to CSV ---
    VCF2TABLE(COMBINE_VARIANTS.out.combined)

}