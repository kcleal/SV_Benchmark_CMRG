#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


// Process to set up conda environment and install tools
process SETUP {
    
    publishDir params.outdir, mode: 'symlink'
    conda "$projectDir/environment.yml"
    
    output:
        path "versions.json", emit: versions
        path "sawfish", emit: sawfish
        path "delly", emit: delly
    
    script:
    """
 
    # Install Sniffles
    wget https://github.com/fritzsedlazeck/Sniffles/archive/refs/tags/v2.4.tar.gz
    tar -xvf v2.4.tar.gz && rm v2.4.tar.gz
    cd Sniffles && python setup.py install && cd ../
    
    # Install Sawfish
    wget https://github.com/PacificBiosciences/sawfish/releases/download/v0.12.7/sawfish-v0.12.7-x86_64-unknown-linux-gnu.tar.gz
    tar -xvf sawfish-v0.12.7-x86_64-unknown-linux-gnu.tar.gz
    mv sawfish-v0.12.7-x86_64-unknown-linux-gnu/bin/sawfish .

    # Install Delly
    wget https://github.com/dellytools/delly/releases/download/v1.3.1/delly_v1.3.1_linux_x86_64bit
    mv delly_v1.3.1_linux_x86_64bit delly && chmod +x delly
    
    # Install NGSEP
    wget https://github.com/NGSEP/NGSEPcore/releases/download/v5.0.0/NGSEPcore_5.0.0.jar
   
    # Capture versions
    dysgu_ver=\$(dysgu --version | cut -f3 -d ' ')
    sniffles_ver=\$(sniffles --version | cut -f3 -d ' ')
    cutesv_ver=\$(cuteSV --version | cut -f2 -d ' ')
    severus_ver=\$(severus --version)
    sawfish_ver=\$(./sawfish --version | cut -f2 -d ' ')
    
    # Write versions to JSON
    cat <<- EOF > versions.json
    {
        "dysgu": "\$dysgu_ver",
        "sniffles": "\$sniffles_ver",
        "cutesv": "\$cutesv_ver",
        "severus": "\$severus_ver",
        "sawfish": "\$sawfish_ver",
        "delly": "1.3.1"
    }
EOF
    """

}


process RUN_CUTESV {
    tag "${data_type}_cutesv" 
    publishDir "${params.outdir}/calls", mode: 'copy'

    input:
        path ref
        path input_bam
        path idx
        val versions
        val data_type

    output:
        path "HG002.${data_type}.cuteSV_*.vcf", emit: vcf

    script:
    def prefix = "HG002.${data_type}"
    """
    mkdir -p wd_cuteSV
    cuteSV -t ${task.cpus} -s 3 --genotype ${input_bam} ${ref} ${prefix}.cuteSV_${versions.cutesv}.vcf wd_cuteSV
    """
}


process RUN_SEVERUS {
    tag "${data_type}_severus"
    publishDir "${params.outdir}/calls", mode: 'copy'

    input:
        path ref
        path input_bam
        path idx
        val versions
        val data_type
        path convert_script

    output:
        path "HG002.${data_type}.severus_*.vcf", emit: vcf

    script:
    def prefix = "HG002.${data_type}"
    """
    curl -O https://raw.githubusercontent.com/KolmogorovLab/Severus/main/vntrs/human_GRCh38_no_alt_analysis_set.trf.bed
    
    severus --target-bam ${input_bam} --out-dir severus_out_${data_type} -t ${task.cpus} \
        --vntr-bed ./human_GRCh38_no_alt_analysis_set.trf.bed
    cat severus_out_${data_type}/all_SVs/severus_all.vcf | python ${convert_script} | \
        bcftools sort -Ov - > ${prefix}.severus_${versions.severus}.vcf
    """
}


process RUN_SNIFFLES {
    tag "${data_type}_sniffles" 
    publishDir "${params.outdir}/calls", mode: 'copy'

    input:
        path ref
        path input_bam
        path idx
        val versions
        val data_type // 'ont' or 'pacbio'

    output:
        path "HG002.${data_type}.sniffles_*.vcf", emit: vcf

    script:
    def prefix = "HG002.${data_type}"
    """
    sniffles -t ${task.cpus} --input ${input_bam} --vcf ${prefix}.sniffles_${versions.sniffles}.vcf
    """
}


process RUN_DYSGU {
    tag "${data_type}_dysgu"
    publishDir "${params.outdir}/calls", mode: 'copy'

    input:
        path ref
        path input_bam
        path idx
        val versions
        val data_type

    output:
        path "HG002.${data_type}.dysgu_*.vcf", emit: vcf

    script:
    def prefix = "HG002.${data_type}"
    def mode = params.caller_params[data_type].dysgu.mode
    """
    if [[ "${input_bam}" == *.cram ]]; then
        dysgu run --mode ${mode} --procs ${task.cpus} -x --clean ${ref} wd ${input_bam} > ${prefix}.dysgu_${versions.dysgu}.vcf
    else
        dysgu call --mode ${mode} --procs ${task.cpus} -x --clean ${ref} wd ${input_bam} > ${prefix}.dysgu_${versions.dysgu}.vcf
    fi 

    """ 
}

process RUN_DELLY {
    publishDir "${params.outdir}/calls", mode: 'copy'

    input:
        path ref
        path input_bam
        path idx
        val versions
        val data_type
        path delly_exe

    output:
        path "HG002.${data_type}.delly_*.vcf", emit: vcf

    script:
    def prefix = "HG002.${data_type}"
    def tech = params.caller_params[data_type].delly.tech
    """
    ./${delly_exe} lr --technology ${tech} -g ${ref} ${input_bam} > ${prefix}.delly_${versions.delly}.vcf
    """
}


process RUN_SAWFISH {
    publishDir "${params.outdir}/calls", mode: 'copy'

    input:
        path ref
        path pacbio_bam
        path idx
        val versions
        path sawfish_exe

    output:
        path "HG002.pacbio.sawfish_*.vcf", emit: vcf

    script:
    """
    ./${sawfish_exe} discover --threads ${task.cpus} --ref ${ref} --bam ${pacbio_bam} --output-dir sawfish_dir
    ./${sawfish_exe} joint-call --threads ${task.cpus} --sample sawfish_dir --output-dir sawfish_call_dir
    gunzip -c sawfish_call_dir/genotyped.sv.vcf.gz > HG002.pacbio.sawfish_${versions.sawfish}.vcf
    """
}


workflow run_ont_callers {

    data_type = 'ont'

    take:
        ref
        input_bam
        idx
        versions
        delly_exe
        convert_severus_script

    main:
        sniffles_res = RUN_SNIFFLES(ref, input_bam, idx, versions, data_type)
        dysgu_res = RUN_DYSGU(ref, input_bam, idx, versions, data_type)
        cutesv_res = RUN_CUTESV(ref, input_bam, idx, versions, data_type)
        severus_res = RUN_SEVERUS(ref, input_bam, idx, versions, data_type, convert_severus_script)
        delly_res = RUN_DELLY(ref, input_bam, idx, versions, data_type, delly_exe)
        
        ont_vcfs = sniffles_res.vcf.mix(
            dysgu_res.vcf,
            cutesv_res.vcf,
            severus_res.vcf,
            delly_res.vcf
        )

    emit:
        ont_vcfs
}


workflow run_pacbio_callers {

    data_type = 'pacbio'

    take:
        ref
        input_bam
        idx
        versions
        delly_exe
        sawfish_exe
        convert_severus_script

    main:
        sniffles_res = RUN_SNIFFLES(ref, input_bam, idx, versions, data_type)
        dysgu_res = RUN_DYSGU(ref, input_bam, idx, versions, data_type)
        cutesv_res = RUN_CUTESV(ref, input_bam, idx, versions, data_type)
        severus_res = RUN_SEVERUS(ref, input_bam, idx, versions, data_type, convert_severus_script)
        delly_res = RUN_DELLY(ref, input_bam, idx, versions, data_type, delly_exe)
        sawfish_res = RUN_SAWFISH(ref, input_bam, idx, versions, sawfish_exe)

        pacbio_vcfs = sniffles_res.vcf.mix(
            dysgu_res.vcf,
            cutesv_res.vcf,
            severus_res.vcf,
            delly_res.vcf,
            sawfish_res.vcf,
        )

    emit:
        pacbio_vcfs
}


process BENCHMARK {
    publishDir params.outdir, mode: 'copy'

    input:
        path ref
        path vcfs
        path cmrg_truth
        path cmrg_idx
        path cmrg_bed
        path giab_truth
        path giab_idx
        path giab_bed
        val versions
        val bench_params

    output:
        path "truvari_*"

    script:
    """
    #!/bin/bash
    echo "Benchmark parameters: ${bench_params}"

    # Iterate over each VCF file
    for vcf_file in ${vcfs}
    do
        base_name=\$(basename \${vcf_file} .vcf)

        tech=\$(echo "\${base_name}" | cut -d'.' -f2)
        caller_name=\$(echo "\${base_name}" | cut -d'.' -f3-)
        caller_name=\${caller_name%.vcf}
        
        # Skip ont for sawfish
        if [[ "\${caller_name}" == *"sawfish"* && "\${tech}" == "ont" ]]; then
            continue
        fi
           
        bgzip -f \${vcf_file}
        tabix -f \${vcf_file}.gz
        
        rm -rf truvari_CMRG_\${tech}_\${caller_name} 

        truvari bench -b ${cmrg_truth} \
                     --includebed ${cmrg_bed} \
                     -c \${vcf_file}.gz \
                     ${bench_params} \
                    -o truvari_CMRG_\${tech}_\${caller_name}

        truvari bench -b ${giab_truth} \
                    --includebed ${giab_bed} \
                    -c \${vcf_file}.gz \
                    ${bench_params} \
                    -o truvari_GIAB_\${tech}_\${caller_name}

        truvari refine -a wfa -R -U -f ${ref} \
                    --regions truvari_GIAB_\${tech}_\${caller_name}/candidate.refine.bed \
                    truvari_GIAB_\${tech}_\${caller_name}/
   
    done
    """
}


process PLOT_RESULTS {

    publishDir params.outdir, mode: 'copy'

    input:
        path plot_script
        path results
    
    output:
        path "*.png"
        path "*.md"

    script:
    """
    python $plot_script
    """
}


workflow {
    SETUP()

    delly_exe = SETUP.out.delly
    sawfish_exe = SETUP.out.sawfish
    convert_severus_script = Channel.fromPath('convert_severus.py')
    plot_script = Channel.fromPath('plot_benchmark.py')

    def versions = SETUP.out.versions.map { versions_file ->
        new groovy.json.JsonSlurper().parseText(versions_file.text)
    }

    // Define channels for each data file
    if (params.test) {
        ont_data = Channel.fromPath('chr1.ont_test.bam')
        ont_data_idx = Channel.fromPath('chr1.ont_test.bam.bai') 
        pacbio_data = Channel.fromPath('chr1.pacbio_test.bam')
        pacbio_data_idx = Channel.fromPath('chr1.pacbio_test.bam.bai') 
    } else {
        ont_data = Channel.fromPath('PAO89685.pass.cram')
        ont_data_idx = Channel.fromPath('PAO89685.pass.cram.crai') 
        pacbio_data = Channel.fromPath('m21009_241011_231051.GRCh38.haplotagged.bam')
        pacbio_data_idx = Channel.fromPath('m21009_241011_231051.GRCh38.haplotagged.bam.bai') 
    }

    ref = Channel.fromPath('GCA_000001405.15_GRCh38_no_alt_analysis_set.fna')
    cmrg_truth = Channel.fromPath('HG002_GRCh38_difficult_medical_gene_SV_benchmark_v0.01.vcf.gz')
    cmrg_idx = Channel.fromPath('HG002_GRCh38_difficult_medical_gene_SV_benchmark_v0.01.vcf.gz.tbi') 
    cmrg_bed = Channel.fromPath('HG002_GRCh38_difficult_medical_gene_SV_benchmark_v0.01.bed')
    giab_truth = Channel.fromPath('GRCh38_HG2-T2TQ100-V1.1.vcf.gz')
    giab_idx = Channel.fromPath('GRCh38_HG2-T2TQ100-V1.1.vcf.gz.tbi') 
    giab_bed = Channel.fromPath('GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed')
    
    ont_vcfs = run_ont_callers(ref, ont_data, ont_data_idx, versions, delly_exe, convert_severus_script)
    pacbio_vcfs = run_pacbio_callers(ref, pacbio_data, pacbio_data_idx, versions, delly_exe, sawfish_exe, convert_severus_script)
    all_vcfs = ont_vcfs.mix(pacbio_vcfs)


    results = BENCHMARK(ref, 
              all_vcfs.collect(),
              cmrg_truth, cmrg_idx, cmrg_bed,
              giab_truth, giab_idx, giab_bed,
              versions,
              params.bench_params)
    
    PLOT_RESULTS(plot_script, results)

}
