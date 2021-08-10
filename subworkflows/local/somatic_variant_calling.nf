/*
========================================================================================
    SOMATIC VARIANT CALLING (PAIR and TUMOR ONLY)
========================================================================================
*/

params.manta_options                  = [:]
params.msisensorpro_msi_options       = [:]
params.strelka_options                = [:]
params.strelka_bp_options             = [:]
params.mutect2_somatic_options        = [:]
params.mutect2_filter_options         = [:]

include { MANTA_SOMATIC as MANTA }                       from '../../modules/nf-core/software/manta/somatic/main'                        addParams(options: params.manta_options)
include { MSISENSORPRO_MSI }                             from '../../modules/nf-core/software/msisensorpro/msi/main'                     addParams(options: params.msisensorpro_msi_options)
include { STRELKA_SOMATIC as STRELKA }                   from '../../modules/nf-core/software/strelka/somatic/main'                      addParams(options: params.strelka_options)
include { STRELKA_SOMATIC_BEST_PRACTICES as STRELKA_BP } from '../../modules/nf-core/software/strelka/somaticbp/main'                    addParams(options: params.strelka_bp_options)
include { MUTECT2_SOMATIC }                              from '../../modules/nf-core/software/gatk4/mutect2/somatic/main'                addParams(options: params.mutect2_somatic_options)
include { MUTECT2_MERGE_STATS }                          from '../../modules/nf-core/software/gatk4/mutect2/merge_stats/main'   
include { MUTECT2_PILEUP_SUMMARIES }                     from '../../modules/nf-core/software/gatk4/mutect2/merge_pileup_summaries/main'   
include { MUTECT2_MERGE_VCF }                            from '../../modules/nf-core/software/gatk4/mutect2/merge_vcf/main'   
include { MUTECT2_PILEUP_SUMMARIES }                     from '../../modules/nf-core/software/gatk4/mutect2/merge_pileup_summaries/main'   
include { MUTECT2_CONTAMINATION }                        from '../../modules/nf-core/software/gatk4/mutect2/contamination/main'   
include { MUTECT2_FILTER }                               from '../../modules/nf-core/software/gatk4/mutect2/filter/main'                 addParams(options: params.mutect2_filter_options)

workflow SOMATIC_VARIANT_CALLING {
    take:
        tools
        cram                     // channel: [mandatory] cram
        dbsnp                    // channel: [mandatory] dbsnp
        dbsnp_tbi                // channel: [mandatory] dbsnp_tbi
        dict                     // channel: [mandatory] dict
        fai                      // channel: [mandatory] fai
        fasta                    // channel: [mandatory] fasta
        intervals                // channel: [mandatory] intervals
        somatic_force_tumor_only // boolean
        msisensorpro_scan        // channel: [optional]  msisensorpro_scan
        target_bed               // channel: [optional]  target_bed
        target_bed_gz_tbi        // channel: [optional]  target_bed_gz_tbi
        germline_resource        // channel: [optional]  germline_resource
        germline_resource_tbi    // channel: [optional]  germline_resource_tbi
        panel_of_normals         // channel: [optional]  panel_of_normals
        panel_of_normals_tbi     // channel: [optional]  panel_of_normals_tbi

    main:

        // split the input into tumor and normal samples
        cram.map { meta, cram, crai ->
            patient = meta.patient
            sample  = meta.sample
            gender  = meta.gender
            status  = meta.status
            [patient, sample, gender, status, cram, crai]
        }.branch {
            normal: it[3] == 0
            tumor:  it[3] == 1
        }.set{ cram_to_cross }

        // create a channel with the paired tumor-normal samples
        cram_pair = cram_to_cross.normal.cross( cram_to_cross.tumor ).map { normal, tumor ->
            def meta = [:]
            meta.patient = normal[0]
            meta.normal  = normal[1]
            meta.tumor   = tumor[1]
            meta.gender  = normal[2]
            meta.id      = "${meta.tumor}_vs_${meta.normal}".toString()
            [meta, normal[4], normal[5], tumor[4], tumor[5]]
        }

        // create a channel with the tumor samples (all of them or just the ones without a paired normal)
        // TODO duplicated code here
        //cram_tumor = Channel.empty()
        //if (somatic_force_tumor_only) {
        //    cram_tumor = cram_to_cross.tumor.map { tumor -> 
        //        def meta = [:]
        //        meta.patient = tumor[0]
        //        meta.tumor   = tumor[1]
        //        meta.gender  = tumor[2]
        //        meta.id      = "${meta.tumor}".toString()
        //        [meta, tumor[4], tumor[5]]
        //    }
        //} else {
        //    cram_tumor = cram_to_cross.tumor.cross( cram_to_cross.normal, remainder=True ).filter {
        //        it[1] == null
        //    }.map { tumor, normal ->
        //        def meta = [:]
        //        meta.patient = tumor[0]
        //        meta.tumor   = tumor[1]
        //        meta.gender  = tumor[2]
        //        meta.id      = "${meta.tumor}".toString()
        //        [meta, tumor[4], tumor[5]]
        //    }
        //}

        manta_vcf                 = Channel.empty()
        msipro_results            = Channel.empty()
        strelka_vcf               = Channel.empty()
        mutect2_vcf               = Channel.empty()
        //msipro_tumor_only_results = Channel.empty()
        //strelka_tumor_only_vcf    = Channel.empty()
        //mutect2_tumor_only_vcf    = Channel.empty()

        if ('manta' in tools) {
            MANTA(
                cram_pair,
                fasta,
                fai,
                target_bed_gz_tbi
            )

            manta_candidate_small_indels_vcf = MANTA.out.candidate_small_indels_vcf
            manta_candidate_sv_vcf           = MANTA.out.candidate_sv_vcf
            manta_diploid_sv_vcf             = MANTA.out.diploid_sv_vcf
            manta_somatic_sv_vcf             = MANTA.out.somatic_sv_vcf
            manta_csi_for_strelka_bp         = MANTA.out.manta_csi_for_strelka_bp
            
            manta_vcf = manta_vcf.mix(manta_candidate_small_indels_vcf, manta_candidate_sv_vcf, manta_diploid_sv_vcf, manta_somatic_sv_vcf)
        }

        if ('msisensorpro' in tools) {
            MSISENSORPRO_MSI(
                cram_pair,
                msisensorpro_scan
            )
            
            msi_dis_list      = MSISENSORPRO_MSI.out.dis_list
            msi_germline_list = MSISENSORPRO_MSI.out.germline_list
            msi_somatic_list  = MSISENSORPRO_MSI.out.somatic_list

            msipro_results = msipro_results.mix(msi_dis_list, msi_germline_list, msi_somatic_list)

            // TODO add tumor-only

        }

        if ('strelka' in tools) {
            STRELKA(
                cram_pair,
                fasta,
                fai,
                target_bed_gz_tbi
            )

            strelka_indels_vcf = STRELKA.out.indels_vcf
            strelka_snvs_vcf   = STRELKA.out.snvs_vcf

            strelka_vcf = strelka_vcf.mix(strelka_indels_vcf, strelka_snvs_vcf)

            // TODO add tumor-only
        }

        if ('mutect2' in tools) {
            panel_of_normals.dump()
            no_intervals = intervals == [] 

            // add intervals to the channels for Mutect2
            cram_pair.combine(intervals).map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals ->
                    new_meta = meta.clone()
                    new_meta.id = meta.id + "_" + intervals.baseName
                    [new_meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals]
            }.set{cram_pair_intervals}

            // TODO make the Mutect2 variant calling a subworkflow

            // 1. Mutect2 variant calling 
            MUTECT2(
                cram_pair_intervals,
                panel_of_normals,
                panel_of_normals_tbi,
                dict,
                fasta,
                fai,
                no_intervals,
                germline_resource,
                germline_resource_tbi
            )

            // 2. Merge stats from Mutect2 variant calling
            MUTECT2_MERGE_STATS(
                MUTECT2.out.stats, 
                dict,
                fasta,
                fai,
                germline_resource,
                germline_resource_tbi,
                intervals
            )

            // 3. Merge VCFs from Mutect2 variant calling
            MUTECT2_MERGE_VCF(
                MUTECT2.out.vcf,
                fai,
                target_bed
            )

            // 4. Pileup summaries for Mutect2 
            MUTECT2_PILEUP_SUMMARIES(
                cram_pair_intervals,
                germline_resource,
                germline_resource_tbi
            )

            pileupSummaries = MUTECT2_PILEUP_SUMMARIES.out.pileupSummaries.groupTuple(by:[0,1])

            // 5. Merge pileup summaries
            MUTECT2_MERGE_PILEUP_SUMMARIES(
                pileupSummaries,
                dict
            )

            // 6. Calculate Contamination
            MUTECT2_CONTAMINATION(
                MUTECT2_MERGE_PILEUP_SUMMARIES.out.mergedPileupFile
            )

            mutect2CallsToFilter = MUTECT2_MERGE_VCF.out.vcfConcatenatedForFilter.map{
                variantCaller, idPatient, idSample, vcf, tbi ->
                [idPatient, idSample, vcf, tbi]
            }.join(MUTECT2_MERGE_STATS.out.mergedStatsFile, by:[0,1]).join(MUTECT2_CONTAMINATION.out.contaminationTable, by:[0,1])

            // 7. Filter Mutect2 calls 
            MUTECT2_FILTER(
                mutect2CallsToFilter
            ) 

            mutect2_vcf = MUTECT2_FILTER.filtered_vcf

        }

    emit:
        manta_vcf                 = manta_vcf
        msipro_results            = msipro_results
        strelka_vcf               = strelka_vcf
        mutect2_vcf               = mutect2_vcf
        //strelka_tumor_only_vcf    = strelka_tumor_only_vcf
        //msipro_tumor_only_results = msipro_tumor_only_results
        //mutect2_tumor_only_vcf    = mutect2_tumor_only_vcf

}
