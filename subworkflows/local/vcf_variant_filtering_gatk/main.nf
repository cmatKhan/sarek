include { GATK4_VARIANTFILTRATION                                       } from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { GATK4_SELECTVARIANTS                                          } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_CNNSCOREVARIANTS      as CNNSCOREVARIANTS               } from '../../../modules/nf-core/gatk4/cnnscorevariants/main'
include { GATK4_FILTERVARIANTTRANCHES as FILTERVARIANTTRANCHES          } from '../../../modules/nf-core/gatk4/filtervarianttranches/main'

workflow VCF_VARIANT_FILTERING_GATK {

    take:
    vcf             // meta, vcf, tbi, intervals
    fasta
    fasta_fai
    dict
    intervals_bed_combined
    known_sites
    known_sites_tbi

    main:

    versions = Channel.empty()

    if (params.gatk_hard_filter){

        GATK4_VARIANTFILTRATION(
            vcf,
            fasta.map{ it -> [ [:], it]},
            fasta_fai.map{ it -> [ [:], it]},
            dict.map{ it -> [ [:], it]},)

        GATK4_SELECTVARIANTS(
            GATK4_VARIANTFILTRATION.out.vcf
                .join(GATK4_VARIANTFILTRATION.out.tbi, failOnDuplicate: true, failOnMismatch: true)
                .map{ meta, vcf, tbi -> [ meta, vcf, tbi, [] ] })

        versions = versions.mix(GATK4_SELECTVARIANTS.out.versions)

        filtered_vcf = GATK4_SELECTVARIANTS.out.vcf
            // remove no longer necessary field: num_intervals
            .map{ meta, vcf -> [ meta - meta.subMap('num_intervals'), vcf ] }

    } else {

        // Don't scatter/gather by intervals, because especially for small regions (targeted or WGS), it easily fails with 0 SNPS in region
        cnn_in = vcf.combine(intervals_bed_combined).map{ meta, vcf, tbi, intervals -> [ meta, vcf, tbi, [], intervals ] }

        CNNSCOREVARIANTS(cnn_in, fasta, fasta_fai, dict, [], [])

        FILTERVARIANTTRANCHES(CNNSCOREVARIANTS.out.vcf.join(CNNSCOREVARIANTS.out.tbi, failOnDuplicate: true, failOnMismatch: true).combine(intervals_bed_combined), known_sites, known_sites_tbi, fasta, fasta_fai, dict)

        filtered_vcf = FILTERVARIANTTRANCHES.out.vcf
            // remove no longer necessary field: num_intervals
            .map{ meta, vcf -> [ meta - meta.subMap('num_intervals'), vcf ] }

        versions = versions.mix(CNNSCOREVARIANTS.out.versions)
        versions = versions.mix(FILTERVARIANTTRANCHES.out.versions)

    }


    emit:
    filtered_vcf

    versions
}
