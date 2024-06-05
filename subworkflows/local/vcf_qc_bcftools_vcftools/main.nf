include { VCFTOOLS as VARIANT_FILTER     } from '../../../modules/nf-core/vcftools/main'
include { BCFTOOLS_STATS                  } from '../../../modules/nf-core/bcftools/stats/main'
include { VCFTOOLS as VCFTOOLS_SUMMARY    } from '../../../modules/nf-core/vcftools/main'
include { VCFTOOLS as VCFTOOLS_TSTV_COUNT } from '../../../modules/nf-core/vcftools/main'
include { VCFTOOLS as VCFTOOLS_TSTV_QUAL  } from '../../../modules/nf-core/vcftools/main'

workflow VCF_QC_BCFTOOLS_VCFTOOLS {
    take:
    vcf
    target_bed
    exclude_variant_intervals

    main:

    versions = Channel.empty()
    VARIANT_FILTER(vcf, exclude_variant_intervals, [])
    BCFTOOLS_STATS(VARIANT_FILTER.out.vcf.map{ meta, vcf -> [ meta, vcf, [] ] }, [[:],[]], [[:],[]], [[:],[]], [[:],[]], [[:],[]])
    VCFTOOLS_TSTV_COUNT(VARIANT_FILTER.out.vcf, target_bed, [])
    VCFTOOLS_TSTV_QUAL(VARIANT_FILTER.out.vcf, target_bed, [])
    VCFTOOLS_SUMMARY(VARIANT_FILTER.out.vcf, target_bed, [])

    versions = versions.mix(VARIANT_FILTER.out.versions)
    versions = versions.mix(BCFTOOLS_STATS.out.versions)
    versions = versions.mix(VCFTOOLS_TSTV_COUNT.out.versions)

    emit:
    variant_filter_vcf      = VARIANT_FILTER.out.vcf
    bcftools_stats          = BCFTOOLS_STATS.out.stats
    vcftools_tstv_counts    = VCFTOOLS_TSTV_COUNT.out.tstv_count
    vcftools_tstv_qual      = VCFTOOLS_TSTV_QUAL.out.tstv_qual
    variant_filter_summary  = VCFTOOLS_SUMMARY.out.filter_summary

    versions
}
