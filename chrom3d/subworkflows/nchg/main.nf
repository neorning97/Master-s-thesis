// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

// Workaround for optional input files: https://github.com/nextflow-io/nextflow/issues/1694
def make_optional_input(path) {
  if (path?.trim()) {
    return [file(path)]
  }
  return []
}

workflow NCHG {

  take:
    sample_sheet

    // masking
    mad_max
    bad_bin_fraction
    cytoband
    gaps

    // cis
    use_cis_interactions
    fdr_cis
    log_ratio_cis

    // trans
    use_trans_interactions
    fdr_trans
    log_ratio_trans

    // plotting
    skip_expected_plots
    skip_sign_interaction_plots
    hic_tgt_resolution_plots
    plot_sig_interactions_cmap_lb
    plot_sig_interactions_cmap_ub
    plot_format

    // misc
    zstd_compression_lvl

  main:

    interaction_types = []
    if (use_cis_interactions) {
      interaction_types.add("cis")
    }
    if (use_trans_interactions) {
      interaction_types.add("trans")
    }

    sample_sheet
      .splitCsv(sep: "\t", header: true)
      .map { row -> tuple(row.sample,
                          file(row.hic_file, checkIfExists: true),
                          row.resolution)
      }
      .set { hic_files }

    sample_sheet
      .splitCsv(sep: "\t", header: true)
      .map { row -> tuple(row.sample, make_optional_input(row.domains))
      }
      .set { domains }

    sample_sheet
      .splitCsv(sep: "\t", header: true)
      .map { row -> tuple(row.sample, make_optional_input(row.mask))
      }
      .set { bin_masks }

    GENERATE_MASK(
      bin_masks,
      make_optional_input(cytoband),
      make_optional_input(gaps),
      zstd_compression_lvl
    )

    hic_files
      .combine(GENERATE_MASK.out.bed, by: 0)
      .set { nchg_expected_tasks }

    EXPECTED(
      nchg_expected_tasks,
      interaction_types,
      mad_max
    )

    DUMP_CHROM_SIZES(
      hic_files
    )

    PREPROCESS_DOMAINS(
      domains
    )

    CARTESIAN_PRODUCT(
      PREPROCESS_DOMAINS.out.domains
        .combine(DUMP_CHROM_SIZES.out.tsv, by: 0)
    )

    CARTESIAN_PRODUCT.out.domains
      .map { sample, file ->
        if (file.extension == 'skip') {
          file = []
        }
        tuple(sample, file)
      }
      .set { valid_domains }

    GENERATE_CHROMOSOME_PAIRS(
      hic_files
        .combine(valid_domains, by: 0)
        .combine(interaction_types)
    )

    GENERATE_CHROMOSOME_PAIRS.out
      .splitCsv(sep: "\t", header: ["sample", "chrom1", "chrom2"])
      .map { tuple(it.sample, it.chrom1, it.chrom2) }
      .combine(hic_files, by: 0)
      .combine(valid_domains, by: 0)
      .combine(EXPECTED.out.h5, by: 0)
      .set { nchg_compute_tasks }

    COMPUTE(
      nchg_compute_tasks,
      bad_bin_fraction,
    )

    COMPUTE.out.parquet
      .branch {
        cis: it[1] == it[2]
        trans: true
      }
      .set { significant_interactions }

    significant_interactions.cis
      .map { tuple(it[0], it[3]) }
      .groupTuple()
      .combine(["cis"])
      .join(DUMP_CHROM_SIZES.out.tsv)
      .set { nchg_merge_cis_tasks }

    significant_interactions.trans
      .map { tuple(it[0], it[3]) }
      .groupTuple()
      .combine(["trans"])
      .join(DUMP_CHROM_SIZES.out.tsv)
      .set { nchg_merge_trans_tasks }

    nchg_merge_cis_tasks
      .mix(nchg_merge_trans_tasks)
      .set { nchg_merge_tasks }

    MERGE(
      nchg_merge_tasks
    )


    MERGE.out.parquet
      .branch {
        cis: it[1] == "cis"
        trans: true
      }
      .set { nchg_merge_output }

    nchg_merge_output.cis
      .map {
        tuple(it[0],
              it[1],
              it[2],
              fdr_cis,
              log_ratio_cis
        )
      }
      .set { nchg_filter_cis_tasks }

    nchg_merge_output.trans
      .map {
        tuple(it[0],
              it[1],
              it[2],
              fdr_trans,
              log_ratio_trans
        )
      }
      .set { nchg_filter_trans_tasks }

    nchg_filter_cis_tasks
      .mix(nchg_filter_trans_tasks)
      .set { nchg_filter_tasks }

    FILTER(
      nchg_filter_tasks
    )

    MERGE.out.parquet
      .groupTuple()
      .map { tuple(it[0], 'unfiltered', it[2]) }
      .set { nchg_unfiltered_concat_tasks }

    FILTER.out.parquet
      .groupTuple()
      .map { tuple(it[0], 'filtered', it[2]) }
      .set { nchg_filtered_concat_tasks }

    nchg_unfiltered_concat_tasks
      .mix(nchg_filtered_concat_tasks)
      .set { nchg_concat_tasks }

    CONCAT(
      nchg_concat_tasks
    )

    CONCAT.out.parquet
      .filter { it[1] == 'filtered' }
      .map { tuple(it[0], it[2]) }
      .set { nchg_view_tasks }

    VIEW(
      nchg_view_tasks
    )

    if (!skip_expected_plots) {
      PLOT_EXPECTED(
        EXPECTED.out.h5,
        interaction_types,
        plot_format
      )
    }

    if (!skip_sign_interaction_plots) {
      GET_HIC_PLOT_RESOLUTION(
        hic_files,
        hic_tgt_resolution_plots
      )

      GET_HIC_PLOT_RESOLUTION.out.tsv
        .splitCsv(header: ["sample", "resolution"],
                  sep: "\t")
        .map { tuple(it.sample, it.resolution) }
        .set { plotting_resolutions }

      hic_files
        .map { tuple(it[0], it[1]) }
        .join(plotting_resolutions)
        .join(VIEW.out.tsv)
        .set { plotting_tasks }

      plot_sig_interactions_cmap_lb = plot_sig_interactions_cmap_lb
      if (!plot_sig_interactions_cmap_lb) {
        plot_sig_interactions_cmap_lb = Math.min(log_ratio_cis,
                                                 log_ratio_trans)
      }

      plot_sig_interactions_cmap_ub = Math.max(plot_sig_interactions_cmap_lb,
                                               plot_sig_interactions_cmap_ub)

      PLOT_SIGNIFICANT(
        plotting_tasks,
        plot_sig_interactions_cmap_lb,
        plot_sig_interactions_cmap_ub,
        plot_format
      )

      PLOT_EXPECTED.out.plots
        .mix(PLOT_SIGNIFICANT.out.plots)
        .set { PLOTS }
    } else {
      Channel.empty().set { PLOTS }
    }

  emit:
    expected = EXPECTED.out.h5
    parquets = CONCAT.out.parquet
    plots = PLOTS
    tsv = VIEW.out.tsv

}

process GENERATE_MASK {
  label 'duration_very_short'

  tag "$sample"

  cpus 1

  input:
    tuple val(sample),
          path(mask)
    path cytoband
    path gaps

    val zstd_compression_lvl

  output:
    tuple val(sample),
          path("*.bed.gz"),
    emit: bed

  script:
    opts=[]
    if (!mask.toString().isEmpty()) {
      opts.push(mask)
    }
    if (!gaps.toString().isEmpty()) {
      opts.push(gaps)
    }
    skip=opts.size() == 0
    if (!cytoband.toString().isEmpty()) {
      opts.push("--cytoband='${cytoband}'")
    }

    opts=opts.join(" ")

    if (skip) {
      """
      echo "" | gzip -9 > __'$sample'.mask.bed.gz
      """
    } else {
      """
      set -o pipefail

      nchg_nf_generate_bin_mask.py $opts | gzip -9 > __'$sample'.mask.bed.gz
      """
    }
}

process DUMP_CHROM_SIZES {
  label 'process_very_short'
  tag "$sample"

  input:
    tuple val(sample),
          path(hic),
          val(resolution)

  output:
    tuple val(sample),
          path("*.chrom.sizes"),
    emit: tsv

  script:
    outname="${sample}.chrom.sizes"
    """
    hictk dump -t chroms '$hic' --resolution '$resolution' > '$outname'
    """
}

process PREPROCESS_DOMAINS {
  label 'process_short'
  tag "$sample"

  input:
    tuple val(sample),
          path(domains)

  output:
    tuple val(sample),
          path("*.{zst,skip}"),
    emit: domains

  script:
    outprefix="${sample}.domains"
    """
    if [ '$domains' = '' ]; then
      touch '$outprefix'.skip
      exit 0
    fi

    tmpfile='$outprefix'.tmp

    set -o pipefail
    nchg_nf_preprocess_domains.py '$domains' | zstd -13 -o "\$tmpfile"
    set +o pipefail

    num_cols="\$(zstdcat "\$tmpfile" | head -n 1 | wc -w)"

    if [ "\$num_cols" -eq 6 ]; then
      mv "\$tmpfile" '$outprefix'.bedpe.zst
    else
      mv "\$tmpfile" '$outprefix'.bed.zst
    fi
    """
}

process CARTESIAN_PRODUCT {
  label 'process_short'
  tag "$sample"

  input:
    tuple val(sample),
          path(domains),
          path(chrom_sizes)

  output:
    tuple val(sample),
          path("*.{bedpe.zst,skip}"),
    emit: domains

  script:
    outprefix="${sample}.domains.ok"
    """
    set -o pipefail

    if [[ '$domains' == *.skip ]]; then
      touch '$outprefix'.skip
      exit 0
    fi

    NCHG cartesian-product \\
      --chrom-sizes '$chrom_sizes' \\
      <(zstd -dcf '$domains') |
      zstd -13 -o '$outprefix'.bedpe.zst
    """
}

process GENERATE_CHROMOSOME_PAIRS {
  label 'process_very_short'
  tag "$sample ($interaction_type)"

  input:
    tuple val(sample),
          path(hic),
          val(resolution),
          path(domains),
          val(interaction_type)

  output:
    stdout emit: tsv

  script:
    opts=[
      "--resolution='${resolution}'",
      "--interaction-type='${interaction_type}'"
    ]

    if (domains.size() != 0) {
      opts.push("--domains='${domains}'")
    }

    opts=opts.join(" ")
    """
    nchg_nf_generate_chromosome_pairs.py \\
      '$sample' \\
      '$hic' \\
      $opts
    """
}

// TODO optimize: trans expected values can be computed in parallel
process EXPECTED {
  tag "$sample"

  input:
    tuple val(sample),
          path(hic),
          val(resolution),
          path(bin_mask)

    val interaction_types
    val mad_max

  output:
    tuple val(sample),
          path(outname),
    emit: h5

  script:
    opts=[
      "--resolution='${resolution}'",
      "--mad-max='${mad_max}'"
    ]

    suffix=""
    if (interaction_types.size() == 1) {
      if ("cis" in interaction_types) {
        opts.push("--cis-only")
        suffix=".cis"
      } else if ("trans" in interaction_types) {
        opts.push("--trans-only")
        suffix=".trans"
      }
    }

    outname="expected_values_${sample}${suffix}.h5"
    opts=opts.join(" ")

    """
    trap 'rm -rf ./tmp/' EXIT

    mkdir tmp/
    bin_mask_plain="\$(mktemp -p ./tmp)"

    zcat -f '$bin_mask' > "\$bin_mask_plain"

    NCHG expected \\
      '$hic' \\
      --output='$outname' \\
      --bin-mask="\$bin_mask_plain" \\
      $opts
    """
}

process COMPUTE{
  tag "$sample [$chrom1:$chrom2]"


  input:
    tuple val(sample),
          val(chrom1),
          val(chrom2),
          path(hic),
          val(resolution),
          path(domains),
          path(expected_values)

    val bad_bin_fraction

  output:
    tuple val(sample),
          val(chrom1),
          val(chrom2),
          path("*.parquet", optional: true),
    emit: parquet

  script:
    outname="${sample}.${chrom1}.${chrom2}.parquet"

    opts=[]

    if (!domains.toString().isEmpty()) {
      opts.push("--domains='domains.bed'")
    }

    opts=opts.join(" ")
    """
    set -o pipefail

    if [ -n '$domains' ]; then
      zstdcat -f '$domains' > domains.bed
    fi

    NCHG compute \\
      '$hic' \\
      out/'$outname' \\
      --resolution='$resolution' \\
      --chrom1='$chrom1' \\
      --chrom2='$chrom2' \\
      --expected-values='$expected_values' \\
      --bad-bin-fraction='$bad_bin_fraction' \\
      $opts

    mv out/'$outname' '$outname'
    """
}

process MERGE {
  tag "$sample ($interaction_type)"

  cpus 2

  input:
    tuple val(sample),
          path(parquets),
          val(interaction_type),
          path(chrom_sizes)

  output:
    tuple val(sample),
          val(interaction_type),
          path(outname),
    emit: parquet

  script:
    input_prefix="${sample}"
    outname="${sample}.${interaction_type}.parquet"
    """
    NCHG merge --input-prefix='$input_prefix' \\
      --output '$outname' \\
      --ignore-report-file \\
      --threads='${task.cpus}'
    """
}

process FILTER {
  tag "$sample ($interaction_type)"

  cpus 2

  input:
    tuple val(sample),
          val(interaction_type),
          path(parquet),
          val(fdr),
          val(log_ratio)

  output:
    tuple val(sample),
          val(interaction_type),
          path(outname),
    emit: parquet

  script:
    outname="${sample}.${interaction_type}.filtered.parquet"
    """
    NCHG filter \\
      '$parquet' \\
      '$outname' \\
      --fdr='$fdr' \\
      --log-ratio='$log_ratio' \\
      --threads='${task.cpus}'
    """
}

process CONCAT {
  label 'process_medium'
  tag "$sample"

  input:
    tuple val(sample),
          val(type),
          path(parquets)

  output:
    tuple val(sample),
          val(type),
          path("*.parquet"),
    emit: parquet

  script:
    outname=(type == 'filtered') ?
        "${sample}.filtered.parquet" :
        "${sample}.parquet"
    """
    mkdir tmp
    trap 'rm -rf ./tmp/' EXIT

    TMPDIR=./tmp/ \\
    NCHG merge \\
      --input-files *.parquet \\
      --output '$outname' \\
      --ignore-report-file \\
      --threads='${task.cpus}'
    """
}

process VIEW {
  label 'process_medium'
  tag "$sample"

  input:
    tuple val(sample),
          path(parquet)

  output:
    tuple val(sample),
          path("*.tsv.gz"),
    emit: tsv

  script:
    tsv="${parquet.baseName}.tsv.gz"
    """
    NCHG view '$parquet' |
      pigz -9 -p ${task.cpus} > '$tsv'
    """
}

process PLOT_EXPECTED {
  tag "$sample"

  input:
    tuple val(sample),
          path(h5)

    val interaction_types
    val plot_format

  output:
    tuple val(sample),
          path("*.${plot_format}"),
    emit: plots

  script:
    plot_cis="cis" in interaction_types
    plot_trans="trans" in interaction_types
    outname_cis="${sample}_cis.${plot_format}"
    outname_trans="${sample}_trans.${plot_format}"
    """
    if [[ '$plot_cis' == true ]]; then
      nchg_nf_plot_expected_values.py \\
        '$h5' \\
        '$outname_cis' \\
        --yscale-log \\
        --plot-cis
    fi

    if [[ '$plot_trans' == true ]]; then
      nchg_nf_plot_expected_values.py \\
        '$h5' \\
        '$outname_trans' \\
        --plot-trans
    fi
    """
}

process GET_HIC_PLOT_RESOLUTION {
  label 'process_very_short'
  tag "$sample"

  input:
    tuple val(sample),
          path(hic),
          val(resolution)

    val tgt_resolution

  output:
    stdout emit: tsv

  script:
    """
    #!/usr/bin/env python3

    import hictkpy

    best_res = int("$resolution")
    tgt_res = int("$tgt_resolution")
    sample = "$sample"

    try:
        resolutions = hictkpy.MultiResFile("$hic").resolutions()

        for res in resolutions:
            delta1 = abs(res - tgt_res)
            delta2 = abs(best_res - tgt_res)

            if delta1 < delta2:
                best_res = res

    except RuntimeError:
        pass
    finally:
        print(f"{sample}\\t{best_res}")
    """
}

process PLOT_SIGNIFICANT {
  label 'process_very_high'
  tag "$sample"

  input:
    tuple val(sample),
          path(hic),
          val(resolution),
          path(tsv)

    val cmap_lb
    val cmap_ub
    val plot_format

  output:
    tuple val(sample),
          path("*.${plot_format}"),
    emit: plots

  script:
    """
    export MPLCONFIGDIR=./mpl
    trap "rm -rf '\$MPLCONFIGDIR'" EXIT
    mkdir "\$MPLCONFIGDIR"

    nchg_nf_plot_significant_interactions.py \\
      '$hic' \\
      '$tsv' \\
      '$sample' \\
      --resolution '$resolution' \\
      --plot-format '$plot_format' \\
      --min-value '$cmap_lb' \\
      --max-value '$cmap_ub' \\
      --nproc '${task.cpus}'
    """
}
