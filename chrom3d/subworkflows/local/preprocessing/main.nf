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

workflow PREPROCESSING {

  take:
    sample_sheet
    sig_interactions_cis
    sig_interactions_trans
    ploidy
    masked_chromosomes

  main:

    sample_sheet
      .splitCsv(sep: "\t", header: true)
      .map { row ->
        tuple(row.sample,
              file(row.hic_file, checkIfExists: true),
              row.resolution)
      }
      .set { hic_files }

    sample_sheet
      .splitCsv(sep: "\t", header: true)
      .map { row ->
        tuple(row.sample, make_optional_input(row.domains))
      }
      .set { domains }

    sample_sheet
      .splitCsv(sep: "\t", header: true)
      .map { row ->
        tuple(row.sample, make_optional_input(row.periphery_constraints))
      }
      .set { periphery_constraints }

    sig_interactions_cis.join(sig_interactions_trans)
      .map { tuple(it[0], tuple(it[1], it[2])) }
      .set { sig_interactions }

    MERGE(
      sig_interactions
    )

    DUMP_CHROM_SIZES(
      hic_files
    )

    DUMP_BINS(
      hic_files
    )

    domains.join(DUMP_BINS.out.bed)
      .map {
        def doms = it[1]
        def bins = it[2]

        if (doms.size() == 0) {
            doms = bins
        } else {
            doms = doms[0]
        }

        return tuple(it[0], file(doms))
      }
      .set { beads }


    DUMP_CHROM_SIZES.out.chrom_sizes
      .join(beads)
      .join(periphery_constraints)
      .join(MERGE.out.tsv)
      .set { make_bead_gtrack_tasks }

    MAKE_BEAD_GTRACK(
      make_bead_gtrack_tasks,
      masked_chromosomes
    )

    CHANGE_PLOIDY(
      MAKE_BEAD_GTRACK.out.gtrack,
      ploidy
    )

  emit:
    gtrack = CHANGE_PLOIDY.out.gtrack
    sig_interactions = MERGE.out.tsv
}

process MERGE {
  label 'duration_very_short'
  tag "${sample}"

  cpus 1

  input:
    tuple val(sample),
          path(interactions)

  output:
    tuple val(sample),
          path(outname),
    emit: tsv

  script:
    outname="${sample}.sig_interactions.tsv.gz"
    """
    #!/usr/bin/env python3

    import pandas as pd

    paths = "$interactions".split(" ")
    dfs = [pd.read_table(p) for p in paths]

    df = pd.concat(dfs).sort_values(["chrom1", "start1", "chrom2", "start2"])

    df.to_csv("$outname", sep="\\t", header=True, index=False, compression="gzip")
    """
}

process DUMP_CHROM_SIZES {
  label 'duration_very_short'
  tag "${sample}"

  cpus 1

  input:
    tuple val(sample),
          path(hic_file),
          val(resolution)

  output:
    tuple val(sample),
          path(outname),
    emit: chrom_sizes

  script:
    outname="${sample}.chrom.sizes"
    """
    hictk dump -t chroms '$hic_file' > '$outname'
    """
}

process DUMP_BINS {
  label 'duration_very_short'
  tag "${sample}"

  cpus 1

  input:
    tuple val(sample),
          path(hic_file),
          val(resolution)

  output:
    tuple val(sample),
          path(outname),
    emit: bed

  script:
    outname="${sample}.bins.bed"
    """
    hictk dump -t bins '$hic_file' --resolution '$resolution' > '$outname'
    """
}

process MAKE_BEAD_GTRACK {
  label 'duration_very_short'
  tag "${sample}"

  cpus 1

  input:
    tuple val(sample),
          path(chrom_sizes),
          path(beads),
          path(periphery_constraints),
          path(sig_interactions)

    val masked_chromosomes

  output:
    tuple val(sample),
          path("*.gtrack"),
    emit: gtrack

  script:
    outprefix="${sample}"
    args = []
    if (!periphery_constraints.toString().isEmpty()) {
        args.push("--periphery-constraints='${periphery_constraints}'")
    }
    if (!masked_chromosomes.isEmpty()) {
        args.push("--masked-chromosomes='${masked_chromosomes}'")
    }
    args=args.join(" ")
    """
    chrom3d_nf_make_bead_file.py \\
        '$sig_interactions' \\
        '$beads' \\
        '$chrom_sizes' \\
        $args \\
        > '$sample'.beads.gtrack
    """
}

process CHANGE_PLOIDY {
  label 'duration_very_short'
  tag "${sample}"

  cpus 1

  input:
    tuple val(sample),
          path(gtrack)

    val ploidy

  output:
    tuple val(sample),
          path("*.gtrack"),
    emit: gtrack

  script:
    outname="${gtrack.baseName}.${ploidy}.gtrack"
    """
    if [ '$ploidy' -eq 1 ]; then
        cp '$gtrack' '$outname'
    else
        chrom3d_nf_change_ploidy_gtrack.py '$gtrack' '$ploidy' > '$outname'
    fi
    """
}
