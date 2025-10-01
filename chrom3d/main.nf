#!/usr/bin/env nextflow

// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

nextflow.enable.dsl=2
nextflow.preview.output=1


include { SAMPLESHEET } from './subworkflows/local/samplesheet'
include { NCHG as NCHG_CIS } from './subworkflows/nchg'
include { NCHG as NCHG_TRANS } from './subworkflows/nchg'
include { PREPROCESSING } from './subworkflows/local/preprocessing'
include { CHROM3D } from './subworkflows/local/chrom3d'


workflow {
  main:
    log.info("-- PARAMETERS")
    log.info("")
    if (params.sample_sheet) {
      log.info("-- sample_sheet: ${params.sample_sheet}")
    } else {
      log.info("-- sample: ${params.sample}")
      log.info("-- hic_file: ${params.hic_file}")
      log.info("-- resolution: ${params.resolution}")
      log.info("-- beads: ${params.beads}")
      log.info("-- periphery_constraints: ${params.periphery_constraints}")
      log.info("-- mask: ${params.mask}")
    }
    log.info("-- publish_dir_mode: ${params.publish_dir_mode}")

    log.info("-- cytoband: ${params.cytoband}")
    log.info("-- assembly_gaps: ${params.assembly_gaps}")
    log.info("-- masked_chromosomes: ${params.masked_chromosomes}")

    log.info("-- chrom3d_args: ${params.chrom3d_args}")
    log.info("-- ploidy: ${params.ploidy}")
    log.info("-- number_of_models: ${params.number_of_models}")
    log.info("-- archive_models: ${params.archive_models}")

    log.info("-- nchg_mad_max: ${params.nchg_mad_max}")
    log.info("-- nchg_bad_bin_fraction: ${params.nchg_bad_bin_fraction}")
    log.info("-- nchg_fdr_cis: ${params.nchg_fdr_cis}")
    log.info("-- nchg_log_ratio_cis: ${params.nchg_log_ratio_cis}")
    log.info("-- nchg_fdr_trans: ${params.nchg_fdr_trans}")
    log.info("-- nchg_log_ratio_trans: ${params.nchg_log_ratio_trans}")

    log.info("-- plot_format: ${params.plot_format}")
    log.info("-- nchg_skip_plots: ${params.nchg_skip_plots}")
    log.info("-- nchg_hic_tgt_resolution_plots: ${params.nchg_hic_tgt_resolution_plots}")
    log.info("-- nchg_plot_sig_interactions_cmap_lb: ${params.nchg_plot_sig_interactions_cmap_lb}")
    log.info("-- nchg_plot_sig_interactions_cmap_ub: ${params.nchg_plot_sig_interactions_cmap_ub}")

    log.info("")

    SAMPLESHEET(
      params.sample_sheet,
      params.sample,
      params.hic_file,
      params.resolution,
      params.beads,
      params.periphery_constraints,
      params.mask
    )

    SAMPLESHEET.out.tsv.set { sample_sheet }
    SAMPLESHEET.out.nchg_cis.set { nchg_cis_sample_sheet }
    SAMPLESHEET.out.nchg_trans.set { nchg_trans_sample_sheet }

    NCHG_CIS(
      nchg_cis_sample_sheet,
      params.nchg_mad_max,
      params.nchg_bad_bin_fraction,
      params.cytoband,
      params.assembly_gaps,
      true,  // use_cis_interactions
      params.nchg_fdr_cis,
      params.nchg_log_ratio_cis,
      false, // use_trans_interactions
      params.nchg_fdr_trans,
      params.nchg_log_ratio_trans,
      params.nchg_skip_plots,  // skip_expected_plots
      params.nchg_skip_plots,  // skip_sign_interaction_plots
      params.nchg_hic_tgt_resolution_plots,
      params.nchg_plot_sig_interactions_cmap_lb,
      params.nchg_plot_sig_interactions_cmap_ub,
      params.plot_format,
      params.zstd_compression_lvl
    )

    NCHG_TRANS(
      nchg_trans_sample_sheet,
      params.nchg_mad_max,
      params.nchg_bad_bin_fraction,
      params.cytoband,
      params.assembly_gaps,
      false, // use_cis_interactions
      params.nchg_fdr_cis,
      params.nchg_log_ratio_cis,
      true,  // use_trans_interactions
      params.nchg_fdr_trans,
      params.nchg_log_ratio_trans,
      params.nchg_skip_plots,  // skip_expected_plots
      params.nchg_skip_plots,  // skip_sign_interaction_plots
      params.nchg_hic_tgt_resolution_plots,
      params.nchg_plot_sig_interactions_cmap_lb,
      params.nchg_plot_sig_interactions_cmap_ub,
      params.plot_format,
      params.zstd_compression_lvl
    )

    NCHG_CIS.out.tsv
      .map {
        def sample = it[0] - ~/_cis$/
        tuple(sample, it[1])
      }
      .set { nchg_sig_interactions_cis }

    NCHG_TRANS.out.tsv
      .map {
        def sample = it[0] - ~/_trans$/
        tuple(sample, it[1])
      }
      .set { nchg_sig_interactions_trans }

    PREPROCESSING(
      sample_sheet,
      nchg_sig_interactions_cis,
      nchg_sig_interactions_trans,
      params.ploidy,
      params.masked_chromosomes
    )

    CHROM3D(
      PREPROCESSING.out.gtrack,
      params.chrom3d_args,
      params.number_of_models,
      params.archive_models
    )

  publish:
    sample_sheet = SAMPLESHEET.out.tsv
    nchg_cis_sample_sheet = SAMPLESHEET.out.nchg_cis
    nchg_trans_sample_sheet = SAMPLESHEET.out.nchg_trans
    nchg_cis_parquet = NCHG_CIS.out.parquets
    nchg_cis_expected = NCHG_CIS.out.expected
    nchg_cis_plots = NCHG_CIS.out.plots
    nchg_trans_parquet = NCHG_TRANS.out.parquets
    nchg_trans_expected = NCHG_TRANS.out.expected
    nchg_trans_plots = NCHG_TRANS.out.plots
    sig_interactions = PREPROCESSING.out.sig_interactions
    gtrack = PREPROCESSING.out.gtrack
    chrom3d_models = CHROM3D.out.models
    chrom3d_tar = CHROM3D.out.tar
}

output {
  // TODO use named outputs
  sample_sheet {
    path '.'
    mode params.publish_dir_mode
  }
  nchg_cis_sample_sheet {
    path '.'
    mode params.publish_dir_mode
  }
  nchg_trans_sample_sheet {
    path '.'
    mode params.publish_dir_mode
  }
  nchg_cis_parquet {
    path {
      def sample = it[0] - ~/_cis$/
      "$sample/nchg/"
    }
    mode params.publish_dir_mode
  }
  nchg_cis_expected {
    path {
      def sample = it[0] - ~/_cis$/
      "$sample/nchg/"
    }
    mode params.publish_dir_mode
  }
  nchg_cis_plots {
    path {
      def sample = it[0] - ~/_cis$/
      "$sample/nchg/plots/"
    }
    mode params.publish_dir_mode
  }
  nchg_trans_parquet {
    path {
      def sample = it[0] - ~/_trans$/
      "$sample/nchg/"
    }
    mode params.publish_dir_mode
  }
  nchg_trans_expected {
    path {
      def sample = it[0] - ~/_trans$/
      "$sample/nchg/"
    }
    mode params.publish_dir_mode
  }
  nchg_trans_plots {
    path {
      def sample = it[0] - ~/_trans$/
      "$sample/nchg/plots/"
    }
    mode params.publish_dir_mode
  }
  sig_interactions {
    path { "${it[0]}/" }
    mode params.publish_dir_mode
  }
  gtrack {
    path { "${it[0]}/" }
    mode params.publish_dir_mode
  }
  chrom3d_models {
    path { "${it[0]}/models/" }
    mode params.publish_dir_mode
  }
  chrom3d_tar {
    path { "${it[0]}/" }
    mode params.publish_dir_mode
  }
}
