#!/usr/bin/env nextflow

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { hash_files }             from './modules/hash_files.nf'
include { concatenate_assemblies } from './modules/plasmid_tree.nf'
include { mafft }                  from './modules/plasmid_tree.nf'
include { iqtree }                 from './modules/plasmid_tree.nf'
include { snp_sites }              from './modules/plasmid_tree.nf'
include { pipeline_provenance }    from './modules/provenance.nf'
include { collect_provenance }     from './modules/provenance.nf'


workflow {

  ch_start_time = Channel.of(workflow.start)
  ch_pipeline_name = Channel.of(workflow.manifest.name)
  ch_pipeline_version = Channel.of(workflow.manifest.version)

  ch_pipeline_provenance = pipeline_provenance(ch_pipeline_name.combine(ch_pipeline_version).combine(ch_start_time))

  ch_tree_id = Channel.of(params.tree_id)

  if (params.samplesheet_input != 'NO_FILE') {
    ch_assemblies = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['ASSEMBLY']] }
  } else {
    ch_assemblies = Channel.fromPath(params.assembly_search_path).map{ it -> [it.getName().split('_')[0], it] }.unique{ it -> it[0] }
  }

  main:
    hash_files(ch_assemblies.combine(Channel.of("assembly-input")))

    concatenate_assemblies(ch_assemblies.collect{ it -> it[1] })

    mafft(ch_tree_id.combine(concatenate_assemblies.out))

    iqtree(mafft.out.aln)

    ch_provenance = mafft.out.provenance
    ch_provenance = ch_provenance.join(iqtree.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    
    ch_provenance = ch_provenance.join(hash_files.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(ch_tree_id.combine(ch_pipeline_provenance)).map{ it -> [it[0], it[1] << it[2]] }

    collect_provenance(ch_provenance)
}
