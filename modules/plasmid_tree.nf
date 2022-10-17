process concatenate_assemblies {

  executor 'local'

  input:
  path(assemblies)

  output:
  path("concatenated_assemblies.fa")

  script:
  """
  cat ${assemblies} > concatenated_assemblies.fa
  """
}

process mafft {

  tag { tree_id }

  publishDir params.versioned_outdir ? "${params.outdir}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}", mode: 'copy', pattern: "${tree_id}.aln.fa"

  input:
  tuple val(tree_id), path(assemblies)

  output:
  tuple val(tree_id), path("${tree_id}.aln.fa"), emit: aln
  tuple val(tree_id), path("${tree_id}_mafft_provenance.yml"), emit: provenance

  script:
  """
  printf -- "- process_name: mafft\\n"                                          >> ${tree_id}_mafft_provenance.yml
  printf -- "  tools:\\n"                                                       >> ${tree_id}_mafft_provenance.yml
  printf -- "    - tool_name: mafft\\n"                                         >> ${tree_id}_mafft_provenance.yml
  printf -- "      tool_version: \$(mafft --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${tree_id}_mafft_provenance.yml

  mafft \
    --thread ${task.cpus} \
    --auto \
    ${assemblies} \
    > ${tree_id}.aln.fa
  """
}

process iqtree {

  tag { tree_id }

  publishDir params.versioned_outdir ? "${params.outdir}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}", mode: 'copy', pattern: "${tree_id}.treefile"

  input:
  tuple val(tree_id), path(alignment)

  output:
  tuple val(tree_id), path("${tree_id}.treefile"), emit: tree
  tuple val(tree_id), path("${tree_id}_iqtree_provenance.yml"), emit: provenance

  script:
  """
  printf -- "- process_name: iqtree\\n"                                          >> ${tree_id}_iqtree_provenance.yml
  printf -- "  tools:\\n"                                                        >> ${tree_id}_iqtree_provenance.yml
  printf -- "    - tool_name: iqtree\\n"                                         >> ${tree_id}_iqtree_provenance.yml
  printf -- "      tool_version: \$(iqtree --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${tree_id}_iqtree_provenance.yml

  iqtree \
    -nt ${task.cpus} \
    -m 'GTR+G' \
    -s ${alignment} \
    --prefix ${tree_id}
  """
}

process snp_sites {

  publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", mode: 'copy', pattern: "mafft.aln.fa"

  input:
  tuple val(tree_id), path(alignment)

  output:
  tuple val(tree_id), path("mafft.aln.fa"), emit: aln
  tuple val(tree_id), path("snp_sites_provenance.yml"), emit: provenance

  script:
  """
  printf -- "- process_name: snp_sites\\n"                                          >> snp_sites_provenance.yml
  printf -- "  tools:\\n"                                                           >> snp_sites_provenance.yml
  printf -- "    - tool_name: snp-sites\\n"                                         >> snp_sites_provenance.yml
  printf -- "      tool_version: \$(snp-sites --version 2>&1 | cut -d ' ' -f 2)\\n" >> snp_sites_provenance.yml

  snp-sites \
    ${alignment} \
    > mafft.aln.fa
  """
}

