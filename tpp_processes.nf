/////////////////////////
// Process definitions //
/////////////////////////


process comet_search {
    // Search an mzXML file with Comet
    tag "$mzXML"
    cpus "$comet_threads"
    // Estimated based on Human, 3 variable mods, semi, 2 missed cleavages
    // and some margin for safety
    memory 30.GB
    
    publishDir 'Results/Comet', mode: 'link'

    input:
    file mzXML
    file comet_params
    file protein_db
    val comet_threads

    output:
    path '*.pep.xml', emit: cometPepOutRaw
    path mzXML, emit: cometMzXMLOut


    """
    # Set proteins DB
    sed -i s,db_path,$protein_db, $comet_params
    sed -i 's,num_threads = 0,num_threads = ${params.comet_threads},' $comet_params

    comet -P$comet_params $mzXML
    sed -ri 's|/tmp/nxf.{11}|${workflow.launchDir}/Results/Comet/|'g ${mzXML.simpleName}.pep.xml
    sed -ri 's|<search_database local_path="|<search_database local_path="${workflow.launchDir}/Results/Comet/|'g ${mzXML.simpleName}.pep.xml
    """
}


process xinteract {
    tag "$pepxml"
    cpus "$xinteract_threads"
    publishDir 'Results/Comet', mode: 'link'
    
    input:
    file pepxml		// If doing a polled search collect() files before calling process
    file protein_db
    file mzXML
    val tpp
    val decoy
    val no_pool
    val xinteract_threads
    
    output:
    path '*.pep.xml', emit: tppPepOutRaw
    path '*.ipro.pep.xml', emit: tppPepIproOutRaw, optional: true
    path '*.ptm.ipro.pep.xml', emit: tppPepPtmOutRaw, optional: true
    
    path '*.prot.xml', emit: tppProtOutRaw, optional: true
    path '*.ipro.prot.xml', emit: tppProtIproOutRaw, optional: true
    path '*.ptm.ipro.prot.xml', emit: tppProtPtmOutRaw, optional: true
    
    path '*.pep-MODELS.html', emit: tppPepModelOut
    path '*.ipro.pep-MODELS.html', emit: tppPepModelIproOut, optional: true
    path '*.ptm.ipro.pep-MODELS.html', emit: tppPepModelPtmOut, optional: true
    
    path '*.pep.xml.index'	// Catchall
    path '*.pep.xml.pIstats'

    path '*.prot-MODELS.html', emit: tppProtModelOut, optional: true
    path '*.ipro.prot-MODELS.html', emit: tppProtModelIproOut, optional: true
    path '*.ptm.ipro.prot-MODELS.html', emit: tppProtModelPtmOut, optional: true

    file(protein_db) // Required for ProteinProphet visualization

    script:
    if(!no_pool) {
	// xinteract all files together
	"""
	xinteract $tpp -THREADS=$xinteract_threads -d$decoy -Ncomet_merged.pep.xml $pepxml
	"""
    }
    else {
	// xinteract each file individually
	"""
	xinteract $tpp -THREADS=$xinteract_threads -d$decoy -N${pepxml}_sep.pep.xml $pepxml
	"""
    }
}


process proteinprophet {
    tag "$pepxml"
    publishDir 'Results/Comet', mode: 'link'

    input:
    file pepxml

    output:
    path 'comet_merged.prot.xml', emit: tppProtOutRaw

    script:
    """
    ProteinProphet $pepxml comet_merged.prot.xml
    """
}


process refactor_pepxml {
    // pepXml uses absolute links. Fix them after moving them to the
    // publishDir
    tag "$pepxml"
    publishDir 'Results/Comet', mode: 'link'

    input:
    file pepxml

    output:
    file pepxml

    """
    sed -ri 's|/tmp/nxf.{11}|${workflow.launchDir}/Results/Comet/|'g $pepxml
    """
}


process refactor_protxml {
    // protXml uses absolute links. Fix them after moving them to the
    // publishDir
    tag "$protxml"
    publishDir 'Results/Comet', mode: 'link'

    input:
    file protxml

    output:
    file protxml

    """
    sed -ri 's|/work/.{2}/.{30}|/Results/Comet|'g $protxml 
    sed -ri 's|/tmp/nxf.{11}|${workflow.launchDir}/Results/Comet/|'g $protxml
    """
}


process st_peter {
    tag "$protxml"
    
    input:
    file pepxml
    file protxml
    file mzXML

    output:
    file protxml
    
    script:
    """
    StPeter $protxml
    """
}


process st_peter2matrix {
    tag "$protxmls"
    publishDir 'Results/TppExport', mode: 'link'
    
    input:
    file merged_protxml
    file protxmls
    
    output:
    file 'merged_stpeter.tsv'

    script:
    """
    StPeter2Matrix -g $merged_protxml $protxmls merged_stpeter.tsv
    """
}


process mayu {
    tag "$pepxml"
    publishDir 'Results/Comet', mode: 'link'

    errorStrategy 'ignore'
    
    input:
    file pepxml
    file comet_params
    file protein_db
    val decoy
    
    output:
    file("mayu_*")
    
    """
    Mayu.pl -A $pepxml \
    -C $protein_db \
    -E $decoy \
    -M mayu_$pepxml \
    -P pep FDR=0.01:1
    """
}


process tpp_stat {
    // For each TPP analysis run calctppstat
    tag "$pepxml"
    publishDir 'Results/Comet', mode: 'link'
    
    input:
    file pepxml

    output:
    file '*.summary.txt'
    
    """
    /usr/local/tpp/cgi-bin/calctppstat.pl -i $pepxml -d $params.decoy --full > ${pepxml}.summary.txt
    """
}


process pepxml2tsv {
    tag "$pepXml"
    
    publishDir 'Results/TppExport', mode: 'link'

    input:
    file pepXml
    file pepXmlModels
    file pepxsl
    
    output:
    file '*.tsv'

    """
    PROB=\$(get_prophet_prob.py -i $pepXmlModels)
    xsltproc --param p_threshold \$PROB $pepxsl $pepXml > ${pepXml}.tsv
    """
}


process protxml2tsv {
    tag "$protXml"
    
    publishDir 'Results/TppExport', mode: 'link'

    input:
    file protXml
    file protXmlModels
    file protxsl
    
    output:
    file '*.tsv'

    """
    PROB=\$(get_prophet_prob.py -i $protXmlModels)
    xsltproc --param p_threshold \$PROB $protxsl $protXml > ${protXml}.tsv
    """
}


// Filters the peptide list created by pepxml2tsv by removing peptides
// not assigned to a protein included in the protxml2tsv generated
// protein list.
process filter_pep_tsv {
    tag "$pepTsv - $protTsv"
    
    publishDir 'Results/TppExport', mode: 'link'

    input:
    file pepTsv
    file protTsv

    output:
    file '*.tsv'
    
    script:
    """
    filter_pep_pro.py -p $pepTsv -P $protTsv -o ${pepTsv.baseName}_filtered.tsv
    """
}
