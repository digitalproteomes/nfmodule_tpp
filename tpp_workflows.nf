//////////////////////////
// Workflow definitions //
//////////////////////////


include {comet_search;
	 xinteract;
	 refactor_pepxml;
	 refactor_pepxml as refactor_pepxml_advanced;
	 refactor_protxml;
	 refactor_protxml as refactor_protxml_advanced;
	 st_peter;
	 st_peter2matrix;
	 mayu;
	 tpp_stat;
	 pepxml2tsv;
	 protxml2tsv;
	 filter_pep_tsv;
	 proteinprophet} from './tpp_processes.nf'


workflow tpp_main {
    // Comet search
    // xinteract either individually or pooled
    
    take:
    dda_folder
    comet_params
    protein_db
    comet_threads
    tpp
    decoy
    no_pool
    xinteract_threads

    main:
    dda_files = channel.fromPath("$dda_folder/*.{mzXML,mzML}")
    comet_search(dda_files,
		 file(comet_params),
		 file(protein_db),
		 comet_threads)
    if(!no_pool) {
	// Merge all search results into a single xinteract analysis
	pepxmls = comet_search.out.cometPepOutRaw.collect()
	mzxmls = comet_search.out.cometMzXMLOut.collect()
	run_xinteract(pepxmls,
		      protein_db,
		      mzxmls,
		      tpp,
		      decoy,
		      no_pool,
		      xinteract_threads)
    }
    else {
	run_xinteract(comet_search.out.cometPepOutRaw,
		      protein_db,
		      comet_search.out.cometMzXMLOut,
		      tpp,
		      decoy,
		      no_pool,
		      xinteract_threads)	
    }

    emit:
    pepxmls = run_xinteract.out.pepxmls_rf
    pepxmlmodels = run_xinteract.out.pepxmlmodels
    protxmls = run_xinteract.out.protxmls_rf
    protxmlmodels = run_xinteract.out.protxmlmodels
}


workflow run_xinteract {
    // Runs xinteract and refactors links in pepxml and protxml
    take:
    pepxmls
    protein_db
    mzXML
    tpp
    decoy
    no_pool
    xinteract_threads
    
    main:
    xinteract(pepxmls,
	      file(protein_db),
	      mzXML,
	      tpp,
	      decoy,
	      no_pool,
	      xinteract_threads)
    if(params.tpp.indexOf("-M") != -1){
	// Refactor everything but emit PTMProphet output only
	pepxmls_rf = refactor_pepxml_advanced(xinteract.out.tppPepPtmOutRaw)
	refactor_pepxml(xinteract.out.tppPepOutRaw)
	protxmls_rf = refactor_protxml_advanced(xinteract.out.tppProtPtmOutRaw)
	refactor_protxml(xinteract.out.tppProtOutRaw)
	pepxmlmodels = xinteract.out.tppPepModelPtmOut
	protxmlmodels = xinteract.out.tppProtModelPtmOut

    }
    else if(params.tpp.indexOf("-i") != -1) {
	// Refactor everything but emit InterProphet output only
	pepxmls_rf = refactor_pepxml_advanced(xinteract.out.tppPepIproOutRaw)
	refactor_pepxml(xinteract.out.tppPepOutRaw)
	protxmls_rf = refactor_protxml_advanced(xinteract.out.tppProtIproOutRaw)
	refactor_protxml(xinteract.out.tppProtOutRaw)
	pepxmlmodels = xinteract.out.tppPepModelIproOut
	protxmlmodels = xinteract.out.tppProtModelIproOut

    }
    else {
	pepxmls_rf = refactor_pepxml(xinteract.out.tppPepOutRaw)
	protxmls_rf = refactor_protxml(xinteract.out.tppProtOutRaw)
	pepxmlmodels = xinteract.out.tppPepModelOut
	protxmlmodels = xinteract.out.tppProtModelOut

    }
    
    emit:
    pepxmls_rf
    pepxmlmodels
    protxmls_rf
    protxmlmodels
}


workflow tpp_summaries {
    // Create tpp_stats and Mayu summaries of TPP analysis
    take:
    pepxml
    comet_params
    protein_db
    decoy

    main:
    mayu(pepxml,
	 file(comet_params),
	 file(protein_db),
	 decoy)

    tpp_stat(pepxml)
    
    emit:
    mayu.out
    tpp_stat.out
}


workflow tpp_exports {
    // Export pep and prot.xmls to tsv
    take:
    pepxml
    pepxmlmodels
    protxml
    protxmlmodels

    main:
    pepxml2tsv(pepxml,
	       pepxmlmodels,
	       file("$baseDir/Xslt/pepxml2tsv.xsl"))

    protxml2tsv(protxml,
		 protxmlmodels,
		 file("$baseDir/Xslt/protxml2tsv.xsl"))

    filter_pep_tsv(pepxml2tsv.out,
		   protxml2tsv.out)
}


workflow run_proteinprophet {
    take:
    pepxml

    main:
    proteinprophet(pepxml)
    protxmls_rf = refactor_protxml(proteinprophet.out.tppProtOutRaw)

    emit:
    protxmls_rf
}


workflow tpp_peter {
    take:
    pepxml
    protxml
    mzxml

    main:
    merged_protxml = run_proteinprophet(pepxml.collect())

    st_peter(pepxml,
	     protxml,
	     mzxml)

    st_peter2matrix(merged_protxml,
		    protxml.collect())

    emit:
    st_peter2matrix.out
}
