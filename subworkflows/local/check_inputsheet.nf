workflow CHECK_INPUTSHEET {
    take:
    ch_input_path // file: /path/to/input_sheet.csv

    main:
    /*
     * Check input_sheet is valid
     */

    ch_input_path.splitCsv(header: true, sep: "\t", strip: true)
              .filter { it ->
                def sample_id = it.Sample_ID
                def path = it.Path
                def genome = it.Genome
                def speciesValue = it.Species
                if (sample_id == null || sample_id == "" ) {
                    exit 1, "Specify \"Sample_ID\" in the input_sheet.tsv."
                }
                if ( ! file(path).exists() ) {
                    exit 1, "The file in \"Path\" does not exist."
                }
                if ( ! file(genome).exists() ) {
                    exit 1, "The file in \"Genome\" does not exist."
                }
                if (speciesValue == null || speciesValue == "" ) {
                    exit 1, "Specify \"Species\" in the input_sheet.tsv."
                }
                return true
              }
              .map { it ->
                [it.Group_ID, it.Sample_ID, it.Analyte_Type, it.Pore_Type, file(it.Path), file(it.Genome), it.Species, it.Sample_Rate ]
              }
              .set{ ch_input }

    emit:
    ch_input
}