package org.molgenis.calibratecadd;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.List;


public class Step8_FullBenchmark
{
	public static void main(String[] args) throws Exception
	{
		//paper submitted: r0.1 (using clinvar: 5 nov 2015). after revision: r0.2, using the same clinvar data but snpeff-annotated ExAC
		//r0.3 (using clinvar: 1 sep 2016)
		//latest: r0.4 (gnomad 2.0.1, clinvar may 11th 2015)
		String version = "r0.4";
		String outFile = "/Users/joeri/github/gavin/data/other/step9_panels_out_"+version+".R";
		String path = "/Users/joeri/github/gavin/data/predictions";
		List<String> datasets = Arrays.asList(new String[]{"Renal", "Pulmonary", "Ophthalmologic", "Oncologic", "Obstetric", "NotInCGD", "Neurologic", "Musculoskeletal", "Hematologic", "Genitourinary", "Gastrointestinal", "Endocrine", "Dermatologic", "Dental", "Craniofacial", "Cardiovascular", "Biochemical", "Audiologic_Otolaryngologic", "Allergy_Immunology_Infectious"});
	//	List<String> tools = Arrays.asList(new String[]{"GAVIN2"});
		List<String> tools = Arrays.asList(new String[]{"GAVIN", "GAVIN2", "CADD_Thr15", "CADD_Thr20", "CADD_Thr25", "MSC_ClinVar95CI", "MSC_HGMD99CI", "PROVEAN", "SIFT", "PolyPhen2", "Condel", "PONP2", "PredictSNP2", "FATHMM", "GWAVA", "FunSeq", "DANN"});

		if(new File(outFile).exists())
		{
			throw new Exception("output file already exists: " + outFile);
		}
		new File(outFile).createNewFile();
		Files.write(Paths.get(outFile), "df <- data.frame()\n".getBytes(), StandardOpenOption.APPEND);

		for(String dataset: datasets)
		{
			for(String tool: tools)
			{
				new Benchmark(path, "/Users/joeri/github/gavin/data/goldstandards/cgdpanels/" + dataset, tool, outFile, version, null);
			}
		}

	}
	
}
