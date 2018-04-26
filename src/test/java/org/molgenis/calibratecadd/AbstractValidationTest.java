package org.molgenis.calibratecadd;

import org.apache.commons.io.FileUtils;

import org.springframework.util.FileCopyUtils;

import java.io.*;
import java.util.Arrays;
import java.util.List;

import static org.testng.Assert.assertEquals;

public class AbstractValidationTest
{
	private String tool;
	public AbstractValidationTest(String tool) throws Exception {
		this.tool = tool;
		beforeClass();
	}

	protected File expected;

	public void beforeClass() throws IOException
	{
		InputStream expectedOutputInpStr = AbstractValidationTest.class.getResourceAsStream("/"+tool+"-test-expected.R");
		expected = new File(FileUtils.getTempDirectory(), tool + "-test-expected.R");
		FileCopyUtils.copy(expectedOutputInpStr, new FileOutputStream(expected));
		System.out.println("Created expected results file at: " + expected.getAbsolutePath());
	}

	public void testPredictionTool() throws Exception
	{
		File observed = new File(FileUtils.getTempDirectory(), tool + "-test-observed.R");
		if(!observed.createNewFile()){
			FileUtils.write(observed, "");
		}
		System.out.println("Created observed results file at: " + observed.getAbsolutePath());

		File predictionsDir = new File("data/predictions");
		System.out.println("Getting predictions from: " + predictionsDir.getAbsolutePath());

		File panelsDir = new File("data/goldstandards/cgdpanels");
		System.out.println("Getting panels from: " + panelsDir.getAbsolutePath());

		List<String> datasets = Arrays.asList(new String[]{"Renal", "Pulmonary", "Ophthalmologic", "Oncologic", "Obstetric", "NotInCGD", "Neurologic", "Musculoskeletal", "Hematologic", "Genitourinary", "Gastrointestinal", "Endocrine", "Dermatologic", "Dental", "Craniofacial", "Cardiovascular", "Biochemical", "Audiologic_Otolaryngologic", "Allergy_Immunology_Infectious"});
		List<String> tools = Arrays.asList(new String[]{"GAVIN", "CADD", "MSC_ClinVar95CI", "MSC_HGMD99CI", "PROVEAN", "SIFT", "PolyPhen2", "Condel", "PONP2", "PredictSNP2", "FATHMM", "GWAVA", "FunSeq", "DANN"});

		for(String dataset: datasets)
		{
			new Benchmark(predictionsDir.getAbsolutePath(), new File(panelsDir, dataset).getAbsolutePath(), tool, observed.getAbsolutePath(), "r0.1", null);
		}

		System.out.println("Going to compare files:\n\n" + expected.getAbsolutePath() + "\nvs.\n" + observed.getAbsolutePath());

		assertEquals(FileUtils.readLines(expected), FileUtils.readLines(observed));

		System.out.println("\n--> they are equal.");
	}

}
