package org.molgenis.calibratecadd.support;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.molgenis.calibratecadd.Step9_Validation;
import org.molgenis.calibratecadd.support.BootStrappingVariant.OutCome;
import org.molgenis.data.Entity;
import org.molgenis.data.annotation.entity.impl.SnpEffAnnotator.Impact;
import org.molgenis.data.annotation.joeri282exomes.CCGGUtils;
import org.molgenis.data.annotation.joeri282exomes.Judgment;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Classification;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Method;
import org.molgenis.data.vcf.VcfRepository;

public class BootStrappingAnalysis
{

	public static void main(String[] args) throws Exception
	{

		new BootStrappingAnalysis(args[0], args[1], args[3]);

	}
	
	private List<BootStrappingVariant> variantClsfResults = new ArrayList<BootStrappingVariant>();
	String outputFile;
	
	public void getStatsOnFullSet() throws Exception
	{
		getStatsOnSet(variantClsfResults);
	}
			
	
	public void getStatsOnSet(List<BootStrappingVariant> set) throws Exception
	{
		int TP = 0;
		int TN = 0;
		int FP = 0;
		int FN = 0;
		int VOUS = 0;
		int judgmentsInCalibratedGenes = 0;
		
		for(BootStrappingVariant bv : set)
		{
			switch (bv.getOutcome()) {
				case TP:  TP++; break;
				case TN:  TN++; break;
				case FP:  FP++; break;
				case FN:  FN++; break;
				case VOUS:  VOUS++; break;
				default : throw new Exception("no outcome");
			}
			if(bv.inCalibGene)
			{
				judgmentsInCalibratedGenes++;
			}
		}
		
		System.out.println("TN = " + TN);
		System.out.println("TP = " + TP);
		System.out.println("FP = " + FP);
		System.out.println("FN = " + FN);
		System.out.println("VOUS = " + VOUS);
		double MCC = ProcessJudgedVariantMVLResults.getMCC(TP, TN, FP, FN);
		System.out.println("MCC = " + MCC);
		double cov = (TP+TN+FP+FN)/(double)(TP+TN+FP+FN+VOUS);
		System.out.println("MCCcovadj = " + (cov*MCC));
		double percCalib = judgmentsInCalibratedGenes/(double)set.size();
		System.out.println("% calibrated: " + judgmentsInCalibratedGenes + "/" + set.size() + " = " + percCalib);
		
		String toR = "row <- data.frame(MCCcovadj = "+(cov*MCC)+", percCalib = "+percCalib+"); df <- rbind(df, row)\n";
		Files.write(Paths.get(outputFile), toR.getBytes(), StandardOpenOption.APPEND);

	}
	
	public List<BootStrappingVariant> randomSubset(int sampleSize)
	{
		ArrayList<BootStrappingVariant> result = new ArrayList<BootStrappingVariant>();
		Collections.shuffle(variantClsfResults);
		for (int i = 0; i < sampleSize; i++)
		{
			result.add(variantClsfResults.get(i));
		}
		return result;
	}
	
	public  List<BootStrappingVariant> downSampleToUniformCalibPercDistr(List<BootStrappingVariant> set) throws Exception
	{
		ArrayList<BootStrappingVariant> result = new ArrayList<BootStrappingVariant>();
		int judgmentsInCalibratedGenes = 0;
		for(BootStrappingVariant bv : set)
		{
			if(bv.inCalibGene)
			{
				judgmentsInCalibratedGenes++;
			}
		}
		
		//uniform distribution of the percentage of variants in calibrated genes
		double wantedPercCalib = new Random().nextDouble();

		System.out.println("set.size() = " + set.size());
		System.out.println("judgmentsInCalibratedGenes = " + judgmentsInCalibratedGenes);
		int judgmentsInUncalibratedGenes = set.size()-judgmentsInCalibratedGenes;
		System.out.println("judgmentsInUncalibratedGenes = " + judgmentsInUncalibratedGenes);
		double currentPercCalib = judgmentsInCalibratedGenes/(double)set.size();
		System.out.println("current percCalib = " + currentPercCalib);
		System.out.println("wanted percCalib = " + wantedPercCalib);

		
		int nrOfUncalibVariantsToDelete = 0;
		int nrOfCalibVariantsToDelete = 0;
		if(wantedPercCalib > currentPercCalib)
		{
			int targetAmountOfUncalib = (int) Math.round((judgmentsInCalibratedGenes/wantedPercCalib)-judgmentsInCalibratedGenes);
			System.out.println("target uncalib amount = " + targetAmountOfUncalib);
			nrOfUncalibVariantsToDelete = judgmentsInUncalibratedGenes-targetAmountOfUncalib;
			System.out.println("nrOfUncalibVariantsToDelete = " + nrOfUncalibVariantsToDelete);
		}
		else
		{
			int targetAmountOfCalib = (int) Math.round((judgmentsInUncalibratedGenes/(1-wantedPercCalib))-judgmentsInUncalibratedGenes);
			System.out.println("target calib amount = " + targetAmountOfCalib);
			nrOfCalibVariantsToDelete = judgmentsInCalibratedGenes-targetAmountOfCalib;
			System.out.println("nrOfCalibVariantsToDelete = " + nrOfCalibVariantsToDelete);
		}
		
		for(BootStrappingVariant bv : set)
		{
			if(bv.inCalibGene && nrOfCalibVariantsToDelete > 0)
			{
				//dont add to result but subtract one from the amount to leave out
				nrOfCalibVariantsToDelete--;
				continue;
			}
			
			if(!bv.inCalibGene && nrOfUncalibVariantsToDelete > 0)
			{
				//dont add to result but subtract one from the amount to leave out
				nrOfUncalibVariantsToDelete--;
				continue;
			}

			result.add(bv);
		}
		
			
		return result;
	}
	

	/**
	 * checked to give the exact same behaviour as Step9_Validation
	 * Except now we do random subsets of variants
	 * @param vcfFile
	 * @param gavinFile
	 * @throws Exception
	 */
	public BootStrappingAnalysis(String vcfFile, String gavinFile, String outputFile) throws Exception
	{
		this.outputFile = outputFile;
		File yourFile = new File(outputFile);
		if(!yourFile.exists()) {
		    yourFile.createNewFile();
			Files.write(Paths.get(outputFile), "df <- data.frame()\n".getBytes(), StandardOpenOption.APPEND);
		} 
		
		CCGGUtils gavin = Step9_Validation.loadCCGG(gavinFile);
		
		// file with combined variants has 25,995 variants
		File variantList = new File(vcfFile);

		
		
	//	System.out.println("selected line nrs: " + randomX.toString());

		VcfRepository vcfRepo = new VcfRepository(variantList, "vcf");

		java.util.Iterator<Entity> vcfRepoIter = vcfRepo.iterator();
		
		int lineNr = 0;
		while (vcfRepoIter.hasNext())
		{
			lineNr++;
			
			Entity record = vcfRepoIter.next();
			
			String chr = record.getString("#CHROM");
			String pos = record.getString("POS");
			String ref = record.getString("REF");
			String alt = record.getString("ALT");
			if(alt.contains(","))
			{
				throw new Exception("Did not expect multiple alt alleles! " + record.toString());
			}
			
			Double getMAF = CCGGUtils.getInfoForAllele(record, "EXAC_AF", alt);
			double MAF = getMAF == null ? 0 : getMAF;
			Double CADDscore = CCGGUtils.getInfoForAllele(record, "CADD_SCALED", alt);
			String ann = record.getString("ANN");
			Set<String> genes = CCGGUtils.getGenesFromAnn(ann);
			String id = record.getString("ID");
			
			//for some variants, we have GENE:MUTATION in the ID
			// for example, in the MVL data, we see "GUSB:c.1222C>T", and in the VariBench data: "PAH:c.617A>G"
			//if this is present, we use this to narrow the scope by matching the annotated genes to the gene symbol here
			String[] idSplit = id.split(":", -1);
			boolean hasGeneId = idSplit.length > 1 ? true : false;
			String geneFromId = idSplit[0];
			
			String mvlClassfc = record.getString("CLSF");
			String mvlName = record.getString("MVL");
			
			ArrayList<Judgment> multipleJudgments = new ArrayList<Judgment>();
			
			boolean geneToIdMatchFound = false;

			for(String gene : genes)
			{
				if(!hasGeneId || gene.equals(geneFromId))
				{
					geneToIdMatchFound = true;
					Impact impact = CCGGUtils.getImpact(ann, gene, alt);
					Judgment judgment = gavin.classifyVariant(gene, MAF, impact, CADDscore);
					multipleJudgments.add(judgment);
				}
			}
			if(hasGeneId && !geneToIdMatchFound)
			{
	//			System.out.println("WARNING: bad data for variant " + chr + ":" + pos + " " + ref + "/" + alt + ", no match from ID field gene to snpeff annotations!");
				multipleJudgments.add(new Judgment(Classification.VOUS, Method.calibrated, "Bad data!"));
			}
			
			//if no judgment, add null for this variant
			if(multipleJudgments.size() == 0)
			{
				throw new Exception("No judgments! should not occur.");
			}
			
			//go through the possible classifications and check if any of them are conflicting
			//also, if we have a calibrated judgment, 
			int nrOfBenignClsf = 0;
			int nrOfPathognClsf = 0;
			boolean hasCalibratedJudgment = false;
			for(Judgment judgment : multipleJudgments)
			{
				if(judgment.getClassification().equals(Classification.Benign))
				{
					nrOfBenignClsf++;
				}
				if(judgment.getClassification().equals(Classification.Pathogn))
				{
					nrOfPathognClsf++;
				}
				if(judgment.getConfidence().equals(Method.calibrated))
				{
					hasCalibratedJudgment = true;
				}
			}
			
			/**
			 * Now we can assign the final verdict for this variant
			 */
			
			//check if we have any conflicts
			//TODO could be improved by prioritizing calibrated over genomewide results for our method
			if(nrOfBenignClsf > 0 && nrOfPathognClsf > 0)
			{
				addToFullSetClsfOutcomes(Classification.VOUS, mvlClassfc, hasCalibratedJudgment);
			//	System.out.println("WARNING: conflicting classification! adding no judgment for this variant: " + chr + ":" + pos + " " + ref + "/" + alt + ", judged: " + multipleJudgments.toString() );
			}
			else
			{
				for(Judgment judgment : multipleJudgments)
				{
					//if we know we have calibrated results, wait for it, then add it, and then break
					if(hasCalibratedJudgment && judgment.getConfidence().equals(Method.calibrated))
					{
					//	addToMVLResults(judgment, mvlClassfc, mvlName, record);
					//	judgmentsInCalibratedGenes++;
						addToFullSetClsfOutcomes(judgment.getClassification(), mvlClassfc, true);
						break;
					}
					//if not, might as well add this one and be done
					//TODO: this means there may be multiple verdicts, e.g. 2x BENIGN for context in two genes, but we only add 1 of them, to keep things a bit more simple
					else if(!hasCalibratedJudgment)
					{
					//	addToMVLResults(judgment, mvlClassfc, mvlName, record);
						
						addToFullSetClsfOutcomes(judgment.getClassification(), mvlClassfc, false);
						break;
					}
				}
			}

		}
	
	}
	
	private void addToFullSetClsfOutcomes (Classification observed, String expected, boolean hasCalibratedJudgment)
	{
		if(observed.equals(Classification.Benign) && (expected.equals("B") ||  expected.equals("LB")))
		{
			variantClsfResults.add(new BootStrappingVariant(OutCome.TN, hasCalibratedJudgment));
		}
		
		if(observed.equals(Classification.Benign) && (expected.equals("P") ||  expected.equals("LP")))
		{
			variantClsfResults.add(new BootStrappingVariant(OutCome.FN, hasCalibratedJudgment));
		}
		
		if(observed.equals(Classification.Pathogn) && (expected.equals("B") ||  expected.equals("LB")))
		{
			variantClsfResults.add(new BootStrappingVariant(OutCome.FP, hasCalibratedJudgment));
		}
		
		if(observed.equals(Classification.Pathogn) && (expected.equals("P") ||  expected.equals("LP")))
		{
			variantClsfResults.add(new BootStrappingVariant(OutCome.TP, hasCalibratedJudgment));
		}
		
		if(observed.equals(Classification.VOUS) && (expected.equals("P") ||  expected.equals("LP") || expected.equals("B") ||  expected.equals("LB")))
		{
			variantClsfResults.add(new BootStrappingVariant(OutCome.VOUS, hasCalibratedJudgment));
		}
				
		
	}

}