package org.molgenis.calibratecadd.support;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.*;

import org.molgenis.calibratecadd.Benchmark;
import org.molgenis.calibratecadd.support.BootStrappingVariant.OutCome;
import org.molgenis.calibratecadd.support.BootStrappingVariant.ExpClsf;
import org.molgenis.data.Entity;
import org.molgenis.data.annotation.entity.impl.gavin.GavinAlgorithm;
import org.molgenis.data.annotation.entity.impl.gavin.GavinEntry;
import org.molgenis.data.annotation.entity.impl.gavin.GavinEntry.Category;
import org.molgenis.data.annotation.entity.impl.snpEff.Impact;

import org.molgenis.data.annotation.entity.impl.gavin.Judgment;
import org.molgenis.data.annotation.entity.impl.gavin.Judgment.Classification;
import org.molgenis.data.annotation.entity.impl.gavin.Judgment.Method;
import org.molgenis.data.vcf.VcfRepository;


public class BootStrappingAnalysis
{

	private List<BootStrappingVariant> variantClsfResults = new ArrayList<BootStrappingVariant>();
	private String outputFile;
	private BootStrappingVariant.Label label;
	private Benchmark.ToolNames toolName;

	public static void main(String[] args) throws Exception
	{
		new BootStrappingAnalysis(args[0], args[1], args[3], Benchmark.ToolNames.valueOf(args[4]));
	}

	/**
	 * checked to give the exact same behaviour as Benchmark
	 * Except now we do random subsets of variants
	 * @param vcfFile
	 * @param gavinFile
	 * @throws Exception
	 */
	public BootStrappingAnalysis(String vcfFile, String gavinFile, String outputFile, Benchmark.ToolNames toolName) throws Exception
	{
		this.outputFile = outputFile;
		this.toolName = toolName;
		Files.write(Paths.get(outputFile), "Label\tCalib\tTool\tAcc\n".getBytes(), StandardOpenOption.APPEND);

		HashMap<String, GavinEntry> gavinData = Benchmark.loadCCGG(gavinFile).getGeneToEntry();
		GavinAlgorithm gavin = new GavinAlgorithm();

		// file with combined variants has 25,995 variants
		File variantList = new File(vcfFile);

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

			Double getMAF = GavinUtils.getInfoForAllele(record, "EXAC_AF", alt);
			double MAF = getMAF == null ? 0 : getMAF;
			Double CADDscore = GavinUtils.getInfoForAllele(record, "CADD_SCALED", alt);
			String ann = record.getString("ANN");
			Set<String> genes = GavinUtils.getGenesFromAnn(ann);
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

			for(String gene : genes) {
				if (!hasGeneId || gene.equals(geneFromId)) {
					Impact impact = GavinUtils.getImpact(ann, gene, alt);
					geneToIdMatchFound = true;
					Judgment judgment;
					if (toolName.equals(Benchmark.ToolNames.GAVIN))
					{
						judgment = gavin.classifyVariant(impact, CADDscore, MAF, gene, null, gavinData);

					}
					else if (toolName.equals(Benchmark.ToolNames.GAVINnocal))
					{
						judgment = gavin.genomewideClassifyVariant(impact, CADDscore, MAF, gene);
					}
					else
					{
						throw new Exception("unknown tool " + toolName);
					}
					multipleJudgments.add(judgment);
				}
			}
			if(hasGeneId && !geneToIdMatchFound)
			{
				//System.out.println("WARNING: bad data for variant " + chr + ":" + pos + " " + ref + "/" + alt + ", no match from ID field gene to snpeff annotations!");
				multipleJudgments.add(new Judgment(Classification.VOUS, Method.calibrated, geneFromId, "Bad data!"));
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
			String gene = null;

			for(Judgment judgment : multipleJudgments)
			{
				if(judgment.getClassification().equals(Classification.Benign))
				{
					nrOfBenignClsf++;
				}
				if(judgment.getClassification().equals(Classification.Pathogenic))
				{
					nrOfPathognClsf++;
				}

				if(judgment.getConfidence().equals(Method.calibrated))
				{
					//if judgment was calibrated, add this gene
					hasCalibratedJudgment = true;
					gene = judgment.getGene();
				}
				else
				{
					//if judgment was not calibrated, add the first gene
					//this should be the gene with the most severe effect and therefore the most relevant
					gene = judgment.getGene();
				}
			}

			/**
			 * Now we can assign the final verdict for this variant
			 */

			//check if we have any conflicts
			//TODO could be improved by prioritizing calibrated over genomewide results for our method
			if(nrOfBenignClsf > 0 && nrOfPathognClsf > 0)
			{
				addToFullSetClsfOutcomes(Classification.VOUS, mvlClassfc, gavinData, gene);
				//	System.out.println("WARNING: conflicting classification! adding no judgment for this variant: " + chr + ":" + pos + " " + ref + "/" + alt + ", judged: " + multipleJudgments.toString() );
			}
			else
			{
				for(Judgment judgment : multipleJudgments)
				{
					//if we know we have calibrated results, wait for it, then add it, and then break
					if(hasCalibratedJudgment && judgment.getConfidence().equals(Method.calibrated))
					{
						addToFullSetClsfOutcomes(judgment.getClassification(), mvlClassfc, gavinData, gene);
						break;
					}
					//if not, might as well add this one and be done
					//TODO: this means there may be multiple verdicts, e.g. 2x BENIGN for context in two genes, but we only add 1 of them, to keep things a bit more simple
					else if(!hasCalibratedJudgment)
					{
						addToFullSetClsfOutcomes(judgment.getClassification(), mvlClassfc, gavinData, gene);
						break;
					}
				}
			}

		}

		//count proportion in full set
		int nrExpPatho = 0;
		int nrExpBenign = 0;
		for(BootStrappingVariant bv : variantClsfResults)
		{
			switch (bv.getExpClsf()) {
				case B:  nrExpBenign++; break;
				case P:  nrExpPatho++; break;
				default : throw new Exception("no expclsf");
			}
		}

		System.out.println("nrExpBenign = " + nrExpBenign);
		System.out.println("nrExpPatho = " + nrExpPatho);

	}

	private void addToFullSetClsfOutcomes (Classification observed, String expected, HashMap<String, GavinEntry> gavinData, String gene)
	{

		BootStrappingVariant.Label label = null;

		if(gavinData.get(gene) != null)
		{
			if(gavinData.get(gene).category.equals(GavinEntry.Category.C1) || gavinData.get(gene).category.equals(Category.C2))
			{
				label = BootStrappingVariant.Label.C1_C2;
			}
			if(gavinData.get(gene).category.equals(Category.C3))
			{
				label = BootStrappingVariant.Label.C3;
			}
			if(gavinData.get(gene).category.equals(Category.C4))
			{
				label = BootStrappingVariant.Label.C4;
			}
		}

		if(observed.equals(Classification.Benign) && (expected.equals("B") ||  expected.equals("LB")))
		{
			variantClsfResults.add(new BootStrappingVariant(OutCome.TN, label, ExpClsf.B));
		}

		if(observed.equals(Classification.Benign) && (expected.equals("P") ||  expected.equals("LP")))
		{
			variantClsfResults.add(new BootStrappingVariant(OutCome.FN, label, ExpClsf.P));
		}

		if(observed.equals(Classification.Pathogenic) && (expected.equals("B") ||  expected.equals("LB")))
		{
			variantClsfResults.add(new BootStrappingVariant(OutCome.FP, label, ExpClsf.B));
		}

		if(observed.equals(Classification.Pathogenic) && (expected.equals("P") ||  expected.equals("LP")))
		{
			variantClsfResults.add(new BootStrappingVariant(OutCome.TP, label, ExpClsf.P));
		}

		if(observed.equals(Classification.VOUS) && (expected.equals("P") ||  expected.equals("LP")))
		{
			variantClsfResults.add(new BootStrappingVariant(OutCome.VOUS, label, ExpClsf.P));
		}

		if(observed.equals(Classification.VOUS) && (expected.equals("B") ||  expected.equals("LB")))
		{
			variantClsfResults.add(new BootStrappingVariant(OutCome.VOUS, label, ExpClsf.B));
		}

	}

	public void getStatsOnSet(List<BootStrappingVariant> set) throws Exception
	{
		if(set == null)
		{
			return;
		}

		int TP = 0;
		int TN = 0;
		int FP = 0;
		int FN = 0;
		int VOUS = 0;
		int nrExpPatho = 0;
		int nrExpBenign = 0;

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
			switch (bv.getExpClsf()) {
				case B:  nrExpBenign++; break;
				case P:  nrExpPatho++; break;
				default : throw new Exception("no expclsf");
			}
		}

		System.out.println("TN = " + TN);
		System.out.println("TP = " + TP);
		System.out.println("FP = " + FP);
		System.out.println("FN = " + FN);
		System.out.println("VOUS = " + VOUS);
		double MCC = ProcessJudgedVariantMVLResults.getMCC(TP, TN, FP, FN);
		double acc = (double)(TN+TP)/(double)(TN+TP+FP+FN+VOUS);
		System.out.println("MCC = " + MCC);
		double cov = (TP+TN+FP+FN)/(double)(TP+TN+FP+FN+VOUS);
		System.out.println("acc = " + (acc));
		System.out.println("% true-pathogenic: " + (nrExpPatho/(double)(nrExpPatho+nrExpBenign)));

		String toR = label+"_"+toolName + "\t" + label + "\t" + toolName + "\t" + (acc) + "\n";
		Files.write(Paths.get(outputFile), toR.getBytes(), StandardOpenOption.APPEND);

	}

	public List<BootStrappingVariant> randomSubset(int benignSamples, int pathoSamples, BootStrappingVariant.Label label)
	{
		this.label = label;
		ArrayList<BootStrappingVariant> result = new ArrayList<BootStrappingVariant>();
		Collections.shuffle(variantClsfResults);
		for (int i = 0; i < variantClsfResults.size(); i++)
		{
			BootStrappingVariant bv = variantClsfResults.get(i);
			if(label.equals(bv.getLabel()) && bv.getExpClsf().equals(ExpClsf.B) && benignSamples > 0)
			{
				result.add(bv);
				benignSamples--;
			}
			else if(label.equals(bv.getLabel()) && bv.getExpClsf().equals(ExpClsf.P) && pathoSamples > 0)
			{
				result.add(bv);
				pathoSamples--;
			}
			if(benignSamples == 0 && pathoSamples == 0)
			{
				break;
			}
		}
		return result;
	}
}
