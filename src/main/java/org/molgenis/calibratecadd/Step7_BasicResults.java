package org.molgenis.calibratecadd;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.stat.descriptive.rank.Percentile.EstimationType;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;

/**
 * 
 * Read:
 * (e.g. clinvar.patho.fix.snpeff.exac.genes.tsv, output of step 4)
 * 
 * Gene	Category	Chr	Start	End	NrOfPopulationVariants	NrOfPathogenicVariants	NrOfOverlappingVariants	NrOfFilteredPopVariants	PathoMAFThreshold	PopImpactHighPerc	PopImpactModeratePerc	PopImpactLowPerc	PopImpactModifierPerc	PathoImpactHighPerc	PathoImpactModeratePerc	PathoImpactLowPerc	PathoImpactModifierPerc	PopImpactHighEq	PopImpactModerateEq	PopImpactLowEq	PopImpactModifierEq
 * SPINT2	Cx	19	38755433	38782579	269	6	2	6	2.326425E-4	1.92	22.36	16.29	59.42	50.00	33.33	16.67	0.00	50.00	33.33	16.67	0.00
 * NDST1	T1	5	149921113	149925129	100	4	1	0	7.000599999999998E-6	0.00	25.83	20.83	53.33	0.00	100.00	0.00	0.00
 * KIAA1109	N1	4	123128323	123128323	18	1	1	0	0	0.00	50.00	22.22	27.78	100.00	0.00	0.00	0.00
 * 
 * 
 * And:
 * (e.g. clinvar.patho.fix.snpeff.exac.withcadd.tsv, output of step 6)
 * 
 * gene	chr	pos	ref	alt	group	cadd
 * IFT172	2	27680627	A	T	PATHOGENIC	28.0
 * IFT172	2	27700177	A	T	PATHOGENIC	15.99
 * IFT172	2	27693963	C	T	PATHOGENIC	36
 * 
 * 
 * Write out:
 * (e.g. clinvar.patho.fix.snpeff.exac.genesumm.tsv)
 *
 * Gene	Category	Chr	Start	End	NrOfPopulationVariants	NrOfPathogenicVariants	NrOfOverlappingVariants	NrOfFilteredPopVariants	PathoMAFThreshold	PopImpactHighPerc	PopImpactModeratePerc	PopImpactLowPerc	PopImpactModifierPerc	PathoImpactHighPerc	PathoImpactModeratePerc	PathoImpactLowPerc	PathoImpactModifierPerc	PopImpactHighEq	PopImpactModerateEq	PopImpactLowEq	PopImpactModifierEq	NrOfCADDScoredPopulationVars	NrOfCADDScoredPathogenicVars	MeanPopulationCADDScore	MeanPathogenicCADDScore	MeanDifference	UTestPvalue	Sens95thPerCADDThreshold	Spec95thPerCADDThreshold
 * SPINT2	C4	19	38755433	38782579	269	6	2	6	2.326425E-4	1.92	22.36	16.29	59.42	50.00	33.33	16.67	0.00	50.00	33.33	16.67	0.00	6	6	22.05	24.85	2.80	0.33666836761003904	22.30	25.18
 * NDST1	T1	5	149921113	149925129	100	4	1	0	7.000599999999998E-6	0.00	25.83	20.83	53.33	0.00	100.00	0.00	0.00
 * KIAA1109	N1	4	123128323	123128323	18	1	1	0	0	0.00	50.00	22.22	27.78	100.00	0.00	0.00	0.00
 *
 *
 *
 * Example:
 * E:\Data\clinvarcadd\clinvar.patho.fix.snpeff.exac.genes.tsv
 * E:\Data\clinvarcadd\clinvar.patho.fix.snpeff.exac.withcadd.tsv
 * E:\Data\clinvarcadd\clinvar.patho.fix.snpeff.exac.genesumm.tsv
 *
 */
public class Step7_BasicResults
{
	HashMap<String, String> geneToInfo = new HashMap<String, String>();
	HashMap<String, ArrayList<String>> geneToVariantAndCADD = new HashMap<String, ArrayList<String>>();
	NumberFormat f = new DecimalFormat("#0.00");
	static String NA = "";


	public Step7_BasicResults(String geneInfoFile, String variantInfoFile, String outputFile) throws Exception
	{
		System.out.println("starting..");
		loadGeneInfo(geneInfoFile);
		loadVariantInfo(variantInfoFile);
		processAndWriteOutput(outputFile);
		System.out.println("..done");
	}

	public static void main(String[] args) throws Exception
	{
		new Step7_BasicResults(args[0], args[1], args[2]);
	}
	
	public void loadGeneInfo(String geneInfoFile) throws FileNotFoundException
	{
		/**
		 * read gene info and put in map
		 */
		Scanner geneInfoScanner = new Scanner(new File(geneInfoFile));
		geneInfoScanner.nextLine(); //skip header
		String line = null;
		while(geneInfoScanner.hasNextLine())
		{
			line = geneInfoScanner.nextLine();
			String[] split = line.split("\t", -1);
			line = line.substring(line.indexOf('\t')+1, line.length());
			geneToInfo.put(split[0], line);
		}
		geneInfoScanner.close();
	}
	
	public void loadVariantInfo(String variantInfoFile) throws FileNotFoundException
	{
		/**
		 * read variant + cadd data and put in map
		 */
		Scanner variantsWithCADDScanner = new Scanner(new File(variantInfoFile));
		variantsWithCADDScanner.nextLine();
		String line = null;
		while(variantsWithCADDScanner.hasNextLine())
		{
			line = variantsWithCADDScanner.nextLine();
			String gene = line.split("\t", -1)[0];
			if(geneToVariantAndCADD.containsKey(gene))
			{
				geneToVariantAndCADD.get(gene).add(line);
			}
			else
			{
				ArrayList<String> lines = new ArrayList<String>();
				lines.add(line);
				geneToVariantAndCADD.put(gene, lines);
			}
		}
		variantsWithCADDScanner.close();
	}
	
	public void processAndWriteOutput(String outputFile) throws Exception
	{
		/**
		 * process everything and write out
		 */
		PrintWriter pw = new PrintWriter(new File(outputFile));
		pw.println("Gene" + "\t" + "Category" + "\t" + "Chr" + "\t" + "Start" + "\t" + "End" + "\t" + "NrOfPopulationVariants" + "\t" + "NrOfPathogenicVariants" + "\t" + "NrOfOverlappingVariants" + "\t" + "NrOfFilteredPopVariants" + "\t" + "PathoMAFThreshold" + "\t" + "PopImpactHighPerc" + "\t" + "PopImpactModeratePerc" + "\t" + "PopImpactLowPerc" + "\t" + "PopImpactModifierPerc" + "\t" + "PathoImpactHighPerc" + "\t" + "PathoImpactModeratePerc" + "\t" + "PathoImpactLowPerc" + "\t" + "PathoImpactModifierPerc" + "\t" + "PopImpactHighEq" + "\t" + "PopImpactModerateEq" + "\t" + "PopImpactLowEq" + "\t" + "PopImpactModifierEq" + "\t" + "NrOfCADDScoredPopulationVars" + "\t" + "NrOfCADDScoredPathogenicVars" + "\t" + "MeanPopulationCADDScore" + "\t" + "MeanPathogenicCADDScore" + "\t" + "MeanDifference" + "\t" + "UTestPvalue" + "\t" + "Sens95thPerCADDThreshold" +"\t" + "Spec95thPerCADDThreshold");
		
		int nrOfGenesPathGtPopPval_5perc = 0;
		int nrOfGenesPathGtPopPval_1perc = 0;
		
		for(String gene : geneToInfo.keySet())
		{
			if(!geneToVariantAndCADD.containsKey(gene))
			{
				pw.println(gene + "\t" + geneToInfo.get(gene) + StringUtils.repeat("\t" + NA, 8));
			}
			else
			{
				ArrayList<Double> caddPatho = new ArrayList<Double>();
				ArrayList<Double> caddPopul = new ArrayList<Double>();
				for(String lineForGene : geneToVariantAndCADD.get(gene))
				{
					String[] split = lineForGene.split("\t", -1);
					String group = split[5];
					double cadd = Double.parseDouble(split[8]);
					if(group.equals("PATHOGENIC"))
					{
						caddPatho.add(cadd);
					}
					else if(group.equals("POPULATION"))
					{
						caddPopul.add(cadd);
					}
					else
					{
						pw.close();
						throw new Exception("unknown group " + group);
					}
				}
				
				//it can happen that variants for one group did not pass CADD webservice, e.g. for PRRT2 we have only 1 population variant and when that fails, we have cannot assess...
				//replace 'Cx' with 'N3'
				if(caddPatho.size() == 0 || caddPopul.size() == 0)
				{
					pw.println(gene + "\t" + "N3" + geneToInfo.get(gene).substring(2, geneToInfo.get(gene).length()) + "\t" + caddPopul.size() + "\t" + caddPatho.size() + StringUtils.repeat("\t" + NA, 6));
					continue;
				}
				
				double[] caddPathoPrim = new double[caddPatho.size()];
				for(int i = 0; i < caddPatho.size(); i++)
				{
					caddPathoPrim[i] = caddPatho.get(i);
				}
				
				double[] caddPopulPrim = new double[caddPopul.size()];
				for(int i = 0; i < caddPopul.size(); i++)
				{
					caddPopulPrim[i] = caddPopul.get(i);
				}
				
				Mean mean = new Mean();
				double pathoMean = mean.evaluate(caddPathoPrim);
				double populMean = mean.evaluate(caddPopulPrim);
				double meanDiff = pathoMean-populMean;
				
				MannWhitneyUTest utest = new MannWhitneyUTest();
				double pval = utest.mannWhitneyUTest(caddPathoPrim, caddPopulPrim);
				
				//get thresholds for 95% sensitivity and 95% specificity
				//sensitive: we catch 95% of the known pathogenic variants, no matter how many population variants we also include when using this threshold
				//specific: we only allow a 'top 5%' error (=finding population variants) in finding pathogenic variants, no matter how many pathogenic variants we have above this threshold right now
				Percentile perc = new Percentile().withEstimationType(EstimationType.R_7);
				double sensThres = perc.evaluate(caddPathoPrim, 5);
				double specThres = perc.evaluate(caddPopulPrim, 95);
				
				String cat = null;
				//to show some stats in the sysout
				if(pval <= 0.05 && pathoMean > populMean)
				{
					cat = "C2";
					nrOfGenesPathGtPopPval_5perc ++;
					if(pval <= 0.01)
					{
						cat = "C1";
						nrOfGenesPathGtPopPval_1perc++;
					}
				}
				
				if(pval > 0.05)
				{
					if(caddPathoPrim.length < 5 || caddPopulPrim.length < 5)
					{
						cat = "C3";
					}
					else
					{
						cat = "C4";
					}
				}
				
				if(cat == null)
				{
					cat = "C5";
				}
				
				//add info to gene
				geneToInfo.put(gene, geneToInfo.get(gene) + "\t" + caddPopulPrim.length + "\t" + caddPathoPrim.length + "\t" + f.format(populMean) + "\t" + f.format(pathoMean) + "\t" + f.format(meanDiff) + "\t" + pval + "\t" + f.format(sensThres) + "\t" + f.format(specThres));
				
				//write table
				pw.println(gene + "\t" + cat + geneToInfo.get(gene).substring(2, geneToInfo.get(gene).length()));
			}
		}
		
		System.out.println("total nr of genes: " + geneToVariantAndCADD.keySet().size());
		System.out.println("nr of genes where patho > pop, pval < 0.05: " + nrOfGenesPathGtPopPval_5perc);
		System.out.println("nr of genes where patho > pop, pval < 0.01: " + nrOfGenesPathGtPopPval_1perc);
		
		pw.flush();
		pw.close();
	}

}
