package org.molgenis.data.annotation.joeri282exomes;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;

import org.molgenis.calibratecadd.support.VariantClassificationException;
import org.molgenis.calibratecadd.support.JudgedVariant.ExpertClassification;
import org.molgenis.data.Entity;
import org.molgenis.data.annotation.entity.impl.SnpEffAnnotator.Impact;
import org.molgenis.data.annotation.joeri282exomes.CCGGEntry.Category;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Classification;
import org.molgenis.data.annotation.joeri282exomes.Judgment.Method;

public class CCGGUtils
{
	HashMap<String, CCGGEntry> geneToEntry = new HashMap<String, CCGGEntry>();
	
	public CCGGUtils(File ccgg) throws Exception
	{	
		Scanner s = new Scanner(ccgg);
		
		//skip header
		s.nextLine();
		
		String line;
		while(s.hasNextLine())
		{
			line = s.nextLine();
			
			CCGGEntry e = new CCGGEntry(line);
			geneToEntry.put(e.gene, e);
		}
		
	}
	
	public Category getCategory(String gene)
	{
		return geneToEntry.get(gene).category;
	}
	
	public boolean contains(String gene)
	{
		return geneToEntry.containsKey(gene) ? true : false;
	}
	
	public Judgment classifyVariant(String gene, Double MAF, Impact impact, Double CADDscore, ExpertClassification mvlClassfc) throws Exception
	{
		//if we have no data for this gene, immediately fall back to the naive method
		if(!geneToEntry.containsKey(gene))
		{
			return naiveClassifyVariant(gene, MAF, impact, CADDscore, mvlClassfc);
		}

		CCGGEntry entry = geneToEntry.get(gene);
		CCGGEntry.Category category = entry.category;
		
		// MAF based classification, calibrated
		if(MAF > entry.PathoMAFThreshold)
		{
			if(mvlClassfc.equals(ExpertClassification.B) || mvlClassfc.equals(ExpertClassification.LB)){ System.out.println("true_benign_pathomaf"); }
			else if(mvlClassfc.equals(ExpertClassification.P) || mvlClassfc.equals(ExpertClassification.LP)){ System.out.println("false_benign_pathomaf"); }
			else { System.out.println("VOUS_benign_pathomaf"); }
				
			return new Judgment(Classification.Benign, Method.calibrated, "Variant MAF of " + MAF + " is greater than the pathogenic 95th percentile MAF of "+ entry.PathoMAFThreshold + ".");
		}
		
		String mafReason = "the variant MAF of " + MAF + " is lesser than the pathogenic 95th percentile MAF of "+ entry.PathoMAFThreshold + ".";
		
		// Impact based classification, calibrated
		if(impact != null)
		{
			if(category.equals(Category.I1) && impact.equals(Impact.HIGH))
			{
				if(mvlClassfc.equals(ExpertClassification.B) || mvlClassfc.equals(ExpertClassification.LB)){ System.out.println("false_patho_i1"); }
				else if(mvlClassfc.equals(ExpertClassification.P) || mvlClassfc.equals(ExpertClassification.LP)){ System.out.println("true_patho_i1"); }
				else { System.out.println("VOUS_patho_i1"); }
			
				return new Judgment(Judgment.Classification.Pathogn,  Method.calibrated, "Variant is of high impact, while there are no known high impact variants in the population. Also, " + mafReason);
			}
			else if(category.equals(Category.I2) && (impact.equals(Impact.MODERATE) || impact.equals(Impact.HIGH)))
			{
				if(mvlClassfc.equals(ExpertClassification.B) || mvlClassfc.equals(ExpertClassification.LB)){ System.out.println("false_patho_i2"); }
				else if(mvlClassfc.equals(ExpertClassification.P) || mvlClassfc.equals(ExpertClassification.LP)){ System.out.println("true_patho_i2"); }
				else { System.out.println("VOUS_patho_i2"); }
	
				return new Judgment(Judgment.Classification.Pathogn,  Method.calibrated, "Variant is of high/moderate impact, while there are no known high/moderate impact variants in the population. Also, " + mafReason);
			}
			else if(category.equals(Category.I3) && (impact.equals(Impact.LOW) || impact.equals(Impact.MODERATE) || impact.equals(Impact.HIGH)))
			{
				if(mvlClassfc.equals(ExpertClassification.B) || mvlClassfc.equals(ExpertClassification.LB)){ System.out.println("false_patho_i3"); }
				else if(mvlClassfc.equals(ExpertClassification.P) || mvlClassfc.equals(ExpertClassification.LP)){ System.out.println("true_patho_i3"); }
				else { System.out.println("VOUS_patho_i3"); }
	
				return new Judgment(Judgment.Classification.Pathogn,  Method.calibrated, "Variant is of high/moderate/low impact, while there are no known high/moderate/low impact variants in the population. Also, " + mafReason);
			}
			else if(impact.equals(Impact.MODIFIER))
			{
				if(mvlClassfc.equals(ExpertClassification.B) || mvlClassfc.equals(ExpertClassification.LB)){ System.out.println("true_benign_modf"); }
				else if(mvlClassfc.equals(ExpertClassification.P) || mvlClassfc.equals(ExpertClassification.LP)){ System.out.println("false_benign_modf"); }
				else { System.out.println("VOUS_benign_modf"); }
	
				return new Judgment(Judgment.Classification.Benign,  Method.calibrated, "Variant is of 'modifier' impact, and therefore unlikely to be pathogenic. However, " + mafReason);
			}
		}

		// CADD score based classification, calibrated
		if(CADDscore != null)
		{
			if((category.equals(Category.C1) || category.equals(Category.C2)))
			{
				if(CADDscore > entry.MeanPathogenicCADDScore)
				{
					if(mvlClassfc.equals(ExpertClassification.B) || mvlClassfc.equals(ExpertClassification.LB)){ System.out.println("false_patho_c12"); }
					else if(mvlClassfc.equals(ExpertClassification.P) || mvlClassfc.equals(ExpertClassification.LP)){ System.out.println("true_patho_c12"); }
					else { System.out.println("VOUS_patho_c12"); }
		
					return new Judgment(Judgment.Classification.Pathogn,  Method.calibrated, "Variant CADD score of " + CADDscore + " is greater than the mean pathogenic score of " + entry.MeanPathogenicCADDScore + " in a gene for which CADD scores are informative. Also, " + mafReason);
				}
				else if(CADDscore < entry.MeanPopulationCADDScore)
				{
					if(mvlClassfc.equals(ExpertClassification.B) || mvlClassfc.equals(ExpertClassification.LB)){ System.out.println("true_benign_c12"); }
					else if(mvlClassfc.equals(ExpertClassification.P) || mvlClassfc.equals(ExpertClassification.LP)){ System.out.println("false_benign_c12"); }
					else { System.out.println("VOUS_benign_c12"); }
		
					return new Judgment(Judgment.Classification.Benign,  Method.calibrated, "Variant CADD score of " + CADDscore + " is lesser than the mean population score of " + entry.MeanPathogenicCADDScore + " in a gene for which CADD scores are informative, although " + mafReason);
				}
			}
			else if((category.equals(Category.C3) || category.equals(Category.C4) || category.equals(Category.C5)))
			{
				if(CADDscore > entry.Spec95thPerCADDThreshold)
				{
					if(mvlClassfc.equals(ExpertClassification.B) || mvlClassfc.equals(ExpertClassification.LB)){ System.out.println("false_patho_c345"); }
					else if(mvlClassfc.equals(ExpertClassification.P) || mvlClassfc.equals(ExpertClassification.LP)){ System.out.println("true_patho_c345"); }
					else { System.out.println("VOUS_patho_c345"); }
		
					return new Judgment(Judgment.Classification.Pathogn,  Method.calibrated, "Variant CADD score of " + CADDscore + " is greater than the 95% specificity threhold of " + entry.Spec95thPerCADDThreshold + " for this gene. Also, " + mafReason);
				}
				else if(CADDscore < entry.Sens95thPerCADDThreshold)
				{
					if(mvlClassfc.equals(ExpertClassification.B) || mvlClassfc.equals(ExpertClassification.LB)){ System.out.println("true_benign_c345"); }
					else if(mvlClassfc.equals(ExpertClassification.P) || mvlClassfc.equals(ExpertClassification.LP)){ System.out.println("false_benign_c345"); }
					else { System.out.println("VOUS_benign_c345"); }
		
					return new Judgment(Judgment.Classification.Benign,  Method.calibrated, "Variant CADD score of " + CADDscore + " is lesser than the 95% sensitivity threhold of " + entry.MeanPathogenicCADDScore + " for this gene, although " + mafReason);
				}
			}
		}
		
		//if everything so far has failed, we can still fall back to the naive method
		return naiveClassifyVariant(gene, MAF, impact, CADDscore, mvlClassfc);
	}
	
	
	public Judgment naiveClassifyVariant(String gene, Double MAF, Impact impact, Double CADDscore, ExpertClassification mvlClassfc) throws Exception
	{
		if(MAF > 0.00474)
		{
			if(mvlClassfc.equals(ExpertClassification.B) || mvlClassfc.equals(ExpertClassification.LB)){ System.out.println("true_benign_naive_maf"); }
			else if(mvlClassfc.equals(ExpertClassification.P) || mvlClassfc.equals(ExpertClassification.LP)){ System.out.println("false_benign_naive_maf"); }
			else { System.out.println("VOUS_benign_naive_maf"); }

			return new Judgment(Judgment.Classification.Benign, Method.naive, "MAF > 0.00474");
		}
		if(impact.equals(Impact.MODIFIER))
		{
			if(mvlClassfc.equals(ExpertClassification.B) || mvlClassfc.equals(ExpertClassification.LB)){ System.out.println("true_benign_naive_modifier"); }
			else if(mvlClassfc.equals(ExpertClassification.P) || mvlClassfc.equals(ExpertClassification.LP)){ System.out.println("false_benign_naive_modifier"); }
			else { System.out.println("VOUS_benign_naive_modifier"); }
			return new Judgment(Judgment.Classification.Benign, Method.naive, "Impact.MODIFIER");
		}
		else
		{
			if(CADDscore != null && CADDscore > 34)
			{
				if(mvlClassfc.equals(ExpertClassification.B) || mvlClassfc.equals(ExpertClassification.LB)){ System.out.println("false_patho_naive_cadd"); }
				else if(mvlClassfc.equals(ExpertClassification.P) || mvlClassfc.equals(ExpertClassification.LP)){ System.out.println("true_patho_naive_cadd"); }
				else { System.out.println("VOUS_patho_naive_cadd"); }

				return new Judgment(Judgment.Classification.Pathogn, Method.naive, "CADDscore > 34");
			}
			else if(CADDscore != null && CADDscore < 2)
			{
				if(mvlClassfc.equals(ExpertClassification.B) || mvlClassfc.equals(ExpertClassification.LB)){ System.out.println("true_benign_naive_cadd"); }
				else if(mvlClassfc.equals(ExpertClassification.P) || mvlClassfc.equals(ExpertClassification.LP)){ System.out.println("false_benign_naive_cadd"); }
				else { System.out.println("VOUS_benign_naive_cadd"); }
				return new Judgment(Judgment.Classification.Benign, Method.naive, "CADDscore < 2");
			}
			else
			{
				throw new VariantClassificationException("Unable to classify " + gene + "! ");
			}
		}
	}
	
	public static Set<String> getGenesFromAnn(String ann) throws Exception
	{
		Set<String> genes = new HashSet<String>();
		String[] annSplit = ann.split(",", -1);
		for(String oneAnn : annSplit)
		{
			String[] fields = oneAnn.split("\\|", -1);
			String gene = fields[3];
			genes.add(gene);
		}
		if(genes.size() == 0)
		{
			throw new Exception("No genes for " + ann);
		}
		return genes;
	}
	
	public static Double getInfoForAllele(Entity record, String infoField, String altAllele) throws Exception
	{
		String info_STR = record.get(infoField) == null ? null : record.get(infoField).toString();
		if(info_STR == null)
		{
			return null;
		}
		String[] alts = record.getString("ALT").split(",", -1);
		String[] info_split = info_STR.split(",", -1);
	
		if(alts.length != info_split.length)
		{
			throw new Exception("length of alts not equal to length of info field for " + record);
		}
		
		for (int i = 0; i < alts.length; i++)
		{
			if(alts[i].equals(altAllele))
			{
				return  (info_split[i] != null && !info_split[i].equals(".")) ? Double.parseDouble(info_split[i]) : null;
			}
		}
		return null;
	}
	
	public static Impact getImpact(String ann, String gene, String allele) throws Exception
	{
		String findAnn = getAnn(ann, gene, allele);
		if(findAnn == null)
		{
			System.out.println("FAILED TO GET IMPACT FOR " + ann);
			return null;
		}
		else
		{
			String[] fields = findAnn.split("\\|", -1);
			String impact = fields[2];
			return Impact.valueOf(impact);
		}
	}
	
	public static String getAnn(String ann, String gene, String allele) throws Exception
	{
		String[] annSplit = ann.split(",", -1);
		for(String oneAnn : annSplit)
		{
			String[] fields = oneAnn.split("\\|", -1);
			String geneFromAnn = fields[3];
			if(!gene.equals(geneFromAnn))
			{
				continue;
			}
			String alleleFromAnn = fields[0];
			if(!allele.equals(alleleFromAnn))
			{
				continue;
			}
			return oneAnn;
		}
		System.out.println("warning: annotation could not be found for " + gene + ", allele=" + allele + ", ann=" + ann);
		return null;
	}

	public static void main(String[] args) throws Exception
	{
		File ccgg = new File(args[0]);
		new CCGGUtils(ccgg);

	}

}
