package org.molgenis.calibratecadd;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;
import java.util.zip.GZIPInputStream;

import org.molgenis.calibratecadd.support.ClinVarVariant;

public class Step1_GetClinVarPathogenic
{

	/**
	 * Program arguments, example:
	 * E:\Data\clinvarcadd\variant_summary.txt
	 * E:\Data\clinvarcadd\clinvar.patho.vcf
	 */
	public static String CLINVAR_INFO = "CLINVAR";
	public static String VKGL_INFO = "VKGL";


	// download @ ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
	public static void main(String[] args) throws Exception
	{
		HashMap<String, List<ClinVarVariant>> cvv = Step1_GetClinVarPathogenic.getAsMap(new File(args[0]));
		PrintWriter pw = new PrintWriter(new File(args[1]));
		
		pw.println("##fileformat=VCFv4.1");
		pw.println("##reference=hg19.fasta");
		for(String chr : new String[]{"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"})
		{
			pw.println("##contig=<ID="+chr+">");
		}
		pw.println("##INFO=<ID="+CLINVAR_INFO+",Number=1,Type=String,Description=\"Variant classification according to ClinVar database\">");
		pw.println("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO");
		
		for(String gene : cvv.keySet())
		{
			for(ClinVarVariant v : cvv.get(gene))
			{
				pw.println(v.chrom + "\t" + v.pos + "\t" + v.id + "\t" + v.ref + "\t" + v.alt + "\t" + "." + "\t" + "." + "\t" + CLINVAR_INFO + "=" + v.clinvarInfoToString());
			}
		}
		pw.flush();
		pw.close();
		
	}
	
	
	public static HashMap<String, List<ClinVarVariant>> getAsMap(File clinvarFile) throws Exception
	{
		System.out.println("loading clinvar..");
		HashMap<String, List<ClinVarVariant>> res = new HashMap<String, List<ClinVarVariant>>();

	//	Scanner s = new Scanner(clinvarFile);

		GZIPInputStream gzip = new GZIPInputStream(new FileInputStream(clinvarFile));
		BufferedReader s = new BufferedReader(new InputStreamReader(gzip));

		String line = null;

		int lost = 0;
		int totalvariants = 0;

		//skip header
		line = s.readLine();

		while (true)
		{
			line = s.readLine();
			if(line == null ) { break; }

			String[] lineSplit = line.split("\t", -1);

			// needs to be GRCh37, also sanity check other options
			String genomeBuild = lineSplit[16];
			if (!genomeBuild.equals("GRCh37") && !genomeBuild.equals("GRCh38") && !genomeBuild.equals("NCBI36") && !genomeBuild.equals(""))
			{
				throw new Exception("bad genome build: " + genomeBuild);
			}
			if (!genomeBuild.equals("GRCh37"))
			{
				continue;
			}

			// needs to contain 'pathogenic'
			String clinsig = lineSplit[6];
		//	System.out.println(clinsig);
			if (!clinsig.toLowerCase().contains("pathogenic"))
			{
				continue;
			}

			//but not 'conflicting' or 'GRCh38/hg38'
			if (clinsig.toLowerCase().contains("conflicting") || clinsig.contains("GRCh38/hg38"))
			{
				continue;
			}

			// eg. NM_005343.2(HRAS):c.37G>C (p.Gly13Arg)
			String name = lineSplit[2];
			// eg. HRAS
			String gene = lineSplit[4];

			// System.out.println(line);

			String geneFromName = null;
			if (name.contains("("))
			{
				geneFromName = name.substring(name.indexOf('(') + 1, name.indexOf(')'));
			}

			if (gene.equals("-") || gene.equals(""))
			{
				if(gene.equals("")) { System.out.println("empty gene, override: '"+gene+"' to '"+geneFromName+"'"); }
				gene = geneFromName;
			}

			if (gene == null)
			{
				lost++;
				System.out.println("lost due to no gene symbol: " + line);
				continue;
			}

			String chrom = lineSplit[18];
			String pos = lineSplit[19];
			String rsid = lineSplit[9];
			String ref = lineSplit[21];
			String alt = lineSplit[22];

			// no RS id and no ref/alt notation means we have no way to fix this variant
			if(rsid.equals("-1") && ref.equals("na") && alt.equals("na"))
			{
				lost++;
			//	System.out.println("lost due to -1 na na: " + line);
				continue;
			}

			//replace 'na' by '-' to make it easier later
			if(ref.equals("na"))
			{
				ref = "-";
			}
			if(alt.equals("na"))
			{
				alt = "-";
			}
			ClinVarVariant cvv = new ClinVarVariant(chrom, pos, rsid, ref, alt, name, gene, clinsig);
			
			if(res.containsKey(gene))
			{
				res.get(gene).add(cvv);
			}
			else
			{
				List<ClinVarVariant> cvvList = new ArrayList<ClinVarVariant>();
				cvvList.add(cvv);
				res.put(gene, cvvList);
			}
			
			totalvariants++;

		}
		
		s.close();

		System.out.println("..done, put " + totalvariants + " 'pathogenic' variants in " + res.size() + " genes, lost due to no gene symbol/bad data:" + lost);
		
		return res;

	}

}
