package org.molgenis.calibratecadd;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Scanner;
import java.util.concurrent.ExecutionException;

import org.molgenis.calibratecadd.support.ChrPosRefAltUniqueVariants;
import org.molgenis.data.Entity;
import org.molgenis.data.annotator.tabix.TabixVcfRepository;
import org.molgenis.gavin2.util.FixNotation;
import org.molgenis.gavin2.util.FixVcfNotation;

public class Step2_FixIndelNotation
{




	/**
	 * TODO: after all curation, there are about 156 duplicate lines - must keep track and don't write them out
	 *
	 * Uses:
	 * [0] file produced in step 1
	 * [1] ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
	 * [2] output file
	 *
	 * Example:
	 * E:\Data\clinvarcadd\clinvar.patho.vcf
	 * E:\Data\clinvarcadd\clinvar.vcf.gz
	 * E:\Data\clinvarcadd\clinvar.patho.fix.vcf
	 *
	 * try to fix 'na' and '-'
	 *
	 * note: http://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/
	 * "At present, ClinVar's VCF file is limited to records that have been assigned rs#."
	 *
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception
	{
		Step2_FixIndelNotation step2 = new Step2_FixIndelNotation(args[0], args[1], args[2]);
	}

	public Step2_FixIndelNotation(String args0, String args1, String args2) throws Exception
	{
		Scanner s = new Scanner(new File(args0));
		PrintWriter pw = new PrintWriter(new File(args2));
		PrintWriter pwlost = new PrintWriter(new File(args2+"failed"));


		TabixVcfRepository clinvarVcf = new TabixVcfRepository(new File(args1), "clinvar");

		int fixes = 0;
		int lost = 0;
		int passed = 0;

		//failed fixed due to missing the target: 627 @ 10, 605 @ 20, 604 @ 25+
		int fixWindowSearchSize = 150;

		String line;
		while(s.hasNextLine())
		{
			line = s.nextLine();
			if(line.startsWith("#"))
			{
				//header, just print and continue
				pw.println(line);
				continue;
			}

			if(passed % 100 == 0)
			{
				System.out.println("passed: " + passed + ", fixes: " + fixes + ", lost: " + lost);
			}

			String[] lineSplit = line.split("\t", -1);

			String chrom = lineSplit[0];
			int pos = Integer.parseInt(lineSplit[1]);
			String rsid = lineSplit[2];
			String ref = lineSplit[3];
			String alt = lineSplit[4];
			String info = lineSplit[7].replace(";", ","); //prevent INFO field termination

			//if both alleles are present (not '-'), only check if trimming is needed
			if( !(ref.equals("-") || alt.equals("-")))
			{
				String trimmed = ref + "\t" + alt;
				if(ref.length()>1 && alt.length()>1)
				{
					trimmed = FixNotation.backTrimRefAlt(ref, alt, "\t");
					fixes++;
					//System.out.println("needs trimming: " + ref + " " + alt + " to " + trimmed);

				}
			//	System.out.println("alleles OK: " + trimmed);
				pw.println(chrom + "\t" + pos + "\t" + rsid + "\t" + trimmed + "\t" + "." + "\t" + "." + "\t" + info);
				passed++;
				continue;
			}

			// these are all cases that need fixing
			if(ref.equals("-") || alt.equals("-"))
			{

				//there is an RS id, we have a quick chance to fix it
				if (!rsid.equals("-1"))
				{
					String fix = fixByClinVarVCFRSid(chrom, pos, rsid, clinvarVcf, fixWindowSearchSize);
					if(fix != null)
					{
						pw.println(fix + "\t" + "." + "\t" + "." + "\t" + info);
						fixes++;
						passed++;
						continue;
					}
				}

				//if there is no RS id at all, or the lookup in ClinVar VCF failed, we can try to fix notation via UCSC
				String fix2 = FixVcfNotation.process(chrom, pos, ref, alt, 100, "-");

				if(fix2 != null)
				{
					pw.println(fix2 + "\t" + "." + "\t" + "." + "\t" + info);
					fixes++;
					passed++;
					continue;
				}
			}

			System.out.println("unable to fix variant: " + line);
			pwlost.println(line);
			lost++;
		}

		pw.flush();
		pw.close();
		pwlost.flush();
		pwlost.close();
		s.close();
		clinvarVcf.close();

		System.out.println("passed: " + passed);
		System.out.println("fixes needed: " + fixes);
		System.out.println("lost: " + lost);
	}


	public String fixByClinVarVCFRSid(String chrom, int pos, String rsid, TabixVcfRepository clinvarVcf, int fixWindowSearchSize)
	{
		//NOTE: deletions will be found some position before!
		//e.g. "47705505 T -" will be in the VCF as "47705504 AT A" !
		// "47635557 AG -" will be at "47635554	rs63749848 CAG C"
		//tricky? yes... but we can grab a window and match by RS id
		List<Entity> records = clinvarVcf.query(chrom, pos-fixWindowSearchSize, pos+fixWindowSearchSize);
		for(Entity e : records)
		{
			String idFromClinVarVCF = e.getString("RS");
			if(idFromClinVarVCF == null)
			{
				continue;
			}
			if(idFromClinVarVCF.equals(rsid))
			{
				return e.getString("#CHROM") + "\t" + e.getString("POS") + "\t" + rsid + "\t" + e.getString("REF") + "\t" + e.getString("ALT");
			}
		}
		return null;
	}
}
