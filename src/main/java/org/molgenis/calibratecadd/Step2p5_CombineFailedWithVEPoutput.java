package org.molgenis.calibratecadd;

import org.molgenis.data.Entity;
import org.molgenis.data.annotator.tabix.TabixVcfRepository;
import org.molgenis.gavin2.util.FixNotation;
import org.molgenis.gavin2.util.FixVcfNotation;

import java.io.File;
import java.io.PrintWriter;
import java.util.*;

public class Step2p5_CombineFailedWithVEPoutput
{


	/**
	 * args[0] : VCF downloaded output from Ensembl VEP (http://grch37.ensembl.org/Homo_sapiens/Tools/VEP) on HGVS notation of the variants that failed (manually taken from clinvar.patho.fix.vcffailed of Step2)
	 * args[1] : the variants from Step 2 that failed (e.g. clinvar.patho.fix.vcffailed, original file)
	 * @param args
	 * @throws Exception
     */
	public static void main(String[] args) throws Exception
	{
		Step2p5_CombineFailedWithVEPoutput step2 = new Step2p5_CombineFailedWithVEPoutput(args[0], args[1]);
	}

	public Step2p5_CombineFailedWithVEPoutput(String args0, String args1) throws Exception {
		Scanner s = new Scanner(new File(args0));

		HashMap<String,String> chromPosToRefAlt = new HashMap<String,String>();
		Set<String> dupChromPos = new HashSet<String>();

		String line;
		while (s.hasNextLine()) {
			line = s.nextLine();
			if (line.startsWith("#")) {
				continue;
			}

			String[] lineSplit = line.split("\t", -1);

			String chrom = lineSplit[0];
			String pos = lineSplit[1];
			String ref = lineSplit[3];
			String alt = lineSplit[4];

			String key = chrom+"_"+pos;
			if(!chromPosToRefAlt.containsKey(key) && !dupChromPos.contains(key))
			{
				chromPosToRefAlt.put(key, ref+"\t"+alt);
			}
			else{
				System.out.println("BEWARE: duplicate chrom/pos found, you need to manually fix the output here!! " + line);
				dupChromPos.add(key);
				chromPosToRefAlt.remove(key);
			}
		}

		System.out.println("chromPosToRefAlt size " +  chromPosToRefAlt.size());
		PrintWriter pw = new PrintWriter(new File(args0) + ".rescued");

		// majority is problems is with deletions that are now annotated at position -1, so assign this first
		s = new Scanner(new File(args1));
		while (s.hasNextLine()) {
			line = s.nextLine();
			String[] lineSplit = line.split("\t", -1);
			String chrom = lineSplit[0]; int pos = Integer.parseInt(lineSplit[1]); String rsid = lineSplit[2]; String info = lineSplit[7];

		//	System.out.println("line" + line);

			String key = "";
			if(info.contains("del"))
			{
				key = chrom+"_"+(pos-1);
			}
			else
			{
			//	System.out.println("");
			}

			if(chromPosToRefAlt.containsKey(key)){
				pw.println(chrom + "\t" + pos + "\t" + rsid + "\t" + chromPosToRefAlt.get(key) + "\t" + "." + "\t" + "." + "\t" + info);
				chromPosToRefAlt.remove(key);
				//System.out.println("match: " + line);
			}
		}
		System.out.println("there are " + chromPosToRefAlt.size() + " left");
		System.out.println("unmapped: " + chromPosToRefAlt.keySet());


		pw.flush();
		pw.close();
	}
}
