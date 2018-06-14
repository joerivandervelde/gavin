package org.molgenis.calibratecadd.support;

import java.io.*;
import java.util.HashMap;
import java.util.Scanner;
import java.util.zip.GZIPInputStream;

import org.apache.commons.lang3.StringUtils;

public class LoadCADDWebserviceOutput
{

	/**
	 * "chr_pos_ref_alt" to CADD PHRED score
	 * @param caddFile
	 * @return
	 * @throws FileNotFoundException
	 */
	public static HashMap<String, Double> load(File caddFile) throws Exception
	{
		Scanner cadd = new Scanner(caddFile);

		HashMap<String, Double> caddScores = new HashMap<String, Double>();

		String line = null;
		while(cadd.hasNextLine())
		{
			line = cadd.nextLine();
			if(line.startsWith("#"))
			{
				continue;
			}
			String[] split = line.split("\t", -1);
			if(split.length != 6)
			{
				throw new Exception("Expected 6 columns in CADD webservice output file, found " + split.length);
			}
			caddScores.put(split[0] + "_" + split[1] + "_" + split[2] + "_" + split[3], Double.parseDouble(split[5]));
		}
		cadd.close();
		return caddScores;
	}


	/**
	 * "chr_pos_ref_alt" to full line of annotations (116 elements)
	 * from GZ file
	 * @param caddFile
	 * @return
	 * @throws FileNotFoundException
	 */
	public static HashMap<String, String> loadAnnot(File caddFile) throws Exception
	{
		InputStream fileStream = new FileInputStream(caddFile);
		InputStream gzipStream = new GZIPInputStream(fileStream);
		Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
		BufferedReader cadd = new BufferedReader(decoder);

		HashMap<String, String> caddAnnot = new HashMap<String, String>();

		String line = null;
		while(true)
		{
			line = cadd.readLine();
			if(line == null)
			{
				break;
			}
			if(line.startsWith("#"))
			{
				continue;
			}
			String[] split = line.split("\t", -1);

			// sanity checks
			if(split.length != 116)
			{
				throw new Exception("Expected 116 columns in annotated CADD webservice output file, found " + split.length);
			}

			/*
			column 2 is ref, 3 is "ancestral", 4 is alt
			1	949503	AC	AC	A
			1	949523	C	C	T
			 */
			if(split[2].equals(split[4]))
			{
				throw new Exception("Expected column 2 to be different from 4, found: " + split[3] + " vs " + split[4]);
			}

			// chrom, pos, ref, alt
			String key = split[0] + "_" + split[1] + "_" + split[2] + "_" + split[4];

			// add only once, although there are usually 2-3 lines for this key
			// the most severe/complete effect is on line 1
			if(!caddAnnot.containsKey(key))
			{
				caddAnnot.put(key, line);
			}
		}
		cadd.close();
		return caddAnnot;
	}


}
