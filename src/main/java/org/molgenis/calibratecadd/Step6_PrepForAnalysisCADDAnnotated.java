package org.molgenis.calibratecadd;

import org.molgenis.calibratecadd.support.LoadCADDWebserviceOutput;
import org.molgenis.data.vcf.utils.FixVcfAlleleNotation;

import java.io.File;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Scanner;

/**
 * Example usage:
 * E:\Data\clinvarcadd\cadd_output.txt
 * E:\Data\clinvarcadd\clinvar.patho.fix.snpeff.exac.vcf.info
 * E:\Data\clinvarcadd\clinvar.patho.fix.snpeff.exac.withcadd.tsv
 * 
 * 
 * Takes the CADD annotated output (Step 5), which looks like this:
 *
 * ## CADD v1.3 (c) University of Washington and Hudson-Alpha Institute for Biotechnology 2013-2015. All rights reserved.
 * #Chrom	Pos	Ref	Anc	Alt	Type	Length	isTv	isDerived	AnnoType	Consequence	ConsScore	ConsDetail	GC	CpG	mapAbility20bp	mapAbility35bp	scoreSegDup	priPhCons	mamPhCons	verPhCons	priPhyloP	mamPhyloP	verPhyloP	GerpN	GerpS	GerpRS	GerpRSpval	bStatistic	mutIndex	dnaHelT	dnaMGW	dnaProT	dnaRoll	mirSVR-Score	mirSVR-E	mirSVR-Aln	targetScan	fitCons	cHmmTssA	cHmmTssAFlnk	cHmmTxFlnk	cHmmTx	cHmmTxWk	cHmmEnhG	cHmmEnh	cHmmZnfRpts	cHmmHet	cHmmTssBiv	cHmmBivFlnk	cHmmEnhBiv	cHmmReprPC	cHmmReprPCWk	cHmmQuies	EncExp	EncH3K27Ac	EncH3K4Me1	EncH3K4Me3	EncNucleo	EncOCC	EncOCCombPVal	EncOCDNasePVal	EncOCFairePVal	EncOCpolIIPVal	EncOCctcfPVal	EncOCmycPVal	EncOCDNaseSig	EncOCFaireSig	EncOCpolIISig	EncOCctcfSig	EncOCmycSig	Segway	tOverlapMotifs	motifDist	motifECount	motifEName	motifEHIPos	motifEScoreChng	TFBS	TFBSPeaks	TFBSPeaksMax	isKnownVariant	ESP_AF	ESP_AFR	ESP_EUR	TG_AF	TG_ASN	TG_AMR	TG_AFR	TG_EUR	minDistTSS	minDistTSE	GeneID	FeatureID	CCDS	GeneName	cDNApos	relcDNApos	CDSpos	relCDSpos	protPos	relProtPos	Domain	Dst2Splice	Dst2SplType	Exon	Intron	oAA	nAA	Grantham	PolyPhenCat	PolyPhenVal	SIFTcat	SIFTval	RawScore	PHRED
 * 1	949490	C	C	G	SNV	0	TRUE	TRUE	CodingTranscript	NON_SYNONYMOUS	7	missense	0.68	0.13	1	1	NA	0.306	0.944	0.937	0.557	0.543	0.599	4.25	2.34	51.8	8.07672e-10	934	117	2.43	0.13	0.86	-1.42	NA	NA	NA	NA	0.550478	0.126	0.622	0.031	0.000	0.016	0.000	0.110	0.000	0.000	0.016	0.055	0.008	0.000	0.008	0.008	3826.84	53.16	24.56	77.64	0.30	1	9.350	7.100	3.660	6.550	1.310	0.000	0.19	0.03	0.12	0.11	0.02	GS	NA	NA	NA	NA	NA	NA	2	2	37.0320	FALSE	NA	NA	NA	NA	NA	NA	NA	NA	687	431	ENSG00000187608	ENST00000379389	CCDS6.1	ISG15	281	0.40	130	0.26	44	0.27	ndomain	NA	NA	2/2	NA	R	G	125	probably_damaging	0.952	deleterious	0.03	5.837416	27.2
 * 1	949490	C	C	G	SNV	0	TRUE	TRUE	Intergenic	UPSTREAM	1	upstream	0.68	0.13	1	1	NA	0.306	0.944	0.937	0.557	0.543	0.599	4.25	2.34	51.8	8.07672e-10	934	117	2.43	0.13	0.86	-1.42	NA	NA	NA	NA	0.550478	0.126	0.622	0.031	0.000	0.016	0.000	0.110	0.000	0.000	0.016	0.055	0.008	0.000	0.008	0.008	3826.84	53.16	24.56	77.64	0.30	1	9.350	7.100	3.660	6.550	1.310	0.000	0.19	0.03	0.12	0.11	0.02	GS	NA	NA	NA	NA	NA	NA	2	2	37.0320	FALSE	NA	NA	NA	NA	NA	NA	NA	NA	687	431	ENSG00000224969	ENST00000458555	NA	RP11-54O7.11	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	5.837416	27.2
 * 1	949490	C	C	G	SNV	0	TRUE	TRUE	RegulatoryFeature	REGULATORY	4	regulatory	0.68	0.13	1	1	NA	0.306	0.944	0.937	0.557	0.543	0.599	4.25	2.34	51.8	8.07672e-10	934	117	2.43	0.13	0.86	-1.42	NA	NA	NA	NA	0.550478	0.126	0.622	0.031	0.000	0.016	0.000	0.110	0.000	0.000	0.016	0.055	0.008	0.000	0.008	0.008	3826.84	53.16	24.56	77.64	0.30	1	9.350	7.100	3.660	6.550	1.310	0.000	0.19	0.03	0.12	0.11	0.02	GS	NA	NA	NA	NA	NA	NA	2	2	37.0320	FALSE	NA	NA	NA	NA	NA	NA	NA	NA	687	431	NA	ENSR00000528877	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	5.837416	27.2
 * etc
 *
 *
 * Plus all accompanying information outputted directly by Step 4:
 * 
 * gene	chr	pos	ref	alt	group
 * IFT172	2	27680627	A	T	PATHOGENIC	splice_region_variant&intron_variant	LOW
 * IFT172	2	27700177	A	T	PATHOGENIC	missense_variant	MODERATE
 * IFT172	2	27693963	C	T	PATHOGENIC	splice_acceptor_variant&intron_variant	HIGH
 * etc
 * 
 * 
 * Combine / clean this into:
 * 
 * gene	chr	pos	ref	alt	group	effect	impact	cadd	#Chrom	Pos	Ref	Anc	Alt	Type	Length	isTv	isDerived	AnnoType	Consequence	ConsScore	ConsDetail	GC	CpG	mapAbility20bp	mapAbility35bp	scoreSegDup	priPhCons	mamPhCons	verPhCons	priPhyloP	mamPhyloP	verPhyloP	GerpN	GerpS	GerpRS	GerpRSpval	bStatistic	mutIndex	dnaHelT	dnaMGW	dnaProT	dnaRoll	mirSVR-Score	mirSVR-E	mirSVR-Aln	targetScan	fitCons	cHmmTssA	cHmmTssAFlnk	cHmmTxFlnk	cHmmTx	cHmmTxWk	cHmmEnhG	cHmmEnh	cHmmZnfRpts	cHmmHet	cHmmTssBiv	cHmmBivFlnk	cHmmEnhBiv	cHmmReprPC	cHmmReprPCWk	cHmmQuies	EncExp	EncH3K27Ac	EncH3K4Me1	EncH3K4Me3	EncNucleo	EncOCC	EncOCCombPVal	EncOCDNasePVal	EncOCFairePVal	EncOCpolIIPVal	EncOCctcfPVal	EncOCmycPVal	EncOCDNaseSig	EncOCFaireSig	EncOCpolIISig	EncOCctcfSig	EncOCmycSig	Segway	tOverlapMotifs	motifDist	motifECount	motifEName	motifEHIPos	motifEScoreChng	TFBS	TFBSPeaks	TFBSPeaksMax	isKnownVariant	ESP_AF	ESP_AFR	ESP_EUR	TG_AF	TG_ASN	TG_AMR	TG_AFR	TG_EUR	minDistTSS	minDistTSE	GeneID	FeatureID	CCDS	GeneName	cDNApos	relcDNApos	CDSpos	relCDSpos	protPos	relProtPos	Domain	Dst2Splice	Dst2SplType	Exon	Intron	oAA	nAA	Grantham	PolyPhenCat	PolyPhenVal	SIFTcat	SIFTval	RawScore	PHRED
 * IFT172	2	27680627	A	T	PATHOGENIC	splice_region_variant	LOW	11.29	2	27680627	A	A	T	SNV	0	TRUE	TRUE	Intergenic	DOWNSTREAM	1	downstream	0.60	0.03	1	1	NA	0.956	0.978	0.123	0.530	0.724	1.017	6.08	4.12	1648.6	2.19966e-177	45	-26	-0.39	1.13	1.19	7.27	NA	NA	NA	NA	0.106103	0.000	0.000	0.000	0.709	0.252	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.008	0.031	5.38	4.72	5.00	3.64	1.70	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	TF1	NA	NA	NA	NA	NA	NA	NA	NA	NA	FALSE	NA	NA	NA	NA	NA	NA	NA	NA	1354	134	ENSG00000235267	ENST00000417130	NA	AC074117.13	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	1.113829	11.29
 * IFT172	2	27672182	G	GG	PATHOGENIC	frameshift_variant	HIGH	35	2	27672182	G	G	GG	INS	1	NA	TRUE	CodingTranscript	FRAME_SHIFT	7	frameshift,feature_elongation	0.50	0.00	1	1	NA	0.996	1.000	1.000	0.550	2.504	5.698	5.37	5.37	1794	1.36987e-94	45	16	-0.74	0.45	-1.87	6.76	NA	NA	NA	NA	0.723164	0.000	0.000	0.000	0.701	0.260	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.039	62.92	4.40	4.00	2.24	1.50	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	TF1	NA	NA	NA	NA	NA	NA	NA	NA	NA	FALSE	NA	NA	NA	NA	NA	NA	NA	NA	4772	400	ENSG00000138002	ENST00000260570	CCDS1755.1	IFT172	4393	0.81	4289	0.82	1430	0.82	NA	NA	NA	39/48	NA	NA	NA	NA	NA	NA	NA	NA	8.966829	35
 * IFT172	2	27703928	A	G	PATHOGENIC	missense_variant	MODERATE	29.1	2	27703928	A	A	G	SNV	0	FALSE	TRUE	CodingTranscript	NON_SYNONYMOUS	7	missense	0.50	0.04	1	1	NA	0.957	1.000	0.997	0.452	2.225	4.900	5.85	5.85	900.8	8.48132e-146	53	-17	0.46	0.26	1.85	2.57	NA	NA	NA	NA	0.666978	0.000	0.000	0.000	0.472	0.457	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.008	0.063	25.25	5.72	4.00	3.00	1.70	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	R5	NA	NA	NA	NA	NA	NA	NA	NA	NA	FALSE	NA	NA	NA	NA	NA	NA	NA	NA	8618	8718	ENSG00000138002	ENST00000260570	CCDS1755.1	IFT172	874	0.16	770	0.15	257	0.15	ndomain	16	DONOR	8/48	NA	L	P	98	possibly_damaging	0.847	deleterious	0.01	6.303687	29.1
 * etc
 * 
 *
 */
public class Step6_PrepForAnalysisCADDAnnotated
{

	public static void main(String[] args) throws Exception
	{
		HashMap<String, String> caddScores = LoadCADDWebserviceOutput.loadAnnot(new File(args[0]));
		Scanner variants = new Scanner(new File(args[1]));
		PrintWriter pw = new PrintWriter(new File(args[2]));
		
		String line;
		
		//write header of output
		// from "#Chrom" on this is from Annotated CADD output
		pw.println("gene" + "\t" + "chr" + "\t" + "pos" + "\t" + "ref" + "\t" + "alt" + "\t" + "group" + "\t" + "effect" + "\t" + "impact" + "\t" + "cadd" + "\t" + "#Chrom" + "\t" + "Pos" + "\t" + "Ref" + "\t" + "Anc" + "\t" + "Alt" + "\t" + "Type" + "\t" + "Length" + "\t" + "isTv" + "\t" + "isDerived" + "\t" + "AnnoType" + "\t" + "Consequence" + "\t" + "ConsScore" + "\t" + "ConsDetail" + "\t" + "GC" + "\t" + "CpG" + "\t" + "mapAbility20bp" + "\t" + "mapAbility35bp" + "\t" + "scoreSegDup" + "\t" + "priPhCons" + "\t" + "mamPhCons" + "\t" + "verPhCons" + "\t" + "priPhyloP" + "\t" + "mamPhyloP" + "\t" + "verPhyloP" + "\t" + "GerpN" + "\t" + "GerpS" + "\t" + "GerpRS" + "\t" + "GerpRSpval" + "\t" + "bStatistic" + "\t" + "mutIndex" + "\t" + "dnaHelT" + "\t" + "dnaMGW" + "\t" + "dnaProT" + "\t" + "dnaRoll" + "\t" + "mirSVR-Score" + "\t" + "mirSVR-E" + "\t" + "mirSVR-Aln" + "\t" + "targetScan" + "\t" + "fitCons" + "\t" + "cHmmTssA" + "\t" + "cHmmTssAFlnk" + "\t" + "cHmmTxFlnk" + "\t" + "cHmmTx" + "\t" + "cHmmTxWk" + "\t" + "cHmmEnhG" + "\t" + "cHmmEnh" + "\t" + "cHmmZnfRpts" + "\t" + "cHmmHet" + "\t" + "cHmmTssBiv" + "\t" + "cHmmBivFlnk" + "\t" + "cHmmEnhBiv" + "\t" + "cHmmReprPC" + "\t" + "cHmmReprPCWk" + "\t" + "cHmmQuies" + "\t" + "EncExp" + "\t" + "EncH3K27Ac" + "\t" + "EncH3K4Me1" + "\t" + "EncH3K4Me3" + "\t" + "EncNucleo" + "\t" + "EncOCC" + "\t" + "EncOCCombPVal" + "\t" + "EncOCDNasePVal" + "\t" + "EncOCFairePVal" + "\t" + "EncOCpolIIPVal" + "\t" + "EncOCctcfPVal" + "\t" + "EncOCmycPVal" + "\t" + "EncOCDNaseSig" + "\t" + "EncOCFaireSig" + "\t" + "EncOCpolIISig" + "\t" + "EncOCctcfSig" + "\t" + "EncOCmycSig" + "\t" + "Segway" + "\t" + "tOverlapMotifs" + "\t" + "motifDist" + "\t" + "motifECount" + "\t" + "motifEName" + "\t" + "motifEHIPos" + "\t" + "motifEScoreChng" + "\t" + "TFBS" + "\t" + "TFBSPeaks" + "\t" + "TFBSPeaksMax" + "\t" + "isKnownVariant" + "\t" + "ESP_AF" + "\t" + "ESP_AFR" + "\t" + "ESP_EUR" + "\t" + "TG_AF" + "\t" + "TG_ASN" + "\t" + "TG_AMR" + "\t" + "TG_AFR" + "\t" + "TG_EUR" + "\t" + "minDistTSS" + "\t" + "minDistTSE" + "\t" + "GeneID" + "\t" + "FeatureID" + "\t" + "CCDS" + "\t" + "GeneName" + "\t" + "cDNApos" + "\t" + "relcDNApos" + "\t" + "CDSpos" + "\t" + "relCDSpos" + "\t" + "protPos" + "\t" + "relProtPos" + "\t" + "Domain" + "\t" + "Dst2Splice" + "\t" + "Dst2SplType" + "\t" + "Exon" + "\t" + "Intron" + "\t" + "oAA" + "\t" + "nAA" + "\t" + "Grantham" + "\t" + "PolyPhenCat" + "\t" + "PolyPhenVal" + "\t" + "SIFTcat" + "\t" + "SIFTval" + "\t" + "RawScore" + "\t" + "PHRED");
		
		//skip header of input
		variants.nextLine();
		
		while(variants.hasNextLine())
		{
			line = variants.nextLine();
			String[] split = line.split("\t", -1);
			String gene = split[0];
			String chr = split[1];
			String pos = split[2];
			String ref = split[3];
			String alt = split[4];
			String group = split[5];
			String effect = split[6];
			String impact = split[7];

			effect = effect.substring(0, (effect.contains("&") ? effect.indexOf("&") : effect.length()));
			String printMe = (gene + "\t" + chr + "\t" + pos + "\t" + ref + "\t" + alt + "\t" + group + "\t" + effect + "\t" + impact).replace("/", "_");

			String key = chr + "_" + pos + "_" + ref + "_" + alt;
			if(caddScores.containsKey(key))
			{
				//FIXME: need to replace '/' to prevent problems in R later on when writing plots based on gene names..
				String[] caddAnnotSplit = caddScores.get(key).split("\t");
				pw.println(printMe + "\t" + caddAnnotSplit[115] + "\t" + caddScores.get(key));
			}
			else
			{
				String backTrimmedRefAlt = FixVcfAlleleNotation.backTrimRefAlt(ref, alt, "_");
				String btkey = chr + "_" + pos + "_" + backTrimmedRefAlt;

				String frontTrimmedRefAlt = FixVcfAlleleNotation.frontTrimRefAlt(ref, alt, "_");
				String ftkey = chr + "_" + pos + "_" + frontTrimmedRefAlt;

				if(caddScores.containsKey(btkey))
				{
					System.out.println("RESOLVED by back-trimming ref alt " + ref + "_" + alt + " to " + backTrimmedRefAlt);
					String[] caddAnnotSplit = caddScores.get(btkey).split("\t");
					pw.println(printMe + "\t" + caddAnnotSplit[115] + "\t" + caddScores.get(key));
				}
				else if(caddScores.containsKey(ftkey))
				{
					System.out.println("RESOLVED by front-trimming ref alt " + ref + "_" + alt + " to " + frontTrimmedRefAlt);
					String[] caddAnnotSplit = caddScores.get(ftkey).split("\t");
					pw.println(printMe + "\t" + caddAnnotSplit[115] + "\t" + caddScores.get(key));
				}
				else
				{
					System.out.println("WARNING: could not get CADD score for " + key + " using either " + backTrimmedRefAlt + " or " + frontTrimmedRefAlt);
				}


			}
		}

		variants.close();
		pw.flush();
		pw.close();
		
		
	}
}
