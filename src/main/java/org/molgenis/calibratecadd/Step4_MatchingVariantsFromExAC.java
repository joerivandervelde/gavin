package org.molgenis.calibratecadd;

import org.molgenis.calibratecadd.structs.GeneCalibResult;
import org.molgenis.calibratecadd.support.*;
import org.molgenis.data.Entity;
import org.molgenis.data.annotation.entity.impl.gavin.GavinEntry;
import org.molgenis.data.annotator.tabix.TabixVcfRepository;
import org.molgenis.data.vcf.VcfRepository;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class Step4_MatchingVariantsFromExAC {

	/**
	 * Uses:
	 * [0] file produced in step 3
	 * [1] ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/ExAC.r0.3.sites.vep.vcf.gz (SnpEff annotated + in the same folder ExAC.r0.3.sites.vep.vcf.gz.tbi )
	 * [2] output file (produce 3 files: output.cadd, output.genes, output.variants)
	 * <p>
	 * Example:
	 * E:\Data\clinvarcadd\clinvar.patho.fix.snpeff.vcf
	 * E:\Data\clinvarcadd\ExAC.r0.3.sites.vep.vcf.gz
	 * E:\Data\clinvarcadd\clinvar.patho.fix.snpeff.exac.vcf
	 */
	public static void main(String[] args) throws Exception {

		File clinvarPathoLoc = new File(args[0]);
		File exacFile = new File(args[1]);
		File outFile = new File(args[2]);
		Step4_MatchingVariantsFromExAC step4 = new Step4_MatchingVariantsFromExAC(clinvarPathoLoc, exacFile, outFile);
		step4.go();

	}

	public void go() throws Exception {
		HashMap<String, List<Entity>> clinvarPatho = loadClinvarPatho(clinvarPathoLoc);
		System.out.println("loaded LP/P variants for " + clinvarPatho.size() + " genes");
		createMatchingExACsets(exacFile, clinvarPatho, outFile);

		//	printVariantsToFile(outFile, clinvarPatho, calibResults);


	}

	File clinvarPathoLoc;
	File exacFile;
	File outFile;

	public Step4_MatchingVariantsFromExAC(File clinvarPathoLoc, File exacFile, File outFile) throws Exception {
		this.clinvarPathoLoc = clinvarPathoLoc;
		this.exacFile = exacFile;
		this.outFile = outFile;
	}


	private void createMatchingExACsets(File exacFile, HashMap<String, List<Entity>> clinvarPatho, File outFile) throws Exception {
		PrintWriter pw_variantInfo = new PrintWriter(outFile + ".variants.tsv");
		PrintWriter pw_forCADD = new PrintWriter(outFile + ".cadd.tsv");
		PrintWriter pw_geneInfo = new PrintWriter(outFile + ".genes.tsv");

		pw_variantInfo.println( "gene" + "\t" + "chr" + "\t" + "pos" + "\t" + "ref" + "\t" + "alt" + "\t" + "group");
		pw_geneInfo.println( "Gene" + "\t" + "Category" + "\t" + "Chr" + "\t" + "Start" + "\t" + "End" + "\t" + "NrOfPopulationVariants" + "\t" + "NrOfPathogenicVariants" + "\t" + "NrOfOverlappingVariants" + "\t" + "NrOfFilteredPopVariants" + "\t" + "PathoMAFThreshold" + "\t" + "PopImpactHighPerc" + "\t" + "PopImpactModeratePerc" + "\t" + "PopImpactLowPerc" + "\t" + "PopImpactModifierPerc" + "\t" + "PathoImpactHighPerc" + "\t" + "PathoImpactModeratePerc" + "\t" + "PathoImpactLowPerc" + "\t" + "PathoImpactModifierPerc" + "\t" + "PopImpactHighEq" + "\t" + "PopImpactModerateEq" + "\t" + "PopImpactLowEq" + "\t" + "PopImpactModifierEq");


		Step4_Helper st4h = new Step4_Helper();
		System.out.println("loading matching exac variants..");

		int index = 0;

		nextGene:
		for (String gene : clinvarPatho.keySet()) {
			index++;


			/** DEV **/
//			if(!gene.equals("TYROBP"))
//			{
//				continue;
//			}
//			if(index % 100 == 0)
//			{
//				break;
//			}

			String chrom = null;
			long leftMostPos = -1;
			long rightMostPos = -1;

			for (Entity cvv : clinvarPatho.get(gene)) {
				long pos = cvv.getLong("POS");
				String cvvchrom = cvv.getString("#CHROM");
				if (pos > rightMostPos) {
					rightMostPos = pos;
				}

				if (pos < leftMostPos || leftMostPos == -1) {
					leftMostPos = pos;
				}
				if (chrom == null) {
					chrom = cvvchrom;
				}
			}

			// include (biggest) part of exon the variant(s) are in, typical exon is 147 nt
			// http://nar.oxfordjournals.org/content/early/2012/07/11/nar.gks652.full
			leftMostPos = leftMostPos - 100;
			rightMostPos = rightMostPos + 100;

			TabixVcfRepository tr = new TabixVcfRepository(exacFile, "exac");

			List<Entity> exacVariants = new ArrayList<Entity>();

			try {
				exacVariants = tr.query(chrom, leftMostPos, rightMostPos);
			} catch (java.lang.ArrayIndexOutOfBoundsException e) {
				// no chrom in tabix or so
			}

			System.out.println("\n#####\n");
			System.out.println(gene + " (" + index + " of " + clinvarPatho.keySet().size() + ") at " + chrom + ":" + leftMostPos + "-" + rightMostPos + " " + " has " + clinvarPatho.get(gene).size() + " patho variants, and " + exacVariants.size() + " exac sites (possibly multi-allelic)");


			//there are a lot of genes with only 1 pathogenic variant.. a bit silly to consider them seriously for calibration work, drop and report as N1 (N=1)
			if (clinvarPatho.get(gene).size() < 2) {
				System.out.println("DROPPED - Too few clinvar variants");
				GeneCalibResult gcr = new GeneCalibResult(new GavinEntry(gene, GavinEntry.Category.N1, chrom, leftMostPos, rightMostPos, null), null);
				printGeneToFile(gcr, pw_geneInfo);
				continue nextGene;
			}

			// too few exac variants
			if (exacVariants.size() < 2) {
				System.out.println("DROPPED - Too few exac variants");
				GeneCalibResult gcr = new GeneCalibResult(new GavinEntry(gene, GavinEntry.Category.N2, chrom, leftMostPos, rightMostPos, null), null);
				printGeneToFile(gcr, pw_geneInfo);
				continue nextGene;
			}

			// intersect the sets
			VariantIntersectResult vir = st4h.intersectVariants(exacVariants, clinvarPatho.get(gene), gene);

			System.out.println("INFO - VariantIntersectResult for '" + gene + "' unique variants, clinvaronly: " + vir.inClinVarOnly.size() + ", exaconly: " + vir.inExACOnly.size() + ", both: " + vir.inBoth_exac.size());

			//calculate MAF for shared variants, and use them to filter the other ExAC variants
			//this way, we use the overlap to determine a fair cutoff for the 'assumed benign' variants
			//if we have nothing to go on, we will set this to 0 and only select singleton variants
			double pathogenicMAF = st4h.calculatePathogenicMAF(vir.inBoth_exac, vir.inClinVarOnly.size());
			List<EntityPlus> exacFilteredByMAF = st4h.filterExACvariantsByMAF(vir.inExACOnly, pathogenicMAF);

			if (exacFilteredByMAF.size() == 0) {
				System.out.println("DROPPED - Too few exac variants after filtering with pathogenic MAF" + pathogenicMAF);
				GeneCalibResult gcr = new GeneCalibResult(new GavinEntry(gene, GavinEntry.Category.T1, chrom, leftMostPos, rightMostPos, pathogenicMAF), null);
				printGeneToFile(gcr, pw_geneInfo);
				continue nextGene;
			}

			System.out.println("INFO - Using pathogenic MAF, ExAC-only filtered down to " + exacFilteredByMAF.size() + " variants " + pathogenicMAF);


			//calculate impact ratios over all clinvar variants, and use them to 'shape' the remaining ExAC variants
			//they must become a set that looks just like the ClinVar variants, including same distribution of impact types
			ImpactRatios pathoImpactRatio = st4h.calculateImpactRatios(Stream.concat(vir.inClinVarOnly.stream(), vir.inBoth_clinvar.stream()).collect(Collectors.toList()), gene);
			System.out.println("INFO - ExAC-only impact ratio: " + st4h.calculateImpactRatios(vir.inExACOnly, gene));
			List<EntityPlus> exacFilteredByMAFandImpact = st4h.shapeExACvariantsByImpactRatios(exacFilteredByMAF, pathoImpactRatio);

			System.out.println("INFO - Using Pathogenic impact ratio " + pathoImpactRatio + ", ExAC-only filtered down to " + exacFilteredByMAFandImpact.size() + " variants");

			if (exacFilteredByMAFandImpact.size() == 0) {
				System.out.println("DROPPED - No match for impact profile, but we did learn something");
				ImpactRatios MAFfilteredExacImpactRatio = st4h.calculateImpactRatios(exacFilteredByMAF, gene);
				GavinEntry.Category cat = st4h.determineImpactFilterCat(MAFfilteredExacImpactRatio, pathoImpactRatio);
				GeneCalibResult gcr = new GeneCalibResult(new GavinEntry(gene, cat, chrom, leftMostPos, rightMostPos, pathogenicMAF), null);
				printGeneToFile(gcr, pw_geneInfo);
				continue nextGene;
			}


			GeneCalibResult gcr = new GeneCalibResult(new GavinEntry(gene, GavinEntry.Category.Cx, chrom, leftMostPos, rightMostPos, pathogenicMAF), exacFilteredByMAFandImpact);
			System.out.println("INFO - calibrated gene ! impact ratio of population reference: " + st4h.calculateImpactRatios(exacFilteredByMAFandImpact, gene));

			//print gene and variant info
			printGeneToFile(gcr, pw_geneInfo);
			printVariantsToFile(gene, clinvarPatho.get(gene), gcr, pw_variantInfo, pw_forCADD);



			pw_variantInfo.flush();
			pw_geneInfo.flush();
			pw_forCADD.flush();


		}


		pw_variantInfo.flush();
		pw_variantInfo.close();

		pw_geneInfo.flush();
		pw_geneInfo.close();

		pw_forCADD.flush();
		pw_forCADD.close();

		System.out.println();
		System.out.println("#### done ####");
		System.out.println();

	}

	/**
	 * @param calibResult
	 * @param pw_geneInfo
     */
	private void printGeneToFile(GeneCalibResult calibResult, PrintWriter pw_geneInfo)
	{
		pw_geneInfo.println(calibResult.toString());
	}

	private void printVariantsToFile(String gene, List<Entity> clinvarPatho, GeneCalibResult calibResult, PrintWriter pw_variantInfo, PrintWriter pw_forCADD) throws Exception
	{
		for(Entity variant : clinvarPatho)
		{
			pw_forCADD.println(variant.getString("#CHROM") + "\t" + variant.getString("POS") + "\t" + "." + "\t" + variant.getString("REF") + "\t" + variant.getString("ALT"));
			pw_variantInfo.println(gene + "\t" + variant.getString("#CHROM") + "\t" + variant.getString("POS") + "\t" + variant.getString("REF") + "\t" + variant.getString("ALT") + "\t" + "PATHOGENIC" + "\t" + GavinUtils.getEffect(variant.getString("ANN"), gene, variant.getString("ALT")) + "\t" + GavinUtils.getImpact(variant.getString("ANN"), gene, variant.getString("ALT")));
		}
		for(EntityPlus variant : calibResult.matchedVariants)
		{
			pw_forCADD.println(variant.getE().getString("#CHROM") + "\t" + variant.getE().getString("POS") + "\t" + "." + "\t"+ variant.getE().getString("REF") + "\t" + variant.getKeyVal().get("ALT").toString());
			pw_variantInfo.println(gene + "\t" + variant.getE().getString("#CHROM") + "\t" + variant.getE().getString("POS") + "\t" + variant.getE().getString("REF") + "\t" + variant.getKeyVal().get("ALT").toString() + "\t" + "POPULATION" + "\t" + variant.getKeyVal().get("EFFECT").toString() + "\t" + variant.getKeyVal().get("IMPACT").toString());
		}
	}


	private HashMap<String, List<Entity>> loadClinvarPatho(File clinvarPathoLoc) throws Exception {
		System.out.println("loading clinvar pathogenic variants from " + clinvarPathoLoc + " ..");
		VcfRepository vcfRepo = new VcfRepository(clinvarPathoLoc, "clinvar");
		HashMap<String, List<Entity>> clinvarPatho = new HashMap<String, List<Entity>>();
		java.util.Iterator<Entity> vcfRepoIter = vcfRepo.iterator();

		while (vcfRepoIter.hasNext())
		{
			Entity record = vcfRepoIter.next();
			if (record.getString("ANN") == null)
			{
				vcfRepo.close();
				throw new Exception("Please annotate the VCF with a recent snpEff version! ANN field not found for "
						+ record.toString());
			}

			if (record.getString("ALT").contains(","))
			{
				//System.out.println("skipping multi-allelic clinvar variant...");
				continue;
			}

			if (record.getString(Step1_GetClinVarPathogenic.CLINVAR_INFO) == null && record.getString(Step1_GetClinVarPathogenic.VKGL_INFO) == null)
			{
				vcfRepo.close();
				throw new Exception("Did you create this VCF running Step1? " + Step1_GetClinVarPathogenic.CLINVAR_INFO
						+ " or " + Step1_GetClinVarPathogenic.VKGL_INFO + "field not found for " + record.toString());
			}

			String geneAccordingToClinVarOrVKGL = null;

			if (record.getString(Step1_GetClinVarPathogenic.CLINVAR_INFO) != null)
			{
				geneAccordingToClinVarOrVKGL = record.getString(Step1_GetClinVarPathogenic.CLINVAR_INFO).split("\\|", -1)[1];
			}
			if (record.getString(Step1_GetClinVarPathogenic.VKGL_INFO) != null)
			{
				geneAccordingToClinVarOrVKGL = record.getString(Step1_GetClinVarPathogenic.VKGL_INFO).split(":", -1)[0];
			}

			//sanity check
			if(geneAccordingToClinVarOrVKGL.length() == 0)
			{
				vcfRepo.close();
				throw new Exception("geneAccordingToClinVar length 0");
			}

			Set<String> genesAccordingToSnpEff = GavinUtils.getGenesFromAnn(record.getString("ANN"));
			//TODO: require the gene to match clinvar/VKGL annotation?
			for(String snpEffGene : genesAccordingToSnpEff)
			{
				if(!snpEffGene.isEmpty()) {
					if (clinvarPatho.containsKey(snpEffGene)) {
						clinvarPatho.get(snpEffGene).add(record);
					} else {
						List<Entity> variants = new ArrayList<Entity>();
						variants.add(record);
						clinvarPatho.put(snpEffGene, variants);
					}
				}
			}
		}

		vcfRepo.close();

		return clinvarPatho;
	}





}
