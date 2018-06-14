package org.molgenis.calibratecadd;

import org.molgenis.data.annotation.entity.impl.gavin.GavinAlgorithm;
import org.molgenis.data.annotation.entity.impl.gavin.GavinEntry;
import org.molgenis.data.annotation.entity.impl.gavin.GavinEntry.Category;

import java.io.File;
import java.io.PrintWriter;
import java.util.HashMap;

/**
 * Created by joeri on 9/22/16.
 */
public class GavinCreateRuleGuide {

    public static void main(String[] args) throws Exception
    {
        if(args.length != 2)
        {
            throw new Exception("please provide: gavin file (e.g. /github/gavin/data/predictions), and version, (e.g. r0.2)");
        }

        String version = args[1];
        String gavinCalibrationFile = args[0] + File.separatorChar + "GAVIN_calibrations_"+version+".tsv";
        String output = args[0] + File.separatorChar + "GAVIN_ruleguide_"+version+".tsv";
        File outputFile = new File(output);

        //dev
        outputFile.delete();

        outputFile.createNewFile();
        new GavinCreateRuleGuide(gavinCalibrationFile, outputFile).createRuleGuide(version);
    }

    HashMap<String, GavinEntry> gavinData;
    PrintWriter pw;

    public GavinCreateRuleGuide(String gavinCalibrationFile, File outputFile) throws Exception
    {
        this.gavinData = Benchmark.loadGAVIN(gavinCalibrationFile).getGeneToEntry();
        this.pw = new PrintWriter(outputFile);
    }

    public void createRuleGuide(String version) throws Exception {

        String sep = "\t";
        String NA = "n/a";

        double GENOMEWIDE_MAF_THRESHOLD = -1;
        if(version.equals("r0.1")) { GENOMEWIDE_MAF_THRESHOLD = GavinAlgorithm.GENOMEWIDE_MAF_THRESHOLD_r0_1; }
        if(version.equals("r0.2")) { GENOMEWIDE_MAF_THRESHOLD = GavinAlgorithm.GENOMEWIDE_MAF_THRESHOLD_r0_2; }
        if(version.equals("r0.3")) { GENOMEWIDE_MAF_THRESHOLD = GavinAlgorithm.GENOMEWIDE_MAF_THRESHOLD_r0_3; }
        if(version.equals("r0.4")) { GENOMEWIDE_MAF_THRESHOLD = GavinAlgorithm.GENOMEWIDE_MAF_THRESHOLD_r0_4; }
        if(version.equals("r0.5")) { GENOMEWIDE_MAF_THRESHOLD = GavinAlgorithm.GENOMEWIDE_MAF_THRESHOLD_r0_5; }

        pw.println("################################");
        pw.println("## GAVIN Applied Rule Guide " + version);
        pw.println("################################");
        pw.println("## ");
        pw.println("## Variant can be interpreted by following the columns from left to right.");
        pw.println("## This classification scheme was implemented in Java, and used to");
        pw.println("## benchmark GAVIN in the paper (Van der Velde et al., using r0.2), see:");
        pw.println("## https://github.com/molgenis/molgenis and https://github.com/molgenis/gavin");
        pw.println("## ");
        pw.println("## Genome-wide rules are used if the gene-specific rules fail to classify.");
        pw.println("## These rules are applied as follows:");
        pw.println("## 1) If impact equals MODIFIER -> benign,");
        pw.println("## 2) if MAF greater than " + GENOMEWIDE_MAF_THRESHOLD + " -> benign,");
        pw.println("## 3) if CADD greater than "+GavinAlgorithm.GENOMEWIDE_CADD_THRESHOLD+" -> pathogenic,");
        pw.println("## 4) if CADD less than "+GavinAlgorithm.GENOMEWIDE_CADD_THRESHOLD+" -> benign.");
        pw.println("## ");
        pw.println("## Explanation of the gene calibration categories:");
        pw.println("## C1 = CADD scores highly significantly predictive for pathogenicity (pval < 0.01).");
        pw.println("## C2 = CADD scores significantly predictive for pathogenicity (pval < 0.05).");
        pw.println("## C3 = CADD scores may be predictive for pathogenicity (pval > 0.05 but with few samples).");
        pw.println("## C4 = CADD scores less predictive for pathogenicity (pval > 0.05 with enough samples).");
        pw.println("## C5 = CADD scores less predictive for pathogenicity (population CADD > pathogenic CADD).");
        pw.println("## I1 = HIGH impact unique for, thus predictive, for pathogenic variants.");
        pw.println("## I2 = MODERATE or HIGH impact unique, thus predictive, for pathogenic variants.");
        pw.println("## I3 = LOW or MODERATE or HIGH impact unique, thus predictive, for pathogenic variants.");
        pw.println("## T1 = Too few exac variants after filtering with pathogenic 95th percentile MAF.");
        pw.println("## T2 = Too few exac variants after filtering with impact distribution.");
        pw.println("## N1 = Too few ClinVar variants for calibration at this time.");
        pw.println("## N2 = Too few ExAC variants found for calibration.");
        pw.println("## ");
        pw.println("## For C1 and C2, CADD score thresholds are means of stratified benign and pathogenic variants.");
        pw.println("## For C3, C4 and C5, CADD score thresholds are 95th sensitivity/specificity percentiles of stratified benign and pathogenic variants.");
        pw.println("## ");
        pw.println("Gene" + sep + "PathogenicIfCADDScoreGreaterThan" + sep + "BenignIfCADDScoreLessThan" + sep + "BenignIfMAFGreaterThan" + sep + "PathogenicIfImpactEquals" + sep + "CalibrationCategory" );

        for(String gene : gavinData.keySet())
        {
            Category category = gavinData.get(gene).category;
            Double pathoMAFThreshold = gavinData.get(gene).PathoMAFThreshold != null ?
                    gavinData.get(gene).PathoMAFThreshold * GavinAlgorithm.extraSensitivityFactor * 2 : null;
            Double meanPathogenicCADDScore = gavinData.get(gene).MeanPathogenicCADDScore != null ?
                    gavinData.get(gene).MeanPathogenicCADDScore - GavinAlgorithm.extraSensitivityFactor : null;
            Double meanPopulationCADDScore = gavinData.get(gene).MeanPopulationCADDScore != null ?
                    gavinData.get(gene).MeanPopulationCADDScore - GavinAlgorithm.extraSensitivityFactor : null;
            Double spec95thPerCADDThreshold = gavinData.get(gene).Spec95thPerCADDThreshold != null ?
                    gavinData.get(gene).Spec95thPerCADDThreshold - GavinAlgorithm.extraSensitivityFactor : null;
            Double sens95thPerCADDThreshold = gavinData.get(gene).Sens95thPerCADDThreshold != null ?
                    gavinData.get(gene).Sens95thPerCADDThreshold - GavinAlgorithm.extraSensitivityFactor : null;

            pathoMAFThreshold = pathoMAFThreshold != null ? Math.min(pathoMAFThreshold, 1) : null;

            switch (category) {
                case C1: case C2:
                    pw.println(gene + sep + Math.max(meanPathogenicCADDScore, 0) + sep + Math.max(meanPopulationCADDScore, 0) + sep + pathoMAFThreshold + sep + NA + sep + category);
                    break;
                case C3: case C4: case C5:
                    pw.println(gene + sep + Math.max(spec95thPerCADDThreshold, 0) + sep + Math.max(sens95thPerCADDThreshold, 0) + sep + pathoMAFThreshold + sep + NA + sep + category);
                    break;
                case I1: case I2: case I3:
                    pw.println(gene + sep + NA + sep + NA + sep + pathoMAFThreshold + sep + (category == Category.I1 ? "HIGH" : category == Category.I2 ? "HIGH,MODERATE" : category == Category.I3 ? "HIGH,MODERATE,LOW" : null) + sep + category);
                    break;
                case T1: case T2:
                    pw.println(gene + sep + NA + sep + NA + sep + pathoMAFThreshold + sep + NA + sep + category);
                    break;
                case N1: case N2: case N3:
                    pw.println(gene + sep + NA + sep + NA + sep + NA + sep + NA + sep + category);
                    break;

                default:
                    throw new Exception("Invalid category: " + category);
            }
        }

        pw.flush();
        pw.close();
    }


}
