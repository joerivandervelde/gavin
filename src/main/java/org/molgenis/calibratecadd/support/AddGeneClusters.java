package org.molgenis.calibratecadd.support;


import java.io.File;
import java.io.PrintWriter;
import java.util.*;


public class AddGeneClusters {

    /**
     * arg[0] is output from APcluster, ie.
     *
     * Cluster 82, exemplar CENPF:
     CENPF
     Cluster 83, exemplar CEP78:
     AHI1 ALMS1 C2CD3 CENPJ CEP135 CEP152 CEP164 CEP290 CEP41 CEP57 CEP63
     CEP78 MKS1 NDE1 OFD1 PCNT PLK4 SDCCAG8
     *
     * args[1] is GAVIN variants, assumed is a header and Gene in column 1, ie.
     * gene	chr	pos	ref	alt	group	effect	impact	cadd
     NUP107	12	69107589	G	A	PATHOGENIC	splice_dono
     NUP107	12	69090619	G	T	PATHOGENIC	missense_va
     NUP107	12	69084526	G	A	PATHOGENIC	missense_va
     NUP107	12	69129114	A	C	PATHOGENIC	missense_va

     args[2] is output file location

     we add clusters to each row
     with one quirk: singletons are put in a special big cluster of their own!




     * @param args
     * @throws Exception
     */
    public static void main(String[] args) throws Exception
    {
        new AddGeneClusters(new File(args[0]), new File(args[1]), new File(args[2]));
    }

    public AddGeneClusters(File APClustResults, File variants, File out) throws Exception {

        HashMap<String, List<String>> exemplarToClusterMembers = getExemplarToClusterMembers(APClustResults);
        HashMap<String, List<String>> exemplarToClusterMembersST = singletonsToSpecialCluster(exemplarToClusterMembers);
        HashMap<String, String> clusterMembersToExemplar = getClusterMemberToExemplar(exemplarToClusterMembersST);

        PrintWriter pw = new PrintWriter(out);

        Scanner s = new Scanner(variants);

        String header = s.nextLine(); //print header back + cluster
        pw.println("cluster" + "\t" + header);

        HashMap<String, Integer> variantsPerCluster = new  HashMap<String, Integer>();

        while(s.hasNextLine())
        {
            String line = s.nextLine();
            String[] split = line.split("\t", -1);
            String gene = split[0];

            String cluster = "UNKNOWN";
            if(clusterMembersToExemplar.containsKey(gene))
            {
                cluster = clusterMembersToExemplar.get(gene);
            }

            if(variantsPerCluster.containsKey(cluster))
            { variantsPerCluster.put(cluster, variantsPerCluster.get(cluster) + 1); }
            else { variantsPerCluster.put(cluster, 1); }

            pw.println("CP_CLUST_" + cluster + "\t" + line);
        }

        for(String cluster : variantsPerCluster.keySet())
        {
            System.out.println(cluster + " - " + variantsPerCluster.get(cluster));
        }

        s.close();
        pw.flush();
        pw.close();

    }

    /**
     * put all singletons in one new big cluster
     * @param exemplarToClusterMembers
     * @return
     * @throws Exception
     */
    public HashMap<String, List<String>> singletonsToSpecialCluster(HashMap<String, List<String>> exemplarToClusterMembers) throws Exception {

        List<String> singletons = new ArrayList<String>();
        for(String exemplar : exemplarToClusterMembers.keySet()) {
            System.out.println(exemplar + " size = " + exemplarToClusterMembers.get(exemplar).size());
            // for size == 1, exemplar is also the only member
           if(exemplarToClusterMembers.get(exemplar).size() == 1)
           {
               singletons.add(exemplar);
           }
        }

        // now remove these exemplar groups
        exemplarToClusterMembers.keySet().removeAll(singletons);

        exemplarToClusterMembers.put("SINGLETONS", singletons);

        System.out.println("new singtleton cluster: " + exemplarToClusterMembers.get("SINGLETONS").size());
        System.out.println("now there are " + exemplarToClusterMembers.keySet().size() + " clusters");

        return exemplarToClusterMembers;
    }


        /**
         * Note that members also includes the exemplar itself, so no additional lookup is needed
         * so if B -> A, C, D, then A->B, B->B, C->B, D->B
         * @param exemplarToClusterMembers
         * @return
         */
    public HashMap<String, String> getClusterMemberToExemplar(HashMap<String, List<String>> exemplarToClusterMembers) throws Exception {
        HashMap<String, String> clusterMembersToExemplar = new HashMap<String, String>();

        for(String exemplar : exemplarToClusterMembers.keySet())
        {

            for(String member : exemplarToClusterMembers.get(exemplar))
            {
                if(clusterMembersToExemplar.containsKey(member))
                {
                    throw new Exception("member as key already in map: " + member);
                }
                // examplars are also members
                clusterMembersToExemplar.put(member, exemplar);
            }
        }
        return clusterMembersToExemplar;
    }


    public HashMap<String, List<String>> getExemplarToClusterMembers(File APClustResults) throws Exception {
        HashMap<String, List<String>> exemplarToCluster = new HashMap<String, List<String>>();
        Scanner apresReader = new Scanner(APClustResults);

        String currentExemplar = null;
        String line;
        boolean startFound = false;
        while (apresReader.hasNextLine()) {
            line = apresReader.nextLine();
            if (!startFound) {
                if (line.startsWith("Clusters")) {
                    startFound = true;
                    System.out.println("Start found at line: " + line);
                    continue;
                } else {
                    continue;
                }
            }

            line = line.trim();
            if (line.startsWith("Cluster")) {
                String exemplar = line.substring(0, line.length() - 1).replaceAll("Cluster [0-9]+?, exemplar ", "");
                System.out.println("Current exemplar: " + exemplar + ", found in line " + line);
                currentExemplar = exemplar;
                ArrayList<String> cluster = new ArrayList<String>();
                exemplarToCluster.put(exemplar, cluster);
                continue;
            }

            // line with members, e.g. "AHI1 ALMS1 C2CD3 CENPJ CEP135 CEP152 CEP164"
            String[] members = line.split(" ", -1);
            exemplarToCluster.get(currentExemplar).addAll(new ArrayList<String>(Arrays.asList(members)));
        }
        System.out.println("Done, loaded " + exemplarToCluster.size() + " clusters");
        apresReader.close();
        return exemplarToCluster;
    }


}

