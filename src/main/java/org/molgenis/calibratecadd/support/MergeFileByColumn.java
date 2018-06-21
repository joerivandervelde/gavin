package org.molgenis.calibratecadd.support;

import java.io.File;
import java.io.PrintWriter;
import java.util.Scanner;

public class MergeFileByColumn {

    public static void main(String[] args) throws Exception
    {
        File a = new File(args[0]);
        File b = new File(args[1]);
        File m = new File(args[2]);

        Scanner as = new Scanner(a);
        Scanner bs = new Scanner(b);
        PrintWriter mw = new PrintWriter(m);

        while(as.hasNextLine() && bs.hasNextLine()){
            mw.println(as.nextLine() + bs.nextLine());
        }

        if(as.hasNextLine()) {System.out.println("WARNING: file " + a.getName() + " still has lines left!");}
        if(bs.hasNextLine()) {System.out.println("WARNING: file " + b.getName() + " still has lines left!");}

        as.close();
        bs.close();
        mw.flush();
        mw.close();
    }

}
