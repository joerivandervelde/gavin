package org.molgenis.calibratecadd;

import org.testng.annotations.Test;

/**
 * Created by joeri on 9/1/16.
 */
public class CADDValidationTest {

    @Test
    public void testPredictionTool() throws Exception {
        new AbstractValidationTest("CADD_Thr15").testPredictionTool();
    }

}
