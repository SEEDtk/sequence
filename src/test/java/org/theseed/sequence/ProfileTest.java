/**
 *
 */
package org.theseed.sequence;

import static org.hamcrest.MatcherAssert.assertThat;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.theseed.genome.Feature;
import org.theseed.genome.FeatureList;
import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.sequence.blast.BlastDB;
import org.theseed.sequence.blast.BlastHit;
import org.theseed.sequence.blast.BlastParms;
import org.theseed.sequence.blast.DnaBlastDB;
import org.theseed.sequence.blast.ProteinBlastDB;
import org.theseed.sequence.blast.ProteinProfiles;

import junit.framework.TestCase;
import static org.hamcrest.Matchers.*;

/**
 * @author Bruce Parrello
 *
 */
public class ProfileTest extends TestCase {

    private final File tempDir = new File("data", "temp");

    public ProfileTest() {
    }

    /**
     * @param name
     */
    public ProfileTest(String name) {
        super(name);
    }

    public void testPsiBlast() throws IOException, InterruptedException {
        Genome g2 = new Genome(new File("data", "1685.390.gto"));
        BlastDB gdb = DnaBlastDB.create(new File(tempDir, "pblast.fa"), g2);
        File pFile = new File("data", "PhenTrnaSyntAlph.smp");
        Map<String,String> qMap = new HashMap<String, String>();
        qMap.put("PhenTrnaSyntAlph", "Phenylalanyl-tRNA synthetase alpha chain (EC 6.1.1.20)");
        BlastParms parms = new BlastParms().maxE(1e-10);
        List<BlastHit> results = gdb.psiBlast(pFile, parms, qMap);
        assertThat(results.size(), equalTo(1));
        Feature feat = g2.getFeature("fig|1685.390.peg.2038");
        assertTrue(feat.getLocation().contains(results.get(0).getSubjectLoc()));
        gdb = ProteinBlastDB.create(new File(tempDir, "pblast2.fa"), g2);
        results = gdb.psiBlast(pFile, parms, qMap);
        assertThat(results.size(), equalTo(1));
        assertThat(results.get(0).getSubjectId(), equalTo("fig|1685.390.peg.2038"));
    }

    public void testProfiles() throws IOException, InterruptedException {
        Genome g2 = new Genome(new File("data", "1685.390.gto"));
        DnaBlastDB gdb = DnaBlastDB.create(new File(tempDir, "pblast.fa"), g2);
        BlastParms parms = new BlastParms().maxE(1e-10).minQueryBitScore(1.5);
        ProteinProfiles profiler = new ProteinProfiles(new File("data", "Profiles"));
        Map<String, List<BlastHit>> hitMap = profiler.profile(gdb, parms);
        assertThat(hitMap.size(), equalTo(3));
        for (String contigId : hitMap.keySet()) {
            // Get all the features on this contig.
            FeatureList pegs = g2.getContigFeatures(contigId);
            // For each hit, get the features that correspond to it.
            for (BlastHit hit : hitMap.get(contigId)) {
                Location hitLoc = hit.getSubjectLoc();
                Collection<Feature> feats = pegs.inRegion(hitLoc.getLeft(), hitLoc.getRight());
                for (Feature peg : feats) {
                    assertThat(hit.getQueryDef(), equalTo(peg.getFunction()));
                }
            }

        }
    }

}
