/**
 *
 */
package org.theseed.sequence;

import java.io.File;
import java.io.IOException;

import org.apache.commons.io.FileUtils;
import org.theseed.genome.Genome;
import org.theseed.sequence.blast.BlastDB;
import org.theseed.sequence.blast.BlastParms;

import junit.framework.TestCase;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;
import static org.hamcrest.io.FileMatchers.*;

/**
 *
 * Tests for the blast support
 *
 * @author Bruce Parrello
 *
 */
public class BlastTest extends TestCase {

    public void testFasta() throws IOException, InterruptedException {
        File gtoFile = new File("src/test", "1313.7001.gto");
        Genome gto = new Genome(gtoFile);
        File tempDir = new File("src/test", "temp");
        if (! tempDir.isDirectory()) tempDir.mkdirs();
        FileUtils.cleanDirectory(tempDir);
        // Create and test a DNA database.
        File fastaFile = new File(tempDir, "temp.fna");
        gto.saveDna(fastaFile);
        BlastDB newBlastDb = new BlastDB(fastaFile, BlastDB.Type.DNA);
        String[] suffixes = new String[] { ".nhr", ".nin", ".nsq" };
        for (String suffix : suffixes) {
            File testFile = new File(tempDir, "temp.fna" + suffix);
            assertThat(testFile, aReadableFile());
        }
        assertThat(newBlastDb.getType(), equalTo(BlastDB.Type.DNA));
        BlastDB oldBlastDb = new BlastDB(fastaFile);
        assertThat(oldBlastDb.getType(), equalTo(BlastDB.Type.DNA));
        // Create and test a protein database.
        File protFile = new File(tempDir, "temp.faa");
        gto.savePegs(protFile);
        BlastDB protBlastDb = new BlastDB(protFile, BlastDB.Type.PROTEIN);
        assertThat(protBlastDb.getType(), equalTo(BlastDB.Type.PROTEIN));
        suffixes = new String[] { ".phr", ".pin", ".psq" };
        for (String suffix : suffixes) {
            File testFile = new File(tempDir, "temp.faa" + suffix);
            assertThat(testFile, aReadableFile());
        }
        oldBlastDb = new BlastDB(protFile);
        assertThat(oldBlastDb.getType(), equalTo(BlastDB.Type.PROTEIN));
        // Verify that updating the fasta file causes a regen.
        File checkFile = new File(tempDir, "temp.fna.nsq");
        Thread.sleep(1000);
        gto.saveDna(fastaFile);
        assertThat(checkFile.lastModified(), lessThan(fastaFile.lastModified()));
        oldBlastDb = new BlastDB(fastaFile);
        assertThat(checkFile.lastModified(), greaterThanOrEqualTo(fastaFile.lastModified()));
    }

    public void testBlastParms() {
        BlastParms parms = new BlastParms().set("-a").set("-b", 100).db_gen_code(11).maxE(1e-20)
                .maxPerQuery(5).minPercent(50).num_threads(6).query_gen_code(4).pctLenOfQuery(0.5);
        assertThat(parms.getPctLenOfQuery(), equalTo(0.5));
        assertThat(parms.get(), contains("-a", "-b", "100", "-db_gen_code", "11", "-evalue", "1.0E-20",
                "-max_target_seqs", "5","-num_threads", "6", "-perc_identity", "50", "-query_genetic_code", "4"));

    }

}
