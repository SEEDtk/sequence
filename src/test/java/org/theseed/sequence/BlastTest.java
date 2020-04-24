/**
 *
 */
package org.theseed.sequence;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.theseed.genome.Genome;
import org.theseed.sequence.blast.BlastDB;
import org.theseed.sequence.blast.BlastHit;
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

    private final File tempDir = new File("src/test", "temp");

    @Override
    protected void setUp() throws Exception {
        // Set up our temp directory for blast files.
        if (! tempDir.isDirectory())
            FileUtils.forceMkdir(tempDir);
        FileUtils.cleanDirectory(tempDir);
    }

    public void testFasta() throws IOException, InterruptedException {
        File gtoFile = new File("src/test", "1313.7001.gto");
        Genome gto = new Genome(gtoFile);
        FileUtils.cleanDirectory(tempDir);
        // Create and test a DNA database.
        File fastaFile = new File(tempDir, "temp.fna");
        gto.saveDna(fastaFile);
        BlastDB newBlastDb = new BlastDB(fastaFile, 11);
        String[] suffixes = new String[] { ".nhr", ".nin", ".nsq" };
        for (String suffix : suffixes) {
            File testFile = new File(tempDir, "temp.fna" + suffix);
            assertThat(testFile, aReadableFile());
        }
        assertThat(newBlastDb.getType(), equalTo(BlastDB.Type.DNA));
        assertThat(newBlastDb.getGeneticCode(), equalTo(11));
        BlastDB oldBlastDb = new BlastDB(fastaFile);
        assertThat(oldBlastDb.getType(), equalTo(BlastDB.Type.DNA));
        assertThat(oldBlastDb.getGeneticCode(), equalTo(11));
        // Create and test a protein database.
        File protFile = new File(tempDir, "temp.faa");
        gto.savePegs(protFile);
        BlastDB protBlastDb = new BlastDB(protFile, BlastDB.PROTEIN);
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
        BlastParms parms = new BlastParms().set("-a").set("-b", 100).db_gencode(11).maxE(1e-20)
                .maxPerQuery(5).minPercent(50).num_threads(6).query_gen_code(4).pctLenOfQuery(0.5);
        assertThat(parms.getPctLenOfQuery(), equalTo(0.5));
        assertThat(parms.get(), contains("-a", "-b", "100", "-db_gencode", "11", "-evalue", "1.0E-20",
                "-max_target_seqs", "5","-num_threads", "6", "-perc_identity", "50", "-query_genetic_code", "4"));

    }

    public void testActualBlast() throws IOException, InterruptedException, CloneNotSupportedException {
        Genome g2 = new Genome(new File("src/test", "1685.390.gto"));
        File g1Pegs = new File("src/test", "g1.faa");
        File g1dna = new File("src/test", "g1.fna");
        File g3dna = new File("src/test", "g3.fna");
        File g2Pegs = new File(tempDir, "g2.faa");
        File g2Contigs = new File(tempDir, "g2.fna");
        g2.saveDna(g2Contigs);
        g2.savePegs(g2Pegs);
        BlastDB g2ContigBlast = new BlastDB(g2Contigs, 11);
        BlastParms parms = new BlastParms().maxE(1e-10).maxPerQuery(5).pctLenOfQuery(50);
        List<BlastHit> results = g2ContigBlast.blastProteins(FastaInputStream.readAll(g1Pegs), parms);
        assertThat(results.size(), equalTo(3));
        BlastHit hit = results.get(0);
        assertThat(hit.getQueryId(), equalTo("fig|47466.134.peg.336"));
        assertThat(hit.getQueryDef(), equalTo("Thioredoxin reductase (EC 1.8.1.9)"));
        assertThat(hit.getQueryLoc().toString(), equalTo("fig|47466.134.peg.336+[18, 324]"));
        assertThat(hit.getQueryLen(), equalTo(325));
        assertThat(hit.getSubjectId(), equalTo("1685.390.con.0085"));
        assertThat(hit.getSubjectDef(), equalTo("SSOQ01000007.1"));
        assertThat(hit.getSubjectLoc().toString(), equalTo("1685.390.con.0085+[82090, 83118]"));
        assertThat(hit.getSubjectLen(), equalTo(114301));
        assertThat(hit.getBitScore(), closeTo(192, 0.01));
        assertThat(hit.getNumIdentical(), equalTo(117));
        assertThat(hit.getNumGap(), equalTo(38));
        assertThat(hit.getPositives(), equalTo(176));
        assertThat(hit.getEvalue(), closeTo(9.09e-56, 1e-58));
        assertThat(hit.getAlignLen(), equalTo(344));
        assertThat(hit.getQuerySeq().length(), equalTo(344));
        assertThat(hit.getQuerySeq().substring(0, 10), equalTo("RTELNSVK-D"));
        assertThat(hit.getSubjectSeq().length(), equalTo(344));
        assertThat(hit.getSubjectSeq().substring(0, 10), equalTo("RTEGNHMEHN"));
        assertThat(hit.getPercentIdentity(), closeTo(34.0, 0.5));
        assertThat(hit.getPercentSimilarity(), closeTo(85.2, 0.1));
        assertThat(hit.getSubjectPercentMatch(), closeTo(0.76, 0.01));
        assertThat(hit.getQueryPercentMatch(), closeTo(90.1, 0.1));
        for (BlastHit result : results) {
            assertThat(result.getQuerySeq().length(), equalTo(result.getAlignLen()));
            assertThat(result.getSubjectSeq().length(), equalTo(result.getAlignLen()));
            assertThat(result.getNumIdentical(), lessThanOrEqualTo(result.getAlignLen()));
            assertThat(result.getPositives(), lessThanOrEqualTo(result.getAlignLen()));
            assertThat(result.getPositives() + result.getNumIdentical() + result.getNumGap(),
                    lessThanOrEqualTo(result.getAlignLen()));
            assertThat(result.getEvalue(), lessThanOrEqualTo(1e-10));
            assertThat(result.getQueryPercentMatch(), greaterThanOrEqualTo(50.0));
        }
        BlastDB g2PegBlast = new BlastDB(g2Pegs, BlastDB.PROTEIN);
        results = g2PegBlast.blastProteins(FastaInputStream.readAll(g1Pegs), parms);
        assertThat(results.size(), equalTo(3));
        hit = results.get(0);
        assertThat(hit.getQueryId(), equalTo("fig|47466.134.peg.336"));
        assertThat(hit.getQueryDef(), equalTo("Thioredoxin reductase (EC 1.8.1.9)"));
        assertThat(hit.getQueryLoc().toString(), equalTo("fig|47466.134.peg.336+[25, 324]"));
        assertThat(hit.getQueryLen(), equalTo(325));
        assertThat(hit.getSubjectId(), equalTo("fig|1685.390.peg.1900"));
        assertThat(hit.getSubjectDef(), equalTo("Thioredoxin reductase (EC 1.8.1.9)"));
        assertThat(hit.getSubjectLoc().toString(), equalTo("fig|1685.390.peg.1900+[3, 337]"));
        assertThat(hit.getSubjectLen(), equalTo(339));
        assertThat(hit.getBitScore(), closeTo(188, 0.01));
        assertThat(hit.getNumIdentical(), equalTo(113));
        assertThat(hit.getNumGap(), equalTo(37));
        assertThat(hit.getPositives(), equalTo(170));
        assertThat(hit.getEvalue(), closeTo(9.99e-59, 1e-61));
        assertThat(hit.getAlignLen(), equalTo(336));
        assertThat(hit.getQuerySeq().length(), equalTo(336));
        assertThat(hit.getQuerySeq().substring(0, 10), equalTo("KDVIIVGSGP"));
        assertThat(hit.getSubjectSeq().length(), equalTo(336));
        assertThat(hit.getSubjectSeq().substring(0, 10), equalTo("HNVIIIGSGP"));
        assertThat(hit.getPercentIdentity(), closeTo(33.6, 0.1));
        assertThat(hit.getPercentSimilarity(), closeTo(84.2, 0.1));
        assertThat(hit.getSubjectPercentMatch(), closeTo(83.5, 0.1));
        assertThat(hit.getQueryPercentMatch(), closeTo(87.1, 0.1));
        results = g2ContigBlast.blastDna(FastaInputStream.readAll(g3dna), parms);
        assertThat(results.size(), equalTo(4));
        hit = results.get(0);
        assertThat(hit.getQueryId(), equalTo("fig|1685.78.peg.1985|NRBB09_1929|"));
        assertThat(hit.getQueryDef(), equalTo("Thioredoxin reductase (EC 1.8.1.9)   [Bifidobacterium breve strain NRBB09 | 1685.78]"));
        assertThat(hit.getQueryLoc().toString(), equalTo("fig|1685.78.peg.1985|NRBB09_1929|+[1, 1020]"));
        assertThat(hit.getQueryLen(), equalTo(1020));
        assertThat(hit.getSubjectId(), equalTo("1685.390.con.0085"));
        assertThat(hit.getSubjectDef(), equalTo("SSOQ01000007.1"));
        assertThat(hit.getSubjectLoc().toString(), equalTo("1685.390.con.0085+[82108, 83127]"));
        assertThat(hit.getSubjectLen(), equalTo(114301));
        assertThat(hit.getBitScore(), closeTo(1845, 0.01));
        assertThat(hit.getNumIdentical(), equalTo(1013));
        assertThat(hit.getNumGap(), equalTo(0));
        assertThat(hit.getPositives(), equalTo(0));
        assertThat(hit.getEvalue(), closeTo(0, 1e-60));
        assertThat(hit.getAlignLen(), equalTo(1020));
        assertThat(hit.getQuerySeq().length(), equalTo(1020));
        assertThat(hit.getQuerySeq().substring(0, 10), equalTo("ATGGAACATA"));
        assertThat(hit.getSubjectSeq().length(), equalTo(1020));
        assertThat(hit.getSubjectSeq().substring(0, 10), equalTo("ATGGAACATA"));
        assertThat(hit.getPercentIdentity(), closeTo(99.3, 0.1));
        assertThat(hit.getPercentSimilarity(), closeTo(99.3, 0.1));
        assertThat(hit.getSubjectPercentMatch(), closeTo(0.89, 0.01));
        assertThat(hit.getQueryPercentMatch(), closeTo(99.3, 0.1));
        results = g2PegBlast.blastDna(FastaInputStream.readAll(g1dna), parms);
        assertThat(results.size(), equalTo(3));
        hit = results.get(0);
        assertThat(hit.getQueryId(), equalTo("fig|47466.134.peg.336|EZU69_01700|"));
        assertThat(hit.getQueryDef(), equalTo("Thioredoxin reductase (EC 1.8.1.9)   [Borrelia miyamotoi strain Yekat-21 | 47466.134]"));
        assertThat(hit.getQueryLoc().toString(), equalTo("fig|47466.134.peg.336|EZU69_01700|+[73, 972]"));
        assertThat(hit.getQueryLen(), equalTo(978));
        assertThat(hit.getSubjectId(), equalTo("fig|1685.390.peg.1900"));
        assertThat(hit.getSubjectDef(), equalTo("Thioredoxin reductase (EC 1.8.1.9)"));
        assertThat(hit.getSubjectLoc().toString(), equalTo("fig|1685.390.peg.1900+[3, 337]"));
        assertThat(hit.getSubjectLen(), equalTo(339));
        assertThat(hit.getBitScore(), closeTo(188, 0.01));
        assertThat(hit.getNumIdentical(), equalTo(113));
        assertThat(hit.getNumGap(), equalTo(37));
        assertThat(hit.getPositives(), equalTo(170));
        assertThat(hit.getEvalue(), closeTo(1.10e-58, 1e-60));
        assertThat(hit.getAlignLen(), equalTo(336));
        assertThat(hit.getQuerySeq().length(), equalTo(336));
        assertThat(hit.getQuerySeq().substring(0, 10), equalTo("KDVIIVGSGP"));
        assertThat(hit.getSubjectSeq().length(), equalTo(336));
        assertThat(hit.getSubjectSeq().substring(0, 10), equalTo("HNVIIIGSGP"));
        assertThat(hit.getPercentIdentity(), closeTo(33.6, 0.1));
        assertThat(hit.getPercentSimilarity(), closeTo(84.2, 0.1));
        assertThat(hit.getSubjectPercentMatch(), closeTo(83.5, 0.1));
        assertThat(hit.getQueryPercentMatch(), closeTo(86.8, 0.1));
    }

}
