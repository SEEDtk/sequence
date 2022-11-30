/**
 *
 */
package org.theseed.sequence;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FileUtils;
import org.theseed.genome.Genome;
import org.theseed.io.LineReader;
import org.theseed.sequence.blast.BlastDB;
import org.theseed.sequence.blast.DnaBlastDB;
import org.theseed.sequence.blast.ProteinBlastDB;

import org.theseed.sequence.blast.BlastHit;
import org.theseed.sequence.blast.BlastParms;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Test;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

/**
 *
 * Tests for the blast support
 *
 * @author Bruce Parrello
 *
 */
public class BlastTest  {

    private static final File tempDir = new File("data", "temp");

    public BlastTest() throws Exception {
        // Set up our temp directory for blast files.
        if (! tempDir.isDirectory())
            FileUtils.forceMkdir(tempDir);
    }

    @AfterAll
    protected static void tearDown() throws Exception {
        // Clean the temp directory for blast files.
        FileUtils.cleanDirectory(tempDir);
    }

    @Test
    public void testFasta() throws IOException, InterruptedException {
        File gtoFile = new File("data", "1313.7001.gto");
        Genome gto = new Genome(gtoFile);
        // Create and test a DNA database.
        File fastaFile = new File(tempDir, "temp.fna");
        BlastDB newBlastDb = DnaBlastDB.create(fastaFile, gto);
        assertThat(((DnaBlastDB) newBlastDb).getGeneticCode(), equalTo(11));
        String[] suffixes = new String[] { ".nhr", ".nin", ".nsq" };
        for (String suffix : suffixes) {
            File testFile = new File(tempDir, "temp.fna" + suffix);
            assertThat(testFile.canRead(), equalTo(true));
        }
        BlastDB oldBlastDb = BlastDB.load(fastaFile);
        assertThat(oldBlastDb, instanceOf(DnaBlastDB.class));
        assertThat(((DnaBlastDB) oldBlastDb).getGeneticCode(), equalTo(11));
        // Create and test a protein database.
        File protFile = new File(tempDir, "temp.faa");
        newBlastDb = ProteinBlastDB.create(protFile, gto);
        suffixes = new String[] { ".phr", ".pin", ".psq" };
        for (String suffix : suffixes) {
            File testFile = new File(tempDir, "temp.faa" + suffix);
            assertThat(testFile.canRead(), equalTo(true));
        }
        oldBlastDb = BlastDB.load(protFile);
        assertThat(oldBlastDb, instanceOf(ProteinBlastDB.class));
        // Verify that updating the fasta file causes a regen.
        File checkFile = new File(tempDir, "temp.fna.nsq");
        Thread.sleep(1000);
        gto.saveDna(fastaFile);
        assertThat(checkFile.lastModified(), lessThan(fastaFile.lastModified()));
        oldBlastDb = BlastDB.load(fastaFile);
        assertThat(checkFile.lastModified(), greaterThanOrEqualTo(fastaFile.lastModified()));
    }

    @Test
    public void testBlastParms() {
        BlastParms parms = new BlastParms().set("-a").set("-b", 100).db_gencode(11).maxE(1e-20)
                .maxPerQuery(5).minPercent(50).num_threads(6).query_gencode(4).pctLenOfQuery(0.5)
                .minQueryBitScore(0.75).minQueryIdentity(0.35).minLen(100);
        assertThat(parms.getPctLenOfQuery(), equalTo(0.5));
        assertThat(parms.getPctIdentity(), equalTo(50.0));
        assertThat(parms.getMinQueryBitScore(), equalTo(0.75));
        assertThat(parms.getMinQueryIdentity(), equalTo(0.35));
        assertThat(parms.getMinLen(), equalTo(100));
        assertThat(parms.get(), contains("-a", "-b", "100", "-db_gencode", "11", "-evalue", "1.0E-20",
                "-max_target_seqs", "5","-num_threads", "6", "-query_gencode", "4"));
        // Test the blastn stuff.
        parms.dustOff().penalty(-1);
        assertThat(parms.get(), contains("-a", "-b", "100", "-db_gencode", "11", "-evalue", "1.0E-20",
                "-max_target_seqs", "5","-num_threads", "6", "-query_gencode", "4"));
        parms.customizeForBlastn();
        assertThat(parms.get(), contains("-a", "-b", "100", "-db_gencode", "11", "-dust", "no", "-evalue", "1.0E-20",
                "-max_target_seqs", "5","-num_threads", "6", "-penalty", "-1", "-query_gencode", "4"));
    }

    @Test
    public void testActualBlast() throws IOException, InterruptedException, CloneNotSupportedException {
        Genome g2 = new Genome(new File("data", "1685.390.gto"));
        File g1Pegs = new File("data", "g1.faa");
        File g1dna = new File("data", "g1.fna");
        File g3dna = new File("data", "g3.fna");
        File g2Pegs = new File(tempDir, "g2.faa");
        File g2Contigs = new File(tempDir, "g2.fna");
        BlastDB g2ContigBlast = DnaBlastDB.create(g2Contigs, g2);
        BlastParms parms = new BlastParms().maxE(1e-10).maxPerQuery(5).pctLenOfQuery(50);
        List<BlastHit> sortTest = new ArrayList<BlastHit>(50);
        try (ProteinInputStream pStream = new ProteinInputStream(g1Pegs)) {
            List<BlastHit> results = pStream.blast(g2ContigBlast, parms);
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
            assertThat(hit.getPositives(), equalTo(59));
            assertThat(hit.getNumSimilar(), equalTo(176));
            assertThat(hit.getEvalue(), closeTo(9.09e-56, 1e-58));
            assertThat(hit.getAlignLen(), equalTo(344));
            assertThat(hit.getQuerySeq().length(), equalTo(344));
            assertThat(hit.getQuerySeq().substring(0, 10), equalTo("RTELNSVK-D"));
            assertThat(hit.getSubjectSeq().length(), equalTo(344));
            assertThat(hit.getSubjectSeq().substring(0, 10), equalTo("RTEGNHMEHN"));
            assertThat(hit.getPercentIdentity(), closeTo(34.0, 0.5));
            assertThat(hit.getPercentSimilarity(), closeTo(51.2, 0.1));
            assertThat(hit.getSubjectPercentMatch(), closeTo(0.46, 0.01));
            assertThat(hit.getQueryPercentMatch(), closeTo(54.2, 0.1));
            assertThat(hit.getQueryIdentity(), closeTo(0.36, 0.005));
            assertThat(hit.getQueryBitScore(), closeTo(0.59, 0.001));
            for (BlastHit result : results) {
                assertThat(result.getQuerySeq().length(), equalTo(result.getAlignLen()));
                assertThat(result.getSubjectSeq().length(), equalTo(result.getAlignLen()));
                assertThat(result.getNumIdentical(), lessThanOrEqualTo(result.getAlignLen()));
                assertThat(result.getPositives(), lessThanOrEqualTo(result.getAlignLen()));
                assertThat(result.getPositives(), greaterThanOrEqualTo(0));
                assertThat(result.getPositives() + result.getNumIdentical() + result.getNumGap(),
                        lessThanOrEqualTo(result.getAlignLen()));
                assertThat(result.getEvalue(), lessThanOrEqualTo(1e-10));
                assertThat(result.getQueryPercentMatch(), greaterThanOrEqualTo(50.0));
            }
            sortTest.addAll(results);
            assertThat(g2ContigBlast.getBlastType(), equalTo("tblastn"));
        }
        BlastDB g2PegBlast = ProteinBlastDB.create(g2Pegs, g2);
        try (ProteinInputStream pStream = new ProteinInputStream(g1Pegs)) {
            List<BlastHit> results = pStream.blast(g2PegBlast, parms);
            if (BlastDB.haveDiamond()) {
                assertThat(results.size(), equalTo(2));
                BlastHit hit = results.get(0);
                assertThat(hit.getQueryId(), equalTo("fig|47466.134.peg.336"));
                assertThat(hit.getQueryDef(), equalTo("Thioredoxin reductase (EC 1.8.1.9)"));
                assertThat(hit.getQueryLoc().toString(), equalTo("fig|47466.134.peg.336+[19, 323]"));
                assertThat(hit.getQueryLen(), equalTo(325));
                assertThat(hit.getSubjectId(), equalTo("fig|1685.390.peg.443"));
                assertThat(hit.getSubjectDef(), equalTo("Alkyl hydroperoxide reductase protein F"));
                assertThat(hit.getSubjectLoc().toString(), equalTo("fig|1685.390.peg.443+[2, 310]"));
                assertThat(hit.getSubjectLen(), equalTo(637));
                assertThat(hit.getBitScore(), closeTo(187, 0.01));
                assertThat(hit.getNumIdentical(), equalTo(108));
                assertThat(hit.getNumGap(), equalTo(8));
                assertThat(hit.getPositives(), equalTo(66));
                assertThat(hit.getNumSimilar(), equalTo(174));
                assertThat(hit.getEvalue(), closeTo(8.99e-56, 1e-58));
                assertThat(hit.getAlignLen(), equalTo(311));
                assertThat(hit.getQuerySeq().length(), equalTo(305));
                assertThat(hit.getQuerySeq().substring(0, 10), equalTo("TELNSVKDVI"));
                assertThat(hit.getSubjectSeq().length(), equalTo(309));
                assertThat(hit.getSubjectSeq().substring(0, 10), equalTo("TQTNNLYDVV"));
                assertThat(hit.getPercentIdentity(), closeTo(34.7, 0.1));
                assertThat(hit.getPercentSimilarity(), closeTo(55.9, 0.1));
                assertThat(hit.getSubjectPercentMatch(), closeTo(27.3, 0.1));
                assertThat(hit.getQueryPercentMatch(), closeTo(53.5, 0.1));
                sortTest.addAll(results);
                assertThat(g2PegBlast.getBlastType(), equalTo("blastp"));
            } else {
                assertThat(results.size(), equalTo(3));
                BlastHit hit = results.get(0);
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
                assertThat(hit.getPositives(), equalTo(57));
                assertThat(hit.getNumSimilar(), equalTo(170));
                assertThat(hit.getEvalue(), closeTo(9.99e-59, 1e-61));
                assertThat(hit.getAlignLen(), equalTo(336));
                assertThat(hit.getQuerySeq().length(), equalTo(336));
                assertThat(hit.getQuerySeq().substring(0, 10), equalTo("KDVIIVGSGP"));
                assertThat(hit.getSubjectSeq().length(), equalTo(336));
                assertThat(hit.getSubjectSeq().substring(0, 10), equalTo("HNVIIIGSGP"));
                assertThat(hit.getPercentIdentity(), closeTo(33.6, 0.1));
                assertThat(hit.getPercentSimilarity(), closeTo(50.6, 0.1));
                assertThat(hit.getSubjectPercentMatch(), closeTo(50.1, 0.1));
                assertThat(hit.getQueryPercentMatch(), closeTo(52.3, 0.1));
                sortTest.addAll(results);
                assertThat(g2PegBlast.getBlastType(), equalTo("blastp"));
            }
        }
        try (DnaInputStream dStream = new DnaInputStream(g3dna, 11)) {
            List<BlastHit> results = dStream.blast(g2ContigBlast, parms);
            assertThat(results.size(), equalTo(4));
            BlastHit hit = results.get(0);
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
            assertThat(hit.getNumSimilar(), equalTo(1013));
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
            sortTest.addAll(results);
            assertThat(g2ContigBlast.getBlastType(), equalTo("blastn"));
        }
        try (DnaInputStream dStream = new DnaInputStream(g1dna, 11)) {
            List<BlastHit> results = dStream.blast(g2PegBlast, parms);
            if (BlastDB.haveDiamond()) {
                assertThat(results.size(), equalTo(2));
                BlastHit hit = results.get(0);
                assertThat(hit.getQueryId(), equalTo("fig|47466.134.peg.336|EZU69_01700|"));
                assertThat(hit.getQueryDef(), equalTo("Thioredoxin reductase (EC 1.8.1.9)   [Borrelia miyamotoi strain Yekat-21 | 47466.134]"));
                assertThat(hit.getQueryLoc().toString(), equalTo("fig|47466.134.peg.336|EZU69_01700|+[55, 969]"));
                assertThat(hit.getQueryLen(), equalTo(978));
                assertThat(hit.getSubjectId(), equalTo("fig|1685.390.peg.443"));
                assertThat(hit.getSubjectDef(), equalTo("Alkyl hydroperoxide reductase protein F"));
                assertThat(hit.getSubjectLoc().toString(), equalTo("fig|1685.390.peg.443+[2, 310]"));
                assertThat(hit.getSubjectLen(), equalTo(637));
                assertThat(hit.getBitScore(), closeTo(187, 0.01));
                assertThat(hit.getNumIdentical(), equalTo(108));
                assertThat(hit.getNumGap(), equalTo(8));
                assertThat(hit.getPositives(), equalTo(66));
                assertThat(hit.getNumSimilar(), equalTo(174));
                assertThat(hit.getEvalue(), closeTo(9.27e-56, 1e-58));
                assertThat(hit.getAlignLen(), equalTo(311));
                assertThat(hit.getQuerySeq().length(), equalTo(915));
                assertThat(hit.getQuerySeq().substring(0, 10), equalTo("ACGGAGCTAA"));
                assertThat(hit.getSubjectSeq().length(), equalTo(309));
                assertThat(hit.getSubjectSeq().substring(0, 10), equalTo("TQTNNLYDVV"));
                assertThat(hit.getPercentIdentity(), closeTo(34.7, 0.1));
                assertThat(hit.getPercentSimilarity(), closeTo(55.9, 0.1));
                assertThat(hit.getSubjectPercentMatch(), closeTo(27.3, 0.1));
                assertThat(hit.getQueryPercentMatch(), closeTo(53.4, 0.1));
                sortTest.addAll(results);
                assertThat(g2PegBlast.getBlastType(), equalTo("blastx"));
            } else {
                assertThat(results.size(), equalTo(3));
                BlastHit hit = results.get(0);
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
                assertThat(hit.getPositives(), equalTo(57));
                assertThat(hit.getNumSimilar(), equalTo(170));
                assertThat(hit.getEvalue(), closeTo(1.10e-58, 1e-60));
                assertThat(hit.getAlignLen(), equalTo(336));
                assertThat(hit.getQuerySeq().length(), equalTo(336));
                assertThat(hit.getQuerySeq().substring(0, 10), equalTo("KDVIIVGSGP"));
                assertThat(hit.getSubjectSeq().length(), equalTo(336));
                assertThat(hit.getSubjectSeq().substring(0, 10), equalTo("HNVIIIGSGP"));
                assertThat(hit.getPercentIdentity(), closeTo(33.6, 0.1));
                assertThat(hit.getPercentSimilarity(), closeTo(50.6, 0.1));
                assertThat(hit.getSubjectPercentMatch(), closeTo(50.1, 0.1));
                assertThat(hit.getQueryPercentMatch(), closeTo(52.1, 0.1));
                sortTest.addAll(results);
                assertThat(g2PegBlast.getBlastType(), equalTo("blastx"));
            }
        }
        // Test sorting by query location.
        Comparator<BlastHit> compare = new BlastHit.ByLoc(BlastHit.QUERY);
        sortTest.sort(compare);
        for (int i = 1; i < sortTest.size(); i++) {
            BlastHit prev = sortTest.get(i-1);
            BlastHit curr = sortTest.get(i);
            assertThat(compare.compare(prev, curr), not(equalTo(0)));
            assertThat(prev.getQueryLoc(), lessThanOrEqualTo(curr.getQueryLoc()));
        }
    }

    /**
     * test the result sort
     * @throws IOException
     */
    @Test
    public void testResultSort() throws IOException {
        List<BlastHit> results0 = new ArrayList<BlastHit>(15);
        Map<String, String> qMap = new HashMap<String, String>();
        qMap.put("q1", "qtitle 1");
        qMap.put("q2", "qtitle 2");
        try (LineReader testStream = new LineReader(new File("data", "results.txt"))) {
            for (String line : testStream) {
                BlastHit result = new BlastHit(line, qMap, true, true);
                results0.add(result);
            }
        }
        Map<String, List<BlastHit>> sortMap = BlastHit.sort(results0);
        List<BlastHit> results = sortMap.get("q1");
        assertThat(results.size(), equalTo(4));
        assertThat(results.get(0).getSubjectId(), equalTo("s2"));
        assertThat(results.get(1).getSubjectId(), equalTo("s3"));
        assertThat(results.get(2).getSubjectId(), equalTo("s1"));
        assertThat(results.get(2).getEvalue(), closeTo(4e-13, 1e-15));
        assertThat(results.get(3).getSubjectId(), equalTo("s4"));
        results = sortMap.get("q2");
        assertThat(results.size(), equalTo(2));
        assertThat(results.get(0).getSubjectId(), equalTo("s1"));
        assertThat(results.get(1).getSubjectId(), equalTo("s5"));
        for (Map.Entry<String, List<BlastHit>> entry : sortMap.entrySet()) {
            List<BlastHit> list = entry.getValue();
            for (BlastHit result : list) {
                assertThat(result.getQueryId(), equalTo(entry.getKey()));
                assertThat(results0.contains(result), equalTo(true));
            }

        }
    }

}
