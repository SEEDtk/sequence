/**
 *
 */
package org.theseed.sequence;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.util.SortedSet;

/**
 * Unit test for minhash stuff.
 *
 * @author Bruce Parrello
 *
 */
public class SeqTest extends TestCase {

    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public SeqTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( SeqTest.class );
    }

    /** jumble arrays for mutations */
    private static final int[] JUMBLE1 = new int[] { 67, 61, 96, 65, 56, 95, 2, 75, 17, 33, 87, 41, 16, 12, 15, 94, 6, 68, 0, 49, 38, 7, 34, 51, 92, 5, 24, 13, 45, 44, 93, 14, 11, 27, 74, 29, 88, 76, 19, 90, 8, 66, 80, 82, 50, 4, 25, 46, 70, 91};
    private static final int[] JUMBLE2 = new int[] { 85, 60, 9, 58, 52, 72, 54, 39, 20, 81, 30, 86, 10, 69, 77, 26, 47, 35, 73, 89, 98, 53, 28, 36, 62, 84, 32, 42, 57, 97, 21, 63, 83, 59, 37, 79, 40, 48, 78, 31, 55, 3, 22, 43, 18, 71, 64, 1, 23, 99};
    private static final int[] JUMBLE3 = new int[] { 3, 73, 47, 25, 57, 14, 31, 45, 26, 91, 13, 65, 7, 68, 61, 42, 98, 78, 37, 39, 79, 70, 54, 74, 40, 80, 6, 51, 52, 94, 86, 66, 33, 85, 81, 56, 0, 49, 44, 88, 12, 72, 16, 77, 35, 64, 59, 55, 9, 92};
    private static final int[] JUMBLE4 = new int[] { 1, 67, 15, 58, 11, 97, 46, 43, 63, 20, 75, 89, 10, 96, 60, 41, 53, 87, 71, 90, 93, 48, 29, 22, 76, 82, 84, 32, 38, 4, 24, 19, 69, 50, 5, 62, 17, 18, 21, 8, 34, 27, 95, 83, 2, 28, 30, 23, 36, 99};
    private static final int[] JUMBLE5 = new int[] { 63, 2, 36, 47, 52, 61, 37, 80, 5, 98, 46, 39, 49, 57, 88, 17, 86, 32, 74, 25, 87, 48, 34, 95, 97, 27, 77, 13, 78, 66, 28, 53, 9, 3, 70, 60, 21, 11, 7, 24, 59, 31, 69, 73, 94, 76, 91, 43, 20, 89};
    private static final int[][] JUMBLES = new int[][] { JUMBLE1, JUMBLE2, JUMBLE3, JUMBLE4, JUMBLE5 };

    public void testSeqHash() {
        ProteinKmers.setKmerSize(8);
        // We will test using a set of protein kmers.  There will be five very different starting proteins
        // of length 100. For each, we will generate a set of mutations of 10, 20, 30, 40, and 50 amino acids.
        // These are put into the hash, and then we ask for the closest N proteins of each starting protein.
        // This is a crude test, but is a good place to start.
        String p1 = "MDIQITHQVTEFDKEELLAGLRSYNAQFVDFSKNGQLGVYCRNESGEMVGGLIADRKGPWLCIDYLWVSESARNCGLGSKLMAMAEKEGLRKGCAHGLVD";
        String p2 = "MSKDEISYQILYRYSLEKLYSTLTRRVDNVLSFALIFLGVGVTINVGSPFILGPGIVGIAILKRVLRFGTRSAQADRQSRAWLKLFNTQHRFPSDKTLFL";
        String p3 = "MELQLMLNHFFERVRKDANFNAFLIDLEYNNIAYYIYFVATGNVKIITHAGHFISIKSNRKLIKVNSTPNTQLIKLTSDKHFSGEHSYEKYCTDLATAGV";
        String p4 = "MTNITLSTQHYRIHRSDVEPVKEKTTEKDIFAKSITAVRNSFISLSTSLSDRFSLHQQTDIPTTHFHRGSASEGRAVLTSKTVKDFMLQKLNSLDIKGNA";
        String p5 = "MSKDPAYARQTCEAILSAVYSNNKDHCCKLLISKGVSITPFLKEIGEAAQNAGLPGEIKNGVFTPGGAGANPFVVPLIASASIKYPHMFINHNQQVSFKA";
        // Test the distances of the original 5.
        ProteinKmers k1 = new ProteinKmers(p1);
        ProteinKmers k2 = new ProteinKmers(p2);
        ProteinKmers k3 = new ProteinKmers(p3);
        ProteinKmers k4 = new ProteinKmers(p4);
        ProteinKmers k5 = new ProteinKmers(p5);
        ProteinKmers[] kArray = new ProteinKmers[] { k1, k2, k3, k4, k5 };
        for (int i = 0; i < 4; i++)
            for (int j = i + 1; j < 5; j++)
                assertThat(kArray[i].distance(kArray[j]), greaterThan(0.99));
        String[] prots = new String[] { p1, p2, p3, p4, p5 };
        // Create the hash.
        LSHSeqHash hash = new LSHSeqHash(200, 15, 20, 142857);
        // Now create the mutations.
        for (int i = 0; i < 5; i++) {
            String p = prots[i];
            int[] jumble = JUMBLES[i];
            // This will be the starting point of the current mutation block.
            int mutations = 0;
            while (mutations <= 46) {
                int limit = mutations + 4;
                for (int k = mutations; k < limit; k++)
                    p = mutate(p, jumble[k]);
                mutations = limit;
                String label = String.format("p%d.%04d", i, mutations);
                ProteinKmers kmers = new ProteinKmers(p);
                double dist = kmers.distance(kArray[i]);
                double skDist = hash.testSketches(kmers, kArray[i]);
                assertThat(label, Math.abs(dist - skDist), lessThan(0.1));
                hash.add(kmers, label);
            }
        }
        // Get the closest 5 to each original protein.
        for (int i = 0; i < 5; i++) {
            LSHSeqHash.Result[] r1 = new LSHSeqHash.Result[8];
            String prefix = "p" + Integer.toString(i) + ".";
            SortedSet<LSHSeqHash.Result> rSet = hash.getClosest(kArray[i], 8, 0.99);
            assertThat(prefix, rSet.size(), lessThanOrEqualTo(8));
            assertThat(prefix, rSet.size(), greaterThanOrEqualTo(3));
            r1 = rSet.toArray(r1);
            for (int k = 0; k < r1.length - 1; k++) {
                LSHSeqHash.Result r = r1[k];
                if (r != null) {
                    LSHSeqHash.Result r0 = r1[k+1];
                    if (r0 != null)
                        assertThat(r.getDistance(), lessThanOrEqualTo(r0.getDistance()));
                    assertThat(r.getDistance(), lessThanOrEqualTo(0.99));
                    assertThat(r.getTarget().substring(0, 3), equalTo(prefix));
                }
            }
        }
    }


    private static final String AA_LIST = "GPAVLIMCFYWHKRQNEDSTG";

    /**
     * Mutate a single amino acid in a protein string.
     *
     * @param prot	protein string to mutate
     * @param pos	position of the amino acid to mutate
     *
     * @return the mutated protein
     */
    private String mutate(String prot, int pos) {
        char[] acids = prot.toCharArray();
        char aa = acids[pos];
        acids[pos] = AA_LIST.charAt(AA_LIST.indexOf(aa) + 1);
        String retVal = String.valueOf(acids);
        return retVal;
    }

}
