/**
 *
 */
package org.theseed.sequence;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import org.junit.Test;
import org.theseed.genome.Contig;

/**
 * Test tetramer profiles
 *
 * @author Bruce Parrello
 *
 */
public class TetraTest {

    @Test
    public void test() {
        // Test all the tetramers to make sure they revcmp properly.
        StringBuffer dnaBuffer = new StringBuffer(256 * 4);
        String[] nucleotides = new String[] { "a", "c", "g", "t" };
        int oldT = -1;
        for (int i1 = 0; i1 < 4; i1++)
            for (int i2 = 0; i2 < 4; i2++)
                for (int i3 = 0; i3 < 4; i3++)
                    for (int i4 = 0; i4 < 4; i4++) {
                        String seq = nucleotides[i1] + nucleotides[i2] + nucleotides[i3] + nucleotides[i4];
                        String rSeq = Contig.reverse(seq);
                        int t1 = TetramerProfile.tetraNum(seq, 1);
                        int rT1 = TetramerProfile.tetraNum(rSeq, 1);
                        assertThat(seq + "/" + rSeq, rT1, equalTo(TetramerProfile.revTetraNum(t1)));
                        assertThat(t1, not(equalTo(oldT)));
                        oldT = t1;
                        dnaBuffer.append(seq);
                    }
        TetramerProfile testProfile = new TetramerProfile(dnaBuffer.toString());
        assertThat(testProfile, not(nullValue()));
    }

    @Test
    public void profiles() {
        String allA = "aaaaaaaaaaaaaaaa";
        TetramerProfile profile = new TetramerProfile(allA);
        assertThat(profile.toString(), startsWith("1.0\t0.0"));
        String test = "taaccgttctccaaaaattaaacgagttttagctgttgaaattcctgaaacaaggtcatttatgcatttggacacggttttcaccatggttaattttgcg";
        profile = new TetramerProfile(test, 1, 100);
        assertThat(profile, not(nullValue()));
    }

}
