/**
 *
 */
package org.theseed.proteins.kmer.hash;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.util.HashSet;
import java.util.Set;

import org.junit.jupiter.api.Test;

/**
 * @author Bruce Parrello
 *
 */
class ProteinEncoderTest {

    private static final String[] KMERS = new String[] {
            "VTPAAPPAPVKATAPVP",
            "MQSQSRIK",
            "DLKRKPGQ",
            "GANLYNLS",
            "LGGLIIIG",
            "GLVWSLPG",
            "TTGNLEIZ",
            "QSAPAPXX",
            "QXXPAPAS",
            "LCWSLPGS",
            "CWSLPGSL",
            "WSLPGSLC",
            "VSLPGSLC",
            "XSLPGSLC",
            "WSLPGTLC",
            "VQSTPDSKPAAPRLPAP",
            "VTPAAPPAPVKATAPV",
            "VTPAAPPAPVKATAP",
            "VTPAAPPAPVKATA",
            "WTPAAPPAPVKAT",
            "VTPAAPPAPVKAT",
            "VTPAAPPAPVKA",
            "VTPAAPPAPVK",
            "VTPAAPPAPV",
            "VTPAAPPAP",
            "VTPAAPPA",
            "VTPAAPP",
            "VTPAAP",
            "VTPAA",
            "VTPA",
            "VTP",
            "VT",
            "V",
            ""
    };

    @Test
    void testEncoding() {
        Set<ProteinEncoding> checker = new HashSet<ProteinEncoding>(KMERS.length * 4 / 3 + 1);
        for (String kmer : KMERS) {
            ProteinEncoding coding = new ProteinEncoding(kmer);
            assertThat(kmer, coding, not(in(checker)));
            checker.add(coding);
        }
    }

}
