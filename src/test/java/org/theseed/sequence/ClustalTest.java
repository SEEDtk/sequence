package org.theseed.sequence;

import org.junit.jupiter.api.Test;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

import org.theseed.sequence.clustal.ClustalPipeline;

public class ClustalTest  {

    @Test
    public void testClustal() throws IOException, InterruptedException {
        ClustalPipeline tester = new ClustalPipeline(new File("data", "seedProt.fa"));
        List<Sequence> output = tester.run();
        FastaInputStream testStream = new FastaInputStream(new File("data", "alignment.fa"));
        int i = 0;
        PrintWriter savedFile = new PrintWriter(new File("data", "align.html"));
        savedFile.println("<div id=\"Aligned\"><table>");
        for (Sequence testSeq : testStream) {
            Sequence seq = output.get(i);
            assertThat(seq, equalTo(testSeq));
            savedFile.format("<tr><th>%s</th><td>%s</td></tr>%n", seq.getLabel(), seq.getSequence());
            i++;
        }
        savedFile.println("</div></table>");
        savedFile.close();
        assertThat(i, equalTo(output.size()));
        testStream.close();
    }
}
