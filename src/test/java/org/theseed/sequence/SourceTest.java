/**
 *
 */
package org.theseed.sequence;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

import org.theseed.genome.Contig;
import org.theseed.genome.Genome;
import org.theseed.sequence.blast.BlastDB;
import org.theseed.sequence.blast.DnaBlastDB;
import org.theseed.sequence.blast.ProteinBlastDB;
import org.theseed.sequence.blast.Source;

import org.junit.jupiter.api.Test;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

/**
 * Test BLAST sources
 *
 * @author Bruce Parrello
 *
 */
public class SourceTest {

    @Test
    public void testBlastDB() throws IOException, InterruptedException {
        File tempDir = new File("data", "temp");
        File gFile = new File("data", "1313.7001.gto");
        File dnaFile = new File("data", "g1.fna");
        File pFile = new File("data", "g1.faa");
        BlastDB db = Source.contigs.subject(tempDir, gFile, 4, false);
        assertThat(db.isProtein(), equalTo(false));
        assertThat(db, instanceOf(DnaBlastDB.class));
        BlastDB db2 = Source.db.subject(tempDir, db.getFile(), 0, true);
        assertThat(db2.isProtein(), equalTo(false));
        assertThat(db2.getFile(), equalTo(db.getFile()));
        assertThat(db2, instanceOf(DnaBlastDB.class));
        assertThat(((DnaBlastDB) db2).getGeneticCode(), equalTo(11));
        db = Source.dna.subject(tempDir, dnaFile, 1, false);
        assertThat(db.isProtein(), equalTo(false));
        assertThat(((DnaBlastDB) db).getGeneticCode(), equalTo(1));
        assertThat(db.getFile(), equalTo(dnaFile));
        db = Source.features.subject(tempDir, gFile, 0, false);
        assertThat(db.isProtein(), equalTo(false));
        assertThat(db, instanceOf(DnaBlastDB.class));
        assertThat(((DnaBlastDB) db).getGeneticCode(), equalTo(11));
        db = Source.pegs.subject(tempDir, gFile, 0, false);
        assertThat(db.isProtein(), equalTo(true));
        assertThat(db, instanceOf(ProteinBlastDB.class));
        db = Source.prot.subject(tempDir, pFile, 0, false);
        assertThat(db.isProtein(), equalTo(true));
        assertThat(db.getFile(), equalTo(pFile));
        assertThat(db, instanceOf(ProteinBlastDB.class));
    }

    @Test
    public void testStream() throws IOException {
        File tempDir = new File("data", "temp");
        File gFile = new File("data", "1313.7001.gto");
        File dnaFile = new File("data", "g1.fna");
        File pFile = new File("data", "g1.faa");
        File dbFile = new File("data", "target.fa");
        SequenceStream s = Source.contigs.query(tempDir, gFile, 0);
        assertThat(s.isProtein(), equalTo(false));
        assertThat(s, instanceOf(DnaStream.class));
        assertThat(((DnaStream) s).getGeneticCode(), equalTo(11));
        s = Source.db.query(tempDir, dbFile, 0);
        assertThat(s.isProtein(), equalTo(false));
        assertThat(s, instanceOf(DnaStream.class));
        assertThat(((DnaStream) s).getGeneticCode(), equalTo(4));
        s = Source.dna.query(tempDir, dnaFile, 1);
        assertThat(s.isProtein(), equalTo(false));
        assertThat(s, instanceOf(DnaStream.class));
        assertThat(((DnaStream) s).getGeneticCode(), equalTo(1));
        s = Source.features.query(tempDir, gFile, 0);
        assertThat(s.isProtein(), equalTo(false));
        assertThat(s, instanceOf(DnaStream.class));
        assertThat(((DnaStream) s).getGeneticCode(), equalTo(11));
        s = Source.pegs.query(tempDir, gFile, 0);
        assertThat(s.isProtein(), equalTo(true));
        assertThat(s, instanceOf(ProteinStream.class));
        s = Source.prot.query(tempDir, pFile, 0);
        assertThat(s.isProtein(), equalTo(true));
        assertThat(s, instanceOf(ProteinStream.class));
    }

    /**
     * Test creating a DNA stream from a genome's contigs.
     * @throws IOException
     */
    @Test
    public void testGenomeDnaStream() throws IOException {
        Genome gto = new Genome(new File("data", "1313.7001.gto"));
        DnaStream stream = new DnaDataStream(gto);
        assertThat(stream.getGeneticCode(), equalTo(gto.getGeneticCode()));
        Iterator<Sequence> iter = stream.iterator();
        for (Contig contig : gto.getContigs()) {
            assertThat(iter.hasNext(), equalTo(true));
            Sequence seq = iter.next();
            assertThat(seq.getLabel(), equalTo(contig.getId()));
            assertThat(seq.getComment(), equalTo(contig.getDescription()));
            assertThat(seq.getSequence(), equalTo(contig.getSequence()));
        }
    }

}
