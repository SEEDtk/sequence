/**
 *
 */
package org.theseed.sequence.blast;

import java.io.File;
import java.io.IOException;

import org.theseed.genome.Genome;
import org.theseed.io.MarkerFile;
import org.theseed.sequence.DnaInputStream;
import org.theseed.sequence.ProteinInputStream;
import org.theseed.sequence.SequenceInputStream;
import org.theseed.utils.IDescribable;

/**
 * This enumeration represents a source for sequence information.  It can be used
 * to connect to a BLAST database or produce a sequence stream.
 *
 * @author Bruce Parrello
 */
public enum Source implements IDescribable {
    /** existing blast database */
    db(false, false, "Existing Blast Database"),
    /** DNA FASTA */
    dna(false, false, "DNA FASTA file"),
    /** protein FASTA */
    prot(false, false, "Protein FASTA file"),
    /** genome contigs */
    contigs(true, true, "Contigs in a GTO"),
    /** genome feature proteins */
    pegs(true, true, "Proteins in a GTO"),
    /** genome feature DNA */
    features(true, true, "Feature DNA in a GTO"),
    /** genome protein DNA */
    pegs_dna(true, true, "DNA of PEGs in a GTO"),
    /** genome RNA */
    rna(true, true, "RNA features in a GTO");

    // FIELDS
    private boolean needsTempFiles;
    private boolean needsGTO;
    private String description;

    private Source(boolean needsFile, boolean needsGTO, String description) {
        this.needsTempFiles = needsFile;
        this.description = description;
        this.needsGTO = needsGTO;
    }

    /**
     * @return a BLAST database for the specified source
     *
     * @param dbFile	source file containing sequence data
     * @param tempDir	directory for temporary files (contigs, pegs, features)
     * @param gc		genetic code (dna)
     * @param keep		TRUE to keep temporary files (contigs, pegs, features)
     *
     * @throws InterruptedException
     * @throws IOException
     */
    public BlastDB subject(File tempDir, File dbFile, int gc, boolean keep) throws IOException, InterruptedException {
        BlastDB retVal = null;
        // Insure we have a storage file if we are extracting sequences from a genome.  We also update the genetic
        // code to the genome's code.
        File tempFile = null;
        Genome genome = null;
        if (this.needsTempFiles) {
            genome = new Genome(dbFile);
            gc = genome.getGeneticCode();
            if (keep)
                tempFile = new File(tempDir, genome.getId() + "." + this.name() + ".fa");
            else
                tempFile = File.createTempFile("blast", ".fa", tempDir);
        }
        // Connect to the BLAST database.
        switch (this) {
        case db:
            retVal = BlastDB.load(dbFile);
            break;
        case dna:
            retVal = DnaBlastDB.create(dbFile, gc);
            break;
        case prot:
            retVal = ProteinBlastDB.create(dbFile);
            break;
        case contigs:
            genome.saveDna(tempFile);
            retVal = DnaBlastDB.create(tempFile, gc);
            retVal.setName("Contigs in " + genome.toString());
            break;
        case pegs:
            genome.savePegs(tempFile);
            retVal = ProteinBlastDB.create(tempFile);
            retVal.setName("Proteins in " + genome.toString());
            break;
        case features:
            genome.saveFeatures(tempFile);
            retVal = DnaBlastDB.create(tempFile, gc);
            retVal.setName("Features in " + genome.toString());
            break;
        case pegs_dna:
            genome.saveFeatures(tempFile, "CDS");
            retVal = DnaBlastDB.create(tempFile, gc);
            retVal.setName("Peg DNA in " + genome.toString());
            break;
        case rna:
            genome.saveFeatures(tempFile, "rna");
            retVal = DnaBlastDB.create(tempFile, gc);
            retVal.setName("RNAs in " + genome.toString());
            break;
        }
        // Insure we free up temporary files.
        if (! keep) {
            if (this.needsTempFiles)
                retVal.deleteOnExit();
            else if (this != db)
                retVal.cleanOnExit();
        }
        return retVal;
    }

    /**
     * @return a sequence data stream for the specified source
     *
     * @param qFile		file containing sequence data
     * @param tempDir	directory for working files (contigs, pegs, features)
     * @param gc		genetic code (dna)
     *
     * @throws IOException
     */
    public SequenceInputStream query(File tempDir, File qFile, int gc) throws IOException {
        SequenceInputStream retVal = null;
        // Insure we have a temporary storage file if we need one.
        File tempFile = null;
        Genome genome = null;
        if (this.needsTempFiles) {
            tempFile = File.createTempFile("query", ".fa", tempDir);
            genome = new Genome(qFile);
        }
        // Connect to the stream.
        switch (this) {
        case db:
            File checkFile = new File(qFile.getPath() + ".psq");
            if (checkFile.exists())
                retVal = new ProteinInputStream(qFile);
            else {
                int gCode = MarkerFile.readInt(new File(qFile.getPath() + ".gc"));
                retVal = new DnaInputStream(qFile, gCode);
            }
            break;
        case dna:
            retVal = new DnaInputStream(qFile, gc);
            break;
        case prot:
            retVal = new ProteinInputStream(qFile);
            break;
        case contigs:
            genome.saveDna(tempFile);
            retVal = new DnaInputStream(tempFile, genome.getGeneticCode());
            break;
        case pegs:
            genome.savePegs(tempFile);
            retVal = new ProteinInputStream(tempFile);
            break;
        case features:
            genome.saveFeatures(tempFile);
            retVal = new DnaInputStream(tempFile, genome.getGeneticCode());
            break;
        case pegs_dna:
            genome.saveFeatures(tempFile, "CDS");
            retVal = new DnaInputStream(tempFile, genome.getGeneticCode());
            break;
        default:
            break;
        }
        // Insure we free up the temporary file.
        if (this.needsTempFiles)
            tempFile.deleteOnExit();
        return retVal;
    }

    /**
     * @return TRUE if the specified file is compatible with this source type
     *
     * @param sourceFile	file to check
     */
    public boolean checkFile(File sourceFile) {
        return (this.needsGTO == sourceFile.getName().endsWith(".gto"));
    }

    @Override
    public String getDescription() {
        return this.description;
    }

}
