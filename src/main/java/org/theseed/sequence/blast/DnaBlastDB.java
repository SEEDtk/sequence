/**
 *
 */
package org.theseed.sequence.blast;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.theseed.genome.Genome;
import org.theseed.io.MarkerFile;
import org.theseed.sequence.DnaStream;
import org.theseed.sequence.ProteinStream;

/**
 * This is a DNA blast database.  Unlike a protein blast database, it knows its genetic code.
 *
 * @author Bruce Parrello
 *
 */
public class DnaBlastDB extends BlastDB {

    // FIELDS
    /** genetic code of this DNA blast database */
    private int geneticCode;
    /** array of file suffixes */
    private static final String[] SUFFIXES = new String[] { ".nin", ".nsq", ".nhr" };

    /**
     * Load an existing DNA blast database.
     *
     * @param fastaFile		name of the FASTA file containing the DNA
     */
    protected DnaBlastDB(File fastaFile) {
        this.setFile(fastaFile);
        File marker = markerFile();
        this.geneticCode = MarkerFile.readInt(marker);
        log.info("Using DNA blast database {} with genetic code {}.", fastaFile, this.geneticCode);
    }

    /**
     * @return the marker file for this database
     */
    protected File markerFile() {
        return new File(this.getFile().getPath() + ".gc");
    }

    public DnaBlastDB() { }

    /**
     * @return a new DNA blast database.
     *
     * @param fastaFile		name of the FASTA file containing the DNA
     * @param int			genetic code of the sequences
     *
     * @throws InterruptedException
     * @throws IOException
     */
    public static DnaBlastDB create(File fastaFile, int geneticCode) throws IOException, InterruptedException {
        DnaBlastDB retVal = new DnaBlastDB();
        retVal.setFile(fastaFile);
        retVal.geneticCode = geneticCode;
        MarkerFile.write(retVal.markerFile(), geneticCode);
        retVal.createDb();
        return retVal;
    }

    /**
     * @return a new DNA blast database created from a genome's contigs
     *
     * @param fastaFile		name of the FASTA file to contain the database
     *
     * @throws InterruptedException
     * @throws IOException
     */
    public static DnaBlastDB create(File fastaFile, Genome genome) throws IOException, InterruptedException {
        genome.saveDna(fastaFile);
        return DnaBlastDB.create(fastaFile, genome.getGeneticCode());
    }

    /**
     * @return a new DNA blast database created from a genome's features
     *
     * @param fastaFile		name of the FASTA file to contain the database
     *
     * @throws InterruptedException
     * @throws IOException
     */
    public static DnaBlastDB createFromFeatures(File fastaFile, Genome genome) throws IOException, InterruptedException {
        genome.saveFeatures(fastaFile);
        return DnaBlastDB.create(fastaFile, genome.getGeneticCode());
    }

    @Override
    protected String getDbParm() {
        return "nucl";
    }

    @Override
    public boolean isProtein() {
        return false;
    }

    @Override
    public List<BlastHit> blast(ProteinStream proteins, BlastParms parms)
            throws IOException, InterruptedException {
        BlastParms myParms = parms.clone().db_gencode(this.geneticCode);
        List<BlastHit> retVal = this.runBlast("tblastn", proteins, myParms);
        return retVal;
    }

    @Override
    public List<BlastHit> blast(DnaStream contigs, BlastParms parms) throws IOException, InterruptedException {
        BlastParms myParms = parms.clone();
        List<BlastHit> retVal = this.runBlast("blastn", contigs, myParms);
        return retVal;
    }

    /**
     * @return the genetic code of this DNA database
     */
    public int getGeneticCode() {
        return geneticCode;
    }

    @Override
    protected String[] getSuffixes() {
        return SUFFIXES;
    }


}
