/**
 *
 */
package org.theseed.sequence.blast;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

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
    private static final String[] SUFFIXES = new String[] { ".nin", ".nsq", ".nhr", ".gc" };

    /**
     * Load an existing DNA blast database.
     *
     * @param fastaFile		name of the FASTA file containing the DNA
     */
    protected DnaBlastDB(File fastaFile) {
        this.setFile(fastaFile);
        getGC(this);
        log.info("Using DNA blast database {} with genetic code {}.", fastaFile, this.geneticCode);
    }

    /**
     * Determine the genetic code for an existing DNA BLAST database.
     *
     * @param db			database of interest
     */
    private static void getGC(DnaBlastDB db) {
        File marker = db.markerFile();
        db.geneticCode = MarkerFile.readInt(marker);
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
     * @param geneticCode	genetic code of the sequences
     *
     * @throws InterruptedException
     * @throws IOException
     */
    public static DnaBlastDB create(File fastaFile, int geneticCode) throws IOException, InterruptedException {
        DnaBlastDB retVal = new DnaBlastDB();
        retVal.setFile(fastaFile);
        setGC(geneticCode, retVal);
        retVal.createDb();
        return retVal;
    }

    /**
     * Specify the genetic code for a DNA database.
     *
     * @param geneticCode	genetic code to store
     * @param db			DNA database to store it in
     */
    private static void setGC(int geneticCode, DnaBlastDB db) {
        db.geneticCode = geneticCode;
        MarkerFile.write(db.markerFile(), geneticCode);
    }

    /**
     * @return a new or existing DNA blast database
     *
     * @param fastaFile		name of the fasta file containing the DNA sequences
     * @param geneticCode	genetic code of the sequences
     *
     * @throws IOException
     * @throws InterruptedException
     */
    public static DnaBlastDB createOrLoad(File fastaFile, int geneticCode) throws IOException, InterruptedException {
        DnaBlastDB retVal = new DnaBlastDB();
        retVal.setFile(fastaFile);
        File testFile = new File(fastaFile.getPath() + ".nsq");
        if (! testFile.exists() || testFile.lastModified() < fastaFile.lastModified()) {
            setGC(geneticCode, retVal);
            retVal.createDb();
        } else
            getGC(retVal);
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
        DnaBlastDB retVal = DnaBlastDB.create(fastaFile, genome.getGeneticCode());
        return retVal;
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
        DnaBlastDB retVal = DnaBlastDB.create(fastaFile, genome.getGeneticCode());
        return retVal;
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
    public List<BlastHit> blast(ProteinStream proteins, BlastParms parms) {
        BlastParms myParms = parms.clone().db_gencode(this.geneticCode);
        List<BlastHit> retVal = this.runBlast("tblastn", proteins, myParms);
        return retVal;
    }

    @Override
    public List<BlastHit> blast(DnaStream contigs, BlastParms parms) {
        BlastParms myParms = parms.clone();
        myParms.customizeForBlastn();
        List<BlastHit> retVal = this.runBlast("blastn", contigs, myParms);
        return retVal;
    }

    /**
     * Run a protein PSI-BLAST against this database.
     *
     * @param pssmFile	protein alignment scoring file for input
     * @param parms		BLAST parameters
     * @param qMap		map of query sequence IDs to descriptions
     *
     * @return the list of blast hits found
     *
     */
    public List<BlastHit> psiBlast(File pssmFile, BlastParms parms, Map<String, String> qMap) {
        this.saveCommand("tblastn", parms);
        BlastParms myParms = parms.clone().set("-in_pssm", pssmFile.getPath());
        List<BlastHit> retVal = this.processBlast("tblastn", myParms, qMap, true);
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
