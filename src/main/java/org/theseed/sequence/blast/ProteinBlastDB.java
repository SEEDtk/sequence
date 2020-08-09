/**
 *
 */
package org.theseed.sequence.blast;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import org.theseed.genome.Genome;
import org.theseed.sequence.DnaStream;
import org.theseed.sequence.ProteinStream;

/**
 * This manages a protein BLAST database.  The database
 * @author parrello
 *
 */
public class ProteinBlastDB extends BlastDB {

    // FIELDS
    /** array of file suffixes */
    private static final String[] SUFFIXES = new String[] { ".pin", ".psq", ".phr" };

    /**
     * Load an existing protein blast database.
     *
     * @param fastaFile		name of the fasta file containing the protein sequences
     */
    protected ProteinBlastDB(File fastaFile) {
        this.setFile(fastaFile);
        log.info("Using protein blast database {}.", fastaFile);
    }

    private ProteinBlastDB() { }

    /**
     * @return a new protein blast database
     *
     * @param fastaFile		name of the fasta file containing the protein sequences
     *
     * @throws IOException
     * @throws InterruptedException
     */
    public static ProteinBlastDB create(File fastaFile) throws IOException, InterruptedException {
        ProteinBlastDB retVal = new ProteinBlastDB();
        retVal.setFile(fastaFile);
        retVal.createDb();
        return retVal;
    }

    /**
     * @return a new or existing protein blast database
     *
     * @param fastaFile		name of the fasta file containing the protein sequences
     *
     * @throws IOException
     * @throws InterruptedException
     */
    public static ProteinBlastDB createOrLoad(File fastaFile) throws IOException, InterruptedException {
        ProteinBlastDB retVal = new ProteinBlastDB();
        retVal.setFile(fastaFile);
        File testFile = new File(fastaFile.getPath() + ".psq");
        if (! testFile.exists() || testFile.lastModified() < fastaFile.lastModified())
            retVal.createDb();
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
     * @throws IOException
     * @throws InterruptedException
     */
    public List<BlastHit> psiBlast(File pssmFile, BlastParms parms, Map<String, String> qMap) {
        this.saveCommand("psiblast", parms);
        BlastParms myParms = parms.clone().set("-in_pssm", pssmFile.getPath());
        List<BlastHit> retVal = this.processBlast("psiblast", myParms, qMap, true);
        return retVal;
    }

    /**
     * @return a new protein blast database created from a genome's PEGs
     *
     * @param fastaFile		name of the fasta file to contain the protein sequences
     * @param genome		genome containing the PEGs of interest
     *
     * @throws IOException
     * @throws InterruptedException
     */
    public static ProteinBlastDB create(File fastaFile, Genome genome) throws IOException, InterruptedException {
        genome.savePegs(fastaFile);
        return ProteinBlastDB.create(fastaFile);
    }

    @Override
    protected String getDbParm() {
        return "prot";
    }

    @Override
    public boolean isProtein() {
        return true;
    }

    @Override
    public List<BlastHit> blast(ProteinStream proteins, BlastParms parms) {
        BlastParms myParms = parms.clone();
        List<BlastHit> retVal = this.runBlast("blastp", proteins, myParms);
        return retVal;
    }

    @Override
    public List<BlastHit> blast(DnaStream contigs, BlastParms parms) {
        BlastParms myParms = parms.clone().query_gencode(contigs.getGeneticCode());
        List<BlastHit> retVal = this.runBlast("blastx", contigs, myParms);
        return retVal;
    }

    @Override
    protected String[] getSuffixes() {
        return SUFFIXES;
    }

}
