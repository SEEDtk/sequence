/**
 *
 */
package org.theseed.sequence.blast;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;
import org.theseed.genome.Genome;
import org.theseed.io.LineReader;
import org.theseed.sequence.DnaStream;
import org.theseed.sequence.ProteinStream;

/**
 * This manages a protein BLAST database.
 *
 */
public class ProteinBlastDB extends BlastDB {

    // FIELDS
    /** array of file suffixes */
    public static final String[] SUFFIXES = new String[] { ".pin", ".psq", ".phr", ".db.dmnd" };
    /** set of programs supported by Diamond */
    private static final Set<String> DIAMOND_TYPES = Set.of("blastp", "blastx");

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
        if (DIAMOND_PATH != null)
            retVal.createDiamondDb();
        return retVal;
    }

    /**
     * Create the DIAMOND helper database for our FASTA file.
     *
     * @throws IOException
     * @throws InterruptedException
     */
    private void createDiamondDb() throws IOException, InterruptedException {
        // Create the command.
        String fileName = this.getFile().getAbsoluteFile().toString();
        ProcessBuilder makeProcess = new ProcessBuilder(new File(DIAMOND_PATH, "diamond").getPath(),
                "makedb", "--in", fileName, "--db", fileName + ".db");
        // Merge STDERR into STDOUT, then start the process.
        makeProcess.redirectErrorStream(true);
        log.info("Running diamond database make for {}.", this.getFile());
        Process process = makeProcess.start();
        // Get a reader for the process output.
        try (LineReader reader = new LineReader(process.getInputStream())) {
            for (String line : reader)
                log.debug("*   {}", line);
            int exitCode = process.waitFor();
            if (exitCode != 0) {
                String commandString = StringUtils.join(makeProcess.command(), " ");
                log.error("Failing diamond command with exit code {} was {}.", exitCode, commandString);
                throw new RuntimeException("Exit code " + exitCode + " from Diamond DB creation.");
            }
            log.info("Diamond database created.");
        }
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
        if (DIAMOND_PATH != null) {
            testFile = new File(fastaFile.getPath() + ".db.dmnd");
            if (! testFile.exists() || testFile.lastModified() < fastaFile.lastModified())
                retVal.createDiamondDb();
        }
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
        ProteinBlastDB retVal = ProteinBlastDB.create(fastaFile);
        return retVal;
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

    @Override
    protected List<String> computeBlastCommand(String blastProgram, BlastParms myParms) {
        List<String> retVal;
        if (DIAMOND_PATH != null && DIAMOND_TYPES.contains(blastProgram)) {
            // Diamond requires its name followed by the blast command.
            String commandString = new File(DIAMOND_PATH, "diamond").toString();
            retVal = myParms.getForDiamond(commandString, blastProgram);
        } else
            retVal = super.computeBlastCommand(blastProgram, myParms);
        return retVal;
    }

}
