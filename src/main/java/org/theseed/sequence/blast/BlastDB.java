/**
 *
 */
package org.theseed.sequence.blast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.ErrorQueue;
import org.theseed.io.LineReader;
import org.theseed.reports.BlastInfo;
import org.theseed.reports.Color;
import org.theseed.reports.IBlastReporter;
import org.theseed.sequence.DnaStream;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.ProteinStream;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.SequenceDataStream;
import org.theseed.sequence.SequenceInputStream;
import org.theseed.sequence.SequenceStream;
import org.theseed.utils.IDescribable;
import org.theseed.utils.ProcessUtils;

import ch.qos.logback.classic.Level;
import ch.qos.logback.classic.LoggerContext;

/**
 * This object represents a blast database stored in a file.  It can be used to reference an
 * existing blast database, or create one from a FASTA file or a genome object.
 *
 * @author Bruce Parrello
 *
 */
public abstract class BlastDB {

    // FIELDS

    /** log stream */
    protected static Logger log = LoggerFactory.getLogger(BlastDB.class);
    /** path to BLAST+ software */
    protected static String BLAST_PATH = System.getenv("BLAST_PATH");
    /** file name of the blast database */
    private File dbFile;
    /** type of last BLAST run */
    private String blastType;
    /** parm string of last BLAST run */
    private String blastParms;
    /** buffer file for BLAST input */
    private File tempFile;
    /** run counter for blasting */
    private static int runCount = 0;
    /** display name for database */
    private String name;

    /**
     * Construct a blank BLAST database.
     */
    protected BlastDB() {
        try {
            // Create the temporary buffer file.
            this.tempFile = File.createTempFile("temp", "fasta");
            this.tempFile.deleteOnExit();
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }

    }

    /**
     * Create a blast database for a FASTA file.  This should only be used for blast databases
     * that already exist.
     *
     * @param fastaFile		the FASTA file to be used
     */
    public static BlastDB load(File fastaFile) throws IOException, InterruptedException {
        BlastDB retVal = null;
        File testFile = new File(fastaFile.getPath() + ".psq");
        if (testFile.exists()) {
            retVal = new ProteinBlastDB(fastaFile);
        } else {
            testFile = new File(fastaFile.getPath() + ".nsq");
            if (testFile.exists()) {
                retVal = new DnaBlastDB(fastaFile);
            } else {
                throw new FileNotFoundException("No support files found for " + fastaFile + ".");
            }
        }
        if (testFile.lastModified() < fastaFile.lastModified()) {
            log.info("BLAST database expired-- recreating.");
            retVal.createDb();
        }
        return retVal;
    }

    /**
     * Specify the logging level relating to BLAST invocations.
     *
     * @param level		desired logging level
     */
    public static void configureLogging(String level) {
        LoggerContext loggerContext = (LoggerContext) LoggerFactory.getILoggerFactory();
        ch.qos.logback.classic.Logger logger = loggerContext.getLogger("org.theseed.sequence.blast");
        logger.setLevel(Level.toLevel(level));
        log.info("BLAST logging set to {}.", level);

    }


    /**
     * Build the database support files using makeblastdb.
     *
     * @throws InterruptedException
     * @throws IOException
     */
    protected void createDb() throws IOException, InterruptedException {
        // Create the command.
        ProcessBuilder makeProcess = new ProcessBuilder(new File(BLAST_PATH, "makeblastdb").getPath(),
                "-in", this.dbFile.getAbsoluteFile().toString(), "-dbtype", this.getDbParm());
        // Merge STDERR into STDOUT, then start the process.
        makeProcess.redirectErrorStream(true);
        log.info("Running blast database make for {} and type {}.", this.dbFile, this.getDbParm());
        Process process = makeProcess.start();
        // Get a reader for the process output.
        try (LineReader reader = new LineReader(process.getInputStream())) {
            for (String line : reader)
                log.debug("*   {}", line);
            int exitCode = process.waitFor();
            if (exitCode != 0) {
                String commandString = StringUtils.join(makeProcess.command(), " ");
                log.error("Failing makeBlastDB command with exit code {} was {}.", exitCode, commandString);
                throw new RuntimeException("Exit code " + exitCode + " from BlastDB creation.");
            }
            log.info("Blast database created.");
        }
    }

    /**
     * @return the database creation type parameter
     */
    protected abstract String getDbParm();

    /**
     * @return the database type
     */
    public abstract boolean isProtein();

    /**
     * BLAST protein sequences against this database.
     *
     * @param proteins	list of protein sequences
     * @param parms		BLAST parameters
     *
     * @return a list of BLAST matches
     *
     * @throws InterruptedException
     * @throws IOException
     */
    public abstract List<BlastHit> blast(ProteinStream proteins, BlastParms parms);

    /**
     * BLAST DNA sequences against this database.
     *
     * @param contigs	list of DNA sequences
     * @param parms		BLAST parameters
     *
     * @return a list of BLAST matches
     *
     * @throws InterruptedException
     * @throws IOException
     */
    public abstract List<BlastHit> blast(DnaStream contigs, BlastParms parms);

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
    public abstract List<BlastHit> psiBlast(File pssmFile, BlastParms parms, Map<String, String> qMap);

    /**
     * Run BLAST against the specified sequences with the specified parameters.
     *
     * @param blastProgram	name of the BLAST program to run
     * @param seqs			list of query sequences
     * @param myParms		current blast parameters
     * @param qProt			TRUE if the query sequence is protein
     *
     * @return a list of BlastHit objects representing the hits
     */
    protected List<BlastHit> runBlast(String blastProgram, SequenceStream seqs, BlastParms myParms) {
        List<BlastHit> retVal = null;
        // Save the command specs.
        saveCommand(blastProgram, myParms);
        // Create a hash to map sequence labels to comments.
        Map<String, String> qMap = new ConcurrentHashMap<String, String>();
        // Copy the input to a temporary buffer file.  Unfortunately, this is required to
        // get BLAST working in all environments.
        try (FastaOutputStream queryStream = new FastaOutputStream(this.tempFile)) {
            // Write the sequences to the temp file and save a map of IDs to comments.
            for (Sequence seq : seqs) {
                // Store this sequence's comment in the map.
                qMap.put(seq.getLabel(), seq.getComment());
                // Write the sequence to the BLAST process.
                queryStream.write(seq);
            }
            queryStream.close();
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
        log.debug("Submitting {} query sequences.", qMap.size());
        // Attach the query file as input.
        myParms.set("-query", this.tempFile.getPath());
        // Execute the BLAST and create the list of hits.
        retVal = processBlast(blastProgram, myParms, qMap, seqs.isProtein());
        // Return the result.
        return retVal;
    }

    /**
     * @param blastProgram
     * @param parms
     */
    protected void saveCommand(String blastProgram, BlastParms parms) {
        this.blastType = blastProgram;
        this.blastParms = parms.toString();
    }

    /**
     * This executes a BLAST.  It is a reduced form that assumes a file is being used as the input,
     * though the caller must know the sequence descriptions already.  It can be used for more
     * exotic types of BLASTing that don't involve sequences.  The query file must already be
     * specified in the parameter object.
     *
     * @param blastProgram		name of the BLAST program to run
     * @param myParms			incoming parameters (will be modified)
     * @param qMap				map of query sequence IDs to query sequence descriptions
     * @param protsIn			TRUE if the input is protein sequences, FALSE if it is DNA
     *
     * @return a list of blast hits
     *
     */
    protected List<BlastHit> processBlast(String blastProgram, BlastParms myParms, Map<String, String> qMap,
            boolean protsIn) {
        // This will be the return value.
        List<BlastHit> retVal = new ArrayList<BlastHit>();
        // Set up the parameters and form the command.
        myParms.set("-outfmt", BlastHit.OUT_FORMAT);
        myParms.set("-db", this.dbFile.getAbsolutePath());
        List<String> command = myParms.get(new File(BLAST_PATH, blastProgram).getPath());
        ProcessBuilder blastCommand = new ProcessBuilder(command);
        // Start the BLAST.
        try {
            log.debug("Running BLAST command {}.", blastProgram);
            Process blastProcess = blastCommand.start();
            // Open the streams.
            try (LineReader blastReader = new LineReader(blastProcess.getInputStream());
                    LineReader logReader = new LineReader(blastProcess.getErrorStream())) {
                // Create a thread to read the blast output.
                HitConsumer hitReader = new HitConsumer(blastReader, protsIn, qMap,
                        myParms, retVal);
                hitReader.start();
                // Create a thread to save error log messages.  Generally there are none.
                List<String> messages = new ArrayList<String>(30);
                ErrorQueue errorReader = new ErrorQueue(logReader, messages);
                errorReader.start();
                // Clean up the process.
                hitReader.join();
                errorReader.join();
                ProcessUtils.finishProcess("BLAST", blastProcess, messages);
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        } catch (InterruptedException e) {
            throw new RuntimeException("Interrupted BLAST execution: " + e.getMessage(), e);
        }

        return retVal;
    }

    /**
     * This class reads hits and creates the list of BlastHit objects.
     */
    private class HitConsumer extends Thread {

        // FIELDS
        private LineReader blastOutput;
        private List<BlastHit> hitList;
        private boolean proteinFlag;
        private Map<String, String> qMap;
        private BlastParms parms;

        /**
         * Construct a thread to create BLAST hits from the output stream.
         *
         * @param blastOutput	line reader for BLAST output
         * @param isProtein		TRUE if the queries are protein sequences, else FALSE
         * @param qMap			map of query sequence IDs to descriptions
         * @param parms			current BLAST parameters
         * @param hitList		output list for blast hit objects
         */
        protected HitConsumer(LineReader blastOutput, boolean isProtein, Map<String, String> qMap,
                BlastParms parms, List<BlastHit> hitList) {
            this.blastOutput = blastOutput;
            this.hitList = hitList;
            this.proteinFlag = isProtein;
            this.qMap = qMap;
            this.parms = parms;
        }

        @Override
        public void run() {
            int count = 0;
            for (String line : blastOutput) {
                BlastHit result = new BlastHit(line, qMap, proteinFlag, BlastDB.this.isProtein());
                count++;
                if (this.parms.acceptable(result))
                    hitList.add(result);
            }
            runCount++;
            log.info("{} results from BLAST #{}, {} returned.", count, runCount, hitList.size());
        }
    }

    /**
     * @return the file name of this database
     */
    public File getFile() {
        return this.dbFile;
    }

    /**
     * Specify the FASTA file for this database.
     *
     * @param file	file containing the main FASTA sequences
     */
    protected void setFile(File file) {
        this.dbFile = file;
    }

    /**
     * Specify that this FASTA database be deleted on exit.
     */
    public void deleteOnExit() {
        this.dbFile.deleteOnExit();
        cleanOnExit();
    }


    /**
     * Specify that this FASTA database's support files be deleted on exit.
     */
    public void cleanOnExit() {
        String[] suffixes = this.getSuffixes();
        String path = this.dbFile.getPath();
        for (String suffix : suffixes) {
            File serviceFile = new File(path + suffix);
            serviceFile.deleteOnExit();
        }
    }


    /**
     * @return a list of the suffixes
     */
    protected abstract String[] getSuffixes();


    /**
     * @return the type of the last BLAST run
     */
    public String getBlastType() {
        return blastType;
    }


    /**
     * @return the parameter string for the last BLAST run
     */
    public String getBlastParms() {
        return blastParms;
    }

    /**
     * @return the name of this blast database
     */
    public String getName() {
        String retVal = this.name;
        if (retVal == null)
            retVal = this.dbFile.getName();
        return retVal;
    }

    /**
     * Specify the name of this database.
     *
     * @param newName	new name to use
     */
    public void setName(String newName) {
        this.name = newName;
    }

    /**
     * Execute a full BLAST and output the results.
     *
     * @param query			query stream
     * @param subject		subject database
     * @param batchSize		query batch size
     * @param parms			blast parameters
     * @param reporter		output object
     *
     * @throws IOException
     * @throws InterruptedException
     */
    public void runBlast(SequenceInputStream query, int batchSize,
            BlastParms parms, IBlastReporter reporter) throws IOException, InterruptedException {
        // Loop through the batches in the query stream.
        Iterator<SequenceDataStream> batcher = query.batchIterator(batchSize);
        int batchCount = 0;
        int seqCount = 0;
        int hitCount = 0;
        int missCount = 0;
        while (batcher.hasNext()) {
            SequenceDataStream batch = batcher.next();
            // This set is used to track queries without hits.
            Set<String> queryIds = batch.stream().map(x -> x.getLabel()).collect(Collectors.toSet());
            // Record the batch.
            batchCount++;
            seqCount += batch.size();
            log.info("Processing input batch {}, {} sequences read.", batchCount, seqCount);
            // BLAST the query against the subject.
            List<BlastHit> results = batch.blast(this, parms);
            for (BlastHit hit : results) {
                queryIds.remove(hit.getQueryId());
                reporter.recordHit(hit);
                hitCount++;
            }
            missCount += queryIds.size();
            if (log.isDebugEnabled()) {
                for (String queryId : queryIds)
                    log.debug("{} had no hits.", queryId);
            }
        }
        log.info("BLAST found {} hits, {} queries had no hits.", hitCount, missCount);
        // Save our parameters and statistics.
        BlastInfo blastInfo = new BlastInfo(this.getBlastParms(), seqCount, missCount, hitCount);
        // Write the report.
        log.info("Writing report.  {} total sequences in {} batches were processed.", seqCount, batchCount);
        reporter.writeReport(this.getBlastType().toUpperCase()
                + " run against " + this.getName(), blastInfo);
    }

    /**
     * This enum indicates how the color of a hit is determined.
     */
    public enum ColorType implements IDescribable {
        /** percent identity */
        ident("percent identity"),
        /** percent similarity */
        sim("percent similarity"),
        /** percent coverage of target sequence */
        covg("percent coverage");

        // FIELDS
        private String name;

        private ColorType(String name) {
            this.name = name;
        }

        /**
         * @return the color for a blast hit by the specified target sequence.
         *
         * @param target	target sequence forming the denominator
         * @param hit		blast hit whose color is needed
         */
        public Color computeColor(BlastHit.SeqData target, BlastHit hit) {
            double fraction = 1.0;
            switch (this) {
            case ident:
                fraction = ((double) hit.getNumIdentical()) / hit.getAlignLen();
                break;
            case sim:
                fraction = ((double) hit.getNumSimilar()) / hit.getAlignLen();
                break;
            case covg:
                fraction = ((double) hit.getNumSimilar()) / target.getLen();
                break;
            }
            Color retVal;
            if (fraction >= 1.0)
                retVal = Color.BLUE;
            else if (fraction >= 0.9)
                retVal = Color.DARK_GREEN.brighten((1.0 - fraction)*5);
            else if (fraction >= 0.7)
                retVal = Color.ORANGE.brighten((0.9 - fraction)*2.5);
            else if (fraction >= 0.5)
                retVal = Color.RED.brighten((0.7 - fraction)*2.5);
            else
                retVal = Color.DARK_GRAY.brighten(0.5 - fraction);
            return retVal;
        }

        public String getDescription() {
            return this.name;
        }

    }

    /**
     * This enum indicates the high-level sort for output.
     */
    public enum SortType implements IDescribable {
        /** sort by query, list subjects within query */
        QUERY(BlastHit.QUERY, "queries", "Sort by Query Sequence"),
        /** sort by subject, list queries within subject */
        SUBJECT(BlastHit.SUBJECT, "subject sequences", "Sort by Subject Sequence");

        // FIELDS
        private int sortIdx;
        private int otherIdx;
        private String plural;
        private String description;

        private SortType(int idx, String plural, String description) {
            this.sortIdx = idx;
            this.otherIdx = 1 - idx;
            this.plural = plural;
            this.description = description;
        }

        /**
         * @return the ID of the sequence being sorted on
         *
         * @param hit	blast hit whose sort ID is desired
         */
        public String idOf(BlastHit hit) {
            return hit.getData(sortIdx).getId();
        }

        /**
         * @return the ID of the sequence being hit
         *
         * @param hit	blast hit whose sequence ID is desired
         */
        public String targetOf(BlastHit hit) {
            return hit.getData(otherIdx).getId();
        }

        /**
         * @return TRUE if the hit has the specified target ID
         *
         * @param hit		blast hit to check
         * @param otherId	ID of the desired target
         */
        public boolean targetEquals(BlastHit hit, String otherId) {
            return otherId.contentEquals(this.targetOf(hit));
        }

        /**
         * @return the sorting sequence data for the specified hit
         */
        public BlastHit.SeqData data(BlastHit hit) {
            return hit.getData(this.sortIdx);
        }

        /**
         * @return the target sequence data for the specified hit
         */
        public BlastHit.SeqData target(BlastHit hit) {
            return hit.getData(this.otherIdx);
        }

        /**
         * @return the plural phrase for the things being sorted
         */
        public String getPlural() {
            return this.plural;
        }

        @Override
        public String getDescription() {
            return this.description;
        }

    }

}
