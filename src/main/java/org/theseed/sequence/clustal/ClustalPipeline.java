/**
 *
 */
package org.theseed.sequence.clustal;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

import org.slf4j.Logger;

import org.slf4j.LoggerFactory;
import org.theseed.io.ErrorQueue;
import org.theseed.io.LineReader;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;
import org.theseed.utils.Parms;
import org.theseed.utils.ProcessUtils;

/**
 * This class reads a sequence stream and produces an output stream that contains the original
 * sequences properly aligned.  It uses Clustal Omega to transform the sequences and presents
 * the output as an iterable.
 *
 * The CLUSTAL_PATH environment variable must be set to the location of the Clustal application.
 *
 * Currently no optional parameters are supported.  This may change. We take input from a
 * file and produce the alignment on the standard output, which is read into memory as
 * Sequence objects.  The list of sequences is presented as an iterable.  Calling the
 * "iterator" method starts the clustal process in the background.  Given the nature
 * of the program, it will chug invisibly for a while, then produce a bunch of output.
 *
 * The output sequences will be ordered by sequence ID.  Duplicate sequence IDs will generate
 * a warning.
 *
 * @author Bruce Parrello
 *
 */
public class ClustalPipeline {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ClustalPipeline.class);
    /** path to CLUSTAL software */
    protected static String CLUSTAL_PATH = System.getenv("CLUSTAL_PATH");
    /** parameter object */
    private Parms clustalParms;
    /** output sequence hash */
    private ConcurrentHashMap<String, Sequence> alignment;

    /**
     * Construct a new clustal pipeline.
     *
     * @param inputFile		file of sequences to align
     */
    public ClustalPipeline(File inputFile) {
        this.clustalParms = new Parms().set("-i", inputFile.getAbsolutePath());
    }

    public List<Sequence> run() throws IOException, InterruptedException {
        // Set up the command line.
        String program = new File(CLUSTAL_PATH, "clustalo").getAbsolutePath();
        List<String> command = this.clustalParms.get(program);
        ProcessBuilder builder = new ProcessBuilder(command);
        // Create the alignment storage buffer.
        this.alignment = new ConcurrentHashMap<String, Sequence>(100);
        // Start Clustal.
        Process process = builder.start();
        // Create the output consumers.
        try (LineReader errorStream = new LineReader(process.getErrorStream());
                FastaInputStream outputStream = new FastaInputStream(process.getInputStream())) {
            // Queue up the error messages.  There will be nothing unless something bad happens.
            // The error message list is owned by the consumer and we CANNOT modify it.
            List<String> errorBuffer = new ArrayList<String>();
            ErrorQueue errorThread = new ErrorQueue(errorStream, errorBuffer);
            errorThread.start();
            // Set up a consumer to process the sequence output.
            SeqConsumer outputThread = this.new SeqConsumer(outputStream);
            outputThread.start();
            // Clean up the process.
            errorThread.join();
            outputThread.join();
            // Wait for the process to finish.
            int exitCode = ProcessUtils.finishProcess("CLUSTALO", process, errorBuffer);
            if (exitCode != 0)
                throw new RuntimeException("CLUSTALO call failed.");
        }
        // Build the list of sequences from the hash map and arrange them in order.
        List<String> seqKeys = new ArrayList<String>(this.alignment.keySet());
        Collections.sort(seqKeys);
        // Produce the output list.
        List<Sequence> retVal = seqKeys.stream().map(x -> this.alignment.get(x)).collect(Collectors.toList());
        return retVal;
    }

    /**
     * This is aested class to read aligned sequences and save them.  The tricky part is we have
     * to resolve duplicate IDs in the hash map.
     */
    protected class SeqConsumer extends Thread {

        // FIELDS
        /** input stream */
        private FastaInputStream sequencesIn;
        /** output hash */
        private Map<String, Sequence> sequencesOut;

        /**
         * Construct the sequence consumer.  Note that the output from the process comes in as an input stream.
         *
         * @param outpuStream		output from the CLUSTAL process
         */
        public SeqConsumer(FastaInputStream outputStream) {
            this.sequencesIn = outputStream;
            this.sequencesOut = ClustalPipeline.this.alignment;
        }

        @Override
        public void run() {
            // Loop through the sequences, adding them to the hash.  Resolve duplicates by adding a spaced number to the ID.
            // The ID is stripped before we send anything back to the client anyway.
            for (Sequence seq : this.sequencesIn) {
                String seqId = seq.getLabel();
                if (this.sequencesOut.containsKey(seqId)) {
                    int count = 1;
                    String newSeqId = seqId + " dup";
                    while (this.sequencesOut.containsKey(newSeqId)) {
                        count++;
                        newSeqId = String.format("%s dup %d", seqId, count);
                    }
                    seqId = newSeqId;
                }
                this.sequencesOut.put(seqId, seq);
            }
        }

    }

}
