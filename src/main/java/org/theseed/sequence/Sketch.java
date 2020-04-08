/**
 *
 */
package org.theseed.sequence;

import java.io.Serializable;

/**
 * This object represents a kmer signature.  Each signature is associated with an identifying name; the name serves
 * as a key for accessing the named object-- it could be a protein family, a genome, a feature, or a role.
 *
 * @author Bruce Parrello
 *
 */
public class Sketch implements Serializable {

    // FIELDS

    private static final long serialVersionUID = -1880805874119392439L;

    /** signature for this entry's sequence */
    private int[] signature;
    /** target object associated with the sequence */
    private String name;

    /**
     * Construct a hash entry.
     *
     * @param seq		kmer object to use for comparison
     * @param string	string associated with the sequence that generated the kmers
     */
    public Sketch(int[] sk, String string) {
        this.signature = sk;
        this.name = string;
    }

    /**
     * @return the distance between another sequence and this entry's sequence
     *
     * @param otherSig	signature of the other sequence
     */
    public double distance(int[] otherSig) {
        return SequenceKmers.signatureDistance(this.signature, otherSig);
    }

    /**
     * @return the sequence target
     */
    public String getName() {
        return name;
    }

    /**
     * @return the signature for this sketch
     */
    public int[] getSignature() {
        return this.signature;
    }

    /**
     * @return the distance between this sketch and another sketch
     *
     * @param sketch	other sketch to compare
     */
    public double distance(Sketch sketch) {
        return this.distance(sketch.signature);
    }

}
