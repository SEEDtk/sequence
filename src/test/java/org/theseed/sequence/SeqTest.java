/**
 *
 */
package org.theseed.sequence;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

/**
 * Unit test for minhash stuff.
 *
 * @author Bruce Parrello
 *
 */
public class SeqTest extends TestCase {

    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public SeqTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( SeqTest.class );
    }

    public void testObjectivism() {
        assertThat("A", equalTo("A"));
    }

}
