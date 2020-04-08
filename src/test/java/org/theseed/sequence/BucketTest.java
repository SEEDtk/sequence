/**
 *
 */
package org.theseed.sequence;

import junit.framework.TestCase;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

/**
 * @author Bruce Parrello
 *
 */
public class BucketTest extends TestCase {

    public void testObjectivism() {
        assertThat("A", equalTo("A"));
    }

}
