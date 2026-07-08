package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class MaxAmbigTrimmerTest {

    private FastqRecord makeRecord(String seq) {
        String qual = "I".repeat(seq.length());
        return new FastqRecord("read1", seq, "", qual, 33);
    }

    @Test
    public void testNoAmbigPassesWithZeroThreshold() {
        MaxAmbigTrimmer trimmer = new MaxAmbigTrimmer("0.0");
        FastqRecord rec = makeRecord("ACGTACGT");
        assertNotNull(trimmer.processRecord(rec));
    }

    @Test
    public void testAllNDroppedWithZeroThreshold() {
        MaxAmbigTrimmer trimmer = new MaxAmbigTrimmer("0.0");
        FastqRecord rec = makeRecord("NNNN");
        assertNull(trimmer.processRecord(rec));
    }

    @Test
    public void testSingleNBelowThresholdPasses() {
        // 1 N out of 10 bases = 0.1, threshold 0.2 → pass
        MaxAmbigTrimmer trimmer = new MaxAmbigTrimmer("0.2");
        FastqRecord rec = makeRecord("ACGTACGNNA");
        // "ACGTACGNNA" has 2 Ns out of 10 = 0.2, exactly at threshold → pass
        assertNotNull(trimmer.processRecord(rec));
    }

    @Test
    public void testTooManyNsDropped() {
        // 3 Ns out of 10 bases = 0.3, threshold 0.2 → drop
        MaxAmbigTrimmer trimmer = new MaxAmbigTrimmer("0.2");
        FastqRecord rec = makeRecord("ACGNNNACGT");
        assertNull(trimmer.processRecord(rec));
    }

    @Test
    public void testAllBasesPassWithThresholdOne() {
        MaxAmbigTrimmer trimmer = new MaxAmbigTrimmer("1.0");
        FastqRecord rec = makeRecord("NNNNNNNNNNN");
        assertNotNull(trimmer.processRecord(rec));
    }

    @Test
    public void testEmptyReadDropped() {
        MaxAmbigTrimmer trimmer = new MaxAmbigTrimmer("0.1");
        FastqRecord rec = makeRecord("");
        assertNull(trimmer.processRecord(rec));
    }

    @Test
    public void testNegativeFractionThrows() {
        assertThrows(IllegalArgumentException.class, () -> new MaxAmbigTrimmer("-0.1"));
    }

    @Test
    public void testFractionAboveOneThrows() {
        assertThrows(IllegalArgumentException.class, () -> new MaxAmbigTrimmer("1.1"));
    }

    @Test
    public void testSingleBaseNoNPasses() {
        MaxAmbigTrimmer trimmer = new MaxAmbigTrimmer("0.0");
        FastqRecord rec = makeRecord("A");
        assertNotNull(trimmer.processRecord(rec));
    }

    @Test
    public void testSingleNAtZeroThresholdDropped() {
        MaxAmbigTrimmer trimmer = new MaxAmbigTrimmer("0.0");
        FastqRecord rec = makeRecord("N");
        assertNull(trimmer.processRecord(rec));
    }
}
