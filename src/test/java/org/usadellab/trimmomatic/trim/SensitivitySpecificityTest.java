package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class SensitivitySpecificityTest {

    private String makeSequence(int len) {
        StringBuilder sb = new StringBuilder(len);
        String pattern = "ACGT";
        for (int i = 0; i < len; i++)
            sb.append(pattern.charAt(i % 4));
        return sb.toString();
    }

    private String makeQuality(char c, int len) {
        StringBuilder sb = new StringBuilder(len);
        for (int i = 0; i < len; i++)
            sb.append(c);
        return sb.toString();
    }

    @Test
    public void testSlidingWindowNoRecovery() {
        // Sensitivity Test: Ensure Sliding Window cuts at the FIRST drop and doesn't keep subsequent good data.
        // Window 4, Threshold 20.
        SlidingWindowTrimmer trimmer = new SlidingWindowTrimmer("4:20");

        String high = makeQuality('?', 20); // Q30
        String low = makeQuality('#', 4);   // Q2
        String highAgain = makeQuality('?', 20); // Q30

        String seq = makeSequence(44);
        String qual = high + low + highAgain;

        FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
        FastqRecord res = trimmer.processRecord(rec);

        assertNotNull(res);
        // Should cut at 20. The subsequent high quality bases are lost.
        assertEquals(20, res.getSequence().length());
    }

    @Test
    public void testLeadingHighQualityNs() {
        // Specificity Test: LEADING is strictly quality-based.
        // It should NOT remove bases if they have high quality.
        // Note: Trimmomatic treats 'N' as quality 0 by default, so we use 'A'.
        
        LeadingTrimmer trimmer = new LeadingTrimmer(20);
        
        // 5 As with Q30, followed by normal seq
        String seq = "AAAAA" + makeSequence(20);
        String qual = makeQuality('?', 25); // All Q30
        
        FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
        FastqRecord res = trimmer.processRecord(rec);
        
        // Should keep the Ns because their quality is high
        assertEquals(25, res.getSequence().length());
        assertEquals(seq, res.getSequence());
    }
    
    @Test
    public void testIlluminaClipSensitivityMaxMismatches() {
        // Sensitivity Test: Verify adapter detection at the limit of allowed mismatches.
        Logger logger = mock(Logger.class);
        // seedMaxMiss = 2.
        IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
        
        String adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"; // 33bp
        trimmer.addClippingSeq(trimmer.new IlluminaLongClippingSeq(adapter), true, false);
        
        // Create read with adapter having exactly 2 mismatches in the seed region (first 16bp).
        // Adapter: AGATCGGAAGAGCACA...
        
        char[] mutAdapterArr = adapter.toCharArray();
        mutAdapterArr[0] = 'C'; // mismatch 1
        mutAdapterArr[1] = 'C'; // mismatch 2
        String mutAdapter = new String(mutAdapterArr);
        
        String dna = "TTCCTTGCCGACGGGCGGTGTGTACAAAGGGCAGGGACTTAATCAACGCA"; // 50bp
        String seq = dna + mutAdapter;
        
        FastqRecord rec = new FastqRecord("head", seq, "", makeQuality('?', seq.length()), 33);
        FastqRecord[] res = trimmer.processRecords(new FastqRecord[]{rec});
        
        // Should clip (2 mismatches <= 2 allowed)
        assertEquals(50, res[0].getSequence().length());
        
        // Now try with a sequence that has no chance of matching the adapter, to test specificity.
        // The original adapter starts with 'AGATCGGAA...'. A string of 'T's should not match.
        String junkAdapter = makeQuality('T', 33);
        String seq_junk = dna + junkAdapter;
        
        FastqRecord rec_junk = new FastqRecord("head", seq_junk, "", makeQuality('?', seq_junk.length()), 33);
        FastqRecord[] res_junk = trimmer.processRecords(new FastqRecord[]{rec_junk});
        
        // Should NOT clip (seed check fails)
        assertEquals(seq_junk.length(), res_junk[0].getSequence().length());
    }
}