package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerZeroThresholdTest {

	private String makeQuality(char c, int len) {
		StringBuilder sb = new StringBuilder(len);
		for (int i = 0; i < len; i++)
			sb.append(c);
		return sb.toString();
	}

	@Test
	public void testLowThreshold() {
		Logger logger = mock(Logger.class);
		// minSequenceLikelihood = 1. Requires ~2bp match (score 1.2).
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 1, 8, false);
		
		String adapter = "AGATCGGAAGAGCACACGTC"; // 20bp non-repetitive
		trimmer.addClippingSeq(trimmer.new IlluminaMediumClippingSeq(adapter), true, false);
		
		// Case 1: Perfect 2bp match at end.
		// Read ends with 'AG'. Use Ts for the rest to avoid matches.
		String seq1 = makeQuality('T', 148) + "AG";
		FastqRecord rec1 = new FastqRecord("r1", seq1, "", makeQuality('I', 150), 33);
		
		FastqRecord[] res1 = trimmer.processRecords(new FastqRecord[] { rec1 });
		// Score ~1.2 > 1. Should clip 2bp.
		assertEquals(148, res1[0].getSequence().length());
		
		// Case 2: Mismatch.
		// Read ends with 'TT'. Adapter is 'AG'.
		// No matches anywhere (T vs A/G). Score 0 < 1. Should NOT clip.
		String seq2 = makeQuality('T', 150);
		FastqRecord rec2 = new FastqRecord("r2", seq2, "", makeQuality('I', 150), 33);
		
		FastqRecord[] res2 = trimmer.processRecords(new FastqRecord[] { rec2 });
		assertEquals(150, res2[0].getSequence().length());
	}
}