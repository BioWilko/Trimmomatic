package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerWildcardTest {

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
	public void testNIsNotWildcard() {
		Logger logger = mock(Logger.class);
		// seedMaxMiss = 0 (Strict seed check)
		// minSequenceLikelihood = 30.
		// Adapter 50bp.
		// If N is neutral (0): 49 * 0.602 = 29.5 < 30. No Clip.
		// If N is match (0.6): 50 * 0.602 = 30.1 > 30. Clip.
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 0, 30, 30, 8, false);
		
		// Adapter with 'N'. 
		// If N were a wildcard, it would match 'A'.
		// If N is a mismatch, it will fail the seed check (since seedMaxMiss=0).
		StringBuilder sb = new StringBuilder();
		sb.append('N');
		for(int i=0; i<49; i++) sb.append('A');
		String adapter = sb.toString(); // 50bp
		
		trimmer.addClippingSeq(trimmer.new IlluminaLongClippingSeq(adapter), true, false);
		
		// Read matches the adapter, but has 'A' where adapter has 'N'.
		String seq = adapter.replace('N', 'A') + makeSequence(100);
		FastqRecord rec = new FastqRecord("head", seq, "", makeQuality('I', 150), 33);
		
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });
		
		assertNotNull(result);
		assertEquals(150, result[0].getSequence().length());
		
		// Now allow 1 mismatch
		trimmer = new IlluminaClippingTrimmer(logger, 1, 30, 10, 8, false);
		trimmer.addClippingSeq(trimmer.new IlluminaLongClippingSeq(adapter), true, false);
		
		FastqRecord[] result2 = trimmer.processRecords(new FastqRecord[] { rec });
		
		// Should clip now because 1 mismatch is allowed.
		assertNull(result2[0]);
	}
}