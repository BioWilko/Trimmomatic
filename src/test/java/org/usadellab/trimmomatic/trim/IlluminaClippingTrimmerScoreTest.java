package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerScoreTest {

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
	public void testLikelihoodThresholds() {
		Logger logger = mock(Logger.class);
		// minSequenceLikelihood = 10
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
		
		// Adapter: 16bp. 
		// Match score: 16 * 0.602 = 9.632.
		// Threshold is 10. Should NOT clip.
		
		String adapter16 = "AGATCGGAAGAGCACA";
		trimmer.addClippingSeq(trimmer.new IlluminaMediumClippingSeq(adapter16), true, false);
		
		String seq = makeSequence(100) + adapter16 + makeSequence(34); // 150bp
		String qual = makeQuality('I', 150);
		
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[]{rec});
		
		// Should NOT clip because 9.632 < 10
		assertEquals(150, result[0].getSequence().length());
		
		// Now try with 17bp adapter -> 17 * 0.602 = 10.234 score > 10. Should clip.
		String adapter17 = "AGATCGGAAGAGCACAC";
		trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
		trimmer.addClippingSeq(trimmer.new IlluminaMediumClippingSeq(adapter17), true, false);
		
		String seq2 = makeSequence(100) + adapter17 + makeSequence(33);
		FastqRecord rec2 = new FastqRecord("head", seq2, "", qual, 33);
		result = trimmer.processRecords(new FastqRecord[]{rec2});
		
		assertEquals(100, result[0].getSequence().length());
	}
	
	@Test
	public void testQualityImpactOnLikelihood() {
		Logger logger = mock(Logger.class);
		// minSequenceLikelihood = 7
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 7, 8, false);
		
		// Adapter 15bp. Score 15 * 0.602 = 9.03.
		// 1 mismatch reduces score by Q/10.
		
		// We need a mismatch in the middle to split the local alignment.
		// 7bp match (score 4.2) - Mismatch - 7bp match (score 4.2).
		// Total with low penalty (1): 4.2 - 1 + 4.2 = 7.4 > 7. Clip.
		// Total with high penalty (4): 4.2 - 4 + 4.2 = 4.4 < 7. No clip.
		String adapter = "AAAAAAACAAAAAAA"; // 15bp
		trimmer.addClippingSeq(trimmer.new IlluminaShortClippingSeq(adapter), true, false);
		
		String dna = makeSequence(100);
		String rest = makeSequence(35);
		String adapterMut = "AAAAAAATAAAAAAA"; // Mismatch at index 7: C -> T
		String seq = dna + adapterMut + rest;
		
		// Case 1: High quality mismatch (Q40). Penalty = 4. Score = 9.03 - 4 = 5.03 < 7. No clip.
		FastqRecord recHigh = new FastqRecord("head", seq, "", makeQuality('I', 150), 33);
		assertEquals(150, trimmer.processRecords(new FastqRecord[]{recHigh})[0].getSequence().length());
		
		// Case 2: Low quality mismatch (Q10). Penalty = 1. Score = 9.03 - 1 = 8.03 > 7. Clip.
		String qualLow = makeQuality('I', 107) + "+" + makeQuality('I', 42); // Q10 at mismatch
		FastqRecord recLow = new FastqRecord("head", seq, "", qualLow, 33);
		assertEquals(100, trimmer.processRecords(new FastqRecord[]{recLow})[0].getSequence().length());
	}
}