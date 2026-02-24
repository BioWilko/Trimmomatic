package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerMergeTest {

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
	public void testMergeLowQualityMismatch() {
		// This test targets the calculateMaximumRange logic.
		// We want a sequence of matches interrupted by a low-quality mismatch.
		// If merging works, the score should be high enough to trigger clipping.
		
		Logger logger = mock(Logger.class);
		// minSequenceLikelihood = 10. 
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
		
		String adapter = "AGATCGGAAGAGCACACGTC"; // 20bp non-repetitive
		trimmer.addClippingSeq(trimmer.new IlluminaMediumClippingSeq(adapter), true, false);
		
		// Read: 100bp DNA + Adapter + Rest
		// Adapter in read has a mismatch in the middle.
		// Adapter: AGATCGGAAG A GCACACGTC
		// Mutated: AGATCGGAAG T GCACACGTC
		
		String dna = makeSequence(100);
		String adapterMut = "AGATCGGAAGTGCACACGTC"; 
		String rest = makeSequence(30);
		
		String seq = dna + adapterMut + rest;
		
		// Quality: High everywhere except at the mismatch.
		// Mismatch at index 100 + 10 = 110.
		String qualPrefix = makeQuality('I', 110);
		String qualMismatch = "#"; // Phred 35 (#) - 33 = 2. Penalty = 0.2.
		// Match score 0.6. 0.6 > 0.2, so it should merge.
		String qualSuffix = makeQuality('I', 39);
		
		String qual = qualPrefix + qualMismatch + qualSuffix;
		
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });
		
		assertNotNull(result);
		// Should clip at 100 because the low quality mismatch is bridged
		assertEquals(100, result[0].getSequence().length());
	}
	
	@Test
	public void testNoMergeHighQualityMismatch() {
		// Same setup, but high quality mismatch.
		// Penalty = 4.0 (Q40). 
		// Match 0.6. 0.6 < 4.0. Should NOT merge.
		// The region will be split into two matches of ~10bp.
		// 10 * 0.6 = 6.
		// Threshold is 10.
		// Neither part is sufficient. Should NOT clip.
		
		Logger logger = mock(Logger.class);
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
		
		String adapter = "AGATCGGAAGAGCACACGTC"; 
		trimmer.addClippingSeq(trimmer.new IlluminaMediumClippingSeq(adapter), true, false);
		
		String dna = makeSequence(100);
		String adapterMut = "AGATCGGAAGTGCACACGTC"; 
		String rest = makeSequence(30);
		
		String seq = dna + adapterMut + rest;
		
		// High quality everywhere
		String qual = makeQuality('I', 150);
		
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });
		
		assertNotNull(result);
		// Should NOT clip because the mismatch breaks the alignment score
		assertEquals(150, result[0].getSequence().length());
	}
}