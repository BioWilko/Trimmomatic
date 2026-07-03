package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerMinPrefixSimpleTest {

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
	public void testMinPrefixIgnoredForSimpleClip() {
		Logger logger = mock(Logger.class);
		// Set minPrefix to 20 (very high).
		// This should affect Palindrome clipping, but NOT Simple clipping.
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 5, 20, false);
		
		// Short adapter (10bp)
		String adapter = "AGATCGGAAG";
		// Add as simple clipping sequence (forward only)
		trimmer.addClippingSeq(trimmer.new IlluminaShortClippingSeq(adapter), true, false);
		
		// Read: 100bp DNA + Adapter + Rest
		String dna = makeSequence(100);
		String rest = makeSequence(40);
		String seq = dna + adapter + rest;
		
		FastqRecord rec = new FastqRecord("head", seq, "", makeQuality('I', 150), 33);
		
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });
		
		assertNotNull(result);
		// Should clip despite minPrefix=20, because simple clipping uses minSequenceLikelihood/Overlap
		// 10bp match score > 5 (threshold).
		assertEquals(100, result[0].getSequence().length());
		assertEquals(dna, result[0].getSequence());
	}
}