package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerOverlapCapTest {

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
	public void testOverlapCapStrictness() {
		Logger logger = mock(Logger.class);
		// Set likelihood to 15.
		// Required overlap approx: 15 / 0.602 = 24.9 -> 25 bases.
		// The internal minSequenceOverlap is capped at 15.
		// This means the trimmer will SCAN for overlaps >= 15.
		// But it should REJECT overlaps that don't meet the score of 15 (i.e., < 25 bases).
		
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 15, 8, false);
		
		String adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTC"; // 32bp
		trimmer.addClippingSeq(trimmer.new IlluminaLongClippingSeq(adapter), true, false);
		
		// Case 1: 20bp overlap.
		// 20 > 15 (cap), so it is scanned.
		// Score: 20 * 0.602 = 12.04.
		// 12.04 < 15 (threshold). Should NOT clip.
		String seq20 = makeSequence(100) + adapter.substring(0, 20) + makeSequence(30);
		FastqRecord rec20 = new FastqRecord("r1", seq20, "", makeQuality('I', 150), 33);
		
		FastqRecord[] res20 = trimmer.processRecords(new FastqRecord[] { rec20 });
		assertEquals(150, res20[0].getSequence().length());
		
		// Case 2: 26bp overlap.
		// 26 > 15. Scanned.
		// Score: 26 * 0.602 = 15.65.
		// 15.65 > 15. Should clip.
		String seq26 = makeSequence(100) + adapter.substring(0, 26) + makeSequence(24);
		FastqRecord rec26 = new FastqRecord("r2", seq26, "", makeQuality('I', 150), 33);
		
		FastqRecord[] res26 = trimmer.processRecords(new FastqRecord[] { rec26 });
		assertEquals(100, res26[0].getSequence().length());
	}
}