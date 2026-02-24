package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerOverlapTest {

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
	public void testOverlapBoundary() {
		Logger logger = mock(Logger.class);
		// minSequenceLikelihood = 10.
		// Log10(4) ~ 0.60206.
		// 10 / 0.60206 = 16.6.
		// So we need at least 17 bases of perfect match to exceed the threshold.
		// 16 bases * 0.60206 = 9.63 < 10.
		// 17 bases * 0.60206 = 10.23 > 10.
		
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
		
		String adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"; // 33bp
		trimmer.addClippingSeq(trimmer.new IlluminaLongClippingSeq(adapter), true, false);
		
		// Case 1: 16bp overlap (should NOT clip)
		// Let's put 16bp match at the end of the read.
		String read1 = makeSequence(134) + adapter.substring(0, 16);
		FastqRecord rec1 = new FastqRecord("r1", read1, "", makeQuality('I', 150), 33);
		
		FastqRecord[] res1 = trimmer.processRecords(new FastqRecord[] { rec1 });
		assertEquals(150, res1[0].getSequence().length());

		// Case 2: 17bp overlap (should clip)
		String read2 = makeSequence(133) + adapter.substring(0, 17);
		FastqRecord rec2 = new FastqRecord("r2", read2, "", makeQuality('I', 150), 33);
		
		FastqRecord[] res2 = trimmer.processRecords(new FastqRecord[] { rec2 });
		// Should clip the adapter part (17bp)
		assertEquals(133, res2[0].getSequence().length());
	}
}