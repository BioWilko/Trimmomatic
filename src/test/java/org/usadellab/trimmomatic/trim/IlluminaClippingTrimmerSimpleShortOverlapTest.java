package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerSimpleShortOverlapTest {

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
	public void testShortOverlapClipping() {
		Logger logger = mock(Logger.class);
		// minSequenceLikelihood = 4.
		// Log10(4) ~ 0.6.
		// 4 / 0.6 = 6.6.
		// Need ~7 bases match.
		
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 4, 8, false);
		
		String adapter = "AGATCGGAAGAGCACACGTC"; // 20bp non-repetitive
		trimmer.addClippingSeq(trimmer.new IlluminaMediumClippingSeq(adapter), true, false);
		
		// Read: 143bp DNA + 7bp Adapter prefix
		String dna = makeSequence(143);
		String adapterPart = adapter.substring(0, 7); // 7bp
		String seq = dna + adapterPart;
		
		FastqRecord rec = new FastqRecord("head", seq, "", makeQuality('I', 150), 33);
		
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });
		
		assertNotNull(result);
		// Should clip the 7bp adapter part
		assertEquals(143, result[0].getSequence().length());
	}
}