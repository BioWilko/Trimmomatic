package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerTinyAdapterTest {

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
	public void testTinyAdapter() {
		Logger logger = mock(Logger.class);
		// Set minSequenceLikelihood to 0 to allow 1bp match (score ~0.6)
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 0, 8, false);
		
		String adapter = "A"; // 1bp adapter
		trimmer.addClippingSeq(trimmer.new IlluminaShortClippingSeq(adapter), true, false);
		
		// Read starts with 'A'
		String seq = "A" + makeSequence(149);
		FastqRecord rec = new FastqRecord("head", seq, "", makeQuality('I', 150), 33);
		
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });
		
		assertNotNull(result);
		// Should clip the 'A' at the start.
		// Since adapter is at start, it might drop the read or return empty?
		// If adapter is found at offset 0, it clips from 0. Result length 0.
		assertNull(result[0]);
	}
}