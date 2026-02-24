package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerShortReadTest {

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
	public void testVeryShortRead() {
		Logger logger = mock(Logger.class);
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
		
		String adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
		trimmer.addClippingSeq(trimmer.new IlluminaLongClippingSeq(adapter), true, false);
		
		// Read length 10 (shorter than 16bp packing window)
		String seq = makeSequence(10);
		FastqRecord rec = new FastqRecord("head", seq, "", makeQuality('I', 10), 33);
		
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });
		
		assertNotNull(result);
		assertEquals(10, result[0].getSequence().length());
	}
}