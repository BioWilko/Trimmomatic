package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerCaseTest {

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
	public void testUppercaseAdapterLowercaseRead() {
		Logger logger = mock(Logger.class);
		// seedMaxMiss = 2
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
		
		String adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"; // Uppercase
		trimmer.addClippingSeq(trimmer.new IlluminaLongClippingSeq(adapter), true, false);
		
		// Read: 100bp DNA + Adapter + Rest
		String dna = makeSequence(100);
		String rest = makeSequence(17);
		String seqUpper = dna + adapter + rest;
		
		// Lowercase version of the read
		String seqLower = seqUpper.toLowerCase();
		
		FastqRecord recLower = new FastqRecord("head", seqLower, "", makeQuality('I', 150), 33);
		
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { recLower });
		
		assertNotNull(result);
		// Should NOT clip because 'a' (packed as 0) != 'A' (packed as 1).
		// The bitwise XOR will show mismatches everywhere.
		assertEquals(150, result[0].getSequence().length());
		
		// Verify Uppercase works (Control)
		FastqRecord recUpper = new FastqRecord("head", seqUpper, "", makeQuality('I', 150), 33);
		FastqRecord[] resultUpper = trimmer.processRecords(new FastqRecord[] { recUpper });
		
		assertNotNull(resultUpper);
		// Should clip
		assertEquals(100, resultUpper[0].getSequence().length());
	}
}