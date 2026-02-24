package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class MaximumInformationTrimmerLongReadTest {

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
	public void testReadTooLong() {
		// MAXINFO has a hardcoded limit of 1000bp.
		// Testing that it fails as expected for longer reads.
		
		MaximumInformationTrimmer trimmer = new MaximumInformationTrimmer("100:0.5");
		
		// 1001 bases
		String seq = makeSequence(1001);
		String qual = makeQuality('I', 1001);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		
		assertThrows(ArrayIndexOutOfBoundsException.class, () -> {
			trimmer.processRecord(rec);
		});
	}
}