package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class MaxLenTrimmerTest {

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
	public void testMaxLenKeepShort() {
		// Max 150. Read 100.
		MaxLenTrimmer trimmer = new MaxLenTrimmer(150);
		String seq = makeSequence(100);
		String qual = makeQuality('?', 100);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);

		assertNotNull(trimmer.processRecord(rec));
	}

	@Test
	public void testMaxLenKeepExact() {
		// Max 150. Read 150.
		MaxLenTrimmer trimmer = new MaxLenTrimmer(150);
		String seq = makeSequence(150);
		String qual = makeQuality('?', 150);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);

		assertNotNull(trimmer.processRecord(rec));
	}

	@Test
	public void testMaxLenDropLong() {
		// Max 150. Read 151.
		MaxLenTrimmer trimmer = new MaxLenTrimmer(150);
		String seq = makeSequence(151);
		String qual = makeQuality('?', 151);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);

		assertNull(trimmer.processRecord(rec));
	}
}