package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class MinLenTrimmerTest {

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
	public void testMinLenKeepLong() {
		// Min 50. Read 150.
		MinLenTrimmer trimmer = new MinLenTrimmer(50);
		String seq = makeSequence(150);
		String qual = makeQuality('?', 150);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);

		assertNotNull(trimmer.processRecord(rec));
	}

	@Test
	public void testMinLenKeepExact() {
		// Min 150. Read 150.
		MinLenTrimmer trimmer = new MinLenTrimmer(150);
		String seq = makeSequence(150);
		String qual = makeQuality('?', 150);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);

		assertNotNull(trimmer.processRecord(rec));
	}

	@Test
	public void testMinLenDropShort() {
		// Min 150. Read 149.
		MinLenTrimmer trimmer = new MinLenTrimmer(150);
		String seq = makeSequence(149);
		String qual = makeQuality('?', 149);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);

		assertNull(trimmer.processRecord(rec));
	}

	@Test
	public void testMinLenDropEmpty() {
		// Min 1. Read 0.
		MinLenTrimmer trimmer = new MinLenTrimmer(1);
		FastqRecord rec = new FastqRecord("head", "", "", "", 33);

		assertNull(trimmer.processRecord(rec));
	}
}