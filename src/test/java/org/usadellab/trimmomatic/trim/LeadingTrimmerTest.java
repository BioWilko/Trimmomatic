package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class LeadingTrimmerTest {

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
	public void testLeadingNoTrim() {
		// Threshold 20. Read all Q30 ('?').
		LeadingTrimmer trimmer = new LeadingTrimmer(20);
		String seq = makeSequence(150);
		String qual = makeQuality('?', 150);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);

		FastqRecord result = trimmer.processRecord(rec);
		assertNotNull(result);
		assertEquals(150, result.getSequence().length());
		assertEquals(seq, result.getSequence());
	}

	@Test
	public void testLeadingTrimStart() {
		// Threshold 20. First 10 bases Q10 ('+'), rest Q30 ('?').
		LeadingTrimmer trimmer = new LeadingTrimmer(20);
		String seq = makeSequence(150);
		String qual = makeQuality('+', 10) + makeQuality('?', 140);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);

		FastqRecord result = trimmer.processRecord(rec);
		assertNotNull(result);
		assertEquals(140, result.getSequence().length());
		assertEquals(seq.substring(10), result.getSequence());
	}

	@Test
	public void testLeadingTrimAllLow() {
		// Threshold 20. All Q10 ('+'). Should drop.
		LeadingTrimmer trimmer = new LeadingTrimmer(20);
		String seq = makeSequence(150);
		String qual = makeQuality('+', 150);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);

		assertNull(trimmer.processRecord(rec));
	}
}