package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class LeadingTrailingTrimmerEdgeTest {

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
	public void testLeadingZeroVsOne() {
		// 150bp read. Starts with 5 Ns.
		// Ns have quality 0 in Trimmomatic's internal representation (when zeroNs=true).
		String seq = "NNNNN" + makeSequence(145);
		// Quality '!' is 0 in Phred33.
		String qual = makeQuality('!', 5) + makeQuality('I', 145);
		
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);

		// LEADING:0 -> Should keep Ns (0 >= 0)
		LeadingTrimmer trimmer0 = new LeadingTrimmer(0);
		FastqRecord res0 = trimmer0.processRecord(rec);
		assertEquals(150, res0.getSequence().length());

		// LEADING:1 -> Should drop Ns (0 < 1)
		LeadingTrimmer trimmer1 = new LeadingTrimmer(1);
		FastqRecord res1 = trimmer1.processRecord(rec);
		assertEquals(145, res1.getSequence().length());
	}

	@Test
	public void testTrailingZeroVsOne() {
		// 150bp read. Ends with 5 Ns.
		String seq = makeSequence(145) + "NNNNN";
		String qual = makeQuality('I', 145) + makeQuality('!', 5);
		
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);

		// TRAILING:0 -> Should keep Ns
		TrailingTrimmer trimmer0 = new TrailingTrimmer(0);
		FastqRecord res0 = trimmer0.processRecord(rec);
		assertEquals(150, res0.getSequence().length());

		// TRAILING:1 -> Should drop Ns
		TrailingTrimmer trimmer1 = new TrailingTrimmer(1);
		FastqRecord res1 = trimmer1.processRecord(rec);
		assertEquals(145, res1.getSequence().length());
	}
}