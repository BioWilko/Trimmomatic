package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class BaseCountTrimmerAdvancedTest {

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
	public void testMultipleBasesCount() {
		// Count 'A' and 'T'. 
		// 150bp read "ACGT..." pattern.
		// A: 38, C: 38, G: 37, T: 37.
		// A+T = 38 + 37 = 75.
		
		// Range 70-80. Should keep.
		BaseCountTrimmer trimmer = new BaseCountTrimmer("AT:70:80");
		FastqRecord rec = new FastqRecord("head", makeSequence(150), "", makeQuality('!', 150), 33);
		assertNotNull(trimmer.processRecord(rec));
		
		// Range 0-70. Should drop.
		trimmer = new BaseCountTrimmer("AT:0:70");
		assertNull(trimmer.processRecord(rec));
	}
	
	@Test
	public void testAllBasesCount() {
		// Count ACGT. Should be 150.
		BaseCountTrimmer trimmer = new BaseCountTrimmer("ACGT:150:150");
		FastqRecord rec = new FastqRecord("head", makeSequence(150), "", makeQuality('!', 150), 33);
		assertNotNull(trimmer.processRecord(rec));
	}
	
	@Test
	public void testFailOnInvalidChars() {
		// If read contains 'N' but we only count ACGT, it shouldn't count N.
		BaseCountTrimmer trimmer = new BaseCountTrimmer("ACGT:150:150");
		String seq = makeSequence(149) + "N";
		FastqRecord rec = new FastqRecord("head", seq, "", makeQuality('!', 150), 33);
		
		// Count is 149. Min is 150. Should drop.
		assertNull(trimmer.processRecord(rec));
	}
}