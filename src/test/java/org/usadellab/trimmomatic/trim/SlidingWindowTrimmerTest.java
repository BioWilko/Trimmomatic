package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class SlidingWindowTrimmerTest {

	private String makeString(char c, int len) {
		StringBuilder sb = new StringBuilder(len);
		for (int i = 0; i < len; i++)
			sb.append(c);
		return sb.toString();
	}

	private String makeSequence(int len) {
		StringBuilder sb = new StringBuilder(len);
		String pattern = "ACGT";
		for (int i = 0; i < len; i++)
			sb.append(pattern.charAt(i % 4));
		return sb.toString();
	}

	@Test
	public void testSlidingWindowCut() {
		// Window size 4, required average quality 20.
		SlidingWindowTrimmer trimmer = new SlidingWindowTrimmer("4:20");

		// 150bp read
		// 0-99: High quality (Q30)
		// 100-103: Low quality (Q10) -> Window average drops here
		// 104-149: High quality
		
		String highQual = makeString('?', 100); // '?' = Q30
		String lowQual = makeString('+', 4);    // '+' = Q10
		String rest = makeString('?', 46);
		
		String seq = makeSequence(150);
		String qual = highQual + lowQual + rest;
		
		FastqRecord rec = new FastqRecord("header", seq, "", qual, 33);

		FastqRecord processed = trimmer.processRecord(rec);

		assertNotNull(processed);
		// Should cut at index 100
		assertEquals(100, processed.getSequence().length());
		assertEquals(highQual, processed.getQuality());
	}

	@Test
	public void testSlidingWindowCutWithBacktrack() {
		// Window size 4, required average quality 20.
		SlidingWindowTrimmer trimmer = new SlidingWindowTrimmer("4:20");

		// 150bp read
		// 0-89: High quality (Q30)
		// 90-99: Borderline/Low individual bases but window average OK until end
		// Let's construct a case where the window fails at 100, but base 99 is also low.
		
		// 0-98: Q30
		// 99: Q10 (Low, but previous window [96,97,98,99] avg is (30+30+30+10)/4 = 25 >= 20. OK)
		// 100-103: Q10 (Window [100,101,102,103] avg 10 < 20. FAIL)
		
		// The trimmer decides to cut at 100.
		// Then it checks base 99. It is Q10 (< 20). It drops it.
		// It checks base 98. It is Q30 (> 20). It stops.
		// Result length should be 99.
		
		String part1 = makeString('?', 99);
		String part2 = "+"; // Index 99, Q10
		String part3 = makeString('+', 50); // The rest is low
		
		String seq = makeSequence(150);
		String qual = part1 + part2 + part3;
		
		FastqRecord rec = new FastqRecord("header", seq, "", qual, 33);
		FastqRecord processed = trimmer.processRecord(rec);

		assertNotNull(processed);
		assertEquals(99, processed.getSequence().length());
	}

	@Test
	public void testSlidingWindowDropAll() {
		SlidingWindowTrimmer trimmer = new SlidingWindowTrimmer("4:20");

		// 150bp Low quality
		String seq = makeSequence(150);
		String qual = makeString('+', 150); // Q10
		FastqRecord rec = new FastqRecord("header", seq, "", qual, 33);

		FastqRecord processed = trimmer.processRecord(rec);

		assertNull(processed);
	}
	
	@Test
	public void testSlidingWindowKeepAll() {
		SlidingWindowTrimmer trimmer = new SlidingWindowTrimmer("4:20");

		// 150bp High quality
		String seq = makeSequence(150);
		String qual = makeString('?', 150); // Q30
		FastqRecord rec = new FastqRecord("header", seq, "", qual, 33);

		FastqRecord processed = trimmer.processRecord(rec);

		assertEquals(150, processed.getSequence().length());
	}
}