package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class TailCropTrimmerTest {

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
	public void testTailCropNormal() {
		// Remove 10 bases from end
		TailCropTrimmer trimmer = new TailCropTrimmer("10");
		
		String seq = makeSequence(150);
		String qual = makeQuality('!', 150);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });
		
		assertNotNull(result);
		assertNotNull(result[0]);
		assertEquals(140, result[0].getSequence().length());
		// Should keep the first 140 bases
		assertEquals(seq.substring(0, 140), result[0].getSequence());
	}

	@Test
	public void testTailCropEntireRead() {
		// Remove 150 bases from end -> Should drop
		TailCropTrimmer trimmer = new TailCropTrimmer("150");
		
		FastqRecord rec = new FastqRecord("head", makeSequence(150), "", makeQuality('!', 150), 33);
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });
		
		assertNotNull(result);
		assertNull(result[0], "Read should be dropped (null) when cropping entire length");
	}

	@Test
	public void testTailCropDifferentRecords() {
		// 6 args: firstBases:firstMax:midBases:midMax:lastBases:lastMax
		// Crop 1 from first, 2 from middle, 3 from last. Max length 150 (no effect)
		TailCropTrimmer trimmer = new TailCropTrimmer("1:150:2:150:3:150");
		
		String seq = makeSequence(150);
		String qual = makeQuality('!', 150);
		
		FastqRecord r1 = new FastqRecord("r1", seq, "", qual, 33);
		FastqRecord r2 = new FastqRecord("r2", seq, "", qual, 33);
		FastqRecord r3 = new FastqRecord("r3", seq, "", qual, 33);
		
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { r1, r2, r3 });
		
		assertNotNull(result);
		assertEquals(3, result.length);
		
		// First record: crop 1
		assertEquals(149, result[0].getSequence().length());
		// Middle record: crop 2
		assertEquals(148, result[1].getSequence().length());
		// Last record: crop 3
		assertEquals(147, result[2].getSequence().length());
	}

	@Test
	public void testTailCropZero() {
		// Crop 0 bases
		TailCropTrimmer trimmer = new TailCropTrimmer("0");
		
		FastqRecord rec = new FastqRecord("head", makeSequence(150), "", makeQuality('!', 150), 33);
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });
		
		assertNotNull(result);
		assertEquals(150, result[0].getSequence().length());
		assertEquals(rec, result[0]); // Should return original object optimization
	}
}