package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class HeadCropTrimmerTest {

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
	public void testHeadCropNormal() {
		// Remove 10 bases from start
		HeadCropTrimmer trimmer = new HeadCropTrimmer("10");
		
		String seq = makeSequence(150);
		String qual = makeQuality('!', 150);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		
		// HeadCropTrimmer implements Trimmer directly, so it takes an array
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });
		
		assertNotNull(result);
		assertNotNull(result[0]);
		assertEquals(140, result[0].getSequence().length());
		assertEquals(seq.substring(10), result[0].getSequence());
	}

	@Test
	public void testHeadCropEntireRead() {
		// Remove 150 bases from 150bp read -> Should drop
		HeadCropTrimmer trimmer = new HeadCropTrimmer("150");
		
		FastqRecord rec = new FastqRecord("head", makeSequence(150), "", makeQuality('!', 150), 33);
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });
		
		assertNotNull(result);
		assertNull(result[0], "Read should be dropped (null) when cropping entire length");
	}

	@Test
	public void testHeadCropDifferentRecords() {
		// 6 args: firstBases:firstMax:midBases:midMax:lastBases:lastMax
		// Crop 1 from first, 2 from middle, 3 from last. Max length 150 (no effect on 150bp read)
		HeadCropTrimmer trimmer = new HeadCropTrimmer("1:150:2:150:3:150");
		
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
	public void testHeadCropZero() {
		// Crop 0 bases
		HeadCropTrimmer trimmer = new HeadCropTrimmer("0");
		
		FastqRecord rec = new FastqRecord("head", makeSequence(150), "", makeQuality('!', 150), 33);
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });
		
		assertNotNull(result);
		assertEquals(150, result[0].getSequence().length());
		assertEquals(rec, result[0]); // Should return original object optimization
	}

	@Test
	public void testHeadCropNegative() {
		HeadCropTrimmer trimmer = new HeadCropTrimmer("-5");
		FastqRecord rec = new FastqRecord("head", makeSequence(150), "", makeQuality('!', 150), 33);
		assertThrows(RuntimeException.class, () -> {
			trimmer.processRecords(new FastqRecord[] { rec });
		});
	}

	@Test
	public void testHeadCropWithMaxLength() {
		// 2 args: bases:maxLength
		// Crop 10 bases. Max length 100.
		// 150bp read.
		// Logic: toTrim = 10.
		// overLen = 150 - 10 - 100 = 40.
		// toTrim += 40 -> 50.
		// Result length = 150 - 50 = 100.
		// It crops from head, so we expect the suffix.
		
		HeadCropTrimmer trimmer = new HeadCropTrimmer("10:100");
		String seq = makeSequence(150);
		FastqRecord rec = new FastqRecord("head", seq, "", makeQuality('!', 150), 33);
		
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });
		
		assertNotNull(result);
		assertEquals(100, result[0].getSequence().length());
		assertEquals(seq.substring(50), result[0].getSequence());
	}

	@Test
	public void testHeadCropPaired() {
		// 4 args: fwdBases:fwdMax:revBases:revMax
		// Fwd: Crop 10, Max 150 (no extra crop)
		// Rev: Crop 20, Max 150 (no extra crop)
		HeadCropTrimmer trimmer = new HeadCropTrimmer("10:150:20:150");
		
		String seq = makeSequence(150);
		FastqRecord r1 = new FastqRecord("r1", seq, "", makeQuality('!', 150), 33);
		FastqRecord r2 = new FastqRecord("r2", seq, "", makeQuality('!', 150), 33);
		
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { r1, r2 });
		
		assertEquals(140, result[0].getSequence().length());
		assertEquals(seq.substring(10), result[0].getSequence());
		
		assertEquals(130, result[1].getSequence().length());
		assertEquals(seq.substring(20), result[1].getSequence());
	}

	@Test
	public void testHeadCropExactLength() {
		// Crop 150 bases from 150bp read -> Should drop
		HeadCropTrimmer trimmer = new HeadCropTrimmer("150");
		
		FastqRecord rec = new FastqRecord("head", makeSequence(150), "", makeQuality('!', 150), 33);
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });
		
		assertNull(result[0]);
	}
}