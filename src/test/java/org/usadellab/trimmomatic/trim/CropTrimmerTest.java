package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class CropTrimmerTest {

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
	public void testCropLongRead() {
		// Crop to 75 bases (150bp read)
		CropTrimmer trimmer = new CropTrimmer(75);
		
		String seq = makeSequence(150);
		String qual = makeQuality('!', 150);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		FastqRecord result = trimmer.processRecord(rec);
		
		assertNotNull(result);
		assertEquals(75, result.getSequence().length());
		assertEquals(seq.substring(0, 75), result.getSequence());
		assertEquals(qual.substring(0, 75), result.getQuality());
	}

	@Test
	public void testCropShortRead() {
		// Crop to 150 bases (read is shorter, 100bp)
		CropTrimmer trimmer = new CropTrimmer(150);
		
		String seq = makeSequence(100);
		String qual = makeQuality('!', 100);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		FastqRecord result = trimmer.processRecord(rec);
		
		assertNotNull(result);
		assertEquals(100, result.getSequence().length());
		assertEquals(seq, result.getSequence());
	}

	@Test
	public void testCropExactLength() {
		// Crop to 150 bases (read is exactly 150bp)
		CropTrimmer trimmer = new CropTrimmer(150);

		String seq = makeSequence(150);
		String qual = makeQuality('!', 150);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		FastqRecord result = trimmer.processRecord(rec);

		assertNotNull(result);
		assertEquals(150, result.getSequence().length());
		assertEquals(seq, result.getSequence());
	}

	@Test
	public void testCropToZero() {
		// Crop to 0 bases
		CropTrimmer trimmer = new CropTrimmer(0);

		String seq = makeSequence(150);
		String qual = makeQuality('!', 150);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		FastqRecord result = trimmer.processRecord(rec);

		assertNotNull(result);
		assertEquals(0, result.getSequence().length());
	}

	@Test
	public void testCropNegativeLength() {
		// CROP:-10
		CropTrimmer trimmer = new CropTrimmer("-10");
		
		FastqRecord rec = new FastqRecord("head", makeSequence(150), "", makeQuality('!', 150), 33);
		
		assertThrows(StringIndexOutOfBoundsException.class, () -> {
			trimmer.processRecord(rec);
		});
	}
}