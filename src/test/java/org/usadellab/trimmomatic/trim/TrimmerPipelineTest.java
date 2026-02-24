package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class TrimmerPipelineTest {

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
	public void testPipeline() {
		// 150bp read
		String seq = makeSequence(150);
		String qual = makeQuality('I', 150);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		
		// 1. HeadCrop 10
		HeadCropTrimmer headCrop = new HeadCropTrimmer("10");
		// 2. TailCrop 10
		TailCropTrimmer tailCrop = new TailCropTrimmer("10");
		// 3. Crop 100
		CropTrimmer crop = new CropTrimmer("100");
		
		// Execute pipeline manually
		FastqRecord[] res1 = headCrop.processRecords(new FastqRecord[] { rec });
		assertNotNull(res1[0]);
		assertEquals(140, res1[0].getSequence().length());
		
		FastqRecord[] res2 = tailCrop.processRecords(res1);
		assertNotNull(res2[0]);
		assertEquals(130, res2[0].getSequence().length());
		
		FastqRecord[] res3 = crop.processRecords(res2);
		assertNotNull(res3[0]);
		assertEquals(100, res3[0].getSequence().length());
		
		// Verify content: Original 10..109
		String expected = seq.substring(10, 110);
		assertEquals(expected, res3[0].getSequence());
	}

	@Test
	public void testPipelineDrop() {
		// 150bp read
		String seq = makeSequence(150);
		String qual = makeQuality('!', 150); // Low quality
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		
		// 1. AvgQual 20 (should drop)
		AvgQualTrimmer avgQual = new AvgQualTrimmer(20);
		// 2. Crop 100
		CropTrimmer crop = new CropTrimmer(100);
		
		FastqRecord[] res1 = avgQual.processRecords(new FastqRecord[] { rec });
		assertNull(res1[0]);
		
		FastqRecord[] res2 = crop.processRecords(res1);
		assertNull(res2[0]);
	}

	@Test
	public void testPipelineWithAdapter() {
		// 1. Adapter Clip -> 2. Crop
		Logger logger = mock(Logger.class);
		IlluminaClippingTrimmer clipper = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
		String adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
		clipper.addClippingSeq(clipper.new IlluminaLongClippingSeq(adapter), true, false);
		
		CropTrimmer crop = new CropTrimmer(50);
		
		// Read: 100bp DNA + Adapter + Rest
		String dna = makeSequence(100);
		String rest = makeSequence(17);
		String seq = dna + adapter + rest; // 100 + 33 + 17 = 150
		String qual = makeQuality('I', 150);
		
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		
		// Step 1: Clip adapter. Should result in 100bp.
		FastqRecord[] res1 = clipper.processRecords(new FastqRecord[] { rec });
		assertNotNull(res1[0]);
		assertEquals(100, res1[0].getSequence().length());
		
		// Step 2: Crop to 50.
		FastqRecord[] res2 = crop.processRecords(res1);
		assertEquals(50, res2[0].getSequence().length());
	}
}