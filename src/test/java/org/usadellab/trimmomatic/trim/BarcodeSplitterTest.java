package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

import java.util.HashMap;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class BarcodeSplitterTest {

	private String makeSequence(String prefix, int len) {
		StringBuilder sb = new StringBuilder(len);
		sb.append(prefix);
		String pattern = "ACGT";
		for (int i = prefix.length(); i < len; i++)
			sb.append(pattern.charAt(i % 4));
		return sb.toString();
	}
	
	private String makeQuality(int len) {
		StringBuilder sb = new StringBuilder(len);
		for (int i = 0; i < len; i++)
			sb.append('I'); // High quality
		return sb.toString();
	}

	@Test
	public void testExactMatchClip() {
		HashMap<String, String> barcodes = new HashMap<>();
		barcodes.put("BC1", "ACGT");
		
		// 0 mismatches, clip = true
		BarcodeSplitter splitter = new BarcodeSplitter(barcodes, 0, true);
		
		String seq = makeSequence("ACGT", 150); // Starts with ACGT
		FastqRecord rec = new FastqRecord("head", seq, "", makeQuality(150), 33);
		
		FastqRecord result = splitter.processRecord(rec);
		
		assertNotNull(result);
		assertEquals("BC1", result.getBarcodeLabel());
		assertEquals(146, result.getSequence().length()); // 150 - 4
		assertEquals(seq.substring(4), result.getSequence());
	}

	@Test
	public void testExactMatchNoClip() {
		HashMap<String, String> barcodes = new HashMap<>();
		barcodes.put("BC1", "ACGT");
		
		// 0 mismatches, clip = false
		BarcodeSplitter splitter = new BarcodeSplitter(barcodes, 0, false);
		
		String seq = makeSequence("ACGT", 150);
		FastqRecord rec = new FastqRecord("head", seq, "", makeQuality(150), 33);
		
		FastqRecord result = splitter.processRecord(rec);
		
		assertNotNull(result);
		assertEquals("BC1", result.getBarcodeLabel());
		assertEquals(150, result.getSequence().length());
		assertEquals(seq, result.getSequence());
	}
	
	@Test
	public void testMismatchMatch() {
		HashMap<String, String> barcodes = new HashMap<>();
		barcodes.put("BC1", "ACGT");
		
		// 1 mismatch allowed
		BarcodeSplitter splitter = new BarcodeSplitter(barcodes, 1, true);
		
		// ACGA (1 mismatch from ACGT)
		String seq = makeSequence("ACGA", 150); 
		FastqRecord rec = new FastqRecord("head", seq, "", makeQuality(150), 33);
		
		FastqRecord result = splitter.processRecord(rec);
		
		assertNotNull(result);
		assertEquals("BC1", result.getBarcodeLabel());
		assertEquals(146, result.getSequence().length());
	}
	
	@Test
	public void testMismatchFail() {
		HashMap<String, String> barcodes = new HashMap<>();
		barcodes.put("BC1", "ACGT");
		
		// 0 mismatches allowed
		BarcodeSplitter splitter = new BarcodeSplitter(barcodes, 0, true);
		
		// ACGA (1 mismatch)
		String seq = makeSequence("ACGA", 150); 
		FastqRecord rec = new FastqRecord("head", seq, "", makeQuality(150), 33);
		
		FastqRecord result = splitter.processRecord(rec);
		
		assertNotNull(result);
		assertEquals("UNKNOWN", result.getBarcodeLabel());
		assertEquals(150, result.getSequence().length()); // Should not clip if unknown
	}
}