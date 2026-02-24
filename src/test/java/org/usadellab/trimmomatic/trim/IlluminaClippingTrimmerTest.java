package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerTest {

	private IlluminaClippingTrimmer trimmer;
	private Logger logger;

	@BeforeEach
	public void setUp() {
		logger = mock(Logger.class);
		// seedMaxMiss=2, minPalindromeLikelihood=30, minSequenceLikelihood=10, minPrefix=8, palindromeKeepBoth=false
		trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
	}

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

	private String reverseComplement(String seq) {
		StringBuilder sb = new StringBuilder(seq.length());
		for (int i = seq.length() - 1; i >= 0; i--) {
			char c = seq.charAt(i);
			switch (c) {
			case 'A':
				sb.append('T');
				break;
			case 'C':
				sb.append('G');
				break;
			case 'G':
				sb.append('C');
				break;
			case 'T':
				sb.append('A');
				break;
			default:
				sb.append(c);
				break;
			}
		}
		return sb.toString();
	}

	@Test
	public void testSimpleClipAdapterAtEnd() {
		// Common Illumina Adapter sequence (33bp)
		String adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
		
		// Add as simple clipping sequence (forward only for this test)
		// We use the inner class IlluminaLongClippingSeq directly
		trimmer.addClippingSeq(trimmer.new IlluminaLongClippingSeq(adapter), true, false);

		// 150bp read: 100bp DNA + Adapter + rest
		String dna = makeSequence(100);
		String rest = makeSequence(17); // 100 + 33 + 17 = 150
		String seq = dna + adapter + rest; 
		String qual = makeString('I', seq.length()); // High quality (Q40)

		FastqRecord rec = new FastqRecord("header", seq, "", qual, 33);
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });

		assertNotNull(result);
		assertEquals(1, result.length);
		// Should be trimmed at 100, removing adapter and everything after
		assertEquals(dna, result[0].getSequence());
	}

	@Test
	public void testSimpleClipAdapterInMiddle() {
		String adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
		trimmer.addClippingSeq(trimmer.new IlluminaLongClippingSeq(adapter), true, false);

		// 50bp DNA + Adapter + rest
		String dna = makeSequence(50);
		String rest = makeSequence(67); // 50 + 33 + 67 = 150
		String seq = dna + adapter + rest;
		String qual = makeString('I', seq.length());

		FastqRecord rec = new FastqRecord("header", seq, "", qual, 33);
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });

		assertNotNull(result);
		assertEquals(dna, result[0].getSequence());
	}

	@Test
	public void testSimpleClipAdapterAtStart() {
		String adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
		trimmer.addClippingSeq(trimmer.new IlluminaLongClippingSeq(adapter), true, false);

		// Adapter + rest
		String rest = makeSequence(117); // 33 + 117 = 150
		String seq = adapter + rest;
		String qual = makeString('I', seq.length());

		FastqRecord rec = new FastqRecord("header", seq, "", qual, 33);
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });

		// If the adapter is at the very start, the entire read is clipped.
		// processRecords returns an array containing null for dropped records 
		// (or filters them depending on implementation, but here it returns null in the slot)
		if (result.length > 0) {
			assertNull(result[0], "Record should be null (dropped) if adapter is at start");
		} else {
			assertEquals(0, result.length);
		}
	}

	@Test
	public void testPalindromeClip() {
		String fwdAdapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
		String revAdapter = "TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"; // RC of fwdAdapter
		
		trimmer.addPrefixPair(fwdAdapter, revAdapter);
		
		// 50bp fragment (shorter than 150bp read)
		String fragment = makeSequence(50);
		String rcFragment = reverseComplement(fragment);
		
		// Read 1: Fragment + FwdAdapter + Garbage
		String seq1 = fragment + fwdAdapter + makeSequence(150 - 50 - fwdAdapter.length());
		// Read 2: RCFragment + RevAdapter + Garbage
		String seq2 = rcFragment + revAdapter + makeSequence(150 - 50 - revAdapter.length());
		
		String qual = makeString('I', 150);
		
		FastqRecord rec1 = new FastqRecord("r1", seq1, "", qual, 33);
		FastqRecord rec2 = new FastqRecord("r2", seq2, "", qual, 33);
		
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec1, rec2 });
		
		assertNotNull(result);
		assertEquals(2, result.length);
		
		// Forward read clipped to fragment
		assertEquals(50, result[0].getSequence().length());
		assertEquals(fragment, result[0].getSequence());
		
		// Reverse read dropped (null) because palindromeKeepBoth is false by default
		assertNull(result[1]);
	}

	@Test
	public void testPalindromeClipKeepBoth() {
		// Re-init trimmer with keepBoth = true
		trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, true);
		
		String fwdAdapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
		String revAdapter = "TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"; // RC of fwdAdapter
		
		trimmer.addPrefixPair(fwdAdapter, revAdapter);
		
		String fragment = makeSequence(50);
		String rcFragment = reverseComplement(fragment);
		
		String seq1 = fragment + fwdAdapter + makeSequence(150 - 50 - fwdAdapter.length());
		String seq2 = rcFragment + revAdapter + makeSequence(150 - 50 - revAdapter.length());
		
		String qual = makeString('I', 150);
		
		FastqRecord rec1 = new FastqRecord("r1", seq1, "", qual, 33);
		FastqRecord rec2 = new FastqRecord("r2", seq2, "", qual, 33);
		
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec1, rec2 });
		
		assertNotNull(result);
		assertEquals(2, result.length);
		
		// Forward read clipped
		assertEquals(50, result[0].getSequence().length());
		assertEquals(fragment, result[0].getSequence());
		
		// Reverse read also kept and clipped
		assertNotNull(result[1]);
		assertEquals(50, result[1].getSequence().length());
		assertEquals(rcFragment, result[1].getSequence());
	}

	@Test
	public void testShortAdapterClip() {
		// Re-init trimmer with lower threshold (5) because 10bp adapter max score is ~6.02 (10 * 0.602)
		// Default threshold of 10 is too high for a 10bp adapter.
		trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 5, 8, false);
		
		// 10bp adapter (< 16) -> Uses IlluminaShortClippingSeq logic
		String adapter = "AGATCGGAAG";
		
		trimmer.addClippingSeq(trimmer.new IlluminaShortClippingSeq(adapter), true, false);

		// 150bp read: 100bp DNA + Adapter + rest
		String dna = makeSequence(100);
		String rest = makeSequence(40); 
		String seq = dna + adapter + rest; 
		String qual = makeString('I', seq.length());

		FastqRecord rec = new FastqRecord("header", seq, "", qual, 33);
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });

		assertNotNull(result);
		assertEquals(1, result.length);
		assertEquals(dna, result[0].getSequence());
	}

	@Test
	public void testMediumAdapterClip() {
		// 20bp adapter (16 <= len < 24) -> Uses IlluminaMediumClippingSeq logic
		String adapter = "AGATCGGAAGAGCACACGTC";
		
		trimmer.addClippingSeq(trimmer.new IlluminaMediumClippingSeq(adapter), true, false);

		// 150bp read: 100bp DNA + Adapter + rest
		String dna = makeSequence(100);
		String rest = makeSequence(30); 
		String seq = dna + adapter + rest; 
		String qual = makeString('I', seq.length());

		FastqRecord rec = new FastqRecord("header", seq, "", qual, 33);
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });

		assertNotNull(result);
		assertEquals(1, result.length);
		assertEquals(dna, result[0].getSequence());
	}

	@Test
	public void testNoClipWhenAdapterAbsent() {
		String adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
		trimmer.addClippingSeq(trimmer.new IlluminaLongClippingSeq(adapter), true, false);

		// 150bp read with NO adapter
		String seq = makeSequence(150);
		String qual = makeString('I', 150);

		FastqRecord rec = new FastqRecord("header", seq, "", qual, 33);
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });

		assertNotNull(result);
		assertEquals(1, result.length);
		assertEquals(150, result[0].getSequence().length());
		assertEquals(seq, result[0].getSequence());
	}

	@Test
	public void testExactAdapterMatch() {
		String adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
		trimmer.addClippingSeq(trimmer.new IlluminaLongClippingSeq(adapter), true, false);

		String seq = adapter;
		String qual = makeString('I', seq.length());

		FastqRecord rec = new FastqRecord("header", seq, "", qual, 33);
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });

		// Should be dropped (null) because the entire read is the adapter
		if (result.length > 0) {
			assertNull(result[0]);
		} else {
			assertEquals(0, result.length);
		}
	}

	@Test
	public void testProcessArrayTooLarge() {
		// IlluminaClippingTrimmer expects 1 or 2 records.
		// If 3 are passed, the current implementation returns an empty array.
		FastqRecord rec = new FastqRecord("head", makeSequence(150), "", makeString('I', 150), 33);
		FastqRecord[] input = new FastqRecord[] { rec, rec, rec };
		
		FastqRecord[] result = trimmer.processRecords(input);
		assertNotNull(result);
		assertEquals(0, result.length);
	}
}