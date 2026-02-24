package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerSeedMismatchTest {

	private String makeSequence(int len) {
		StringBuilder sb = new StringBuilder(len);
		String pattern = "ACGT";
		for (int i = 0; i < len; i++)
			sb.append(pattern.charAt(i % 4));
		return sb.toString();
	}

	@Test
	public void testSeedMismatchBlocking() {
		Logger logger = mock(Logger.class);
		// seedMaxMiss = 2.
		// This means max 2 mismatches allowed in the 16bp seed.
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
		
		// Use 30bp non-repetitive adapter.
		// Seeds (stride 4) at 0, 4, 8, 12.
		// We need to block ALL seeds with > 2 mismatches.
		// Mismatches at 12, 13, 14 will be present in all seeds:
		// Seed 0 (0-15): contains 12,13,14
		// Seed 4 (4-19): contains 12,13,14
		// Seed 8 (8-23): contains 12,13,14
		// Seed 12 (12-27): contains 12,13,14
		String adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAG"; // 30bp
		trimmer.addClippingSeq(trimmer.new IlluminaLongClippingSeq(adapter), true, false);
		
		// Construct read: 100bp DNA + Adapter + Rest
		String dna = makeSequence(100);
		String rest = makeSequence(20);
		
		// Mutate adapter at 12, 13, 14.
		// Original at 12-14: "CAC"
		// Mutate to "GTG"
		String adapterMut = "AGATCGGAAGAGGTGACGTCTGAACTCCAG"; 
		
		String seq = dna + adapterMut + rest;
		
		// Setup qualities
		// High quality everywhere
		char[] quals = new char[seq.length()];
		for(int i=0; i<seq.length(); i++) quals[i] = 'I';
		
		// Set low quality at mismatches (indices 112, 113, 114)
		// Q10 ('+') -> Penalty 1.0
		quals[112] = '+';
		quals[113] = '+';
		quals[114] = '+';
		
		// Score calculation:
		// 30bp length. 3 mismatches. 27 matches.
		// Match score: 27 * 0.6 = 16.2
		// Penalty: 3 * 1.0 = 3.0
		// Total: 13.2 > 10 (threshold).
		// If seed check passes, it would clip.
		// But seed check should fail.
		
		FastqRecord rec = new FastqRecord("head", seq, "", new String(quals), 33);
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });
		
		assertNotNull(result);
		// Should NOT clip because seed check failed (3 mismatches > 2 allowed)
		assertEquals(seq.length(), result[0].getSequence().length());
	}
}