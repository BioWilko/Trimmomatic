package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerOvertrimTest {

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

	private String reverseComplement(String seq) {
		StringBuilder sb = new StringBuilder(seq.length());
		for (int i = seq.length() - 1; i >= 0; i--) {
			char c = seq.charAt(i);
			sb.append(switch (c) {
				case 'A' -> 'T';
				case 'C' -> 'G';
				case 'G' -> 'C';
				case 'T' -> 'A';
				default -> c;
			});
		}
		return sb.toString();
	}

	@Test
	public void testOvertrimPrevention() {
		// This test checks for "overtrimming" by placing a sequence with weak homology
		// to the adapter *before* the actual adapter sequence.
		// A naive or overly-sensitive algorithm might clip at the weak homology site.
		// Trimmomatic should identify the best match and clip at the correct location.

		Logger logger = mock(Logger.class);
		// Use a relatively low simple clip threshold to make a spurious match more likely.
		// minSequenceLikelihood = 5. A perfect 9bp match would score ~5.4.
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 5, 8, false);

		// A real adapter sequence
		String adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"; // 33bp
		trimmer.addClippingSeq(trimmer.new IlluminaLongClippingSeq(adapter), true, false);

		// A sequence with weak homology to the start of the adapter.
		// A 7-base perfect match: "AGATCGG". Score: 7 * 0.6 = 4.2. This is below the threshold of 5.
		String weakHomology = "AGATCGG"; // 7bp perfect match, score < 5

		// Construct the read:
		// [ 50bp DNA ] + [ 7bp weak homology ] + [ 43bp DNA ] + [ 33bp real adapter ] + [ 17bp DNA ]
		// Total length = 50 + 7 + 43 + 33 + 17 = 150
		// The correct clip point is at base 100 (50+7+43).
		String part1 = "TTCCTTGCCGACGGGCGGTGTGTACAAAGGGCAGGGACTTAATCAACGCA"; // 50bp
		String part2 = "AGCCTTGCGACTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"; // 43bp
		String part3 = "TTTTTTTTTTTTTTTTT"; // 17bp

		String seq = part1 + weakHomology + part2 + adapter + part3;
		String qual = makeQuality('I', 150); // High quality throughout

		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });

		assertNotNull(result[0]);

		// The trimmer should ignore the weak match at base 50 and correctly clip at base 100.
		int expectedLength = part1.length() + weakHomology.length() + part2.length();
		assertEquals(expectedLength, result[0].getSequence().length(), "Trimmer should not over-trim at a weak homology site.");
		
		String expectedSequence = part1 + weakHomology + part2;
		assertEquals(expectedSequence, result[0].getSequence());
	}

	@Test
	public void testOvertrimPreventionPaired() {
		Logger logger = mock(Logger.class);
		// Palindrome mode: seedMaxMiss=2, minPalindromeLikelihood=30, minSequenceLikelihood=10, minPrefix=1, keepBoth=false
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);

		String fwdAdapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
		String revAdapter = reverseComplement(fwdAdapter);
		trimmer.addPrefixPair(fwdAdapter, revAdapter);

		// Construct a realistic insert (100bp)
		// We want to place a "decoy" sequence inside the insert that looks like the adapter start.
		// Adapter start: AGATCGGAAG...
		// Decoy: AGATCGG (7bp match)
		
		String prefix = "TTCCTTGCCGACGGGCGGTGTGTACAAAGGGCAGGGACTTAATCAACGCAAGCCTTGCGACTAG"; // 64bp
		String decoy = "AGATCGG"; // 7bp
		String suffix = "CTAGCTAGCTAGCTAGCTAGCTAGCTAGC"; // 29bp
		
		// Total insert = 64 + 7 + 29 = 100bp.
		String insert = prefix + decoy + suffix;
		
		// Read 1: Insert + FwdAdapter + Garbage
		String seq1 = insert + fwdAdapter; // 100 + 33 = 133bp
		
		// Read 2: RC(Insert) + RevAdapter + Garbage
		String rcInsert = reverseComplement(insert);
		String seq2 = rcInsert + revAdapter; // 133bp

		FastqRecord rec1 = new FastqRecord("r1", seq1, "", makeQuality('I', seq1.length()), 33);
		FastqRecord rec2 = new FastqRecord("r2", seq2, "", makeQuality('I', seq2.length()), 33);

		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec1, rec2 });
		
		assertNotNull(result);
		assertEquals(100, result[0].getSequence().length());
		assertEquals(insert, result[0].getSequence());
	}
}