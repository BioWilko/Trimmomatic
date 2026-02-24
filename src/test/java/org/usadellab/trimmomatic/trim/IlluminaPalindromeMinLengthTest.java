package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaPalindromeMinLengthTest {

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
			switch (c) {
			case 'A': sb.append('T'); break;
			case 'C': sb.append('G'); break;
			case 'G': sb.append('C'); break;
			case 'T': sb.append('A'); break;
			default: sb.append(c); break;
			}
		}
		return sb.toString();
	}

	@Test
	public void testMinAdapterLength() {
		Logger logger = mock(Logger.class);
		
		String fwdAdapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
		String revAdapter = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";

		// 150bp read.
		// Construct a case with a LONG insert, meaning a SHORT adapter read-through.
		// Insert: 145bp. Adapter part: 5bp.
		
		String fragment = makeSequence(145);
		String rcFragment = reverseComplement(fragment);
		
		String seq1 = fragment + fwdAdapter.substring(0, 5);
		String seq2 = rcFragment + revAdapter.substring(0, 5);
		
		FastqRecord rec1 = new FastqRecord("r1", seq1, "", makeQuality('I', 150), 33);
		FastqRecord rec2 = new FastqRecord("r2", seq2, "", makeQuality('I', 150), 33);

		// Case 1: Default minPrefix (8).
		// The adapter match is 5bp. 5 < 8. Should NOT clip.
		IlluminaClippingTrimmer trimmerDefault = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
		trimmerDefault.addPrefixPair(fwdAdapter, revAdapter);
		
		FastqRecord[] resDefault = trimmerDefault.processRecords(new FastqRecord[] { rec1, rec2 });
		assertEquals(150, resDefault[0].getSequence().length());
		
		// Case 2: Low minPrefix (1).
		// The adapter match is 5bp. 5 >= 1. Should clip.
		IlluminaClippingTrimmer trimmerLow = new IlluminaClippingTrimmer(logger, 2, 30, 10, 1, false);
		trimmerLow.addPrefixPair(fwdAdapter, revAdapter);
		
		FastqRecord[] resLow = trimmerLow.processRecords(new FastqRecord[] { rec1, rec2 });
		assertEquals(145, resLow[0].getSequence().length());
		
		// Case 3: High minPrefix (20).
		// The adapter match is 5bp. 5 < 20. Should NOT clip.
		IlluminaClippingTrimmer trimmerHigh = new IlluminaClippingTrimmer(logger, 2, 30, 10, 20, false);
		trimmerHigh.addPrefixPair(fwdAdapter, revAdapter);
		
		FastqRecord[] resHigh = trimmerHigh.processRecords(new FastqRecord[] { rec1, rec2 });
		assertEquals(150, resHigh[0].getSequence().length());
	}
}