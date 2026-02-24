package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerAdapterDimerTest {

	private String makeQuality(char c, int len) {
		StringBuilder sb = new StringBuilder(len);
		for (int i = 0; i < len; i++)
			sb.append(c);
		return sb.toString();
	}

	@Test
	public void testAdapterDimer() {
		Logger logger = mock(Logger.class);
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
		
		String fwdAdapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
		String revAdapter = "TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"; // RC of fwdAdapter
		trimmer.addPrefixPair(fwdAdapter, revAdapter);
		
		// Adapter Dimer: Insert size is 0.
		// Read 1 consists of the Forward Adapter (plus potential garbage/polymerase run-off)
		// Read 2 consists of the Reverse Adapter
		
		String seq1 = fwdAdapter + "TTTTTTTT"; // Adapter + garbage
		String seq2 = revAdapter + "AAAAAAAA"; // Adapter + garbage
		
		FastqRecord rec1 = new FastqRecord("r1", seq1, "", makeQuality('I', seq1.length()), 33);
		FastqRecord rec2 = new FastqRecord("r2", seq2, "", makeQuality('I', seq2.length()), 33);
		
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec1, rec2 });
		
		assertNotNull(result);
		// Forward read should be clipped to 0 length (0 insert) -> Dropped (null)
		assertNull(result[0]);
		
		// Reverse read is dropped by default (null)
		assertNull(result[1]);
		
		// If we keep both:
		trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, true);
		trimmer.addPrefixPair(fwdAdapter, revAdapter);
		result = trimmer.processRecords(new FastqRecord[] { rec1, rec2 });
		assertNull(result[0]);
		assertNull(result[1]);
	}
}