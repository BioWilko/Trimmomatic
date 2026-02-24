package org.usadellab.trimmomatic.fastq;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

public class FastqParserEdgeCasesTest {

	@Test
	public void testWindowsLineEndings(@TempDir File tempDir) throws IOException {
		File fastqFile = new File(tempDir, "windows.fastq");
		try (FileWriter writer = new FileWriter(fastqFile)) {
			// Explicitly write \r\n
			writer.write("@seq1\r\nACGT\r\n+\r\nIIII\r\n");
		}

		FastqParser parser = new FastqParser(33);
		parser.open(fastqFile);

		try {
			assertTrue(parser.hasNext());
			FastqRecord rec = parser.next();
			
			// Name should not contain \r
			assertEquals("seq1", rec.getName());
			// Sequence should not contain \r
			assertEquals("ACGT", rec.getSequence());
			// Quality should not contain \r
			assertEquals("IIII", rec.getQuality());
		} finally {
			parser.close();
		}
	}
}