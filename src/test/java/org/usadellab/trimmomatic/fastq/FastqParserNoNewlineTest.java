package org.usadellab.trimmomatic.fastq;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

public class FastqParserNoNewlineTest {

	@Test
	public void testNoTrailingNewline(@TempDir File tempDir) throws IOException {
		File fastqFile = new File(tempDir, "no_newline.fastq");
		try (FileWriter writer = new FileWriter(fastqFile)) {
			// Write record without final \n
			writer.write("@seq1\nACGT\n+\nIIII");
		}

		FastqParser parser = new FastqParser(33);
		parser.open(fastqFile);

		try {
			assertTrue(parser.hasNext());
			FastqRecord rec = parser.next();
			assertEquals("seq1", rec.getName());
			assertEquals("IIII", rec.getQuality());
		} finally {
			parser.close();
		}
	}
}