package org.usadellab.trimmomatic.fastq;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

public class FastqParserWhitespaceTest {

	@Test
	public void testTrailingWhitespace(@TempDir File tempDir) throws IOException {
		File fastqFile = new File(tempDir, "whitespace.fastq");
		
		String seq150 = new String(new char[150]).replace('\0', 'A');
		String qual150 = new String(new char[150]).replace('\0', 'I');
		
		try (FileWriter writer = new FileWriter(fastqFile)) {
			// Sequence has trailing space. Length becomes 151.
			// Quality has trailing space. Length becomes 151.
			// This is technically valid in some parsers but often an error.
			// Trimmomatic's parser usually reads the whole line.
			writer.write("@seq1\n");
			writer.write(seq150 + " \n");
			writer.write("+\n");
			writer.write(qual150 + " \n");
		}

		FastqParser parser = new FastqParser(33);
		parser.open(fastqFile);

		try {
			assertTrue(parser.hasNext());
			FastqRecord rec = parser.next();
			
			// Verify that the space is included in the sequence
			assertEquals(151, rec.getSequence().length());
			assertTrue(rec.getSequence().endsWith(" "));
		} finally {
			parser.close();
		}
	}
}