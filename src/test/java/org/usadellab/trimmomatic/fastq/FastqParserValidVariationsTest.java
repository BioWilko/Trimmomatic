package org.usadellab.trimmomatic.fastq;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

public class FastqParserValidVariationsTest {

	@Test
	public void testCommentLineContent(@TempDir File tempDir) throws IOException {
		File fastqFile = new File(tempDir, "comment.fastq");
		try (FileWriter writer = new FileWriter(fastqFile)) {
			// The '+' line can optionally contain the sequence identifier or other info
			writer.write("@seq1\n");
			writer.write("ACGT\n");
			writer.write("+seq1 description\n");
			writer.write("IIII\n");
		}

		FastqParser parser = new FastqParser(33);
		parser.open(fastqFile);

		try {
			assertTrue(parser.hasNext());
			FastqRecord rec = parser.next();
			assertEquals("seq1", rec.getName());
			assertEquals("seq1 description", rec.getComment());
		} finally {
			parser.close();
		}
	}
}