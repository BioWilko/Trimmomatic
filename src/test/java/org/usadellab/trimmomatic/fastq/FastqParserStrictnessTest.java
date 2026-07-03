package org.usadellab.trimmomatic.fastq;

import static org.junit.jupiter.api.Assertions.assertThrows;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

public class FastqParserStrictnessTest {

	@Test
	public void testBlankLines(@TempDir File tempDir) throws IOException {
		File fastqFile = new File(tempDir, "blank_lines.fastq");
		try (FileWriter writer = new FileWriter(fastqFile)) {
			writer.write("@seq1\nACGT\n+\nIIII\n");
			writer.write("\n"); // Blank line
			writer.write("@seq2\nACGT\n+\nIIII\n");
		}

		FastqParser parser = new FastqParser(33);
		
		// The parser reads line by line.
		// After seq1, it reads the blank line.
		// It expects '@' at char 0.
		// An empty line has length 0, so charAt(0) throws StringIndexOutOfBoundsException.
		// Or if it contains spaces, it fails the '@' check.
		
		try {
			assertThrows(RuntimeException.class, () -> {
				parser.open(fastqFile);
				while(parser.hasNext()) {
					parser.next();
				}
			});
		} finally {
			parser.close();
		}
	}
}