package org.usadellab.trimmomatic.fastq;

import static org.junit.jupiter.api.Assertions.assertThrows;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

public class FastqParserTrailingNewlineTest {

	@Test
	public void testTrailingNewlines(@TempDir File tempDir) throws IOException {
		File fastqFile = new File(tempDir, "trailing_newlines.fastq");
		try (FileWriter writer = new FileWriter(fastqFile)) {
			writer.write("@seq1\nACGT\n+\nIIII\n");
			writer.write("\n"); // Trailing newline
		}

		FastqParser parser = new FastqParser(33);
		
		// The parser expects a valid record start (@) but finds an empty line (or EOF after empty line)
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