package org.usadellab.trimmomatic.fastq;

import static org.junit.jupiter.api.Assertions.assertThrows;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

public class FastqParserTruncatedTest {

	@Test
	public void testTruncatedAtEOF(@TempDir File tempDir) throws IOException {
		File fastqFile = new File(tempDir, "truncated.fastq");
		try (FileWriter writer = new FileWriter(fastqFile)) {
			writer.write("@seq1\nACGT\n+\nIIII\n");
			writer.write("@seq2\nACGT\n"); // Truncated record
		}

		FastqParser parser = new FastqParser(33);
		
		// open() only reads the first record when phredOffset is specified.
		// We need to iterate to hit the truncated second record.
		
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