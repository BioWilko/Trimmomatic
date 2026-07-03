package org.usadellab.trimmomatic.trim;

import static org.mockito.ArgumentMatchers.contains;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.verify;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.UUID;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerDuplicateTest {

	@Test
	public void testDuplicateAdapters() throws IOException {
		// Use a unique filename to avoid file locking issues on Windows if the trimmer doesn't close the file
		String uniqueFileName = "duplicates_" + UUID.randomUUID().toString() + ".fa";
		File adapterFile = new File(uniqueFileName);
		try {
			try (FileWriter writer = new FileWriter(adapterFile)) {
				// Same sequence, different names
				writer.write(">seq1\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCA\n");
				writer.write(">seq2\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCA\n");
			}

			Logger logger = mock(Logger.class);
			// Args: file:seed:palindrome:simple:min:keep
			String args = uniqueFileName + ":2:30:10:8:true";

			IlluminaClippingTrimmer.makeIlluminaClippingTrimmer(logger, args);

			// Verify that the logger warned about the duplicate
			verify(logger).warnln(contains("Skipping duplicate Clipping Sequence"));
		} finally {
			// Attempt to delete. If it fails (due to open handle), force GC and try again, or schedule for exit.
			if (!adapterFile.delete()) {
				System.gc();
				if (!adapterFile.delete()) {
					adapterFile.deleteOnExit();
				}
			}
		}
	}
}