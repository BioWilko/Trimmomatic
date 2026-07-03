package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.util.Logger;

public class TrimmerFactoryFailureTest {

	@Test
	public void testMaxInfoInvalidArgs() {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		// Invalid float for strictness
		assertThrows(NumberFormatException.class, () -> factory.makeTrimmer("MAXINFO:50:abc"));
	}

	@Test
	public void testMinLenInvalidArgs() {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		// Invalid integer
		assertThrows(NumberFormatException.class, () -> factory.makeTrimmer("MINLEN:abc"));
	}
	
	@Test
	public void testCropMissingArgs() {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		// Missing length argument
		assertThrows(NumberFormatException.class, () -> factory.makeTrimmer("CROP:"));
	}

	@Test
	public void testAvgQualInvalidArgs() {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		// Invalid integer
		assertThrows(NumberFormatException.class, () -> factory.makeTrimmer("AVGQUAL:xyz"));
	}
	
	@Test
	public void testBaseCountInvalidArgs() {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		// Invalid min count
		assertThrows(NumberFormatException.class, () -> factory.makeTrimmer("BASECOUNT:ACGT:abc"));
	}
}